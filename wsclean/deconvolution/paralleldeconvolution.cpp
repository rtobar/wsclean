#include "paralleldeconvolution.h"

#include "../multiscale/multiscalealgorithm.h"

#include "../wsclean/imagefilename.h"
#include "../wsclean/primarybeam.h"
#include "../wsclean/wscleansettings.h"

#include "../units/fluxdensity.h"

#include "../image.h"

#include "subdivision.h"

#include "../aocommon/parallelfor.h"

void ParallelDeconvolution::SetAlgorithm(std::unique_ptr<class DeconvolutionAlgorithm> algorithm)
{
	if(_settings.parallelDeconvolutionMaxSize == 0)
	{
		_algorithms.resize(1);
		_algorithms.front() = std::move(algorithm);
	}
	else {
		const size_t
			width = _settings.trimmedImageWidth,
			height = _settings.trimmedImageHeight;
		size_t
			maxSubImageSize = _settings.parallelDeconvolutionMaxSize;
		_horImages = (width+maxSubImageSize-1) / maxSubImageSize,
		_verImages = (height+maxSubImageSize-1) / maxSubImageSize;
		_algorithms.resize(_horImages * _verImages);
		_algorithms.front() = std::move(algorithm);
		size_t threadsPerAlg = (System::ProcessorCount()+_algorithms.size()-1)
			/ _algorithms.size();
		_algorithms.front()->SetThreadCount(threadsPerAlg);
		Logger::Debug << "Parallel cleaning will use " << _algorithms.size() << " subimages.\n";
		for(size_t i=1; i!=_algorithms.size(); ++i)
			_algorithms[i] = _algorithms.front()->Clone();
	}
}

void ParallelDeconvolution::SetRMSFactorImage(Image&& image)
{
	_algorithms.front()->SetRMSFactorImage(std::move(image));
	// TODO cut in pieces
}

void ParallelDeconvolution::SetThreshold(double threshold)
{
	for(auto& alg : _algorithms)
		alg->SetThreshold(threshold);
}

void ParallelDeconvolution::SetAutoMaskMode(bool trackPerScaleMasks, bool usePerScaleMasks)
{
	_trackPerScaleMasks = trackPerScaleMasks;
	_usePerScaleMasks = usePerScaleMasks;
	for(auto& alg : _algorithms)
	{
		class MultiScaleAlgorithm& algorithm =
			static_cast<class MultiScaleAlgorithm&>(*alg);
		algorithm.SetAutoMaskMode(trackPerScaleMasks, usePerScaleMasks);
	}
}

void ParallelDeconvolution::SetCleanMask(const bool* mask)
{
	_mask = mask;
	if(_algorithms.size() == 1)
		_algorithms.front()->SetCleanMask(mask);
}

void ParallelDeconvolution::runSubImage(SubImage& subImg, ImageSet& dataImage, ImageSet& modelImage, const ao::uvector<const double*>& psfImages, double majorIterThreshold, bool findPeakOnly, std::mutex* mutex)
{
	const size_t
		width = _settings.trimmedImageWidth,
		height = _settings.trimmedImageHeight;
	
	std::unique_ptr<ImageSet> subModel, subData;
	{
		std::unique_lock<std::mutex> lock(*mutex);
		subModel = modelImage.Trim(
			subImg.x, subImg.y, subImg.x+subImg.width, subImg.y+subImg.height, width);
		subData = dataImage.Trim(
			subImg.x, subImg.y, subImg.x+subImg.width, subImg.y+subImg.height, width);
	}
	
	// Construct the smaller psfs
	std::vector<Image> subPsfs(psfImages.size());
	ao::uvector<const double*> subPsfVector(psfImages.size());
	for(size_t i=0; i!=psfImages.size(); ++i)
	{
		subPsfs[i] = Image(subImg.width, subImg.height, *_allocator);
		Image::Trim(subPsfs[i].data(), subImg.width, subImg.height, psfImages[i], width, height);
		subPsfVector[i] = subPsfs[i].data();
	}
	_algorithms[subImg.index]->SetCleanMask(subImg.mask.data());
	
	size_t maxNIter = _algorithms[subImg.index]->MaxNIter();
	if(findPeakOnly)
		_algorithms[subImg.index]->SetMaxNIter(0);
	else
		_algorithms[subImg.index]->SetMajorIterThreshold(majorIterThreshold);
	
	if(_usePerScaleMasks || _trackPerScaleMasks)
	{
		std::lock_guard<std::mutex> lock(*mutex);
		MultiScaleAlgorithm& msAlg = static_cast<class MultiScaleAlgorithm&>(*_algorithms[subImg.index]);
		if(!_scaleMasks.empty())
		{
			// During the first iteration, msAlg will not have scales/masks yet and the nr scales has also not
			// been determined yet.
			for(size_t i=0; i!=msAlg.ScaleCount(); ++i)
			{
				ao::uvector<bool>& output = msAlg.GetScaleMask(i);
				output.resize(subImg.width * subImg.height);
				if(i < _scaleMasks.size())
					Image::TrimBox(output.data(), subImg.x, subImg.y, subImg.width, subImg.height, _scaleMasks[i].data(), width, height);
			}
		}
	}
	
	subImg.peak = _algorithms[subImg.index]->ExecuteMajorIteration(
		*subData, *subModel, subPsfVector,
		subImg.width, subImg.height, subImg.reachedMajorThreshold);
	
	if(_trackPerScaleMasks)
	{
		std::lock_guard<std::mutex> lock(*mutex);
		MultiScaleAlgorithm& msAlg = static_cast<class MultiScaleAlgorithm&>(*_algorithms[subImg.index]);
		if(_scaleMasks.empty())
		{
			_scaleMasks.resize(msAlg.ScaleCount());
			for(ao::uvector<bool>& scaleMask : _scaleMasks)
				scaleMask.assign(width * height, false);
		}
		for(size_t i=0; i!=msAlg.ScaleCount(); ++i)
		{
			const ao::uvector<bool>& msMask = msAlg.GetScaleMask(i);
			if(i < _scaleMasks.size())
				Image::CopyMasked(_scaleMasks[i].data(), subImg.x, subImg.y, width, msMask.data(), subImg.width, subImg.height, subImg.mask.data());
		}
	}
	
	if(findPeakOnly)
		_algorithms[subImg.index]->SetMaxNIter(maxNIter);
	
	if(!findPeakOnly)
	{
		std::lock_guard<std::mutex> lock(*mutex);
		dataImage.CopyMasked(*subData, subImg.x, subImg.y, width, subImg.width, subImg.height, subImg.mask.data());
		modelImage.CopyMasked(*subModel, subImg.x, subImg.y, width, subImg.width, subImg.height, subImg.mask.data());
	}
}

void ParallelDeconvolution::ExecuteMajorIteration(class ImageSet& dataImage, class ImageSet& modelImage, const ao::uvector<const double*>& psfImages, bool& reachedMajorThreshold)
{
	const size_t
		width = _settings.trimmedImageWidth,
		height = _settings.trimmedImageHeight;
	if(_algorithms.size() == 1)
	{
		_algorithms.front()->ExecuteMajorIteration(dataImage, modelImage, psfImages, width, height, reachedMajorThreshold);
	}
	else {
		const size_t
			avgHSubImageSize = width / _horImages,
			avgVSubImageSize = height / _verImages;
		
		Image
			image(width, height, *_allocator),
			scratch(width, height, *_allocator),
			dividingLines(width, height, 0.0, *_allocator);
		dataImage.GetLinearIntegrated(image.data());
		
		Subdivision divisor(width, height);
				
		// Make the horizontal subdivisions (i.e. the vertical lines)
		for(size_t divNr=1; divNr!=_horImages; ++divNr)
		{
			Logger::Debug << "Vertical division " << divNr << '\n';
			size_t
				splitStart = width * divNr / _horImages - avgHSubImageSize/4,
				splitEnd = width * divNr / _horImages + avgHSubImageSize/4;
			divisor.DivideVertically(image.data(), scratch.data(), splitStart, splitEnd);
			dividingLines += scratch;
		}
		
		// Make the vertical subdivisions (horizontal lines)
		for(size_t divNr=1; divNr!=_verImages; ++divNr)
		{
			Logger::Debug << "Horizontal division " << divNr << '\n';
			size_t
				splitStart = height * divNr / _verImages - avgVSubImageSize/4,
				splitEnd = height * divNr / _verImages + avgVSubImageSize/4;
			divisor.DivideHorizontally(image.data(), scratch.data(), splitStart, splitEnd);
			dividingLines += scratch;
		}
		
		// Find the bounding boxes and clean masks for each subimage
		ao::uvector<bool>
			mask(width * height),
			visited(width * height,  false);
		std::vector<SubImage> subImages;
		for(size_t y=0; y!=_verImages; ++y)
		{
			for(size_t x=0; x!=_horImages; ++x)
			{
				size_t
					midX = x * width / _horImages + avgHSubImageSize/2,
					midY = y * height / _verImages + avgVSubImageSize/2;
				subImages.emplace_back();
				SubImage& subImage = subImages.back();
				subImage.index = subImages.size()-1;
				divisor.GetBoundingMask(dividingLines.data(), midX, midY, mask.data(), visited.data(), subImage.x, subImage.y, subImage.width, subImage.height);
				Logger::Debug << "Subimage " << subImages.size() << " at (" << subImage.x << "," << subImage.y << ") - (" << subImage.x+subImage.width << "," << subImage.y+subImage.height << ")\n";
				subImage.mask.resize(subImage.width * subImage.height);
				Image::TrimBox(subImage.mask.data(), subImage.x, subImage.y, subImage.width, subImage.height, mask.data(), width, height);
				
				// If a user mask is active, take the union of that mask with the division mask
				// (note that 'mask' is reused as a scratch space)
				if(_mask != nullptr)
				{
					Image::TrimBox(mask.data(), subImage.x, subImage.y, subImage.width, subImage.height, _mask, width, height);
					for(size_t i=0; i!=subImage.mask.size(); ++i)
						subImage.mask[i] = subImage.mask[i] && mask[i];
				}
			}
		}
		std::mutex mutex;
		
		// Find the starting peak over all subimages
		ao::ParallelFor<size_t> loop(System::ProcessorCount());
		loop.Run(0, _algorithms.size(), [&](size_t index, size_t)
		{
			runSubImage(subImages[index], dataImage, modelImage, psfImages, 0.0, true, &mutex);
		});
		double maxValue = 0.0;
		size_t indexOfMax = 0;
		for(SubImage& img : subImages)
		{
			if(img.peak > maxValue)
			{
				maxValue = img.peak;
				indexOfMax = img.index;
			}
		}
		Logger::Info << "Subimage " << (indexOfMax+1) << " has maximum peak of " <<
			FluxDensity::ToNiceString(maxValue) << ".\n";
		double mIterThreshold = maxValue * (1.0-_settings.deconvolutionMGain);
		
		// Run the deconvolution
		loop.Run(0, _algorithms.size(), [&](size_t index, size_t)
		{
			runSubImage(subImages[index], dataImage, modelImage, psfImages, mIterThreshold, false, &mutex);
		});
		
		reachedMajorThreshold = false;
		for(SubImage& img : subImages)
			reachedMajorThreshold = reachedMajorThreshold || img.reachedMajorThreshold;
	}
}

void ParallelDeconvolution::SaveSourceList(CachedImageSet& modelImages, const ImagingTable& table, long double phaseCentreRA, long double phaseCentreDec)
{
	std::string filename = _settings.prefixName + "-sources.txt";
	if(_settings.useMultiscale)
	{
		// TODO correct for parallel case
		MultiScaleAlgorithm& algorithm =
			static_cast<MultiScaleAlgorithm&>(*_algorithms.front());
		algorithm.GetComponentList().Write(filename, algorithm, _settings.pixelScaleX, _settings.pixelScaleY, phaseCentreRA, phaseCentreDec);
	}
	else {
		// TODO correct for parallel case
		_allocator->FreeUnused();
		const size_t
			w = _settings.trimmedImageWidth,
			h = _settings.trimmedImageHeight;
		ImageSet modelSet(&table, *_allocator, _settings.deconvolutionChannelCount,
			_settings.squaredJoins, _settings.linkedPolarizations, w, h);
		modelSet.LoadAndAverage(modelImages);
		ComponentList componentList(w, h, modelSet);
		componentList.WriteSingleScale(filename, *_algorithms.front(), _settings.pixelScaleX, _settings.pixelScaleY, phaseCentreRA, phaseCentreDec);
	}
}

void ParallelDeconvolution::correctChannelForPB(ComponentList& list, const ImagingTableEntry& entry) const
{
	Logger::Debug << "Correcting source list of channel " << entry.outputChannelIndex << " for beam\n";
	ImageFilename filename(entry.outputChannelIndex, entry.outputIntervalIndex);
	filename.SetPolarization(entry.polarization);
	PrimaryBeam pb(_settings);
	PrimaryBeamImageSet beam(_settings.trimmedImageWidth, _settings.trimmedImageHeight, *_allocator);
	pb.Load(beam, filename);
	list.CorrectForBeam(beam, entry.outputChannelIndex);
}

void ParallelDeconvolution::SavePBSourceList(CachedImageSet& modelImages, const ImagingTable& table, long double phaseCentreRA, long double phaseCentreDec) const
{
	// TODO make this work with subimages
	std::unique_ptr<ComponentList> list;
	const size_t
		w = _settings.trimmedImageWidth,
		h = _settings.trimmedImageHeight;
	if(_settings.useMultiscale)
	{
		MultiScaleAlgorithm& algorithm =
			static_cast<MultiScaleAlgorithm&>(*_algorithms.front());
		list.reset(new ComponentList(algorithm.GetComponentList()));
	}
	else {
		_allocator->FreeUnused();
		ImageSet modelSet(&table, *_allocator, _settings.deconvolutionChannelCount, _settings.squaredJoins, _settings.linkedPolarizations, w, h);
		modelSet.LoadAndAverage(modelImages);
		list.reset(new ComponentList(w, h, modelSet));
	}
	
	if(_settings.deconvolutionChannelCount == 0 ||
		_settings.deconvolutionChannelCount == table.SquaredGroupCount())
	{
		// No beam averaging is required
		for(size_t i=0; i!=table.SquaredGroupCount(); ++i)
		{
			const ImagingTableEntry entry = table.GetSquaredGroup(i).Front();
			correctChannelForPB(*list, entry);
		}
	}
	else {
		for(size_t ch=0; ch!=_settings.deconvolutionChannelCount; ++ch)
		{
			PrimaryBeamImageSet beamImages(w, h, *_allocator);
			Logger::Debug << "Correcting source list of channel " << ch << " for averaged beam\n";
			loadAveragePrimaryBeam(beamImages, ch, table);
			list->CorrectForBeam(beamImages, ch);
		}
	}
	
	std::string filename = _settings.prefixName + "-sources-pb.txt";
	if(_settings.useMultiscale)
	{
		MultiScaleAlgorithm& algorithm =
			static_cast<MultiScaleAlgorithm&>(*_algorithms.front());
		list->Write(filename, algorithm, _settings.pixelScaleX, _settings.pixelScaleY, phaseCentreRA, phaseCentreDec);
	}
	else {
		list->WriteSingleScale(filename, *_algorithms.front(), _settings.pixelScaleX, _settings.pixelScaleY, phaseCentreRA, phaseCentreDec);
	}
}

void ParallelDeconvolution::loadAveragePrimaryBeam(PrimaryBeamImageSet& beamImages, size_t imageIndex, const ImagingTable& table) const
{
	Logger::Debug << "Averaging beam for deconvolution channel " << imageIndex << "\n";
	
	beamImages.SetToZero();
	
	ImageBufferAllocator::Ptr scratch;
	_allocator->Allocate(_settings.trimmedImageWidth*_settings.trimmedImageHeight, scratch);
	size_t deconvolutionChannels = _settings.deconvolutionChannelCount;
	
	/// TODO : use real weights of images
	size_t count = 0;
	PrimaryBeam pb(_settings);
	for(size_t sqIndex=0; sqIndex!=table.SquaredGroupCount(); ++sqIndex)
	{
		size_t curImageIndex = (sqIndex*deconvolutionChannels)/table.SquaredGroupCount();
		if(curImageIndex == imageIndex)
		{
			const ImagingTableEntry e = table.GetSquaredGroup(sqIndex).Front();
			Logger::Debug << "Adding beam at " << e.CentralFrequency()*1e-6 << " MHz\n";
			ImageFilename filename(e.outputChannelIndex, e.outputIntervalIndex);
			
			PrimaryBeamImageSet scratch(_settings.trimmedImageWidth, _settings.trimmedImageHeight, *_allocator);
			pb.Load(scratch, filename);
			beamImages += scratch;
			
			count++;
		}
	}
	beamImages *= (1.0 / double(count));
}
