#include "dynamicset.h"
#include "spectralfitter.h"

#include "../wsclean/cachedimageset.h"

void DynamicSet::LoadAndAverage(CachedImageSet& imageSet)
{
	for(size_t i=0; i!=_images.size(); ++i)
		assign(_images[i], 0.0);
	
	ImageBufferAllocator::Ptr scratch;
	_allocator.Allocate(_imageSize, scratch);
	
	ao::uvector<size_t> weights(_images.size(), 0.0);
	size_t imgIndex = 0;
	for(size_t sqIndex=0; sqIndex!=_imagingTable.SquaredGroupCount(); ++sqIndex)
	{
		size_t imgIndexForChannel = imgIndex;
		ImagingTable subTable = _imagingTable.GetSquaredGroup(sqIndex);
		for(size_t eIndex=0; eIndex!=subTable.EntryCount(); ++eIndex)
		{
			const ImagingTableEntry& e = subTable[eIndex];
			for(size_t i=0; i!=e.imageCount; ++i)
			{
				imageSet.Load(scratch.data(), e.polarization, e.outputChannelIndex, i==1);
				add(_images[imgIndex], scratch.data());
				weights[imgIndex]++;
				++imgIndex;
			}
		}
		size_t thisChannelIndex = (sqIndex*_channelsInDeconvolution)/_imagingTable.SquaredGroupCount();
		size_t nextChannelIndex = ((sqIndex+1)*_channelsInDeconvolution)/_imagingTable.SquaredGroupCount();
		if(thisChannelIndex == nextChannelIndex)
			imgIndex = imgIndexForChannel;
	}
	
	for(size_t i=0; i!=_images.size(); ++i)
		multiply(_images[i], 1.0/double(weights[i]));
}

void DynamicSet::LoadAndAveragePSFs(CachedImageSet& psfSet, vector<ao::uvector<double>>& psfImages, PolarizationEnum psfPolarization)
{
	for(size_t chIndex=0; chIndex!=_channelsInDeconvolution; ++chIndex)
		psfImages[chIndex].assign(_imageSize, 0.0);
	
	ImageBufferAllocator::Ptr scratch;
	_allocator.Allocate(_imageSize, scratch);
	
	ao::uvector<size_t> weights(_channelsInDeconvolution, 0.0);
	for(size_t sqIndex=0; sqIndex!=_imagingTable.SquaredGroupCount(); ++sqIndex)
	{
		size_t chIndex = (sqIndex*_channelsInDeconvolution)/_imagingTable.SquaredGroupCount();
		ImagingTable subTable = _imagingTable.GetSquaredGroup(sqIndex);
		const ImagingTableEntry& e = subTable.Front();
		psfSet.Load(scratch.data(), psfPolarization, e.outputChannelIndex, 0);
		add(psfImages[chIndex].data(), scratch.data());
		weights[chIndex]++;
	}
	
	for(size_t chIndex=0; chIndex!=ChannelsInDeconvolution(); ++chIndex)
		multiply(psfImages[chIndex].data(), 1.0/double(weights[chIndex]));
}

void DynamicSet::InterpolateAndStore(CachedImageSet& imageSet, const SpectralFitter& fitter)
{
	if(_channelsInDeconvolution == _imagingTable.SquaredGroupCount())
	{
		directStore(imageSet);
	}
	else {
		std::cout << "Interpolating from " << _channelsInDeconvolution << " to " << _imagingTable.SquaredGroupCount() << " channels...\n";
		
		// TODO SpectralFitter should also do the interpolation of images; here we should
		// just unpack the data structure
		
		// The following loop will make an 'image' with at each pixel
		// the terms of the fit. By doing this first, it is not necessary
		// to have all channel images in memory at the same time.
		// TODO: this assumes that polarizations are not joined!
		size_t nTerms = fitter.NTerms();
		ao::uvector<double> termsImage(_imageSize * nTerms);
		ao::uvector<double> spectralPixel(_channelsInDeconvolution);
		ao::uvector<double> termsPixel(nTerms);
		for(size_t px=0; px!=_imageSize; ++px)
		{
			double isZero = true;
			for(size_t s=0; s!=_images.size(); ++s)
			{
				double value = _images[s][px];
				spectralPixel[s] = value;
				isZero = isZero && (value == 0.0);
			}
			double* termsPtr = &termsImage[px*nTerms];
			// Skip fitting if it is zero; most of model images will be zero, so this can
			// save a lot of time.
			if(isZero) {
				for(double* p=termsPtr; p!=termsPtr+nTerms; ++p)
					*p = 0.0;
			}
			else {
				fitter.Fit(termsPixel, spectralPixel.data());
				for(size_t i=0; i!=nTerms; ++i)
					termsPtr[i] = termsPixel[i];
			}
		}
		
		// Now that we know the fit for each pixel, evaluate the function for each
		// pixel of each output channel.
		ImageBufferAllocator::Ptr scratch;
		_allocator.Allocate(_imageSize, scratch);
		size_t imgIndex = 0;
		for(size_t eIndex=0; eIndex!=_imagingTable.EntryCount(); ++eIndex)
		{
			const ImagingTableEntry& e = _imagingTable[eIndex];
			double freq = e.CentralFrequency();
			for(size_t px=0; px!=_imageSize; ++px)
			{
				const double* termsPtr = &termsImage[px*nTerms];
				for(size_t i=0; i!=nTerms; ++i)
					termsPixel[i] = termsPtr[i];
				scratch[px] = fitter.Evaluate(termsPixel, freq);
			}
			
			imageSet.Store(scratch.data(), e.polarization, e.outputChannelIndex, false);
			++imgIndex;
		}
	}
}

void DynamicSet::AssignAndStore(CachedImageSet& imageSet)
{
	if(_channelsInDeconvolution == _imagingTable.SquaredGroupCount())
	{
		directStore(imageSet);
	}
	else {
		std::cout << "Assigning from " << _channelsInDeconvolution << " to " << _imagingTable.SquaredGroupCount() << " channels...\n";
		size_t imgIndex = 0;
		for(size_t sqIndex=0; sqIndex!=_imagingTable.SquaredGroupCount(); ++sqIndex)
		{
			size_t imgIndexForChannel = imgIndex;
			ImagingTable subTable = _imagingTable.GetSquaredGroup(sqIndex);
			for(size_t eIndex=0; eIndex!=subTable.EntryCount(); ++eIndex)
			{
				const ImagingTableEntry& e = subTable[eIndex];
				for(size_t i=0; i!=e.imageCount; ++i)
				{
					imageSet.Store(_images[imgIndex], e.polarization, e.outputChannelIndex, i==1);
					++imgIndex;
				}
			}
			size_t thisChannelIndex = (sqIndex*_channelsInDeconvolution)/_imagingTable.SquaredGroupCount();
			size_t nextChannelIndex = ((sqIndex+1)*_channelsInDeconvolution)/_imagingTable.SquaredGroupCount();
			if(thisChannelIndex == nextChannelIndex)
				imgIndex = imgIndexForChannel;
		}
	}
}

void DynamicSet::directStore(CachedImageSet& imageSet)
{
	size_t imgIndex = 0;
	for(size_t i=0; i!=_imagingTable.EntryCount(); ++i)
	{
		const ImagingTableEntry& e = _imagingTable[i];
		for(size_t i=0; i!=e.imageCount; ++i)
		{
			imageSet.Store(_images[imgIndex], e.polarization, e.outputChannelIndex, i==1);
			++imgIndex;
		}
	}
}

void DynamicSet::GetSquareIntegrated(double* dest, double* scratch) const
{
	for(size_t sqIndex = 0; sqIndex!=_channelsInDeconvolution; ++sqIndex)
	{
		ImagingTable subTable = _imagingTable.GetSquaredGroup(sqIndex);
		if(subTable.EntryCount() == 1)
		{
			const ImagingTableEntry& entry = subTable[0];
			size_t imageIndex = _tableIndexToImageIndex.find(entry.index)->second;
			assign(scratch, _images[imageIndex]);
		}
		else {
			for(size_t eIndex = 0; eIndex!=subTable.EntryCount(); ++eIndex)
			{
				const ImagingTableEntry& entry = subTable[eIndex];
				size_t imageIndex = _tableIndexToImageIndex.find(entry.index)->second;
				if(eIndex == 0)
				{
					assign(scratch, _images[0]);
					square(scratch);
				}
				else {
					addSquared(scratch, _images[imageIndex]);
				}
			}
			squareRoot(scratch);
		}
		
		if(sqIndex == 0)
			assign(dest, scratch);
		else
			add(dest, scratch);
	}
	if(_channelsInDeconvolution > 0.0)
		multiply(dest, 1.0/_channelsInDeconvolution);
	else
		assign(dest, 0.0);
}

void DynamicSet::GetLinearIntegrated(double* dest) const
{
	size_t addIndex = 0;
	for(size_t sqIndex = 0; sqIndex!=_channelsInDeconvolution; ++sqIndex)
	{
		ImagingTable subTable = _imagingTable.GetSquaredGroup(sqIndex);
		for(size_t eIndex = 0; eIndex!=subTable.EntryCount(); ++eIndex)
		{
			const ImagingTableEntry& entry = subTable[eIndex];
			size_t imageIndex = _tableIndexToImageIndex.find(entry.index)->second;
			if(addIndex == 0)
				assign(dest, _images[imageIndex]);
			else
				add(dest, _images[imageIndex]);
			++addIndex;
		}
	}
	if(_channelsInDeconvolution > 0)
		multiply(dest, 1.0/double(_channelsInDeconvolution));
	else
		assign(dest, 0.0);
}
