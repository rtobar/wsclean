#include "deconvolution.h"

#include "simpleclean.h"
#include "moresane.h"
#include "iuwtdeconvolution.h"
#include "genericclean.h"

#include "../multiscale/multiscalealgorithm.h"

#include "../casamaskreader.h"
#include "../fitsreader.h"
#include "../image.h"
#include "../ndppp.h"
#include "../rmsimage.h"

#include "../units/fluxdensity.h"

#include "../wsclean/imagefilename.h"
#include "../wsclean/imagingtable.h"
#include "../wsclean/primarybeam.h"
#include "../wsclean/wscleansettings.h"

Deconvolution::Deconvolution(const class WSCleanSettings& settings) :
	_settings(settings), _autoMaskIsFinished(false),
	_beamSize(0.0), _pixelScaleX(0.0), _pixelScaleY(0.0)
{
}

Deconvolution::~Deconvolution()
{
	FreeDeconvolutionAlgorithms();
}

void Deconvolution::Perform(const class ImagingTable& groupTable, bool& reachedMajorThreshold, size_t majorIterationNr)
{
	Logger::Info.Flush();
	Logger::Info << " == Cleaning (" << majorIterationNr << ") ==\n";
	
	_imageAllocator->FreeUnused();
	ImageSet
		residualSet(&groupTable, *_imageAllocator, _settings.deconvolutionChannelCount, _settings.squaredJoins, _imgWidth, _imgHeight),
		modelSet(&groupTable, *_imageAllocator, _settings.deconvolutionChannelCount, _settings.squaredJoins, _imgWidth, _imgHeight);
		
	residualSet.LoadAndAverage(*_residualImages);
	modelSet.LoadAndAverage(*_modelImages);
	
	Image integrated(_imgWidth, _imgHeight, *_imageAllocator);
	residualSet.GetLinearIntegrated(integrated.data());
	double stddev = integrated.StdDevFromMAD();
	Logger::Info << "Estimated standard deviation of background noise: " << FluxDensity::ToNiceString(stddev) << '\n';
	if(_settings.autoMask && _autoMaskIsFinished)
	{
		// When we are in the second phase of automasking, don't use
		// the RMS background anymore
		_cleanAlgorithm->SetRMSFactorImage(Image());
	}
	else {
		if(!_settings.localRMSImage.empty())
		{
			Image rmsImage(_imgWidth, _imgHeight, *_imageAllocator);
			FitsReader reader(_settings.localRMSImage);
			reader.Read(rmsImage.data());
			// Normalize the RMS image
			stddev = rmsImage.Min();
			Logger::Info << "Lowest RMS in image: " << FluxDensity::ToNiceString(stddev) << '\n';
			if(stddev <= 0.0)
				throw std::runtime_error("RMS image can only contain values > 0, but contains values <= 0.0");
			for(double& value : rmsImage)
			{
				if(value != 0.0)
					value = stddev / value;
			}
			_cleanAlgorithm->SetRMSFactorImage(std::move(rmsImage));
		}
		else if(_settings.localRMS)
		{
			Image rmsImage;
			// TODO this should use full beam parameters
			switch(_settings.localRMSMethod)
			{
				case WSCleanSettings::RMSWindow:
					RMSImage::Make(rmsImage, integrated, _settings.localRMSWindow, _beamSize, _beamSize, 0.0, _pixelScaleX, _pixelScaleY);
					break;
				case WSCleanSettings::RMSAndMinimumWindow:
					RMSImage::MakeWithNegativityLimit(rmsImage, integrated, _settings.localRMSWindow, _beamSize, _beamSize, 0.0, _pixelScaleX, _pixelScaleY);
					break;
			}
			// Normalize the RMS image relative to the threshold so that Jy remains Jy.
			stddev = rmsImage.Min();
			Logger::Info << "Lowest RMS in image: " << FluxDensity::ToNiceString(stddev) << '\n';
			for(double& value : rmsImage)
			{
				if(value != 0.0)
					value = stddev / value;
			}
			_cleanAlgorithm->SetRMSFactorImage(std::move(rmsImage));
		}
	}
	if(_settings.autoMask && !_autoMaskIsFinished)
		_cleanAlgorithm->SetThreshold(std::max(stddev * _settings.autoMaskSigma, _settings.deconvolutionThreshold));
	else if(_settings.autoDeconvolutionThreshold)
		_cleanAlgorithm->SetThreshold(std::max(stddev * _settings.autoDeconvolutionThresholdSigma, _settings.deconvolutionThreshold));
	integrated.reset();
	
	std::vector<ao::uvector<double>> psfVecs(groupTable.SquaredGroupCount());
	residualSet.LoadAndAveragePSFs(*_psfImages, psfVecs, _psfPolarization);
	
	ao::uvector<const double*> psfs(groupTable.SquaredGroupCount());
	for(size_t i=0; i!=psfVecs.size(); ++i)
		psfs[i] = psfVecs[i].data();
	
	if(_settings.useMultiscale)
	{
		MultiScaleAlgorithm& algorithm = static_cast<MultiScaleAlgorithm&>(*_cleanAlgorithm.get());
		if(_settings.autoMask)
		{
			if(_autoMaskIsFinished)
				algorithm.SetAutoMaskMode(false, true);
			else
				algorithm.SetAutoMaskMode(true, false);
		}
	}
	else {
		if(_settings.autoMask && _autoMaskIsFinished)
		{
			if(_autoMask.empty())
			{
				_autoMask.resize(_imgWidth * _imgHeight);
				for(size_t imgIndex=0; imgIndex!=modelSet.size(); ++imgIndex)
				{
					const double* image = modelSet[imgIndex];
					for(size_t i=0; i!=_imgWidth * _imgHeight; ++i)
					{
						_autoMask[i] = (image[i]==0.0) ? false : true;
					}
				}
			}
			_cleanAlgorithm->SetCleanMask(_autoMask.data());
		}
	}
		
	_cleanAlgorithm->ExecuteMajorIteration(residualSet, modelSet, psfs, _imgWidth, _imgHeight, reachedMajorThreshold);
	
	if(!reachedMajorThreshold && _settings.autoMask && !_autoMaskIsFinished)
	{
		Logger::Info << "Auto-masking threshold reached; continuing next major iteration with deeper threshold and mask.\n";
		_autoMaskIsFinished = true;
		reachedMajorThreshold = true;
	}
	
	if(_settings.majorIterationCount != 0 && majorIterationNr >= _settings.majorIterationCount)
	{
		reachedMajorThreshold = false;
		Logger::Info << "Maximum number of major iterations was reached: not continuing deconvolution.\n";
	}
	
	residualSet.AssignAndStore(*_residualImages);
	modelSet.InterpolateAndStore(*_modelImages, _cleanAlgorithm->Fitter());
}

void Deconvolution::FreeDeconvolutionAlgorithms()
{
	_cleanAlgorithm.reset();
}

void Deconvolution::InitializeDeconvolutionAlgorithm(const ImagingTable& groupTable, PolarizationEnum psfPolarization, class ImageBufferAllocator* imageAllocator, size_t imgWidth, size_t imgHeight, double pixelScaleX, double pixelScaleY, double beamSize, size_t threadCount)
{
	_imageAllocator = imageAllocator;
	_imgWidth = imgWidth;
	_imgHeight = imgHeight;
	_psfPolarization = psfPolarization;
	_beamSize = beamSize;
	_pixelScaleX = pixelScaleX;
	_pixelScaleY = pixelScaleY;
	_autoMaskIsFinished = false;
	FreeDeconvolutionAlgorithms();
	
	_summedCount = groupTable.SquaredGroupCount();
	if(_summedCount == 0)
		throw std::runtime_error("Nothing to clean");
	
	if(!std::isfinite(_beamSize))
	{
		Logger::Warn << "No proper beam size available in deconvolution!\n";
		_beamSize = 0.0;
	}
	
	ImagingTable firstSquaredGroup = groupTable.GetSquaredGroup(0);
	_squaredCount = firstSquaredGroup.EntryCount();
	_polarizations.clear();
	for(size_t p=0; p!=_squaredCount; ++p)
	{
		if(_polarizations.count(firstSquaredGroup[p].polarization) != 0)
			throw std::runtime_error("Two equal polarizations were given to the deconvolution algorithm within a single polarized group");
		else
			_polarizations.insert(firstSquaredGroup[p].polarization);
	}
	
	if(_settings.useMoreSaneDeconvolution)
	{
		_cleanAlgorithm.reset(new MoreSane(_settings.moreSaneLocation, _settings.moreSaneArgs, _settings.moreSaneSigmaLevels, _settings.prefixName, *_imageAllocator));
	}
	else if(_settings.useIUWTDeconvolution)
	{
		IUWTDeconvolution* method = new IUWTDeconvolution();
		_cleanAlgorithm.reset(method);
		method->SetUseSNRTest(_settings.iuwtSNRTest);
	}
	else if(_settings.useMultiscale)
	{
		_cleanAlgorithm.reset(new MultiScaleAlgorithm(*_imageAllocator, beamSize, pixelScaleX, pixelScaleY));
		MultiScaleAlgorithm *algorithm = static_cast<MultiScaleAlgorithm*>(_cleanAlgorithm.get());
		algorithm->SetManualScaleList(_settings.multiscaleScaleList);
		algorithm->SetMultiscaleScaleBias(_settings.multiscaleDeconvolutionScaleBias);
		algorithm->SetMultiscaleNormalizeResponse(_settings.multiscaleNormalizeResponse);
		algorithm->SetMultiscaleGain(_settings.multiscaleGain);
		algorithm->SetShape(_settings.multiscaleShapeFunction);
		algorithm->SetTrackComponents(_settings.saveSourceList);
		algorithm->SetConvolutionPadding(_settings.multiscaleConvolutionPadding);
		algorithm->SetUseFastSubMinorLoop(_settings.multiscaleFastSubMinorLoop);
	}
	else
	{
		_cleanAlgorithm.reset(new GenericClean(*_imageAllocator, _settings.useClarkOptimization));
	}
	
	_cleanAlgorithm->SetMaxNIter(_settings.deconvolutionIterationCount);
	_cleanAlgorithm->SetThreshold(_settings.deconvolutionThreshold);
	_cleanAlgorithm->SetGain(_settings.deconvolutionGain);
	_cleanAlgorithm->SetMGain(_settings.deconvolutionMGain);
	_cleanAlgorithm->SetCleanBorderRatio(_settings.deconvolutionBorderRatio);
	_cleanAlgorithm->SetAllowNegativeComponents(_settings.allowNegativeComponents);
	_cleanAlgorithm->SetStopOnNegativeComponents(_settings.stopOnNegativeComponents);
	_cleanAlgorithm->SetThreadCount(threadCount);
	_cleanAlgorithm->SetSpectralFittingMode(_settings.spectralFittingMode, _settings.spectralFittingTerms);
	
	ao::uvector<double> frequencies, weights;
	calculateDeconvolutionFrequencies(groupTable, frequencies, weights);
	_cleanAlgorithm->InitializeFrequencies(frequencies, weights);
	
	if(!_settings.fitsDeconvolutionMask.empty())
	{
		if(_cleanMask.empty())
		{
			Logger::Info << "Reading mask '" << _settings.fitsDeconvolutionMask << "'...\n";
			FitsReader maskReader(_settings.fitsDeconvolutionMask);
			if(maskReader.ImageWidth() != _imgWidth || maskReader.ImageHeight() != _imgHeight)
				throw std::runtime_error("Specified Fits file mask did not have same dimensions as output image!");
			ao::uvector<float> maskData(_imgWidth*_imgHeight);
			maskReader.Read(maskData.data());
			_cleanMask.assign(_imgWidth*_imgHeight, false);
			for(size_t i=0; i!=_imgWidth*_imgHeight; ++i)
				_cleanMask[i] = maskData[i]!=0.0;
		}
		_cleanAlgorithm->SetCleanMask(_cleanMask.data());
	} else if(!_settings.casaDeconvolutionMask.empty())
	{
		if(_cleanMask.empty())
		{
			Logger::Info << "Reading CASA mask '" << _settings.casaDeconvolutionMask << "'...\n";
			_cleanMask.assign(_imgWidth*_imgHeight, false);
			CasaMaskReader maskReader(_settings.casaDeconvolutionMask);
			if(maskReader.Width() != _imgWidth || maskReader.Height() != _imgHeight)
				throw std::runtime_error("Specified CASA mask did not have same dimensions as output image!");
			maskReader.Read(_cleanMask.data());
		}
		_cleanAlgorithm->SetCleanMask(_cleanMask.data());
	}
}

void Deconvolution::calculateDeconvolutionFrequencies(const ImagingTable& groupTable, ao::uvector<double>& frequencies, ao::uvector<double>& weights)
{
	size_t deconvolutionChannels = _settings.deconvolutionChannelCount;
	if(deconvolutionChannels == 0) deconvolutionChannels = _summedCount;
	frequencies.assign(deconvolutionChannels, 0.0);
	weights.assign(deconvolutionChannels, 0.0);
	ao::uvector<size_t> counts(deconvolutionChannels, 0);
	for(size_t i=0; i!=_summedCount; ++i)
	{
		const ImagingTableEntry& entry = groupTable.GetSquaredGroup(i)[0];
		double freq = entry.CentralFrequency();
		size_t deconvolutionChannel = i * deconvolutionChannels / _summedCount;
		frequencies[deconvolutionChannel] += freq;
		Logger::Debug << "Weight: " << entry.imageWeight << '\n';
		weights[deconvolutionChannel] += entry.imageWeight;
		counts[deconvolutionChannel]++;
	}
	for(size_t i=0; i!=deconvolutionChannels; ++i)
		frequencies[i] /= counts[i];
}

void Deconvolution::SaveSourceList(const class ImagingTable& table, long double phaseCentreRA, long double phaseCentreDec) const
{
	std::string filename = _settings.prefixName + "-sources.txt";
	if(_settings.useMultiscale)
	{
		MultiScaleAlgorithm& algorithm = static_cast<MultiScaleAlgorithm&>(*_cleanAlgorithm);
		algorithm.GetComponentList().Write(filename, algorithm, _pixelScaleX, _pixelScaleY, phaseCentreRA, phaseCentreDec);
	}
	else {
		_imageAllocator->FreeUnused();
		ImageSet modelSet(&table, *_imageAllocator, _settings.deconvolutionChannelCount, _settings.squaredJoins, _imgWidth, _imgHeight);
		modelSet.LoadAndAverage(*_modelImages);
		ComponentList componentList(_imgWidth, _imgHeight, modelSet);
		componentList.WriteSingleScale(filename, *_cleanAlgorithm, _pixelScaleX, _pixelScaleY, phaseCentreRA, phaseCentreDec);
	}
}

void Deconvolution::SavePBSourceList(const class ImagingTable& table, long double phaseCentreRA, long double phaseCentreDec) const
{
	std::unique_ptr<ComponentList> list;
	if(_settings.useMultiscale)
	{
		MultiScaleAlgorithm& algorithm = static_cast<MultiScaleAlgorithm&>(*_cleanAlgorithm);
		list.reset(new ComponentList(algorithm.GetComponentList()));
	}
	else {
		_imageAllocator->FreeUnused();
		ImageSet modelSet(&table, *_imageAllocator, _settings.deconvolutionChannelCount, _settings.squaredJoins, _imgWidth, _imgHeight);
		modelSet.LoadAndAverage(*_modelImages);
		list.reset(new ComponentList(_imgWidth, _imgHeight, modelSet));
	}
	
	for(size_t i=0; i!=table.SquaredGroupCount(); ++i)
	{
		const ImagingTableEntry entry = table.GetSquaredGroup(i).Front();
		Logger::Debug << "Correcting source list of channel " << entry.outputChannelIndex << " for beam\n";
		ImageFilename filename(entry.outputChannelIndex, entry.outputIntervalIndex);
		filename.SetPolarization(entry.polarization);
		PrimaryBeam pb(_settings);
		PrimaryBeamImageSet beam(_imgWidth, _imgHeight, *_imageAllocator);
		pb.Load(beam, filename);
		list->CorrectForBeam(beam, entry.outputChannelIndex);
	}
	
	std::string filename = _settings.prefixName + "-sources-pb.txt";
	if(_settings.useMultiscale)
	{
		MultiScaleAlgorithm& algorithm = static_cast<MultiScaleAlgorithm&>(*_cleanAlgorithm);
		list->Write(filename, algorithm, _pixelScaleX, _pixelScaleY, phaseCentreRA, phaseCentreDec);
	}
	else {
		list->WriteSingleScale(filename, *_cleanAlgorithm, _pixelScaleX, _pixelScaleY, phaseCentreRA, phaseCentreDec);
	}
}
