#include "wsclean.h"

#include "imageweightcache.h"
#include "inversionalgorithm.h"
#include "logger.h"
#include "wsmsgridder.h"
#include "primarybeam.h"
#include "imagefilename.h"

#include "../angle.h"
#include "../areaset.h"
#include "../dftpredictionalgorithm.h"
#include "../fftresampler.h"
#include "../fitswriter.h"
#include "../gaussianfitter.h"
#include "../imageweights.h"
#include "../modelrenderer.h"
#include "../msselection.h"
#include "../msproviders/contiguousms.h"
#include "../progressbar.h"
#include "../uvector.h"

#include "../lofar/lmspredicter.h"

#include "../model/areaparser.h"
#include "../model/model.h"

#include "../deconvolution/deconvolutionalgorithm.h"
#include "../deconvolution/dynamicset.h"

#include <iostream>
#include <memory>

std::string commandLine;

WSClean::WSClean() :
	_globalSelection(),
	_commandLine(),
	_inversionWatch(false), _predictingWatch(false), _deconvolutionWatch(false),
	_isFirstInversion(true), _doReorder(false),
	_currentIntervalIndex(0), _majorIterationNr(0),
	_deconvolution(_settings)
{
}

WSClean::~WSClean()
{
}

void WSClean::initFitsWriter(FitsWriter& writer)
{
	double
		ra = _inversionAlgorithm->PhaseCentreRA(),
		dec = _inversionAlgorithm->PhaseCentreDec(),
		pixelScaleX = _inversionAlgorithm->PixelSizeX(),
		pixelScaleY = _inversionAlgorithm->PixelSizeY(),
		freqHigh = _inversionAlgorithm->HighestFrequencyChannel(),
		freqLow = _inversionAlgorithm->LowestFrequencyChannel(),
		freqCentre = (freqHigh + freqLow) * 0.5,
		bandwidth = _inversionAlgorithm->BandEnd() - _inversionAlgorithm->BandStart(),
		beamSize = _inversionAlgorithm->BeamSize(),
		dateObs = _inversionAlgorithm->StartTime();
		
	writer.SetImageDimensions(_settings.trimmedImageWidth, _settings.trimmedImageHeight, ra, dec, pixelScaleX, pixelScaleY);
	writer.SetFrequency(freqCentre, bandwidth);
	writer.SetDate(dateObs);
	writer.SetPolarization(_inversionAlgorithm->Polarization());
	writer.SetOrigin("WSClean", "W-stacking imager written by Andre Offringa");
	writer.AddHistory(commandLine);
	if(_settings.manualBeamMajorSize != 0.0) {
		writer.SetBeamInfo(_settings.manualBeamMajorSize, _settings.manualBeamMinorSize, _settings.manualBeamPA);
	}
	else {
		writer.SetBeamInfo(beamSize, beamSize, 0.0);
	}
	if(_inversionAlgorithm->HasDenormalPhaseCentre())
		writer.SetPhaseCentreShift(_inversionAlgorithm->PhaseCentreDL(), _inversionAlgorithm->PhaseCentreDM());
	
	writer.SetExtraKeyword("WSCIMGWG", _inversionAlgorithm->ImageWeight());
	writer.SetExtraKeyword("WSCNWLAY", _inversionAlgorithm->WGridSize());
	writer.SetExtraKeyword("WSCDATAC", _inversionAlgorithm->DataColumnName());
	writer.SetExtraKeyword("WSCWEIGH", _inversionAlgorithm->Weighting().ToString());
	writer.SetExtraKeyword("WSCGKRNL", _inversionAlgorithm->AntialiasingKernelSize());
	if(_settings.endChannel!=0)
	{
		writer.SetExtraKeyword("WSCCHANS", _settings.startChannel);
		writer.SetExtraKeyword("WSCCHANE", _settings.endChannel);
	}
	if(_globalSelection.HasInterval())
	{
		writer.SetExtraKeyword("WSCTIMES", _globalSelection.IntervalStart());
		writer.SetExtraKeyword("WSCTIMEE", _globalSelection.IntervalEnd());
	}
	writer.SetExtraKeyword("WSCFIELD", _globalSelection.FieldId());
}

void WSClean::copyWSCleanKeywords(FitsReader& reader, FitsWriter& writer)
{
	const size_t
		N_STRKEYWORDS=2, N_DBLKEYWORDS=17;
	const char* strKeywords[N_STRKEYWORDS] =
		{ "WSCDATAC", "WSCWEIGH" };
	const char* dblKeywords[N_DBLKEYWORDS] =
		{ "WSCIMGWG", "WSCNWLAY", "WSCGKRNL", "WSCCHANS", "WSCCHANE", "WSCTIMES", "WSCTIMEE", "WSCFIELD",
			"WSCNITER", "WSCTHRES", "WSCGAIN", "WSCMGAIN", "WSCNEGCM", "WSCNEGST", "WSCSMPSF",
			"WSCMINOR", "WSCMAJOR"
		};
	for(size_t i=0; i!=N_STRKEYWORDS; ++i)
		writer.CopyStringKeywordIfExists(reader, strKeywords[i]);
	for(size_t i=0; i!=N_DBLKEYWORDS; ++i)
		writer.CopyDoubleKeywordIfExists(reader, dblKeywords[i]);
}

void WSClean::setCleanParameters(FitsWriter& writer)
{
	writer.SetExtraKeyword("WSCNITER", _settings.deconvolutionIterationCount);
	writer.SetExtraKeyword("WSCTHRES", _settings.deconvolutionThreshold);
	writer.SetExtraKeyword("WSCGAIN", _settings.deconvolutionGain);
	writer.SetExtraKeyword("WSCMGAIN", _settings.deconvolutionMGain);
	writer.SetExtraKeyword("WSCNEGCM", _settings.allowNegativeComponents);
	writer.SetExtraKeyword("WSCNEGST", _settings.stopOnNegativeComponents);
}

void WSClean::updateCleanParameters(FitsWriter& writer, size_t minorIterationNr, size_t majorIterationNr)
{
	writer.SetExtraKeyword("WSCMINOR", minorIterationNr);
	writer.SetExtraKeyword("WSCMAJOR", majorIterationNr);
}

void WSClean::multiplyImage(double factor, double* image)
{
	size_t nPix = _settings.trimmedImageWidth * _settings.trimmedImageHeight;
	for(size_t i=0; i!=nPix; ++i)
		image[i] *= factor;
}

void WSClean::imagePSF(size_t currentChannelIndex)
{
	Logger::Info.Flush();
	Logger::Info << " == Constructing PSF ==\n";
	_inversionWatch.Start();
	_inversionAlgorithm->SetDoImagePSF(true);
	_inversionAlgorithm->SetVerbose(_isFirstInversion);
	_inversionAlgorithm->Invert();
	
	size_t centralIndex = _settings.trimmedImageWidth/2 + (_settings.trimmedImageHeight/2) * _settings.trimmedImageWidth;
	if(_settings.normalizeForWeighting)
	{
		const double normFactor = 1.0/_inversionAlgorithm->ImageRealResult()[centralIndex];
		_infoPerChannel[currentChannelIndex].psfNormalizationFactor = normFactor;
		multiplyImage(normFactor, _inversionAlgorithm->ImageRealResult());
		Logger::Debug << "Normalized PSF by factor of " << normFactor << ".\n";
	}
		
	DeconvolutionAlgorithm::RemoveNaNsInPSF(_inversionAlgorithm->ImageRealResult(), _settings.trimmedImageWidth, _settings.trimmedImageHeight);
	initFitsWriter(_fitsWriter);
	_psfImages.SetFitsWriter(_fitsWriter);
	_psfImages.Store(_inversionAlgorithm->ImageRealResult(), *_settings.polarizations.begin(), currentChannelIndex, false);
	_inversionWatch.Pause();
	
	if(_settings.isUVImageSaved)
	{
		saveUVImage(_inversionAlgorithm->ImageRealResult(), *_settings.polarizations.begin(), currentChannelIndex, false, "uvpsf");
	}
	
	_isFirstInversion = false;
	
	double bMaj, bMin, bPA;
	determineBeamSize(bMaj, bMin, bPA, _inversionAlgorithm->ImageRealResult(), _inversionAlgorithm->BeamSize());
	_infoPerChannel[currentChannelIndex].theoreticBeamSize = _inversionAlgorithm->BeamSize();
	_infoPerChannel[currentChannelIndex].beamMaj = bMaj;
	_infoPerChannel[currentChannelIndex].beamMin = bMin;
	_infoPerChannel[currentChannelIndex].beamPA = bPA;
		
	if(std::isfinite(_infoPerChannel[currentChannelIndex].beamMaj))
	{
		_fitsWriter.SetBeamInfo(
			_infoPerChannel[currentChannelIndex].beamMaj,
			_infoPerChannel[currentChannelIndex].beamMin,
			_infoPerChannel[currentChannelIndex].beamPA);
	}
		
	Logger::Info << "Writing psf image... ";
	Logger::Info.Flush();
	const std::string name(ImageFilename::GetPSFPrefix(_settings, currentChannelIndex, _currentIntervalIndex) + "-psf.fits");
	_fitsWriter.Write(name, _inversionAlgorithm->ImageRealResult());
	Logger::Info << "DONE\n";
}

void WSClean::imageGridding()
{
	Logger::Info << "Writing gridding correction image... ";
	Logger::Info.Flush();
	double* gridding = _imageAllocator.Allocate(_inversionAlgorithm->ActualInversionWidth() * _inversionAlgorithm->ActualInversionHeight());
	_inversionAlgorithm->GetGriddingCorrectionImage(&gridding[0]);
	FitsWriter fitsWriter;
	initFitsWriter(fitsWriter);
	fitsWriter.SetImageDimensions(_inversionAlgorithm->ActualInversionWidth(), _inversionAlgorithm->ActualInversionHeight());
	fitsWriter.Write(_settings.prefixName + "-gridding.fits", &gridding[0]);
	_imageAllocator.Free(gridding);
	Logger::Info << "DONE\n";
}

void WSClean::imageMainFirst(PolarizationEnum polarization, size_t joinedChannelIndex)
{
	Logger::Info.Flush();
	Logger::Info << " == Constructing image ==\n";
	_inversionWatch.Start();
	if(_settings.nWLayers != 0)
		_inversionAlgorithm->SetWGridSize(_settings.nWLayers);
	else
		_inversionAlgorithm->SetNoWGridSize();
	_inversionAlgorithm->SetDoImagePSF(false);
	_inversionAlgorithm->SetDoSubtractModel(_settings.subtractModel);
	_inversionAlgorithm->SetVerbose(_isFirstInversion);
	_inversionAlgorithm->Invert();
	_inversionWatch.Pause();
	_inversionAlgorithm->SetVerbose(false);
	
	multiplyImage(_infoPerChannel[joinedChannelIndex].psfNormalizationFactor, _inversionAlgorithm->ImageRealResult());
	storeAndCombineXYandYX(_residualImages, polarization, joinedChannelIndex, false, _inversionAlgorithm->ImageRealResult());
	if(Polarization::IsComplex(polarization))
	{
		multiplyImage(_infoPerChannel[joinedChannelIndex].psfNormalizationFactor, _inversionAlgorithm->ImageImaginaryResult());
		storeAndCombineXYandYX(_residualImages, polarization, joinedChannelIndex, true, _inversionAlgorithm->ImageImaginaryResult());
	}
}

void WSClean::imageMainNonFirst(PolarizationEnum polarization, size_t joinedChannelIndex)
{
	Logger::Info.Flush();
	Logger::Info << " == Constructing image ==\n";
	_inversionWatch.Start();
	_inversionAlgorithm->SetDoSubtractModel(true);
	_inversionAlgorithm->Invert();
	_inversionWatch.Pause();
	
	multiplyImage(_infoPerChannel[joinedChannelIndex].psfNormalizationFactor, _inversionAlgorithm->ImageRealResult());
	storeAndCombineXYandYX(_residualImages, polarization, joinedChannelIndex, false, _inversionAlgorithm->ImageRealResult());
	if(Polarization::IsComplex(polarization))
	{
		multiplyImage(_infoPerChannel[joinedChannelIndex].psfNormalizationFactor, _inversionAlgorithm->ImageImaginaryResult());
		storeAndCombineXYandYX(_residualImages, polarization, joinedChannelIndex, true, _inversionAlgorithm->ImageImaginaryResult());
	}
}

void WSClean::storeAndCombineXYandYX(CachedImageSet& dest, PolarizationEnum polarization, size_t joinedChannelIndex, bool isImaginary, const double* image)
{
	if(polarization == Polarization::YX && _settings.polarizations.count(Polarization::XY)!=0)
	{
		Logger::Info << "Adding XY and YX together...\n";
		double
			*xyImage = _imageAllocator.Allocate(_settings.trimmedImageWidth*_settings.trimmedImageHeight);
		dest.Load(xyImage, Polarization::XY, joinedChannelIndex, isImaginary);
		size_t count = _settings.trimmedImageWidth*_settings.trimmedImageHeight;
		if(isImaginary)
		{
			for(size_t i=0; i!=count; ++i)
				xyImage[i] = (xyImage[i]-image[i])*0.5;
		}
		else {
			for(size_t i=0; i!=count; ++i)
				xyImage[i] = (xyImage[i]+image[i])*0.5;
		}
		dest.Store(xyImage, Polarization::XY, joinedChannelIndex, isImaginary);
		_imageAllocator.Free(xyImage);
	}
	else {
		dest.Store(image, polarization, joinedChannelIndex, isImaginary);
	}
}

void WSClean::predict(PolarizationEnum polarization, size_t joinedChannelIndex)
{
	Logger::Info.Flush();
	Logger::Info << " == Converting model image to visibilities ==\n";
	const size_t size = _settings.trimmedImageWidth*_settings.trimmedImageHeight;
	double
		*modelImageReal = _imageAllocator.Allocate(size),
		*modelImageImaginary = 0;
		
	if(polarization == Polarization::YX)
	{
		_modelImages.Load(modelImageReal, Polarization::XY, joinedChannelIndex, false);
		modelImageImaginary = _imageAllocator.Allocate(size);
		_modelImages.Load(modelImageImaginary, Polarization::XY, joinedChannelIndex, true);
		for(size_t i=0; i!=size; ++i)
			modelImageImaginary[i] = -modelImageImaginary[i];
	}
	else {
		_modelImages.Load(modelImageReal, polarization, joinedChannelIndex, false);
		if(Polarization::IsComplex(polarization))
		{
			modelImageImaginary = _imageAllocator.Allocate(size);
			_modelImages.Load(modelImageImaginary, polarization, joinedChannelIndex, true);
		}
	}
	
	_predictingWatch.Start();
	_inversionAlgorithm->SetAddToModel(false);
	if(Polarization::IsComplex(polarization))
		_inversionAlgorithm->Predict(modelImageReal, modelImageImaginary);
	else
		_inversionAlgorithm->Predict(modelImageReal);
	_predictingWatch.Pause();
	_imageAllocator.Free(modelImageReal);
	_imageAllocator.Free(modelImageImaginary);
}

void WSClean::dftPredict(const ImagingTable& squaredGroup)
{
	Logger::Info.Flush();
	Logger::Info << " == Predicting visibilities ==\n";
	const size_t size = _settings.trimmedImageWidth*_settings.trimmedImageHeight;
	double
		*modelImageReal = _imageAllocator.Allocate(size),
		*modelImageImaginary = 0;
		
	std::unique_ptr<DFTPredictionImage> image(new DFTPredictionImage(_settings.trimmedImageWidth, _settings.trimmedImageHeight, _imageAllocator));
		
	for(size_t i=0; i!=squaredGroup.EntryCount(); ++i)
	{
		const ImagingTableEntry& entry = squaredGroup[i];
		if(entry.polarization == Polarization::YX)
		{
			_modelImages.Load(modelImageReal, Polarization::XY, entry.outputChannelIndex, false);
			modelImageImaginary = _imageAllocator.Allocate(size);
			_modelImages.Load(modelImageImaginary, Polarization::XY, entry.outputChannelIndex, true);
			for(size_t i=0; i!=size; ++i)
				modelImageImaginary[i] = -modelImageImaginary[i];
			image->Add(entry.polarization, modelImageReal, modelImageImaginary);
			_imageAllocator.Free(modelImageReal);
		}
		else {
			_modelImages.Load(modelImageReal, entry.polarization, entry.outputChannelIndex, false);
			if(Polarization::IsComplex(entry.polarization))
			{
				modelImageImaginary = _imageAllocator.Allocate(size);
				_modelImages.Load(modelImageImaginary, entry.polarization, entry.outputChannelIndex, true);
				image->Add(entry.polarization, modelImageReal, modelImageImaginary);
				_imageAllocator.Free(modelImageReal);
			}
			else {
				image->Add(entry.polarization, modelImageReal);
			}
		}
	}
	_imageAllocator.Free(modelImageReal);
	
	casacore::MeasurementSet firstMS(_settings.filenames.front());
	BandData firstBand(firstMS.spectralWindow());
	DFTPredictionInput input;
	image->FindComponents(input, _inversionAlgorithm->PhaseCentreRA(), _inversionAlgorithm->PhaseCentreDec(), _settings.pixelScaleX, _settings.pixelScaleY, _inversionAlgorithm->PhaseCentreDL(), _inversionAlgorithm->PhaseCentreDM(), firstBand.ChannelCount());
	// Free the input model images
	image.reset();
	Logger::Info << "Number of components to be predicted: " << input.ComponentCount() << '\n';
	
	_predictingWatch.Start();
	
	for(size_t filenameIndex=0; filenameIndex!=_settings.filenames.size(); ++filenameIndex)
	{
		for(size_t d=0; d!=_msBands[filenameIndex].DataDescCount(); ++d)
		{
			const std::string& msName = _settings.filenames[filenameIndex];
			
			std::vector<MSProvider*> msProviders(squaredGroup.EntryCount());
			MSSelection selection(_globalSelection);
			if(!selectChannels(selection, filenameIndex, d, squaredGroup.Front()))
				continue;
			for(size_t i=0; i!=squaredGroup.EntryCount(); ++i)
			{
				const ImagingTableEntry& entry = squaredGroup[i];
				msProviders[i] = initializeMSProvider(entry, selection, filenameIndex, d);
			}
			casacore::MeasurementSet ms(msName);
			
			size_t nRow = ms.nrow();
			LMSPredicter predicter(ms, _settings.threadCount);
			predicter.SetApplyBeam(_settings.dftWithBeam);
			predicter.Input() = input;
			
			if(_settings.dftWithBeam)
			{
				Logger::Info << "Converting model to absolute values...\n";
				predicter.Input().ConvertApparentToAbsolute(ms);
			}
			
			Logger::Info << "Creating row mapping...\n";
			std::vector<size_t> msToRowId;
			msProviders[0]->MakeMSRowToRowIdMapping(msToRowId);
			
			ProgressBar progress("Predicting visibilities for " + msName);
			BandData band(ms.spectralWindow());
			predicter.Start();
			LMSPredicter::RowData row;
			ao::uvector<std::complex<float>> buffer[4];
			for(size_t p=0; p!=_settings.polarizations.size(); ++p)
				buffer[p].assign(band.ChannelCount(), std::complex<float>(0.0));
			
			while(predicter.GetNextRow(row))
			{
				// Write to MS provider(s)
				size_t polIndex = 0;
				for(std::set<PolarizationEnum>::iterator pol=_settings.polarizations.begin(); pol!=_settings.polarizations.end(); ++pol)
				{
					for(size_t ch=0; ch!=band.ChannelCount(); ++ch)
					{
						switch(*pol)
						{
							case Polarization::XX:
								buffer[polIndex][ch] = row.modelData[ch][0];
								break;
							case Polarization::XY:
								buffer[polIndex][ch] = row.modelData[ch][1];
								break;
							case Polarization::YX:
								buffer[polIndex][ch] = row.modelData[ch][2];
								break;
							case Polarization::YY:
								buffer[polIndex][ch] = row.modelData[ch][3];
								break;
							case Polarization::StokesI:
								buffer[polIndex][ch] =
									(row.modelData[ch][0] +
									row.modelData[ch][3])*0.5;
								break;
							default:
								throw std::runtime_error("Can't predict for this polarization at the moment");
						}
					}
					++polIndex;
				}
				
				boost::mutex::scoped_lock lock(predicter.IOMutex());
				for(size_t polIndex=0; polIndex!=_settings.polarizations.size(); ++polIndex)
				{
					msProviders[polIndex]->WriteModel(msToRowId[row.rowIndex], buffer[polIndex].data());
				}
				lock.unlock();
					
				predicter.FinishRow(row);
				progress.SetProgress(row.rowIndex+1, nRow);
			}
			for(std::vector<MSProvider*>::iterator provider=msProviders.begin(); provider!=msProviders.end(); ++provider)
				delete *provider;
		}
	}
	
	_predictingWatch.Pause();
}

void WSClean::initializeImageWeights(const ImagingTableEntry& entry)
{
	if(!_settings.mfsWeighting)
	{
		_imageWeightCache->Update(*_inversionAlgorithm, entry.outputChannelIndex, entry.outputTimestepIndex);
		if(_settings.isWeightImageSaved)
			_imageWeightCache->Weights().Save(_settings.prefixName+"-weights.fits");
	}
	_inversionAlgorithm->SetPrecalculatedWeightInfo(&_imageWeightCache->Weights());
}

void WSClean::initializeMFSImageWeights()
{
	Logger::Info << "Precalculating MFS weights for " << _settings.weightMode.ToString() << " weighting...\n";
	_imageWeightCache->ResetWeights();
	if(_doReorder)
	{
		for(size_t sg=0; sg!=_imagingTable.SquaredGroupCount(); ++sg)
		{
			const ImagingTable subTable = _imagingTable.GetSquaredGroup(sg);
			const ImagingTableEntry& entry = subTable.Front();
			for(size_t msIndex=0; msIndex!=_settings.filenames.size(); ++msIndex)
			{
				const ImagingTableEntry::MSInfo& ms = entry.msData[msIndex];
				for(size_t dataDescId=0; dataDescId!=_msBands[msIndex].DataDescCount(); ++dataDescId)
				{
					MSSelection partSelection(_globalSelection);
					partSelection.SetBandId(dataDescId);
					bool hasSelection = selectChannels(partSelection, msIndex, dataDescId, subTable.Front());
					if(hasSelection)
					{
						PartitionedMS msProvider(_partitionedMSHandles[msIndex], ms.bands[dataDescId].partIndex, entry.polarization, dataDescId);
						_imageWeightCache->Weights().Grid(msProvider, partSelection);
					}
				}
			}
		}
	}
	else {
		for(size_t i=0; i!=_settings.filenames.size(); ++i)
		{
			for(size_t d=0; d!=_msBands[i].DataDescCount(); ++d)
			{
				ContiguousMS msProvider(_settings.filenames[i], _settings.dataColumnName, _globalSelection, *_settings.polarizations.begin(), d, _settings.deconvolutionMGain != 1.0);
				_imageWeightCache->Weights().Grid(msProvider,  _globalSelection);
				Logger::Info << '.';
				Logger::Info.Flush();
			}
		}
	}
	_imageWeightCache->Weights().FinishGridding();
	_imageWeightCache->InitializeWeightTapers();
	if(_settings.isWeightImageSaved)
		_imageWeightCache->Weights().Save(_settings.prefixName+"-weights.fits");
}

void WSClean::prepareInversionAlgorithm(PolarizationEnum polarization)
{
	_inversionAlgorithm->SetGridMode(_settings.gridMode);
	_inversionAlgorithm->SetImageWidth(_settings.untrimmedImageWidth);
	_inversionAlgorithm->SetImageHeight(_settings.untrimmedImageHeight);
	_inversionAlgorithm->SetTrimSize(_settings.trimmedImageWidth, _settings.trimmedImageHeight);
	_inversionAlgorithm->SetNWSize(_settings.widthForNWCalculation, _settings.heightForNWCalculation);
	_inversionAlgorithm->SetPixelSizeX(_settings.pixelScaleX);
	_inversionAlgorithm->SetPixelSizeY(_settings.pixelScaleY);
	if(_settings.nWLayers != 0)
		_inversionAlgorithm->SetWGridSize(_settings.nWLayers);
	else
		_inversionAlgorithm->SetNoWGridSize();
	_inversionAlgorithm->SetAntialiasingKernelSize(_settings.antialiasingKernelSize);
	_inversionAlgorithm->SetOverSamplingFactor(_settings.overSamplingFactor);
	_inversionAlgorithm->SetPolarization(polarization);
	_inversionAlgorithm->SetIsComplex(polarization == Polarization::XY || polarization == Polarization::YX);
	_inversionAlgorithm->SetDataColumnName(_settings.dataColumnName);
	_inversionAlgorithm->SetWeighting(_settings.weightMode);
	_inversionAlgorithm->SetWLimit(_settings.wLimit/100.0);
	_inversionAlgorithm->SetSmallInversion(_settings.smallInversion);
	_inversionAlgorithm->SetNormalizeForWeighting(_settings.normalizeForWeighting);
	_inversionAlgorithm->SetVisibilityWeightingMode(_settings.visibilityWeightingMode);
}

void WSClean::validateDimensions()
{
	if(_settings.trimmedImageWidth==0 || _settings.trimmedImageHeight==0)
	{
		_settings.trimmedImageWidth = _settings.untrimmedImageWidth;
		_settings.trimmedImageHeight = _settings.untrimmedImageHeight;
	}
	else if(_settings.trimmedImageWidth > _settings.untrimmedImageWidth || _settings.trimmedImageHeight > _settings.untrimmedImageHeight)
	{
		throw std::runtime_error("Error in specified trim dimensions: at least one dimension of the trimmed image is larger than in the untrimmed image");
	}
}

void WSClean::checkPolarizations()
{
	bool hasXY = _settings.polarizations.count(Polarization::XY)!=0;
	bool hasYX = _settings.polarizations.count(Polarization::YX)!=0;
	if(_settings.joinedPolarizationCleaning)
	{
		if(_settings.polarizations.size() == 1)
			throw std::runtime_error("Joined polarization cleaning requested, but only one polarization is being imaged. Specify multiple polarizatons, or do not request to join the polarizations");
		else if(_settings.polarizations.size() != 2 && _settings.polarizations.size() != 4)
			throw std::runtime_error("Joined polarization cleaning requested, but neither 2 or 4 polarizations are imaged that are suitable for this");
	}
	else {
		if((hasXY || hasYX) && _settings.deconvolutionIterationCount !=0)
			throw std::runtime_error("You are imaging XY and/or YX polarizations and have enabled cleaning (niter!=0). This is not possible -- you have to specify '-joinpolarizations' or disable cleaning.");
	}
	if((hasXY && !hasYX) || (!hasXY && hasYX))
		throw std::runtime_error("You are imaging only one of XY or YX polarizations. This is not possible -- you have to specify both XY and YX polarizations (the output of imaging both polarizations will be the XY and imaginary XY images).");
	if(_settings.IsSpectralFittingEnabled())
	{
		if(_settings.joinedPolarizationCleaning)
			throw std::runtime_error("You have requested spectral fitting, but you are joining multiple polarizations. This is not supported. You probably want to turn off the joining of polarizations (leave out -joinpolarizations).");
		if(!_settings.joinedFrequencyCleaning)
			throw std::runtime_error("You have requested spectral fitting, but you are not joining channels. This is not possible: you probably want to turn channel joining on (add -joinchannels).");
	}
}

void WSClean::performReordering(bool isPredictMode)
{
	_partitionedMSHandles.clear();
	for(size_t i=0; i != _settings.filenames.size(); ++i)
	{
		std::vector<PartitionedMS::ChannelRange> channels;
		std::map<PolarizationEnum, size_t> nextIndex;
		for(size_t j=0; j!=_imagingTable.SquaredGroupCount(); ++j)
		{
			ImagingTable squaredGroup = _imagingTable.GetSquaredGroup(j);
			for(size_t s=0; s!=squaredGroup.EntryCount(); ++s)
			{
				ImagingTableEntry& entry =
					_imagingTable[squaredGroup[s].index];
				for(size_t d=0; d!=_msBands[i].DataDescCount(); ++d)
				{
					MSSelection selection(_globalSelection);
					if(selectChannels(selection, i, d, entry))
					{
						if(entry.polarization == *_settings.polarizations.begin())
						{
							PartitionedMS::ChannelRange r;
							r.dataDescId = d;
							r.start = selection.ChannelRangeStart();
							r.end = selection.ChannelRangeEnd();
							channels.push_back(r);
						}
						entry.msData[i].bands[d].partIndex = nextIndex[entry.polarization];
						++nextIndex[entry.polarization];
					}
				}
			}
		}
		
		bool useModel = _settings.deconvolutionMGain != 1.0 || isPredictMode || _settings.subtractModel;
		bool initialModelRequired = _settings.subtractModel;
		_partitionedMSHandles.push_back(PartitionedMS::Partition(_settings.filenames[i], channels, _globalSelection, _settings.dataColumnName, true, useModel, initialModelRequired, _settings.modelUpdateRequired, _settings.polarizations, _settings.temporaryDirectory));
	}
}

void WSClean::RunClean()
{
	// If no column specified, determine column to use
	if(_settings.dataColumnName.empty())
	{
		casacore::MeasurementSet ms(_settings.filenames.front());
		bool hasCorrected = ms.tableDesc().isColumn("CORRECTED_DATA");
		if(hasCorrected) {
			Logger::Info << "First measurement set has corrected data: tasks will be applied on the corrected data column.\n";
			_settings.dataColumnName = "CORRECTED_DATA";
		} else {
			Logger::Info << "No corrected data in first measurement set: tasks will be applied on the data column.\n";
			_settings.dataColumnName= "DATA";
		}
	}

	validateDimensions();
	
	checkPolarizations();
	
	_settings.GetMSSelection(_globalSelection);
	MSSelection fullSelection = _globalSelection;
	
	for(_currentIntervalIndex=0; _currentIntervalIndex!=_settings.intervalsOut; ++_currentIntervalIndex)
	{
		makeImagingTable();
		
		_globalSelection = selectInterval(fullSelection);
		
		_doReorder = preferReordering();
		
		if(_doReorder) performReordering(false);
		
		_infoPerChannel.assign(_settings.channelsOut, ChannelInfo());
		
		_imageWeightCache.reset(createWeightCache());
		
		if(_settings.mfsWeighting)
			initializeMFSImageWeights();
		
		for(size_t groupIndex=0; groupIndex!=_imagingTable.IndependentGroupCount(); ++groupIndex)
		{
			runIndependentGroup(_imagingTable.GetIndependentGroup(groupIndex));
		}
		
		if(_settings.channelsOut > 1)
		{
			for(std::set<PolarizationEnum>::const_iterator pol=_settings.polarizations.begin(); pol!=_settings.polarizations.end(); ++pol)
			{
				bool psfWasMade = (_settings.deconvolutionIterationCount > 0 || _settings.makePSF || _settings.makePSFOnly) && pol == _settings.polarizations.begin();
				
				if(psfWasMade)
					makeMFSImage("psf.fits", *pol, false, true);
				
				if(!(*pol == Polarization::YX && _settings.polarizations.count(Polarization::XY)!=0) && !_settings.makePSFOnly)
				{
					makeMFSImage("dirty.fits", *pol, false);
					if(_settings.deconvolutionIterationCount == 0)
						makeMFSImage("image.fits", *pol, false);
					else 
					{
						makeMFSImage("residual.fits", *pol, false);
						makeMFSImage("model.fits", *pol, false);
						renderMFSImage(*pol, false);
					}
					if(Polarization::IsComplex(*pol))
					{
						makeMFSImage("dirty.fits", *pol, true);
						if(_settings.deconvolutionIterationCount == 0)
								makeMFSImage("image.fits", *pol, true);
						else
						{
							makeMFSImage("residual.fits", *pol, true);
							makeMFSImage("model.fits", *pol, true);
							renderMFSImage(*pol, true);
						}
					}
				}
			}
		}
	}
}

ImageWeightCache* WSClean::createWeightCache()
{
	ImageWeightCache* cache = new ImageWeightCache(_settings.weightMode, _settings.untrimmedImageWidth, _settings.untrimmedImageHeight, _settings.pixelScaleX, _settings.pixelScaleY, _settings.minUVInLambda, _settings.maxUVInLambda, _settings.rankFilterLevel, _settings.rankFilterSize);
	cache->SetTaperInfo(_settings.gaussianTaperBeamSize, _settings.tukeyTaperInLambda, _settings.tukeyInnerTaperInLambda, _settings.edgeTaperInLambda, _settings.edgeTukeyTaperInLambda);
	return cache;
}

void WSClean::RunPredict()
{
	if(_settings.joinedFrequencyCleaning)
		throw std::runtime_error("Joined frequency cleaning specified for prediction: prediction doesn't clean, parameter invalid");
	if(_settings.joinedPolarizationCleaning)
		throw std::runtime_error("Joined polarization cleaning specified for prediction: prediction doesn't clean, parameter invalid");
	
	_settings.dataColumnName = "DATA";
	
	validateDimensions();
	
	checkPolarizations();
	
	_settings.GetMSSelection(_globalSelection);
	MSSelection fullSelection = _globalSelection;
	
	for(_currentIntervalIndex=0; _currentIntervalIndex!=_settings.intervalsOut; ++_currentIntervalIndex)
	{
		makeImagingTable();
		
		if(_settings.predictionChannels != 0)
		{
			// TODO
		}
		
		_globalSelection = selectInterval(fullSelection);
		
		_doReorder = preferReordering();
		
		if(_doReorder) performReordering(true);
		
		_imageWeightCache.reset(createWeightCache());
		
		for(size_t groupIndex=0; groupIndex!=_imagingTable.SquaredGroupCount(); ++groupIndex)
		{
			predictGroup(_imagingTable.GetSquaredGroup(groupIndex));
		}
	}
}

bool WSClean::selectChannels(MSSelection& selection, size_t msIndex, size_t dataDescId, const ImagingTableEntry& entry)
{
	const BandData& band = _msBands[msIndex][dataDescId];
	double firstCh = band.ChannelFrequency(0), lastCh = band.ChannelFrequency(band.ChannelCount()-1);
	// Some mses have decreasing (i.e. reversed) channel frequencies in them
	bool isReversed = false;
	if(firstCh > lastCh) {
		std::swap(firstCh, lastCh);
		isReversed = true;
	}
	if(band.ChannelCount()!=0 && entry.lowestFrequency <= lastCh && entry.highestFrequency >= firstCh)
	{
		size_t newStart, newEnd;
		if(isReversed)
		{
			BandData::const_reverse_iterator lowPtr, highPtr;
			lowPtr = std::lower_bound(band.rbegin(), band.rend(), entry.lowestFrequency);
			highPtr = std::lower_bound(lowPtr, band.rend(), entry.highestFrequency);
			
			if(highPtr == band.rend())
				--highPtr;
			newStart = band.ChannelCount() - 1 - (lowPtr - band.rbegin());
			newEnd = band.ChannelCount() - (highPtr - band.rbegin());
		}
		else {
			const double *lowPtr, *highPtr;
			lowPtr = std::lower_bound(band.begin(), band.end(), entry.lowestFrequency);
			highPtr = std::lower_bound(lowPtr, band.end(), entry.highestFrequency);
			
			if(highPtr == band.end())
				--highPtr;
			newStart = lowPtr - band.begin();
			newEnd = highPtr - band.begin() + 1;
		}
			
		selection.SetChannelRange(newStart, newEnd);
		return true;
	}
	else {
		return false;
	}
}

void WSClean::runIndependentGroup(const ImagingTable& groupTable)
{
	_inversionAlgorithm.reset(new WSMSGridder(&_imageAllocator, _settings.threadCount, _settings.memFraction, _settings.absMemLimit));
	
	_modelImages.Initialize(_fitsWriter, _settings.polarizations.size(), _settings.channelsOut, _settings.prefixName + "-model", _imageAllocator);
	_residualImages.Initialize(_fitsWriter, _settings.polarizations.size(), _settings.channelsOut, _settings.prefixName + "-residual", _imageAllocator);
	if(groupTable.Front().polarization == *_settings.polarizations.begin())
		_psfImages.Initialize(_fitsWriter, 1, groupTable.SquaredGroupCount(), _settings.prefixName + "-psf", _imageAllocator);
	
	const std::string rootPrefix = _settings.prefixName;
		
	for(size_t joinedIndex=0; joinedIndex!=groupTable.EntryCount(); ++joinedIndex)
	{
		const ImagingTableEntry& entry = groupTable[joinedIndex];
		runFirstInversion(entry);
	}
	
	_deconvolution.InitializeDeconvolutionAlgorithm(groupTable, *_settings.polarizations.begin(), &_imageAllocator, _settings.trimmedImageWidth, _settings.trimmedImageHeight, _settings.pixelScaleX, _settings.pixelScaleY, _settings.channelsOut, _inversionAlgorithm->BeamSize(), _settings.threadCount);

	initFitsWriter(_fitsWriter);
	setCleanParameters(_fitsWriter);
	updateCleanParameters(_fitsWriter, 0, 0);
		
	if(!_settings.makePSFOnly)
	{
		if(_settings.deconvolutionIterationCount > 0)
		{
			// Start major cleaning loop
			_majorIterationNr = 1;
			bool reachedMajorThreshold = false;
			do {
				_deconvolution.InitializeImages(_residualImages, _modelImages, _psfImages);
				_deconvolutionWatch.Start();
				_deconvolution.Perform(groupTable, reachedMajorThreshold, _majorIterationNr);
				_deconvolutionWatch.Pause();
				
				if(_majorIterationNr == 1 && _settings.deconvolutionMGain != 1.0)
					writeFirstResidualImages(groupTable);
		
				if(!reachedMajorThreshold)
					writeModelImages(groupTable);
		
				if(_settings.deconvolutionMGain != 1.0)
				{
					for(size_t sGroupIndex=0; sGroupIndex!=groupTable.SquaredGroupCount(); ++sGroupIndex)
					{
						const ImagingTable sGroupTable = groupTable.GetSquaredGroup(sGroupIndex);
						size_t currentChannelIndex = sGroupTable.Front().outputChannelIndex;
						if(_settings.dftPrediction)
						{
							dftPredict(sGroupTable);
							for(size_t e=0; e!=sGroupTable.EntryCount(); ++e)
							{
								prepareInversionAlgorithm(sGroupTable[e].polarization);
								initializeCurMSProviders(sGroupTable[e]);
								initializeImageWeights(sGroupTable[e]);
			
								imageMainNonFirst(sGroupTable[e].polarization, currentChannelIndex);
								clearCurMSProviders();
							}
						}
						else {
							for(size_t e=0; e!=sGroupTable.EntryCount(); ++e)
							{
								prepareInversionAlgorithm(sGroupTable[e].polarization);
								initializeCurMSProviders(sGroupTable[e]);
								initializeImageWeights(sGroupTable[e]);
			
								predict(sGroupTable[e].polarization, currentChannelIndex);
								
								imageMainNonFirst(sGroupTable[e].polarization, currentChannelIndex);
								clearCurMSProviders();
							} // end of polarization loop
						}
					} // end of joined channels loop
					
					++_majorIterationNr;
				}
				
			} while(reachedMajorThreshold);
			
			Logger::Info << _majorIterationNr << " major iterations were performed.\n";
		}
		
		_inversionAlgorithm->FreeImagingData();
		
		for(size_t joinedIndex=0; joinedIndex!=groupTable.EntryCount(); ++joinedIndex)
			saveRestoredImagesForGroup(groupTable[joinedIndex]);
	}
	
	_imageAllocator.ReportStatistics();
	Logger::Info << "Inversion: " << _inversionWatch.ToString() << ", prediction: " << _predictingWatch.ToString() << ", deconvolution: " << _deconvolutionWatch.ToString() << '\n';
	
	_settings.prefixName = rootPrefix;
	
	// Needs to be destructed before image allocator, or image allocator will report error caused by leaked memory
	_inversionAlgorithm.reset();
}

void WSClean::saveRestoredImagesForGroup(const ImagingTableEntry& tableEntry)
{
	// Restore model to residual and save image
	size_t currentChannelIndex =
		tableEntry.outputChannelIndex;
	
	double
		freqLow = tableEntry.minBandFrequency,
		freqHigh = tableEntry.maxBandFrequency;
	
	PolarizationEnum curPol = tableEntry.polarization;
	for(size_t imageIter=0; imageIter!=tableEntry.imageCount; ++imageIter)
	{
		bool isImaginary = (imageIter == 1);
		double* restoredImage = _imageAllocator.Allocate(_settings.trimmedImageWidth*_settings.trimmedImageHeight);
		_residualImages.Load(restoredImage, curPol, currentChannelIndex, isImaginary);
		
		if(_settings.deconvolutionIterationCount != 0)
			writeFits("residual.fits", restoredImage, curPol, currentChannelIndex, isImaginary);
		
		if(_settings.isUVImageSaved)
			saveUVImage(restoredImage, curPol, currentChannelIndex, isImaginary, "uv");
		
		double* modelImage = _imageAllocator.Allocate(_settings.trimmedImageWidth*_settings.trimmedImageHeight);
		_modelImages.Load(modelImage, curPol, currentChannelIndex, isImaginary);
		ModelRenderer renderer(_fitsWriter.RA(), _fitsWriter.Dec(), _settings.pixelScaleX, _settings.pixelScaleY, _fitsWriter.PhaseCentreDL(), _fitsWriter.PhaseCentreDM());
		double beamMaj = _infoPerChannel[currentChannelIndex].beamMaj;
		double beamMin, beamPA;
		std::string beamStr;
		if(std::isfinite(beamMaj))
		{
			beamMin = _infoPerChannel[currentChannelIndex].beamMin;
			beamPA = _infoPerChannel[currentChannelIndex].beamPA;
			beamStr = "(beam=" + Angle::ToNiceString(beamMin) + "-" +
			Angle::ToNiceString(beamMaj) + ", PA=" +
			Angle::ToNiceString(beamPA) + ")";
		}
		else {
			beamStr = "(beam is neither fitted nor estimated -- using delta scales!)";
			beamMaj = 0.0; beamMin = 0.0; beamPA = 0.0;
		}
		if(_settings.useMultiscale || _settings.useFastMultiscale || _settings.useMoreSaneDeconvolution || _settings.useIUWTDeconvolution)
		{
			Logger::Info << "Rendering sources to restored image " + beamStr + "... ";
			Logger::Info.Flush();
			renderer.Restore(restoredImage, modelImage, _settings.trimmedImageWidth, _settings.trimmedImageHeight, beamMaj, beamMin, beamPA);
			Logger::Info << "DONE\n";
		}
		else {
			Model model;
			// A model cannot hold instrumental pols (xx/xy/yx/yy), hence always use Stokes I here
			DeconvolutionAlgorithm::GetModelFromImage(model, modelImage, _settings.trimmedImageWidth, _settings.trimmedImageHeight, _fitsWriter.RA(), _fitsWriter.Dec(), _settings.pixelScaleX, _settings.pixelScaleY, _fitsWriter.PhaseCentreDL(), _fitsWriter.PhaseCentreDM(), 0.0, _fitsWriter.Frequency(), Polarization::StokesI);
			
			if(beamMaj == beamMin) {
				Logger::Info << "Rendering " << model.SourceCount() << " circular sources to restored image " + beamStr + "... ";
				Logger::Info.Flush();
				renderer.Restore(restoredImage, _settings.trimmedImageWidth, _settings.trimmedImageHeight, model, beamMaj, freqLow, freqHigh, Polarization::StokesI);
			}
			else {
				Logger::Info << "Rendering " << model.SourceCount() << " elliptical sources to restored image " + beamStr + "... ";
				Logger::Info.Flush();
				renderer.Restore(restoredImage, _settings.trimmedImageWidth, _settings.trimmedImageHeight, model, beamMaj, beamMin, beamPA, freqLow, freqHigh, Polarization::StokesI);
			}
			Logger::Info << "DONE\n";
		}
		_imageAllocator.Free(modelImage);
		
		Logger::Info << "Writing restored image... ";
		Logger::Info.Flush();
		writeFits("image.fits", restoredImage, curPol, currentChannelIndex, isImaginary);
		Logger::Info << "DONE\n";
		_imageAllocator.Free(restoredImage);
		
		if(curPol == *_settings.polarizations.rbegin() && _settings.applyPrimaryBeam)
		{
			PrimaryBeam pb(_settings);
			initializeMSProvidersForPB(tableEntry, pb);
			initializeImageWeights(tableEntry);
			pb.SetPhaseCentre(_fitsWriter.RA(), _fitsWriter.Dec(), _fitsWriter.PhaseCentreDL(), _fitsWriter.PhaseCentreDM());
			ImageFilename imageName = ImageFilename(currentChannelIndex, _currentIntervalIndex);
			pb.MakeImage(imageName, tableEntry, currentChannelIndex, _imageWeightCache.get(), _imageAllocator);
			clearCurMSProviders();
		}
	}
}

void WSClean::writeFirstResidualImages(const ImagingTable& groupTable)
{
	Logger::Info << "Writing first iteration image(s)...\n";
	ImageBufferAllocator::Ptr ptr;
	_imageAllocator.Allocate(_settings.trimmedImageWidth*_settings.trimmedImageHeight, ptr);
	for(size_t e=0; e!=groupTable.EntryCount(); ++e)
	{
		const ImagingTableEntry& entry = groupTable[e];
		size_t ch = entry.outputChannelIndex;
		if(entry.polarization == Polarization::YX) {
			_residualImages.Load(ptr.data(), Polarization::XY, ch, true);
			writeFits("first-residual.fits", ptr.data(), Polarization::XY, ch, true);
		}
		else {
			_residualImages.Load(ptr.data(), entry.polarization, ch, false);
			writeFits("first-residual.fits", ptr.data(), entry.polarization, ch, false);
		}
	}
}

void WSClean::writeModelImages(const ImagingTable& groupTable)
{
	Logger::Info << "Writing model image...\n";
	ImageBufferAllocator::Ptr ptr;
	_imageAllocator.Allocate(_settings.trimmedImageWidth*_settings.trimmedImageHeight, ptr);
	for(size_t e=0; e!=groupTable.EntryCount(); ++e)
	{
		const ImagingTableEntry& entry = groupTable[e];
		size_t ch = entry.outputChannelIndex;
		if(entry.polarization == Polarization::YX) {
			_modelImages.Load(ptr.data(), Polarization::XY, ch, true);
			writeFits("model.fits", ptr.data(), Polarization::XY, ch, true);
		}
		else {
			_modelImages.Load(ptr.data(), entry.polarization, ch, false);
			writeFits("model.fits", ptr.data(), entry.polarization, ch, false);
		}
	}
}

void WSClean::predictGroup(const ImagingTable& imagingGroup)
{
	_inversionAlgorithm.reset(new WSMSGridder(&_imageAllocator, _settings.threadCount, _settings.memFraction, _settings.absMemLimit));
	
	_modelImages.Initialize(_fitsWriter, _settings.polarizations.size(), 1, _settings.prefixName + "-model", _imageAllocator);
	
	const std::string rootPrefix = _settings.prefixName;
		
	for(size_t e=0; e!=imagingGroup.EntryCount(); ++e)
	{
		const ImagingTableEntry& entry = imagingGroup[e];
		// load image(s) from disk
		for(size_t i=0; i!=entry.imageCount; ++i)
		{
			std::string prefix = ImageFilename::GetPrefix(_settings, entry.polarization, entry.outputChannelIndex, _currentIntervalIndex, i==1);
			FitsReader reader(prefix + "-model.fits");
			_fitsWriter = FitsWriter(reader);
			_modelImages.SetFitsWriter(_fitsWriter);
			Logger::Info << "Reading " << reader.Filename() << "...\n";
			double* buffer = _imageAllocator.Allocate(_settings.trimmedImageWidth*_settings.trimmedImageHeight);
			if(reader.ImageWidth()!=_settings.trimmedImageWidth || reader.ImageHeight()!=_settings.trimmedImageHeight)
				throw std::runtime_error("Inconsistent image size: input image did not match with specified dimensions.");
			reader.Read(buffer);
			for(size_t j=0; j!=_settings.trimmedImageWidth*_settings.trimmedImageHeight; ++j)
			{
				if(!std::isfinite(buffer[j]))
					throw std::runtime_error("The input image contains non-finite values -- can't predict from an image with non-finite values");
			}
			_modelImages.Store(buffer, entry.polarization, 0, i==1);
			_imageAllocator.Free(buffer);
		}
		
		prepareInversionAlgorithm(entry.polarization);
		initializeCurMSProviders(entry);
		initializeImageWeights(entry);

		predict(entry.polarization, 0);
		
		clearCurMSProviders();
	} // end of polarization loop
	
	_imageAllocator.ReportStatistics();
	Logger::Info << "Inversion: " << _inversionWatch.ToString() << ", prediction: " << _predictingWatch.ToString() << ", cleaning: " << _deconvolutionWatch.ToString() << '\n';
	
	_settings.prefixName = rootPrefix;
	
	// Needs to be destructed before image allocator, or image allocator will report error caused by leaked memory
	_inversionAlgorithm.reset();
}

MSProvider* WSClean::initializeMSProvider(const ImagingTableEntry& entry, const MSSelection& selection, size_t filenameIndex, size_t dataDescId)
{
	if(_doReorder)
		return new PartitionedMS(_partitionedMSHandles[filenameIndex], entry.msData[filenameIndex].bands[dataDescId].partIndex, entry.polarization, dataDescId);
	else
		return new ContiguousMS(_settings.filenames[filenameIndex], _settings.dataColumnName, selection, entry.polarization, dataDescId, _settings.deconvolutionMGain != 1.0);
}

void WSClean::initializeCurMSProviders(const ImagingTableEntry& entry)
{
	_inversionAlgorithm->ClearMeasurementSetList();
	for(size_t i=0; i != _settings.filenames.size(); ++i)
	{
		for(size_t d=0; d!=_msBands[i].DataDescCount(); ++d)
		{
			MSSelection selection(_globalSelection);
			if(selectChannels(selection, i, d, entry))
			{
				MSProvider* msProvider = initializeMSProvider(entry, selection, i, d);
				_inversionAlgorithm->AddMeasurementSet(msProvider, selection);
				_currentPolMSes.push_back(msProvider);
			}
		}
	}
}

void WSClean::initializeMSProvidersForPB(const ImagingTableEntry& entry, PrimaryBeam& pb)
{
	for(size_t i=0; i != _settings.filenames.size(); ++i)
	{
		for(size_t d=0; d!=_msBands[i].DataDescCount(); ++d)
		{
			MSSelection selection(_globalSelection);
			if(selectChannels(selection, i, d, entry))
			{
				MSProvider* msProvider = initializeMSProvider(entry, selection, i, d);
				pb.AddMS(msProvider, selection);
				_currentPolMSes.push_back(msProvider);
			}
		}
	}
}

void WSClean::clearCurMSProviders()
{
	for(std::vector<MSProvider*>::iterator i=_currentPolMSes.begin(); i != _currentPolMSes.end(); ++i)
		delete *i;
	_currentPolMSes.clear();
}

void WSClean::runFirstInversion(const ImagingTableEntry& entry)
{
	initializeCurMSProviders(entry);
	initializeImageWeights(entry);
	
	prepareInversionAlgorithm(entry.polarization);
	
	const bool firstBeforePSF = _isFirstInversion;

	bool isFirstPol = entry.polarization == *_settings.polarizations.begin();
	bool doMakePSF = _settings.deconvolutionIterationCount > 0 || _settings.makePSF || _settings.makePSFOnly;
	if(doMakePSF && isFirstPol)
		imagePSF(entry.outputChannelIndex);
	
	if(!_settings.makePSFOnly)
	{
		initFitsWriter(_fitsWriter);
		_modelImages.SetFitsWriter(_fitsWriter);
		_residualImages.SetFitsWriter(_fitsWriter);
		
		imageMainFirst(entry.polarization, entry.outputChannelIndex);
		
		// If this was the first polarization of this channel, we need to set
		// the info for this channel
		if(isFirstPol)
		{
			_infoPerChannel[entry.outputChannelIndex].weight = _inversionAlgorithm->ImageWeight();
			_infoPerChannel[entry.outputChannelIndex].bandStart = _inversionAlgorithm->BandStart();
			_infoPerChannel[entry.outputChannelIndex].bandEnd = _inversionAlgorithm->BandEnd();
			// If no PSF is made, also set the beam size. If the PSF was made, these would already be set
			// after imaging the PSF.
			if(!doMakePSF)
			{
				if(_settings.theoreticBeam) {
					_infoPerChannel[entry.outputChannelIndex].beamMaj = _inversionAlgorithm->BeamSize();
					_infoPerChannel[entry.outputChannelIndex].beamMin = _inversionAlgorithm->BeamSize();
					_infoPerChannel[entry.outputChannelIndex].beamPA = 0.0;
				}
				else if(_settings.manualBeamMajorSize != 0.0) {
					_infoPerChannel[entry.outputChannelIndex].beamMaj = _settings.manualBeamMajorSize;
					_infoPerChannel[entry.outputChannelIndex].beamMin = _settings.manualBeamMinorSize;
					_infoPerChannel[entry.outputChannelIndex].beamPA = _settings.manualBeamPA;
				}
				else {
					_infoPerChannel[entry.outputChannelIndex].beamMaj = std::numeric_limits<double>::quiet_NaN();
					_infoPerChannel[entry.outputChannelIndex].beamMin = std::numeric_limits<double>::quiet_NaN();
					_infoPerChannel[entry.outputChannelIndex].beamPA = std::numeric_limits<double>::quiet_NaN();
				}
			}
		}
		
		if(_settings.isGriddingImageSaved && firstBeforePSF && _inversionAlgorithm->HasGriddingCorrectionImage())
			imageGridding();
		
		_isFirstInversion = false;
		
		// Set model to zero: already done if this is YX of XY/YX imaging combi
		if(!(entry.polarization == Polarization::YX && _settings.polarizations.count(Polarization::XY)!=0))
		{
			double* modelImage = _imageAllocator.Allocate(_settings.trimmedImageWidth * _settings.trimmedImageHeight);
			memset(modelImage, 0, _settings.trimmedImageWidth * _settings.trimmedImageHeight * sizeof(double));
			_modelImages.Store(modelImage, entry.polarization, entry.outputChannelIndex, false);
			if(Polarization::IsComplex(entry.polarization))
				_modelImages.Store(modelImage, entry.polarization, entry.outputChannelIndex, true);
			_imageAllocator.Free(modelImage);
		}
		
		if(entry.polarization == Polarization::XY && _settings.polarizations.count(Polarization::YX)!=0)
		{ // Skip saving XY of XY/YX combi
		}
		else {
			PolarizationEnum savedPol = entry.polarization;
			if(savedPol == Polarization::YX && _settings.polarizations.count(Polarization::XY)!=0)
				savedPol = Polarization::XY;
			double* dirtyImage = _imageAllocator.Allocate(_settings.trimmedImageWidth * _settings.trimmedImageHeight);
			_residualImages.Load(dirtyImage, savedPol, entry.outputChannelIndex, false);
			Logger::Info << "Writing dirty image...\n";
			writeFits("dirty.fits", dirtyImage, savedPol, entry.outputChannelIndex, false);
			if(Polarization::IsComplex(entry.polarization))
			{
				_residualImages.Load(dirtyImage, savedPol, entry.outputChannelIndex, true);
				writeFits("dirty.fits", dirtyImage, savedPol, entry.outputChannelIndex, true);
			}
			_imageAllocator.Free(dirtyImage);
		}
	}
	
	clearCurMSProviders();
}

void WSClean::makeMFSImage(const string& suffix, PolarizationEnum pol, bool isImaginary, bool isPSF)
{
	double lowestFreq = 0.0, highestFreq = 0.0;
	const size_t size = _settings.trimmedImageWidth * _settings.trimmedImageHeight;
	ao::uvector<double> mfsImage(size, 0.0), addedImage(size), weightImage(size, 0.0);
	double weightSum = 0.0;
	FitsWriter writer;
	for(size_t ch=0; ch!=_settings.channelsOut; ++ch)
	{
		std::string prefixStr = isPSF ?
			ImageFilename::GetPSFPrefix(_settings, ch, _currentIntervalIndex) :
			ImageFilename::GetPrefix(_settings, pol, ch, _currentIntervalIndex, isImaginary);
		const std::string name(prefixStr + '-' + suffix);
		FitsReader reader(name);
		if(ch == 0)
		{
			writer = FitsWriter(reader);
			copyWSCleanKeywords(reader, writer);
			lowestFreq = reader.Frequency() - reader.Bandwidth()*0.5;
			highestFreq = reader.Frequency() + reader.Bandwidth()*0.5;
		}
		else {
			lowestFreq = std::min(lowestFreq, reader.Frequency() - reader.Bandwidth()*0.5);
			highestFreq = std::max(highestFreq, reader.Frequency() + reader.Bandwidth()*0.5);
		}
		double weight;
		if(!reader.ReadDoubleKeyIfExists("WSCIMGWG", weight))
		{
			Logger::Error << "Error: image " << name << " did not have the WSCIMGWG keyword.\n";
			weight = 0.0;
		}
		weightSum += weight;
		reader.Read(addedImage.data());
		for(size_t i=0; i!=size; ++i)
		{
			if(std::isfinite(addedImage[i])) {
				mfsImage[i] += addedImage[i] * weight;
				weightImage[i] += weight;
			}
		}
	}
	for(size_t i=0; i!=size; ++i)
		mfsImage[i] /= weightImage[i];
	
	if(isPSF)
	{
		double bMaj, bMin, bPA;
		determineBeamSize(bMaj, bMin, bPA, mfsImage.data(), _infoPerChannel.back().theoreticBeamSize);
		_infoForMFS.beamMaj = bMaj;
		_infoForMFS.beamMin = bMin;
		_infoForMFS.beamPA = bPA;
	}
	if(std::isfinite(_infoForMFS.beamMaj))
		writer.SetBeamInfo(_infoForMFS.beamMaj, _infoForMFS.beamMin, _infoForMFS.beamPA);
	else
		writer.SetNoBeamInfo();
	
	std::string mfsName(ImageFilename::GetMFSPrefix(_settings, pol, _currentIntervalIndex, isImaginary, isPSF) + '-' + suffix);
	Logger::Info << "Writing " << mfsName << "...\n";
	writer.SetFrequency((lowestFreq+highestFreq)*0.5, highestFreq-lowestFreq);
	writer.SetExtraKeyword("WSCIMGWG", weightSum);
	writer.RemoveExtraKeyword("WSCCHANS");
	writer.RemoveExtraKeyword("WSCCHANE");
	writer.Write(mfsName, mfsImage.data());
}

void WSClean::renderMFSImage(PolarizationEnum pol, bool isImaginary)
{
	const size_t size = _settings.trimmedImageWidth * _settings.trimmedImageHeight;
	
	std::string mfsPrefix(ImageFilename::GetMFSPrefix(_settings, pol, _currentIntervalIndex, isImaginary, false));
	FitsReader residualReader(mfsPrefix + "-residual.fits");
	FitsReader modelReader(mfsPrefix + "-model.fits");
	ao::uvector<double> image(size), modelImage(size);
	residualReader.Read(image.data());
	modelReader.Read(modelImage.data());
	
	ModelRenderer renderer(_fitsWriter.RA(), _fitsWriter.Dec(), _settings.pixelScaleX, _settings.pixelScaleY, _fitsWriter.PhaseCentreDL(), _fitsWriter.PhaseCentreDM());
	double beamMaj = _infoForMFS.beamMaj;
	double beamMin, beamPA;
	std::string beamStr;
	if(std::isfinite(beamMaj))
	{
		beamMin = _infoForMFS.beamMin;
		beamPA = _infoForMFS.beamPA;
		beamStr = "(beam=" + Angle::ToNiceString(beamMin) + "-" +
		Angle::ToNiceString(beamMaj) + ", PA=" +
		Angle::ToNiceString(beamPA) + ")";
	}
	else {
		beamStr = "(beam is neither fitted nor estimated -- using delta scales!)";
		beamMaj = 0.0; beamMin = 0.0; beamPA = 0.0;
	}
	Logger::Info << "Rendering sources to restored image " + beamStr + "... ";
	Logger::Info.Flush();
	renderer.Restore(image.data(), modelImage.data(), _settings.trimmedImageWidth, _settings.trimmedImageHeight, beamMaj, beamMin, beamPA);
	Logger::Info << "DONE\n";
	
	Logger::Info << "Writing " + mfsPrefix + "-image.fits...\n";
	FitsWriter imageWriter(residualReader);
	imageWriter.Write(mfsPrefix + "-image.fits", image.data());
}

void WSClean::writeFits(const string& suffix, const double* image, PolarizationEnum pol, size_t channelIndex, bool isImaginary)
{
	const double
		bandStart = _infoPerChannel[channelIndex].bandStart,
		bandEnd = _infoPerChannel[channelIndex].bandEnd,
		centreFrequency = 0.5*(bandStart+bandEnd),
		bandwidth = bandEnd-bandStart;
	const std::string name(ImageFilename::GetPrefix(_settings, pol, channelIndex, _currentIntervalIndex, isImaginary) + '-' + suffix);
	initFitsWriter(_fitsWriter);
	_fitsWriter.SetPolarization(pol);
	_fitsWriter.SetFrequency(centreFrequency, bandwidth);
	_fitsWriter.SetExtraKeyword("WSCIMGWG", _infoPerChannel[channelIndex].weight);
	_fitsWriter.SetBeamInfo(
		_infoPerChannel[channelIndex].beamMaj,
		_infoPerChannel[channelIndex].beamMin,
		_infoPerChannel[channelIndex].beamPA);
	size_t polIndex;
	if(_settings.joinedPolarizationCleaning)
		polIndex = 0;
	else
		Polarization::TypeToIndex(pol, _settings.polarizations, polIndex);
	setCleanParameters(_fitsWriter);
	if(_deconvolution.IsInitialized())
		updateCleanParameters(_fitsWriter, _deconvolution.GetAlgorithm().IterationNumber(), _majorIterationNr);
	_fitsWriter.Write(name, image);
}

MSSelection WSClean::selectInterval(MSSelection& fullSelection)
{
	if(_settings.intervalsOut == 1)
		return fullSelection;
	else {
		size_t tS, tE;
		if(fullSelection.HasInterval())
		{
			tS = fullSelection.IntervalStart();
			tE = fullSelection.IntervalEnd();
		}
		else {
			casacore::MeasurementSet ms(_settings.filenames[0]);
			Logger::Info << "Counting number of scans... ";
			Logger::Info.Flush();
			casacore::ROScalarColumn<double> timeColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::TIME));
			double time = timeColumn(0);
			size_t timestepIndex = 0;
			for(size_t row = 0; row!=ms.nrow(); ++row)
			{
				if(time != timeColumn(row))
				{
					++timestepIndex;
					time = timeColumn(row);
				}
			}
			Logger::Info << "DONE (" << timestepIndex << ")\n";
			tS = 0;
			tE = timestepIndex;
			// Store the full interval in the selection, so that it doesn't need to be determined again.
			fullSelection.SetInterval(tS, tE);
		}
		MSSelection newSelection(fullSelection);
		newSelection.SetInterval(
			tS + (tE-tS) * _currentIntervalIndex / _settings.intervalsOut,
			tS + (tE-tS) * (_currentIntervalIndex+1) / _settings.intervalsOut
		);
		return newSelection;
	}
}

void WSClean::saveUVImage(const double* image, PolarizationEnum pol, size_t channelIndex, bool isImaginary, const std::string& prefix)
{
	ao::uvector<double>
		realUV(_settings.trimmedImageWidth*_settings.trimmedImageHeight, std::numeric_limits<double>::quiet_NaN()),
		imagUV(_settings.trimmedImageWidth*_settings.trimmedImageHeight, std::numeric_limits<double>::quiet_NaN());
	FFTResampler fft(_settings.trimmedImageWidth, _settings.trimmedImageHeight, _settings.trimmedImageWidth, _settings.trimmedImageHeight, 1, true);
	fft.SingleFT(image, realUV.data(), imagUV.data());
	writeFits(prefix+"-real.fits", realUV.data(), pol, channelIndex, isImaginary);
	writeFits(prefix+"-imag.fits", imagUV.data(), pol, channelIndex, isImaginary);
}

void WSClean::makeImagingTable()
{
	std::set<double> channelSet;
	double highestFreq = 0.0;
	_msBands.assign(_settings.filenames.size(), MultiBandData());
	for(size_t i=0; i!=_settings.filenames.size(); ++i)
	{
		casacore::MeasurementSet ms(_settings.filenames[i]);
		_msBands[i] = MultiBandData(ms.spectralWindow(), ms.dataDescription());
		for(size_t d=0; d!=_msBands[i].DataDescCount(); ++d)
		{
			for(size_t ch=0; ch!=_msBands[i][d].ChannelCount(); ++ch)
			{
				double f = _msBands[i][d].ChannelFrequency(ch);
				channelSet.insert(f);
			}
			if(_msBands[i][d].BandEnd() > highestFreq)
				highestFreq = _msBands[i][d].BandEnd();
		}
	}
	if(channelSet.size() < _settings.channelsOut)
	{
		std::ostringstream str;
		str << "Parameter '-channelsout' was set to an invalid value: " << _settings.channelsOut << " output channels requested, but combined in all specified measurement sets, there are only " << channelSet.size() << " unique channels.";
		throw std::runtime_error(str.str());
	}
	_inputChannelFrequencies = std::vector<double>(channelSet.begin(), channelSet.end());
	
	size_t joinedGroupIndex = 0, squaredGroupIndex = 0;
	_imagingTable.Clear();
	
	//for(size_t interval=0; interval!=_settings.intervalsOut; ++interval)
	//{
		if(_settings.joinedFrequencyCleaning)
		{
			size_t maxLocalJGI = joinedGroupIndex;
			for(size_t outChannelIndex=0; outChannelIndex!=_settings.channelsOut; ++outChannelIndex)
			{
				ImagingTableEntry freqTemplate;
				makeImagingTableEntry(_inputChannelFrequencies, outChannelIndex, freqTemplate);
				
				size_t localJGI = joinedGroupIndex;
				addPolarizationsToImagingTable(localJGI, squaredGroupIndex, outChannelIndex, freqTemplate);
				if(localJGI > maxLocalJGI)
					maxLocalJGI = localJGI;
			}
			joinedGroupIndex = maxLocalJGI;
		}
		else {
			for(size_t outChannelIndex=0; outChannelIndex!=_settings.channelsOut; ++outChannelIndex)
			{
				ImagingTableEntry freqTemplate;
				makeImagingTableEntry(_inputChannelFrequencies, outChannelIndex, freqTemplate);
				
				addPolarizationsToImagingTable(joinedGroupIndex, squaredGroupIndex, outChannelIndex, freqTemplate);
			}
		}
	//}
	_imagingTable.Update();
	_imagingTable.Print();
}

void WSClean::makeImagingTableEntry(const std::vector<double>& channels, size_t outChannelIndex, ImagingTableEntry& entry)
{
	size_t startCh, width;
	if(_settings.endChannel != 0)
	{
		if(_settings.endChannel > channels.size())
			throw std::runtime_error("Bad channel selection -- more channels selected than available");
		startCh = _settings.startChannel;
		width = _settings.endChannel-startCh;
	}
	else {
		startCh = 0;
		width = channels.size();
	}
	
	size_t
		chLowIndex = startCh + outChannelIndex*width/_settings.channelsOut,
		chHighIndex = startCh + (outChannelIndex+1)*width/_settings.channelsOut - 1;
	entry.lowestFrequency = channels[chLowIndex];
	entry.highestFrequency = channels[chHighIndex];
	// TODO this should include the channelwidth
	entry.minBandFrequency = entry.lowestFrequency;
	entry.maxBandFrequency = entry.highestFrequency;
	
	entry.msData.resize(_settings.filenames.size());
	for(size_t msIndex=0; msIndex!=_settings.filenames.size(); ++msIndex)
	{
		entry.msData[msIndex].bands.resize(_msBands[msIndex].BandCount());
	}
}

void WSClean::addPolarizationsToImagingTable(size_t& joinedGroupIndex, size_t& squaredGroupIndex, size_t outChannelIndex, const ImagingTableEntry& templateEntry)
{
	for(std::set<PolarizationEnum>::const_iterator p=_settings.polarizations.begin();
			p!=_settings.polarizations.end(); ++p)
	{
		ImagingTableEntry& entry = _imagingTable.AddEntry();
		entry = templateEntry;
		entry.index = _imagingTable.EntryCount()-1;
		entry.outputChannelIndex = outChannelIndex;
		entry.joinedGroupIndex = joinedGroupIndex;
		entry.squaredDeconvolutionIndex = squaredGroupIndex;
		entry.polarization = *p;
		if(*p == Polarization::XY)
			entry.imageCount = 2;
		else if(*p == Polarization::YX)
			entry.imageCount = 0;
		else
			entry.imageCount = 1;
		
		if(!_settings.joinedPolarizationCleaning)
		{
			++joinedGroupIndex;
			++squaredGroupIndex;
		}
	}
	
	if(_settings.joinedPolarizationCleaning)
	{
		++joinedGroupIndex;
		++squaredGroupIndex;
	}
}

void WSClean::fitBeamSize(double& bMaj, double& bMin, double& bPA, const double* image, double beamEstimate)
{
	GaussianFitter beamFitter;
	Logger::Info << "Fitting beam... ";
	Logger::Info.Flush();
	beamFitter.Fit2DGaussianCentred(
		image,
		_settings.trimmedImageWidth, _settings.trimmedImageHeight,
		beamEstimate,
		bMaj, bMin, bPA);
	if(bMaj < 1.0) bMaj = 1.0;
	if(bMin < 1.0) bMin = 1.0;
	bMaj = bMaj*0.5*(_settings.pixelScaleX+_settings.pixelScaleY);
	bMin = bMin*0.5*(_settings.pixelScaleX+_settings.pixelScaleY);
}

void WSClean::determineBeamSize(double& bMaj, double& bMin, double& bPA, const double* image, double theoreticBeam)
{
	if(_settings.manualBeamMajorSize != 0.0)
	{
		bMaj = _settings.manualBeamMajorSize;
		bMin = _settings.manualBeamMinorSize;
		bPA = _settings.manualBeamPA;
	} else if(_settings.fittedBeam)
	{
		fitBeamSize(bMaj, bMin, bPA, image, theoreticBeam*2.0/(_settings.pixelScaleX+_settings.pixelScaleY));
		Logger::Info << "major=" << Angle::ToNiceString(bMaj) << ", minor=" <<
		Angle::ToNiceString(bMin) << ", PA=" << Angle::ToNiceString(bPA) << ", theoretical=" <<
		Angle::ToNiceString(theoreticBeam)<< ".\n";
		
		if(_settings.circularBeam)
		{
			bMin = bMaj;
			bPA = 0.0;
		}
	}
	else if(_settings.theoreticBeam) {
		bMaj = theoreticBeam;
		bMin = theoreticBeam;
		bPA = 0.0;
		Logger::Info << "Beam size is " << Angle::ToNiceString(theoreticBeam) << '\n';
	} else {
		bMaj = std::numeric_limits<double>::quiet_NaN();
		bMin = std::numeric_limits<double>::quiet_NaN();
		bPA = std::numeric_limits<double>::quiet_NaN();
	}
}
