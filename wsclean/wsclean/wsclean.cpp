#include "wsclean.h"

#include "inversionalgorithm.h"
#include "wsmsgridder.h"

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

#include "imageweightcache.h"

#include <iostream>
#include <memory>

std::string commandLine;

WSClean::WSClean() :
	_untrimmedWidth(2048), _untrimmedHeight(2048),
	_trimWidth(0), _trimHeight(0),
	_channelsOut(1), _intervalCount(1),
	_pixelScaleX(0.01 * M_PI / 180.0), _pixelScaleY(0.01 * M_PI / 180.0),
	_manualBeamMajorSize(0.0), _manualBeamMinorSize(0.0),
	_manualBeamPA(0.0), _fittedBeam(true), _theoreticBeam(false), _circularBeam(false),
	_memFraction(1.0), _absMemLimit(0.0),
	_minUVInLambda(0.0), _maxUVInLambda(0.0), _wLimit(0.0),
	_rankFilterLevel(0.0), _rankFilterSize(16),
	_gaussianTaperBeamSize(0.0),
	_tukeyTaperInLambda(0.0), _tukeyInnerTaperInLambda(0.0),
	_edgeTaperInLambda(0.0), _edgeTukeyTaperInLambda(0.0),
	_nWLayers(0), _antialiasingKernelSize(7), _overSamplingFactor(63),
	_threadCount(sysconf(_SC_NPROCESSORS_ONLN)),
	_startChannel(0), _endChannel(0),
	_joinedPolarizationCleaning(false), _joinedFrequencyCleaning(false),
	_predictionChannels(0),
	_globalSelection(),
	_columnName(),
	_polarizations(),
	_weightMode(WeightMode::UniformWeighted),
	_prefixName("wsclean"),
	_smallInversion(true), _makePSF(false), _isWeightImageSaved(false),
	_isUVImageSaved(false), _isGriddingImageSaved(false),
	_dftPrediction(false), _dftWithBeam(false),
	_temporaryDirectory(),
	_forceReorder(false), _forceNoReorder(false),
	_subtractModel(false),
	_modelUpdateRequired(true),
	_mfsWeighting(false),
	_gridMode(WStackingGridder::KaiserBessel),
	_filenames(),
	_commandLine(),
	_inversionWatch(false), _predictingWatch(false), _deconvolutionWatch(false),
	_isFirstInversion(true), _doReorder(false),
	_currentIntervalIndex(0), _majorIterationNr(0),
	_normalizeForWeighting(true),
	_visibilityWeightingMode(InversionAlgorithm::NormalVisibilityWeighting)
{
	_polarizations.insert(Polarization::StokesI);
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
		
	writer.SetImageDimensions(_trimWidth, _trimHeight, ra, dec, pixelScaleX, pixelScaleY);
	writer.SetFrequency(freqCentre, bandwidth);
	writer.SetDate(dateObs);
	writer.SetPolarization(_inversionAlgorithm->Polarization());
	writer.SetOrigin("WSClean", "W-stacking imager written by Andre Offringa");
	writer.AddHistory(commandLine);
	if(_manualBeamMajorSize != 0.0) {
		writer.SetBeamInfo(_manualBeamMajorSize, _manualBeamMinorSize, _manualBeamPA);
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
	if(_endChannel!=0)
	{
		writer.SetExtraKeyword("WSCCHANS", _startChannel);
		writer.SetExtraKeyword("WSCCHANE", _endChannel);
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
	writer.SetExtraKeyword("WSCNITER", _deconvolution.NIter());
	writer.SetExtraKeyword("WSCTHRES", _deconvolution.Threshold());
	writer.SetExtraKeyword("WSCGAIN", _deconvolution.Gain());
	writer.SetExtraKeyword("WSCMGAIN", _deconvolution.MGain());
	writer.SetExtraKeyword("WSCNEGCM", _deconvolution.AllowNegativeComponents());
	writer.SetExtraKeyword("WSCNEGST", _deconvolution.StopOnNegativeComponents());
	//writer.SetExtraKeyword("WSCSMPSF", clean.ResizePSF());
}

void WSClean::updateCleanParameters(FitsWriter& writer, size_t minorIterationNr, size_t majorIterationNr)
{
	writer.SetExtraKeyword("WSCMINOR", minorIterationNr);
	writer.SetExtraKeyword("WSCMAJOR", majorIterationNr);
}

void WSClean::imagePSF(size_t currentChannelIndex)
{
	std::cout << std::flush << " == Constructing PSF ==\n";
	_inversionWatch.Start();
	_inversionAlgorithm->SetDoImagePSF(true);
	_inversionAlgorithm->SetVerbose(_isFirstInversion);
	_inversionAlgorithm->Invert();
		
	DeconvolutionAlgorithm::RemoveNaNsInPSF(_inversionAlgorithm->ImageRealResult(), _trimWidth, _trimHeight);
	initFitsWriter(_fitsWriter);
	_psfImages.SetFitsWriter(_fitsWriter);
	_psfImages.Store(_inversionAlgorithm->ImageRealResult(), *_polarizations.begin(), currentChannelIndex, false);
	_inversionWatch.Pause();
	
	if(_isUVImageSaved)
	{
		saveUVImage(_inversionAlgorithm->ImageRealResult(), *_polarizations.begin(), currentChannelIndex, false, "uvpsf");
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
		
	std::cout << "Writing psf image... " << std::flush;
	const std::string name(getPSFPrefix(currentChannelIndex) + "-psf.fits");
	_fitsWriter.Write(name, _inversionAlgorithm->ImageRealResult());
	std::cout << "DONE\n";
}

void WSClean::imageGridding()
{
	std::cout << "Writing gridding correction image... " << std::flush;
	double* gridding = _imageAllocator.Allocate(_inversionAlgorithm->ActualInversionWidth() * _inversionAlgorithm->ActualInversionHeight());
	_inversionAlgorithm->GetGriddingCorrectionImage(&gridding[0]);
	FitsWriter fitsWriter;
	initFitsWriter(fitsWriter);
	fitsWriter.SetImageDimensions(_inversionAlgorithm->ActualInversionWidth(), _inversionAlgorithm->ActualInversionHeight());
	fitsWriter.Write(_prefixName + "-gridding.fits", &gridding[0]);
	_imageAllocator.Free(gridding);
	std::cout << "DONE\n";
}

void WSClean::imageMainFirst(PolarizationEnum polarization, size_t joinedChannelIndex)
{
	std::cout << std::flush << " == Constructing image ==\n";
	_inversionWatch.Start();
	if(_nWLayers != 0)
		_inversionAlgorithm->SetWGridSize(_nWLayers);
	else
		_inversionAlgorithm->SetNoWGridSize();
	_inversionAlgorithm->SetDoImagePSF(false);
	_inversionAlgorithm->SetDoSubtractModel(_subtractModel);
	_inversionAlgorithm->SetVerbose(_isFirstInversion);
	_inversionAlgorithm->Invert();
	_inversionWatch.Pause();
	_inversionAlgorithm->SetVerbose(false);
	
	storeAndCombineXYandYX(_residualImages, polarization, joinedChannelIndex, false, _inversionAlgorithm->ImageRealResult());
	if(Polarization::IsComplex(polarization))
		storeAndCombineXYandYX(_residualImages, polarization, joinedChannelIndex, true, _inversionAlgorithm->ImageImaginaryResult());
}

void WSClean::imageMainNonFirst(PolarizationEnum polarization, size_t joinedChannelIndex)
{
	std::cout << std::flush << " == Constructing image ==\n";
	_inversionWatch.Start();
	_inversionAlgorithm->SetDoSubtractModel(true);
	_inversionAlgorithm->Invert();
	_inversionWatch.Pause();
	
	storeAndCombineXYandYX(_residualImages, polarization, joinedChannelIndex, false, _inversionAlgorithm->ImageRealResult());
	if(Polarization::IsComplex(polarization))
		storeAndCombineXYandYX(_residualImages, polarization, joinedChannelIndex, true, _inversionAlgorithm->ImageImaginaryResult());
}

void WSClean::storeAndCombineXYandYX(CachedImageSet& dest, PolarizationEnum polarization, size_t joinedChannelIndex, bool isImaginary, const double* image)
{
	if(polarization == Polarization::YX && _polarizations.count(Polarization::XY)!=0)
	{
		std::cout << "Adding XY and YX together...\n";
		double
			*xyImage = _imageAllocator.Allocate(_trimWidth*_trimHeight);
		dest.Load(xyImage, Polarization::XY, joinedChannelIndex, isImaginary);
		size_t count = _trimWidth*_trimHeight;
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
	std::cout << std::flush << " == Converting model image to visibilities ==\n";
	const size_t size = _trimWidth*_trimHeight;
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
	std::cout << std::flush << " == Predicting visibilities ==\n";
	const size_t size = _trimWidth*_trimHeight;
	double
		*modelImageReal = _imageAllocator.Allocate(size),
		*modelImageImaginary = 0;
		
	std::unique_ptr<DFTPredictionImage> image(new DFTPredictionImage(_trimWidth, _trimHeight, _imageAllocator));
		
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
	
	casacore::MeasurementSet firstMS(_filenames.front());
	BandData firstBand(firstMS.spectralWindow());
	DFTPredictionInput input;
	image->FindComponents(input, _inversionAlgorithm->PhaseCentreRA(), _inversionAlgorithm->PhaseCentreDec(), _pixelScaleX, _pixelScaleY, _inversionAlgorithm->PhaseCentreDL(), _inversionAlgorithm->PhaseCentreDM(), firstBand.ChannelCount());
	// Free the input model images
	image.reset();
	std::cout << "Number of components to be predicted: " << input.ComponentCount() << '\n';
	
	_predictingWatch.Start();
	
	for(size_t filenameIndex=0; filenameIndex!=_filenames.size(); ++filenameIndex)
	{
		for(size_t d=0; d!=_msBands[filenameIndex].DataDescCount(); ++d)
		{
			const std::string& msName = _filenames[filenameIndex];
			
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
			LMSPredicter predicter(ms, _threadCount);
			predicter.SetApplyBeam(_dftWithBeam);
			predicter.Input() = input;
			
			if(_dftWithBeam)
			{
				std::cout << "Converting model to absolute values...\n";
				predicter.Input().ConvertApparentToAbsolute(ms);
			}
			
			std::cout << "Creating row mapping...\n";
			std::vector<size_t> msToRowId;
			msProviders[0]->MakeMSRowToRowIdMapping(msToRowId, _globalSelection);
			
			ProgressBar progress("Predicting visibilities for " + msName);
			BandData band(ms.spectralWindow());
			predicter.Start();
			LMSPredicter::RowData row;
			ao::uvector<std::complex<float>> buffer[4];
			for(size_t p=0; p!=_polarizations.size(); ++p)
				buffer[p].assign(band.ChannelCount(), std::complex<float>(0.0));
			
			while(predicter.GetNextRow(row))
			{
				// Write to MS provider(s)
				size_t polIndex = 0;
				for(std::set<PolarizationEnum>::iterator pol=_polarizations.begin(); pol!=_polarizations.end(); ++pol)
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
				for(size_t polIndex=0; polIndex!=_polarizations.size(); ++polIndex)
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
	if(!_mfsWeighting)
	{
		_imageWeightCache->Update(*_inversionAlgorithm, entry.outputChannelIndex, entry.outputTimestepIndex);
		if(_isWeightImageSaved)
			_imageWeightCache->Weights().Save(_prefixName+"-weights.fits");
	}
	_inversionAlgorithm->SetPrecalculatedWeightInfo(&_imageWeightCache->Weights());
}

void WSClean::initializeMFSImageWeights()
{
	std::cout << "Precalculating MFS weights for " << _weightMode.ToString() << " weighting...\n";
	_imageWeightCache->ResetWeights();
	if(_doReorder)
	{
		for(size_t sg=0; sg!=_imagingTable.SquaredGroupCount(); ++sg)
		{
			const ImagingTable subTable = _imagingTable.GetSquaredGroup(sg);
			const ImagingTableEntry& entry = subTable.Front();
			for(size_t msIndex=0; msIndex!=_filenames.size(); ++msIndex)
			{
				const ImagingTableEntry::MSInfo& ms = entry.msData[msIndex];
				for(size_t d=0; d!=_msBands[msIndex].DataDescCount(); ++d)
				{
					MSSelection partSelection(_globalSelection);
					partSelection.SetBandId(d);
					bool hasSelection = selectChannels(partSelection, msIndex, d, subTable.Front());
					if(hasSelection)
					{
						PartitionedMS msProvider(_partitionedMSHandles[msIndex], ms.bands[d].partIndex, entry.polarization, d);
						_imageWeightCache->Weights().Grid(msProvider, partSelection);
					}
				}
			}
		}
	}
	else {
		for(size_t i=0; i!=_filenames.size(); ++i)
		{
			for(size_t d=0; d!=_msBands[i].DataDescCount(); ++d)
			{
				ContiguousMS msProvider(_filenames[i], _columnName, _globalSelection, *_polarizations.begin(), d, _deconvolution.MGain() != 1.0);
				_imageWeightCache->Weights().Grid(msProvider,  _globalSelection);
				std::cout << '.' << std::flush;
			}
		}
	}
	_imageWeightCache->Weights().FinishGridding();
	_imageWeightCache->InitializeWeightTapers();
	if(_isWeightImageSaved)
		_imageWeightCache->Weights().Save(_prefixName+"-weights.fits");
}

void WSClean::prepareInversionAlgorithm(PolarizationEnum polarization)
{
	static_cast<WSMSGridder&>(*_inversionAlgorithm).SetGridMode(_gridMode);
	_inversionAlgorithm->SetImageWidth(_untrimmedWidth);
	_inversionAlgorithm->SetImageHeight(_untrimmedHeight);
	_inversionAlgorithm->SetTrimWidth(_trimWidth);
	_inversionAlgorithm->SetTrimHeight(_trimHeight);
	_inversionAlgorithm->SetPixelSizeX(_pixelScaleX);
	_inversionAlgorithm->SetPixelSizeY(_pixelScaleY);
	if(_nWLayers != 0)
		_inversionAlgorithm->SetWGridSize(_nWLayers);
	else
		_inversionAlgorithm->SetNoWGridSize();
	_inversionAlgorithm->SetAntialiasingKernelSize(_antialiasingKernelSize);
	_inversionAlgorithm->SetOverSamplingFactor(_overSamplingFactor);
	_inversionAlgorithm->SetPolarization(polarization);
	_inversionAlgorithm->SetIsComplex(polarization == Polarization::XY || polarization == Polarization::YX);
	_inversionAlgorithm->SetDataColumnName(_columnName);
	_inversionAlgorithm->SetWeighting(_weightMode);
	_inversionAlgorithm->SetWLimit(_wLimit/100.0);
	_inversionAlgorithm->SetSmallInversion(_smallInversion);
	_inversionAlgorithm->SetNormalizeForWeighting(_normalizeForWeighting);
	_inversionAlgorithm->SetVisibilityWeightingMode(_visibilityWeightingMode);
}

void WSClean::validateDimensions()
{
	if(_trimWidth==0 || _trimHeight==0)
	{
		_trimWidth = _untrimmedWidth;
		_trimHeight = _untrimmedHeight;
	}
}

void WSClean::checkPolarizations()
{
	bool hasXY = _polarizations.count(Polarization::XY)!=0;
	bool hasYX = _polarizations.count(Polarization::YX)!=0;
	if(JoinPolarizations())
	{
		if(_polarizations.size() == 1)
			throw std::runtime_error("Joined polarization cleaning requested, but only one polarization is being imaged. Specify multiple polarizatons, or do not request to join the polarizations");
		else if(_polarizations.size() != 2 && _polarizations.size() != 4)
			throw std::runtime_error("Joined polarization cleaning requested, but neither 2 or 4 polarizations are imaged that are suitable for this");
	}
	else {
		if((hasXY || hasYX) && _deconvolution.NIter() !=0)
			throw std::runtime_error("You are imaging XY and/or YX polarizations and have enabled cleaning (niter!=0). This is not possible -- you have to specify '-joinpolarizations' or disable cleaning.");
	}
	if((hasXY && !hasYX) || (!hasXY && hasYX))
		throw std::runtime_error("You are imaging only one of XY or YX polarizations. This is not possible -- you have to specify both XY and YX polarizations (the output of imaging both polarizations will be the XY and imaginary XY images).");
	if(_deconvolution.IsSpectralFittingEnabled())
	{
		if(JoinPolarizations())
			throw std::runtime_error("You have requested spectral fitting, but you are joining multiple polarizations. This is not supported. You probably want to turn off the joining of polarizations (leave out -joinpolarizations).");
		if(!JoinChannels())
			throw std::runtime_error("You have requested spectral fitting, but you are not joining channels. This is not possible: you probably want to turn channel joining on (add -joinchannels).");
	}
}

void WSClean::performReordering(bool isPredictMode)
{
	_partitionedMSHandles.clear();
	for(size_t i=0; i != _filenames.size(); ++i)
	{
		std::vector<PartitionedMS::ChannelRange> channels;
		size_t nextIndex = 0;
		for(size_t j=0; j!=_imagingTable.SquaredGroupCount(); ++j)
		{
			ImagingTable squaredGroup = _imagingTable.GetSquaredGroup(j);
			for(size_t s=0; s!=squaredGroup.EntryCount(); ++s)
			{
				ImagingTableEntry& entry =
					_imagingTable[squaredGroup[s].index];
				if(entry.polarization == *_polarizations.begin())
				{
					for(size_t d=0; d!=_msBands[i].DataDescCount(); ++d)
					{
						MSSelection selection(_globalSelection);
						if(selectChannels(selection, i, d, entry))
						{
							PartitionedMS::ChannelRange r;
							r.dataDescId = d;
							r.start = selection.ChannelRangeStart();
							r.end = selection.ChannelRangeEnd();
							entry.msData[i].bands[d].partIndex = nextIndex;
							++nextIndex;
							channels.push_back(r);
						}
					}
				}
				else {
					for(size_t d=0; d!=_msBands[i].DataDescCount(); ++d)
						entry.msData[i].bands[d].partIndex = nextIndex-1;
				}
			}
		}
		
		bool useModel = _deconvolution.MGain() != 1.0 || isPredictMode || _subtractModel;
		bool initialModelRequired = _subtractModel;
		_partitionedMSHandles.push_back(PartitionedMS::Partition(_filenames[i], channels, _globalSelection, _columnName, true, useModel, initialModelRequired, _modelUpdateRequired, _polarizations, _temporaryDirectory));
	}
}

void WSClean::RunClean()
{
	// If no column specified, determine column to use
	if(_columnName.empty())
	{
		casacore::MeasurementSet ms(_filenames.front());
		bool hasCorrected = ms.tableDesc().isColumn("CORRECTED_DATA");
		if(hasCorrected) {
			std::cout << "First measurement set has corrected data: tasks will be applied on the corrected data column.\n";
			_columnName = "CORRECTED_DATA";
		} else {
			std::cout << "No corrected data in first measurement set: tasks will be applied on the data column.\n";
			_columnName= "DATA";
		}
	}

	validateDimensions();
	
	checkPolarizations();
	
	MSSelection fullSelection = _globalSelection;
	
	for(_currentIntervalIndex=0; _currentIntervalIndex!=_intervalCount; ++_currentIntervalIndex)
	{
		makeImagingTable();
		
		_globalSelection = selectInterval(fullSelection);
		
		_doReorder = preferReordering();
		
		if(_doReorder) performReordering(false);
		
		_infoPerChannel.assign(_channelsOut, ChannelInfo());
		
		_imageWeightCache.reset(createWeightCache());
		
		if(_mfsWeighting)
			initializeMFSImageWeights();
		
		for(size_t groupIndex=0; groupIndex!=_imagingTable.IndependentGroupCount(); ++groupIndex)
		{
			runIndependentGroup(_imagingTable.GetIndependentGroup(groupIndex));
		}
		
		if(_channelsOut > 1)
		{
			for(std::set<PolarizationEnum>::const_iterator pol=_polarizations.begin(); pol!=_polarizations.end(); ++pol)
			{
				bool psfWasMade = (_deconvolution.NIter() > 0 || _makePSF) && pol == _polarizations.begin();
				
				if(psfWasMade)
					makeMFSImage("psf.fits", *pol, false, true);
				
				if(!(*pol == Polarization::YX && _polarizations.count(Polarization::XY)!=0))
				{
					makeMFSImage("dirty.fits", *pol, false);
					if(_deconvolution.NIter() == 0)
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
						if(_deconvolution.NIter() == 0)
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
	ImageWeightCache* cache = new ImageWeightCache(_weightMode, _untrimmedWidth, _untrimmedHeight, _pixelScaleX, _pixelScaleY, _minUVInLambda, _maxUVInLambda, _rankFilterLevel, _rankFilterSize);
	cache->SetTaperInfo(_gaussianTaperBeamSize, _tukeyTaperInLambda, _tukeyInnerTaperInLambda, _edgeTaperInLambda, _edgeTukeyTaperInLambda);
	return cache;
}

void WSClean::RunPredict()
{
	if(JoinChannels())
		throw std::runtime_error("Joined frequency cleaning specified for prediction: prediction doesn't clean, parameter invalid");
	if(JoinPolarizations())
		throw std::runtime_error("Joined polarization cleaning specified for prediction: prediction doesn't clean, parameter invalid");
	
	_columnName = "DATA";
	
	validateDimensions();
	
	checkPolarizations();
	
	MSSelection fullSelection = _globalSelection;
	
	for(_currentIntervalIndex=0; _currentIntervalIndex!=_intervalCount; ++_currentIntervalIndex)
	{
		makeImagingTable();
		
		if(_predictionChannels != 0)
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
	double lastCh = band.ChannelFrequency(band.ChannelCount()-1);
	if(band.ChannelCount()!=0 && entry.lowestFrequency <= lastCh && entry.highestFrequency >= band.ChannelFrequency(0))
	{
		const double* lowPtr =
			std::lower_bound(band.begin(), band.end(), entry.lowestFrequency);
		const double* highPtr =
			std::lower_bound(lowPtr, band.end(), entry.highestFrequency);
		if(highPtr == band.end())
			--highPtr;
		size_t
			newStart = lowPtr - band.begin(),
			newEnd = highPtr - band.begin() + 1;

		selection.SetChannelRange(newStart, newEnd);
		return true;
	}
	else {
		return false;
	}
}

void WSClean::runIndependentGroup(const ImagingTable& groupTable)
{
	_inversionAlgorithm.reset(new WSMSGridder(&_imageAllocator, _threadCount, _memFraction, _absMemLimit));
	
	_modelImages.Initialize(_fitsWriter, _polarizations.size(), _channelsOut, _prefixName + "-model", _imageAllocator);
	_residualImages.Initialize(_fitsWriter, _polarizations.size(), _channelsOut, _prefixName + "-residual", _imageAllocator);
	if(groupTable.Front().polarization == *_polarizations.begin())
		_psfImages.Initialize(_fitsWriter, 1, groupTable.SquaredGroupCount(), _prefixName + "-psf", _imageAllocator);
	
	const std::string rootPrefix = _prefixName;
		
	for(size_t joinedIndex=0; joinedIndex!=groupTable.EntryCount(); ++joinedIndex)
	{
		const ImagingTableEntry& entry = groupTable[joinedIndex];
		runFirstInversion(entry);
	}
	
	_deconvolution.InitializeDeconvolutionAlgorithm(groupTable, *_polarizations.begin(), &_imageAllocator, _trimWidth, _trimHeight, _pixelScaleX, _pixelScaleY, _channelsOut, _inversionAlgorithm->BeamSize(), _threadCount);

	initFitsWriter(_fitsWriter);
	setCleanParameters(_fitsWriter);
	updateCleanParameters(_fitsWriter, 0, 0);
		
	if(_deconvolution.NIter() > 0)
	{
		// Start major cleaning loop
		_majorIterationNr = 1;
		bool reachedMajorThreshold = false;
		do {
			_deconvolution.InitializeImages(_residualImages, _modelImages, _psfImages);
			_deconvolutionWatch.Start();
			_deconvolution.Perform(groupTable, reachedMajorThreshold, _majorIterationNr);
			_deconvolutionWatch.Pause();
			
			if(_majorIterationNr == 1 && _deconvolution.MGain() != 1.0)
				writeFirstResidualImages(groupTable);
	
			if(!reachedMajorThreshold)
				writeModelImages(groupTable);
	
			if(_deconvolution.MGain() != 1.0)
			{
				for(size_t sGroupIndex=0; sGroupIndex!=groupTable.SquaredGroupCount(); ++sGroupIndex)
				{
					const ImagingTable sGroupTable = groupTable.GetSquaredGroup(sGroupIndex);
					size_t currentChannelIndex = sGroupTable.Front().outputChannelIndex;
					if(_dftPrediction)
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
		
		std::cout << _majorIterationNr << " major iterations were performed.\n";
	}
	
	_inversionAlgorithm->FreeImagingData();
	
	// Restore model to residuals and save all images
	for(size_t joinedIndex=0; joinedIndex!=groupTable.EntryCount(); ++joinedIndex)
	{
		size_t currentChannelIndex =
			groupTable[joinedIndex].outputChannelIndex;
		
		double
			freqLow = groupTable[joinedIndex].minBandFrequency,
			freqHigh = groupTable[joinedIndex].maxBandFrequency;
		
		PolarizationEnum curPol = groupTable[joinedIndex].polarization;
		for(size_t imageIter=0; imageIter!=groupTable[joinedIndex].imageCount; ++imageIter)
		{
			bool isImaginary = (imageIter == 1);
			double* restoredImage = _imageAllocator.Allocate(_trimWidth*_trimHeight);
			_residualImages.Load(restoredImage, curPol, currentChannelIndex, isImaginary);
			if(_deconvolution.NIter() != 0)
				writeFits("residual.fits", restoredImage, curPol, currentChannelIndex, isImaginary);
			if(_isUVImageSaved)
				saveUVImage(restoredImage, curPol, currentChannelIndex, isImaginary, "uv");
			double* modelImage = _imageAllocator.Allocate(_trimWidth*_trimHeight);
				_modelImages.Load(modelImage, curPol, currentChannelIndex, isImaginary);
			ModelRenderer renderer(_fitsWriter.RA(), _fitsWriter.Dec(), _pixelScaleX, _pixelScaleY, _fitsWriter.PhaseCentreDL(), _fitsWriter.PhaseCentreDM());
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
			if(_deconvolution.MultiScale() || _deconvolution.FastMultiScale() || _deconvolution.UseMoreSane() || _deconvolution.UseIUWT())
			{
				std::cout << "Rendering sources to restored image " + beamStr + "... " << std::flush;
				renderer.Restore(restoredImage, modelImage, _trimWidth, _trimHeight, beamMaj, beamMin, beamPA);
				std::cout << "DONE\n";
			}
			else {
				Model model;
				// A model cannot hold instrumental pols (xx/xy/yx/yy), hence always use Stokes I here
				DeconvolutionAlgorithm::GetModelFromImage(model, modelImage, _trimWidth, _trimHeight, _fitsWriter.RA(), _fitsWriter.Dec(), _pixelScaleX, _pixelScaleY, _fitsWriter.PhaseCentreDL(), _fitsWriter.PhaseCentreDM(), 0.0, _fitsWriter.Frequency(), Polarization::StokesI);
				
				if(beamMaj == beamMin) {
					std::cout << "Rendering " << model.SourceCount() << " circular sources to restored image " + beamStr + "... " << std::flush;
					renderer.Restore(restoredImage, _trimWidth, _trimHeight, model, beamMaj, freqLow, freqHigh, Polarization::StokesI);
				}
				else {
					std::cout << "Rendering " << model.SourceCount() << " elliptical sources to restored image " + beamStr + "... " << std::flush;
					renderer.Restore(restoredImage, _trimWidth, _trimHeight, model, beamMaj, beamMin, beamPA, freqLow, freqHigh, Polarization::StokesI);
				}
				std::cout << "DONE\n";
			}
			std::cout << "Writing restored image... " << std::flush;
			writeFits("image.fits", restoredImage, curPol, currentChannelIndex, isImaginary);
			std::cout << "DONE\n";
			_imageAllocator.Free(restoredImage);
			_imageAllocator.Free(modelImage);
		}
	} // end joined index loop
	
	_imageAllocator.ReportStatistics();
	std::cout << "Inversion: " << _inversionWatch.ToString() << ", prediction: " << _predictingWatch.ToString() << ", deconvolution: " << _deconvolutionWatch.ToString() << '\n';
	
	_prefixName = rootPrefix;
	
	// Needs to be destructed before image allocator, or image allocator will report error caused by leaked memory
	_inversionAlgorithm.reset();
}

void WSClean::writeFirstResidualImages(const ImagingTable& groupTable)
{
	std::cout << "Writing first iteration image(s)...\n";
	ImageBufferAllocator::Ptr ptr;
	_imageAllocator.Allocate(_trimWidth*_trimHeight, ptr);
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
	std::cout << "Writing model image...\n";
	ImageBufferAllocator::Ptr ptr;
	_imageAllocator.Allocate(_trimWidth*_trimHeight, ptr);
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
	_inversionAlgorithm.reset(new WSMSGridder(&_imageAllocator, _threadCount, _memFraction, _absMemLimit));
	
	_modelImages.Initialize(_fitsWriter, _polarizations.size(), 1, _prefixName + "-model", _imageAllocator);
	
	const std::string rootPrefix = _prefixName;
		
	for(size_t e=0; e!=imagingGroup.EntryCount(); ++e)
	{
		const ImagingTableEntry& entry = imagingGroup[e];
		// load image(s) from disk
		for(size_t i=0; i!=entry.imageCount; ++i)
		{
			FitsReader reader(getPrefix(entry.polarization, entry.outputChannelIndex, i==1) + "-model.fits");
			_fitsWriter = FitsWriter(reader);
			_modelImages.SetFitsWriter(_fitsWriter);
			std::cout << "Reading " << reader.Filename() << "...\n";
			double* buffer = _imageAllocator.Allocate(_trimWidth*_trimHeight);
			if(reader.ImageWidth()!=_trimWidth || reader.ImageHeight()!=_trimHeight)
				throw std::runtime_error("Inconsistent image size: input image did not match with specified dimensions.");
			reader.Read(buffer);
			for(size_t j=0; j!=_trimWidth*_trimHeight; ++j)
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
	std::cout << "Inversion: " << _inversionWatch.ToString() << ", prediction: " << _predictingWatch.ToString() << ", cleaning: " << _deconvolutionWatch.ToString() << '\n';
	
	_prefixName = rootPrefix;
	
	// Needs to be destructed before image allocator, or image allocator will report error caused by leaked memory
	_inversionAlgorithm.reset();
}

MSProvider* WSClean::initializeMSProvider(const ImagingTableEntry& entry, const MSSelection& selection, size_t filenameIndex, size_t dataDescId)
{
	if(_doReorder)
		return new PartitionedMS(_partitionedMSHandles[filenameIndex], entry.msData[filenameIndex].bands[dataDescId].partIndex, entry.polarization, dataDescId);
	else
		return new ContiguousMS(_filenames[filenameIndex], _columnName, selection, entry.polarization, dataDescId, _deconvolution.MGain() != 1.0);
}

void WSClean::initializeCurMSProviders(const ImagingTableEntry& entry)
{
	_inversionAlgorithm->ClearMeasurementSetList();
	for(size_t i=0; i != _filenames.size(); ++i)
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

	bool isFirstPol = entry.polarization == *_polarizations.begin();
	bool doMakePSF = _deconvolution.NIter() > 0 || _makePSF;
	if(doMakePSF && isFirstPol)
		imagePSF(entry.outputChannelIndex);
	
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
			if(_theoreticBeam) {
				_infoPerChannel[entry.outputChannelIndex].beamMaj = _inversionAlgorithm->BeamSize();
				_infoPerChannel[entry.outputChannelIndex].beamMin = _inversionAlgorithm->BeamSize();
				_infoPerChannel[entry.outputChannelIndex].beamPA = 0.0;
			}
			else if(_manualBeamMajorSize != 0.0) {
				_infoPerChannel[entry.outputChannelIndex].beamMaj = _manualBeamMajorSize;
				_infoPerChannel[entry.outputChannelIndex].beamMin = _manualBeamMinorSize;
				_infoPerChannel[entry.outputChannelIndex].beamPA = _manualBeamPA;
			}
			else {
				_infoPerChannel[entry.outputChannelIndex].beamMaj = std::numeric_limits<double>::quiet_NaN();
				_infoPerChannel[entry.outputChannelIndex].beamMin = std::numeric_limits<double>::quiet_NaN();
				_infoPerChannel[entry.outputChannelIndex].beamPA = std::numeric_limits<double>::quiet_NaN();
			}
		}
	}
	
	if(_isGriddingImageSaved && firstBeforePSF && _inversionAlgorithm->HasGriddingCorrectionImage())
		imageGridding();
	
	_isFirstInversion = false;
	
	// Set model to zero: already done if this is YX of XY/YX imaging combi
	if(!(entry.polarization == Polarization::YX && _polarizations.count(Polarization::XY)!=0))
	{
		double* modelImage = _imageAllocator.Allocate(_trimWidth * _trimHeight);
		memset(modelImage, 0, _trimWidth * _trimHeight * sizeof(double));
		_modelImages.Store(modelImage, entry.polarization, entry.outputChannelIndex, false);
		if(Polarization::IsComplex(entry.polarization))
			_modelImages.Store(modelImage, entry.polarization, entry.outputChannelIndex, true);
		_imageAllocator.Free(modelImage);
	}
	
	if(entry.polarization == Polarization::XY && _polarizations.count(Polarization::YX)!=0)
	{ // Skip saving XY of XY/YX combi
	}
	else {
		PolarizationEnum savedPol = entry.polarization;
		if(savedPol == Polarization::YX && _polarizations.count(Polarization::XY)!=0)
			savedPol = Polarization::XY;
		double* dirtyImage = _imageAllocator.Allocate(_trimWidth * _trimHeight);
		_residualImages.Load(dirtyImage, savedPol, entry.outputChannelIndex, false);
		std::cout << "Writing dirty image...\n";
		writeFits("dirty.fits", dirtyImage, savedPol, entry.outputChannelIndex, false);
		if(Polarization::IsComplex(entry.polarization))
		{
			_residualImages.Load(dirtyImage, savedPol, entry.outputChannelIndex, true);
			writeFits("dirty.fits", dirtyImage, savedPol, entry.outputChannelIndex, true);
		}
		_imageAllocator.Free(dirtyImage);
	}
	
	clearCurMSProviders();
}

void WSClean::makeMFSImage(const string& suffix, PolarizationEnum pol, bool isImaginary, bool isPSF)
{
	double lowestFreq = 0.0, highestFreq = 0.0;
	const size_t size = _trimWidth * _trimHeight;
	ao::uvector<double> mfsImage(size, 0.0), addedImage(size), weightImage(size, 0.0);
	double weightSum = 0.0;
	FitsWriter writer;
	for(size_t ch=0; ch!=_channelsOut; ++ch)
	{
		std::string prefixStr = isPSF ? getPSFPrefix(ch) : getPrefix(pol, ch, isImaginary);
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
			std::cout << "Error: image " << name << " did not have the WSCIMGWG keyword.\n";
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
	
	std::string mfsName(getMFSPrefix(pol, isImaginary, isPSF) + '-' + suffix);
	std::cout << "Writing " << mfsName << "...\n";
	writer.SetFrequency((lowestFreq+highestFreq)*0.5, highestFreq-lowestFreq);
	writer.SetExtraKeyword("WSCIMGWG", weightSum);
	writer.RemoveExtraKeyword("WSCCHANS");
	writer.RemoveExtraKeyword("WSCCHANE");
	writer.Write(mfsName, mfsImage.data());
}

void WSClean::renderMFSImage(PolarizationEnum pol, bool isImaginary)
{
	const size_t size = _trimWidth * _trimHeight;
	
	std::string mfsPrefix(getMFSPrefix(pol, isImaginary, false));
	FitsReader residualReader(mfsPrefix + "-residual.fits");
	FitsReader modelReader(mfsPrefix + "-model.fits");
	ao::uvector<double> image(size), modelImage(size);
	residualReader.Read(image.data());
	modelReader.Read(modelImage.data());
	
	ModelRenderer renderer(_fitsWriter.RA(), _fitsWriter.Dec(), _pixelScaleX, _pixelScaleY, _fitsWriter.PhaseCentreDL(), _fitsWriter.PhaseCentreDM());
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
	std::cout << "Rendering sources to restored image " + beamStr + "... " << std::flush;
	renderer.Restore(image.data(), modelImage.data(), _trimWidth, _trimHeight, beamMaj, beamMin, beamPA);
	std::cout << "DONE\n";
	
	std::cout << "Writing " + mfsPrefix + "-image.fits...\n";
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
	const std::string name(getPrefix(pol, channelIndex, isImaginary) + '-' + suffix);
	initFitsWriter(_fitsWriter);
	_fitsWriter.SetPolarization(pol);
	_fitsWriter.SetFrequency(centreFrequency, bandwidth);
	_fitsWriter.SetExtraKeyword("WSCIMGWG", _infoPerChannel[channelIndex].weight);
	_fitsWriter.SetBeamInfo(
		_infoPerChannel[channelIndex].beamMaj,
		_infoPerChannel[channelIndex].beamMin,
		_infoPerChannel[channelIndex].beamPA);
	size_t polIndex;
	if(JoinPolarizations())
		polIndex = 0;
	else
		Polarization::TypeToIndex(pol, _polarizations, polIndex);
	setCleanParameters(_fitsWriter);
	if(_deconvolution.IsInitialized())
		updateCleanParameters(_fitsWriter, _deconvolution.GetAlgorithm().IterationNumber(), _majorIterationNr);
	_fitsWriter.Write(name, image);
}

MSSelection WSClean::selectInterval(MSSelection& fullSelection)
{
	if(_intervalCount == 1)
		return fullSelection;
	else {
		size_t tS, tE;
		if(fullSelection.HasInterval())
		{
			tS = fullSelection.IntervalStart();
			tE = fullSelection.IntervalEnd();
		}
		else {
			casacore::MeasurementSet ms(_filenames[0]);
			std::cout << "Counting number of scans... " << std::flush;
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
			std::cout << "DONE (" << timestepIndex << ")\n";
			tS = 0;
			tE = timestepIndex;
			// Store the full interval in the selection, so that it doesn't need to be determined again.
			fullSelection.SetInterval(tS, tE);
		}
		MSSelection newSelection(fullSelection);
		newSelection.SetInterval(
			tS + (tE-tS) * _currentIntervalIndex / _intervalCount,
			tS + (tE-tS) * (_currentIntervalIndex+1) / _intervalCount
		);
		return newSelection;
	}
}

void WSClean::saveUVImage(const double* image, PolarizationEnum pol, size_t channelIndex, bool isImaginary, const std::string& prefix)
{
	ao::uvector<double>
		realUV(_trimWidth*_trimHeight, std::numeric_limits<double>::quiet_NaN()),
		imagUV(_trimWidth*_trimHeight, std::numeric_limits<double>::quiet_NaN());
	FFTResampler fft(_trimWidth, _trimHeight, _trimWidth, _trimHeight, 1, true);
	fft.SingleFT(image, realUV.data(), imagUV.data());
	writeFits(prefix+"-real.fits", realUV.data(), pol, channelIndex, isImaginary);
	writeFits(prefix+"-imag.fits", imagUV.data(), pol, channelIndex, isImaginary);
}

void WSClean::makeImagingTable()
{
	std::set<double> channelSet;
	double highestFreq = 0.0;
	_msBands.assign(_filenames.size(), MultiBandData());
	for(size_t i=0; i!=_filenames.size(); ++i)
	{
		casacore::MeasurementSet ms(_filenames[i]);
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
	if(channelSet.size() < _channelsOut)
	{
		std::ostringstream str;
		str << "Parameter '-channelsout' was set to an invalid value: " << _channelsOut << " output channels requested, but combined in all specified measurement sets, there are only " << channelSet.size() << " unique channels.";
		throw std::runtime_error(str.str());
	}
	_inputChannelFrequencies = std::vector<double>(channelSet.begin(), channelSet.end());
	
	size_t joinedGroupIndex = 0, squaredGroupIndex = 0;
	_imagingTable.Clear();
	
	//for(size_t interval=0; interval!=_intervalCount; ++interval)
	//{
		if(JoinChannels())
		{
			size_t maxLocalJGI = joinedGroupIndex;
			for(size_t outChannelIndex=0; outChannelIndex!=_channelsOut; ++outChannelIndex)
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
			for(size_t outChannelIndex=0; outChannelIndex!=_channelsOut; ++outChannelIndex)
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
	if(_endChannel != 0)
	{
		if(_endChannel > channels.size())
			throw std::runtime_error("Bad channel selection -- more channels selected than available");
		startCh = _startChannel;
		width = _endChannel-startCh;
	}
	else {
		startCh = 0;
		width = channels.size();
	}
	
	size_t
		chLowIndex = startCh + outChannelIndex*width/_channelsOut,
		chHighIndex = startCh + (outChannelIndex+1)*width/_channelsOut - 1;
	entry.lowestFrequency = channels[chLowIndex];
	entry.highestFrequency = channels[chHighIndex];
	// TODO this should include the channelwidth
	entry.minBandFrequency = entry.lowestFrequency;
	entry.maxBandFrequency = entry.highestFrequency;
	
	entry.msData.resize(_filenames.size());
	for(size_t msIndex=0; msIndex!=_filenames.size(); ++msIndex)
	{
		entry.msData[msIndex].bands.resize(_msBands[msIndex].BandCount());
	}
}

void WSClean::addPolarizationsToImagingTable(size_t& joinedGroupIndex, size_t& squaredGroupIndex, size_t outChannelIndex, const ImagingTableEntry& templateEntry)
{
	for(std::set<PolarizationEnum>::const_iterator p=_polarizations.begin();
			p!=_polarizations.end(); ++p)
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
		
		if(!JoinPolarizations())
		{
			++joinedGroupIndex;
			++squaredGroupIndex;
		}
	}
	
	if(JoinPolarizations())
	{
		++joinedGroupIndex;
		++squaredGroupIndex;
	}
}

void WSClean::fitBeamSize(double& bMaj, double& bMin, double& bPA, const double* image, double beamEstimate)
{
	GaussianFitter beamFitter;
	std::cout << "Fitting beam... " << std::flush;
	beamFitter.Fit2DGaussianCentred(
		image,
		_trimWidth, _trimHeight,
		beamEstimate,
		bMaj, bMin, bPA);
	if(bMaj < 1.0) bMaj = 1.0;
	if(bMin < 1.0) bMin = 1.0;
	bMaj = bMaj*0.5*(_pixelScaleX+_pixelScaleY);
	bMin = bMin*0.5*(_pixelScaleX+_pixelScaleY);
}

void WSClean::determineBeamSize(double& bMaj, double& bMin, double& bPA, const double* image, double theoreticBeam)
{
	if(_manualBeamMajorSize != 0.0)
	{
		bMaj = _manualBeamMajorSize;
		bMin = _manualBeamMinorSize;
		bPA = _manualBeamPA;
	} else if(_fittedBeam)
	{
		fitBeamSize(bMaj, bMin, bPA, image, theoreticBeam*2.0/(_pixelScaleX+_pixelScaleY));
		std::cout << "major=" << Angle::ToNiceString(bMaj) << ", minor=" <<
		Angle::ToNiceString(bMin) << ", PA=" << Angle::ToNiceString(bPA) << ", theoretical=" <<
		Angle::ToNiceString(theoreticBeam)<< ".\n";
		
		if(_circularBeam)
		{
			bMin = bMaj;
			bPA = 0.0;
		}
	}
	else if(_theoreticBeam) {
		bMaj = theoreticBeam;
		bMin = theoreticBeam;
		bPA = 0.0;
		std::cout << "Beam size is " << Angle::ToNiceString(theoreticBeam) << '\n';
	} else {
		bMaj = std::numeric_limits<double>::quiet_NaN();
		bMin = std::numeric_limits<double>::quiet_NaN();
		bPA = std::numeric_limits<double>::quiet_NaN();
	}
}
