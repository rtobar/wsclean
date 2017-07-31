#include "wsmsgridder.h"

#include "imagebufferallocator.h"
#include "logger.h"

#include "../imageweights.h"
#include "../buffered_lane.h"
#include "../fftresampler.h"
#include "../image.h"

#include "../msproviders/msprovider.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <iostream>
#include <stdexcept>

WSMSGridder::WSMSGridder(ImageBufferAllocator* imageAllocator, size_t threadCount, double memFraction, double absMemLimit) :
	MSGridderBase(),
	_cpuCount(threadCount),
	_laneBufferSize(std::max<size_t>(_cpuCount*2,1024)),
	_imageBufferAllocator(imageAllocator)
{
	long int pageCount = sysconf(_SC_PHYS_PAGES), pageSize = sysconf(_SC_PAGE_SIZE);
	_memSize = (int64_t) pageCount * (int64_t) pageSize;
	double memSizeInGB = (double) _memSize / (1024.0*1024.0*1024.0);
	if(memFraction == 1.0 && absMemLimit == 0.0) {
		Logger::Info << "Detected " << round(memSizeInGB*10.0)/10.0 << " GB of system memory, usage not limited.\n";
	}
	else {
		double limitInGB = memSizeInGB*memFraction;
		if(absMemLimit!=0.0 && limitInGB > absMemLimit)
			limitInGB = absMemLimit;
		Logger::Info << "Detected " << round(memSizeInGB*10.0)/10.0 << " GB of system memory, usage limited to " << round(limitInGB*10.0)/10.0 << " GB (frac=" << round(memFraction*1000.0)/10.0 << "%, ";
		if(absMemLimit == 0.0)
			Logger::Info << "no limit)\n";
		else
			Logger::Info << "limit=" << round(absMemLimit*10.0)/10.0 << "GB)\n";
		
		_memSize = int64_t((double) pageCount * (double) pageSize * memFraction);
		if(absMemLimit!=0.0 && double(_memSize) > double(1024.0*1024.0*1024.0) * absMemLimit)
			_memSize = int64_t(double(absMemLimit) * double(1024.0*1024.0*1024.0));
	}
}
		
void WSMSGridder::countSamplesPerLayer(MSData& msData)
{
	ao::uvector<size_t> sampleCount(WGridSize(), 0);
	size_t total = 0;
	msData.matchingRows = 0;
	msData.msProvider->Reset();
	while(msData.msProvider->CurrentRowAvailable())
	{
		double uInM, vInM, wInM;
		size_t dataDescId;
		msData.msProvider->ReadMeta(uInM, vInM, wInM, dataDescId);
		const BandData& bandData(msData.bandData[dataDescId]);
		for(size_t ch=msData.startChannel; ch!=msData.endChannel; ++ch)
		{
			double w = wInM / bandData.ChannelWavelength(ch);
			size_t wLayerIndex = _gridder->WToLayer(w);
			if(wLayerIndex < WGridSize())
			{
				++sampleCount[wLayerIndex];
				++total;
			}
		}
		++msData.matchingRows;
		msData.msProvider->NextRow();
	}
	Logger::Debug << "Visibility count per layer: ";
	for(ao::uvector<size_t>::const_iterator i=sampleCount.begin(); i!=sampleCount.end(); ++i)
	{
		Logger::Debug << *i << ' ';
	}
	Logger::Debug << "\nTotal nr. of visibilities to be gridded: " << total << '\n';
}

size_t WSMSGridder::getSuggestedWGridSize() const
{
	size_t wWidth, wHeight;
	if(HasNWSize()) {
		wWidth = NWWidth(); wHeight = NWHeight();
	}
	else {
		wWidth = TrimWidth(); wHeight = TrimHeight();
	}
	double
		maxL = wWidth * PixelSizeX() * 0.5 + fabs(PhaseCentreDL()),
		maxM = wHeight * PixelSizeY() * 0.5 + fabs(PhaseCentreDM()),
		lmSq = maxL * maxL + maxM * maxM;
	double cMinW = IsComplex() ? -_maxW : _minW;
	double radiansForAllLayers;
	if(lmSq < 1.0)
		radiansForAllLayers = 2 * M_PI * (_maxW - cMinW) * (1.0 - sqrt(1.0 - lmSq));
	else
		radiansForAllLayers = 2 * M_PI * (_maxW - cMinW);
	size_t suggestedGridSize = size_t(ceil(radiansForAllLayers));
	if(suggestedGridSize == 0) suggestedGridSize = 1;
	if(suggestedGridSize < _cpuCount)
	{
		// When nwlayers is lower than the nr of cores, we cannot parallellize well. 
		// However, we don't want extra w-layers if we are low on mem, as that might slow down the process
		double memoryRequired = double(_cpuCount) * double(sizeof(double))*double(_actualInversionWidth*_actualInversionHeight);
		if(4.0 * memoryRequired < double(_memSize))
		{
			Logger::Info <<
				"The theoretically suggested number of w-layers (" << suggestedGridSize << ") is less than the number of availables\n"
				"cores (" << _cpuCount << "). Changing suggested number of w-layers to " << _cpuCount << ".\n";
			suggestedGridSize = _cpuCount;
		}
		else {
			Logger::Info <<
				"The theoretically suggested number of w-layers (" << suggestedGridSize << ") is less than the number of availables\n"
				"cores (" << _cpuCount << "), but there is not enough memory available to increase the number of w-layers.\n"
				"Not all cores can be used efficiently.\n";
		}
	}
	if(Verbose())
		Logger::Info << "Suggested number of w-layers: " << ceil(suggestedGridSize) << '\n';
	return suggestedGridSize;
}

void WSMSGridder::gridMeasurementSet(MSData &msData)
{
	const MultiBandData selectedBand(msData.SelectedBand());
	_gridder->PrepareBand(selectedBand);
	ao::uvector<std::complex<float>> modelBuffer(selectedBand.MaxChannels());
	ao::uvector<float> weightBuffer(selectedBand.MaxChannels());
	ao::uvector<bool> isSelected(selectedBand.MaxChannels());
	
	// Samples of the same w-layer are collected in a buffer
	// before they are written into the lane. This is done because writing
	// to a lane is reasonably slow; it requires holding a mutex. Without
	// these buffers, writing the lane was a bottleneck and multithreading
	// did not help. I think.
	std::unique_ptr<lane_write_buffer<InversionWorkSample>[]>
		bufferedLanes(new lane_write_buffer<InversionWorkSample>[_cpuCount]);
	size_t bufferSize = std::max<size_t>(8u, _inversionCPULanes[0].capacity()/8);
	bufferSize = std::min<size_t>(128, std::min(bufferSize, _inversionCPULanes[0].capacity()));
	for(size_t i=0; i!=_cpuCount; ++i)
	{
		bufferedLanes[i].reset(&_inversionCPULanes[i], bufferSize);
	}
	
	InversionRow newItem;
	ao::uvector<std::complex<float>> newItemData(selectedBand.MaxChannels());
	newItem.data = newItemData.data();
			
	size_t rowsRead = 0;
	msData.msProvider->Reset();
	while(msData.msProvider->CurrentRowAvailable())
	{
		size_t dataDescId;
		double uInMeters, vInMeters, wInMeters;
		msData.msProvider->ReadMeta(uInMeters, vInMeters, wInMeters, dataDescId);
		const BandData& curBand(selectedBand[dataDescId]);
		const double
			w1 = wInMeters / curBand.LongestWavelength(),
			w2 = wInMeters / curBand.SmallestWavelength();
		if(_gridder->IsInLayerRange(w1, w2))
		{
			newItem.uvw[0] = uInMeters;
			newItem.uvw[1] = vInMeters;
			newItem.uvw[2] = wInMeters;
			newItem.dataDescId = dataDescId;
			
			// Any visibilities that are not gridded in this pass
			// should not contribute to the weight sum, so set these
			// to have zero weight.
			for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
			{
				double w = newItem.uvw[2] / curBand.ChannelWavelength(ch);
				isSelected[ch] = _gridder->IsInLayerRange(w);
			}
	
			readAndWeightVisibilities<1>(*msData.msProvider, newItem, curBand, weightBuffer.data(), modelBuffer.data(), isSelected.data());
			
			InversionWorkSample sampleData;
			for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
			{
				double wavelength = curBand.ChannelWavelength(ch);
				sampleData.sample = newItem.data[ch];
				sampleData.uInLambda = newItem.uvw[0] / wavelength;
				sampleData.vInLambda = newItem.uvw[1] / wavelength;
				sampleData.wInLambda = newItem.uvw[2] / wavelength;
				size_t cpu = _gridder->WToLayer(sampleData.wInLambda) % _cpuCount;
				bufferedLanes[cpu].write(sampleData);
			}
			
			++rowsRead;
		}
		
		msData.msProvider->NextRow();
	}
	
	for(size_t i=0; i!=_cpuCount; ++i)
		bufferedLanes[i].write_end();
	
	if(Verbose())
		Logger::Info << "Rows that were required: " << rowsRead << '/' << msData.matchingRows << '\n';
	msData.totalRowsProcessed += rowsRead;
}

void WSMSGridder::startInversionWorkThreads(size_t maxChannelCount)
{
	_inversionCPULanes.reset(new ao::lane<InversionWorkSample>[_cpuCount]);
	boost::thread_group group;
	_threadGroup.reset(new boost::thread_group());
	for(size_t i=0; i!=_cpuCount; ++i)
	{
		_inversionCPULanes[i].resize(maxChannelCount * _laneBufferSize);
		set_lane_debug_name(_inversionCPULanes[i], "Work lane (buffered) containing individual visibility samples");
		_threadGroup->add_thread(new boost::thread(&WSMSGridder::workThreadPerSample, this, &_inversionCPULanes[i]));
	}
}

void WSMSGridder::finishInversionWorkThreads()
{
	_threadGroup->join_all();
	_threadGroup.reset();
	_inversionCPULanes.reset();
}

void WSMSGridder::workThreadPerSample(ao::lane<InversionWorkSample>* workLane)
{
	size_t bufferSize = std::max<size_t>(8u, workLane->capacity()/8);
	bufferSize = std::min<size_t>(128,std::min(bufferSize, workLane->capacity()));
	lane_read_buffer<InversionWorkSample> buffer(workLane, bufferSize);
	InversionWorkSample sampleData;
	while(buffer.read(sampleData))
	{
		_gridder->AddDataSample(sampleData.sample, sampleData.uInLambda, sampleData.vInLambda, sampleData.wInLambda);
	}
}

void WSMSGridder::predictMeasurementSet(MSData &msData)
{
	msData.msProvider->ReopenRW();
	const MultiBandData selectedBandData(msData.SelectedBand());
	_gridder->PrepareBand(selectedBandData);
	
	size_t rowsProcessed = 0;
	
	ao::lane<PredictionWorkItem>
		calcLane(_laneBufferSize+_cpuCount),
		writeLane(_laneBufferSize);
	set_lane_debug_name(calcLane, "Prediction calculation lane (buffered) containing full row data");
	set_lane_debug_name(writeLane, "Prediction write lane containing full row data");
	lane_write_buffer<PredictionWorkItem> bufferedCalcLane(&calcLane, _laneBufferSize);
	boost::thread writeThread(&WSMSGridder::predictWriteThread, this, &writeLane, &msData);
	boost::thread_group calcThreads;
	for(size_t i=0; i!=_cpuCount; ++i)
		calcThreads.add_thread(new boost::thread(&WSMSGridder::predictCalcThread, this, &calcLane, &writeLane));

		
	/* Start by reading the u,v,ws in, so we don't need IO access
	 * from this thread during further processing */
	std::vector<double> us, vs, ws;
	std::vector<size_t> rowIds, dataIds;
	msData.msProvider->Reset();
	while(msData.msProvider->CurrentRowAvailable())
	{
		size_t dataDescId;
		double uInMeters, vInMeters, wInMeters;
		msData.msProvider->ReadMeta(uInMeters, vInMeters, wInMeters, dataDescId);
		const BandData& curBand(selectedBandData[dataDescId]);
		const double
			w1 = wInMeters / curBand.LongestWavelength(),
			w2 = wInMeters / curBand.SmallestWavelength();
		if(_gridder->IsInLayerRange(w1, w2))
		{
			us.push_back(uInMeters);
			vs.push_back(vInMeters);
			ws.push_back(wInMeters);
			dataIds.push_back(dataDescId);
			rowIds.push_back(msData.msProvider->RowId());
			++rowsProcessed;
		}
		
		msData.msProvider->NextRow();
	}
	
	for(size_t i=0; i!=us.size(); ++i)
	{
		PredictionWorkItem newItem;
		newItem.u = us[i];
		newItem.v = vs[i];
		newItem.w = ws[i];
		newItem.dataDescId = dataIds[i];
		newItem.data.reset(new std::complex<float>[selectedBandData[dataIds[i]].ChannelCount()]);
		newItem.rowId = rowIds[i];
				
		bufferedCalcLane.write(std::move(newItem));
	}
	if(Verbose())
		Logger::Info << "Rows that were required: " << rowsProcessed << '/' << msData.matchingRows << '\n';
	msData.totalRowsProcessed += rowsProcessed;
	
	bufferedCalcLane.write_end();
	calcThreads.join_all();
	writeLane.write_end();
	writeThread.join();
}

void WSMSGridder::predictCalcThread(ao::lane<PredictionWorkItem>* inputLane, ao::lane<PredictionWorkItem>* outputLane)
{
	lane_write_buffer<PredictionWorkItem> writeBuffer(outputLane, _laneBufferSize);
	
	PredictionWorkItem item;
	while(inputLane->read(item))
	{
		_gridder->SampleData(item.data.get(), item.dataDescId, item.u, item.v, item.w);
		
		writeBuffer.write(std::move(item));
	}
}

void WSMSGridder::predictWriteThread(ao::lane<PredictionWorkItem>* predictionWorkLane, const MSData* msData)
{
	lane_read_buffer<PredictionWorkItem> buffer(predictionWorkLane, std::min(_laneBufferSize, predictionWorkLane->capacity()));
	PredictionWorkItem workItem;
	while(buffer.read(workItem))
	{
		msData->msProvider->WriteModel(workItem.rowId, workItem.data.get());
	}
}

void WSMSGridder::Invert()
{
	std::vector<MSData> msDataVector;
	initializeMSDataVector(msDataVector, 1);
	
	_gridder = std::unique_ptr<WStackingGridder>(new WStackingGridder(_actualInversionWidth, _actualInversionHeight, _actualPixelSizeX, _actualPixelSizeY, _cpuCount, _imageBufferAllocator, AntialiasingKernelSize(), OverSamplingFactor()));
	_gridder->SetGridMode(GridMode());
	if(HasDenormalPhaseCentre())
		_gridder->SetDenormalPhaseCentre(PhaseCentreDL(), PhaseCentreDM());
	_gridder->SetIsComplex(IsComplex());
	//_imager->SetImageConjugatePart(Polarization() == Polarization::YX && IsComplex());
	_gridder->PrepareWLayers(WGridSize(), double(_memSize)*(7.0/10.0), _minW, _maxW);
	
	if(Verbose() && Logger::IsVerbose())
	{
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
			countSamplesPerLayer(msDataVector[i]);
	}
	
	resetVisibilityCounters();
	for(size_t pass=0; pass!=_gridder->NPasses(); ++pass)
	{
		Logger::Info << "Gridding pass " << pass << "... ";
		if(Verbose()) Logger::Info << '\n';
		else Logger::Info.Flush();
		
		//_inversionWorkLane.reset(new ao::lane<InversionWorkItem>(2048));
		//set_lane_debug_name(*_inversionWorkLane, "Inversion work lane containing full row data");
		
		_gridder->StartInversionPass(pass);
		
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
		{
			//_inversionWorkLane->clear();
			
			MSData& msData = msDataVector[i];
			
			const MultiBandData selectedBand(msData.SelectedBand());
			
			startInversionWorkThreads(selectedBand.MaxChannels());
		
			gridMeasurementSet(msData);
			
			//_inversionWorkLane->write_end();
			finishInversionWorkThreads();
		}
		//_inversionWorkLane.reset();
		
		Logger::Info << "Fourier transforms...\n";
		_gridder->FinishInversionPass();
	}
	
	if(Verbose())
	{
		size_t totalRowsRead = 0, totalMatchingRows = 0;
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
		{
			totalRowsRead += msDataVector[i].totalRowsProcessed;
			totalMatchingRows += msDataVector[i].matchingRows;
		}
		
		Logger::Info << "Total rows read: " << totalRowsRead;
		if(totalMatchingRows != 0)
			Logger::Info << " (overhead: " << std::max(0.0, round(totalRowsRead * 100.0 / totalMatchingRows - 100.0)) << "%)";
		Logger::Info << '\n';
	}
	
	if(NormalizeForWeighting())
		_gridder->FinalizeImage(1.0/totalWeight(), false);
	else {
		Logger::Info << "Not dividing by normalization factor of " << totalWeight()/2.0 << ".\n";
		_gridder->FinalizeImage(2.0, true);
	}
	Logger::Info << "Gridded visibility count: " << double(GriddedVisibilityCount());
	if(Weighting().IsNatural())
		Logger::Info << ", effective count after weighting: " << EffectiveGriddedVisibilityCount();
	Logger::Info << '\n';
	
	if(ImageWidth()!=_actualInversionWidth || ImageHeight()!=_actualInversionHeight)
	{
		// Interpolate the image
		// The input is of size _actualInversionWidth x _actualInversionHeight
		FFTResampler resampler(_actualInversionWidth, _actualInversionHeight, ImageWidth(), ImageHeight(), _cpuCount);
		
		if(IsComplex())
		{
			double *resizedReal = _imageBufferAllocator->Allocate(ImageWidth() * ImageHeight());
			double *resizedImag = _imageBufferAllocator->Allocate(ImageWidth() * ImageHeight());
			resampler.Start();
			resampler.AddTask(_gridder->RealImage(), resizedReal);
			resampler.AddTask(_gridder->ImaginaryImage(), resizedImag);
			resampler.Finish();
			_gridder->ReplaceRealImageBuffer(resizedReal);
			_gridder->ReplaceImaginaryImageBuffer(resizedImag);
		}
		else {
			double *resized = _imageBufferAllocator->Allocate(ImageWidth() * ImageHeight());
			resampler.RunSingle(_gridder->RealImage(), resized);
			_gridder->ReplaceRealImageBuffer(resized);
		}
	}
	
	if(TrimWidth() != ImageWidth() || TrimHeight() != ImageHeight())
	{
		Logger::Info << "Trimming " << ImageWidth() << " x " << ImageHeight() << " -> " << TrimWidth() << " x " << TrimHeight() << '\n';
		// Perform trimming
		
		double *trimmed = _imageBufferAllocator->Allocate(TrimWidth() * TrimHeight());
		Image::Trim(trimmed, TrimWidth(), TrimHeight(), _gridder->RealImage(), ImageWidth(), ImageHeight());
		_gridder->ReplaceRealImageBuffer(trimmed);
		
		if(IsComplex())
		{
			double *trimmedImag = _imageBufferAllocator->Allocate(TrimWidth() * TrimHeight());
			Image::Trim(trimmedImag, TrimWidth(), TrimHeight(), _gridder->ImaginaryImage(), ImageWidth(), ImageHeight());
			_gridder->ReplaceImaginaryImageBuffer(trimmedImag);
		}
	}
}

void WSMSGridder::Predict(double* real, double* imaginary)
{
	if(imaginary==0 && IsComplex())
		throw std::runtime_error("Missing imaginary in complex prediction");
	if(imaginary!=0 && !IsComplex())
		throw std::runtime_error("Imaginary specified in non-complex prediction");
	
	std::vector<MSData> msDataVector;
	initializeMSDataVector(msDataVector, 1);
	
	_gridder = std::unique_ptr<WStackingGridder>(new WStackingGridder(_actualInversionWidth, _actualInversionHeight, _actualPixelSizeX, _actualPixelSizeY, _cpuCount, _imageBufferAllocator, AntialiasingKernelSize(), OverSamplingFactor()));
	_gridder->SetGridMode(GridMode());
	if(HasDenormalPhaseCentre())
		_gridder->SetDenormalPhaseCentre(PhaseCentreDL(), PhaseCentreDM());
	_gridder->SetIsComplex(IsComplex());
	//_imager->SetImageConjugatePart(Polarization() == Polarization::YX && IsComplex());
	_gridder->PrepareWLayers(WGridSize(), double(_memSize)*(7.0/10.0), _minW, _maxW);
	
	if(Verbose())
	{
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
			countSamplesPerLayer(msDataVector[i]);
	}
	
	ImageBufferAllocator::Ptr untrimmedReal, untrimmedImag;
	if(TrimWidth() != ImageWidth() || TrimHeight() != ImageHeight())
	{
		Logger::Info << "Untrimming " << TrimWidth() << " x " << TrimHeight() << " -> " << ImageWidth() << " x " << ImageHeight() << '\n';
		// Undo trimming (i.e., extend with zeros)
		// The input is of size TrimWidth() x TrimHeight()
		// This will make the model image of size ImageWidth() x ImageHeight()
		_imageBufferAllocator->Allocate(ImageWidth() * ImageHeight(), untrimmedReal);
		Image::Untrim(untrimmedReal.data(), ImageWidth(), ImageHeight(), real, TrimWidth(), TrimHeight());
		real = untrimmedReal.data();
		
		if(IsComplex())
		{
			_imageBufferAllocator->Allocate(ImageWidth() * ImageHeight(), untrimmedImag);
			Image::Untrim(untrimmedImag.data(), ImageWidth(), ImageHeight(), imaginary, TrimWidth(), TrimHeight());
			imaginary = untrimmedImag.data();
		}
	}
	
	ImageBufferAllocator::Ptr resampledReal, resampledImag;
	if(ImageWidth()!=_actualInversionWidth || ImageHeight()!=_actualInversionHeight)
	{
		// Decimate the image
		// Input is ImageWidth() x ImageHeight()
		FFTResampler resampler(ImageWidth(), ImageHeight(), _actualInversionWidth, _actualInversionHeight, _cpuCount);
		
		_imageBufferAllocator->Allocate(ImageWidth() * ImageHeight(), resampledReal);
		if(imaginary == 0)
		{
			resampler.RunSingle(real, resampledReal.data());
		}
		else {
			_imageBufferAllocator->Allocate(ImageWidth() * ImageHeight(), resampledImag);
			resampler.Start();
			resampler.AddTask(real, resampledReal.data());
			resampler.AddTask(imaginary, resampledImag.data());
			resampler.Finish();
			imaginary = resampledImag.data();
		}
		real = resampledReal.data();
	}
	
	for(size_t pass=0; pass!=_gridder->NPasses(); ++pass)
	{
		Logger::Info << "Fourier transforms for pass " << pass << "... ";
		if(Verbose()) Logger::Info << '\n';
		else Logger::Info.Flush();
		if(imaginary == 0)
			_gridder->InitializePrediction(real);
		else
			_gridder->InitializePrediction(real, imaginary);
		
		_gridder->StartPredictionPass(pass);
		
		Logger::Info << "Predicting...\n";
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
			predictMeasurementSet(msDataVector[i]);
	}
	
	resampledReal.reset();
	resampledImag.reset();
	untrimmedReal.reset();
	untrimmedImag.reset();
	
	size_t totalRowsWritten = 0, totalMatchingRows = 0;
	for(size_t i=0; i!=MeasurementSetCount(); ++i)
	{
		totalRowsWritten += msDataVector[i].totalRowsProcessed;
		totalMatchingRows += msDataVector[i].matchingRows;
	}
	
	Logger::Info << "Total rows written: " << totalRowsWritten;
	if(totalMatchingRows != 0)
		Logger::Info << " (overhead: " << std::max(0.0, round(totalRowsWritten * 100.0 / totalMatchingRows - 100.0)) << "%)";
	Logger::Info << '\n';
}
