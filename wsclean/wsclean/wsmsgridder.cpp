#include "wsmsgridder.h"

#include "imagebufferallocator.h"

#include "../imageweights.h"
#include "../buffered_lane.h"
#include "../fftresampler.h"
#include "../imageoperations.h"
#include "../angle.h"

#include "../msproviders/msprovider.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>

#include <casacore/tables/Tables/ArrColDesc.h>

#include <iostream>
#include <stdexcept>

WSMSGridder::MSData::MSData() : matchingRows(0), totalRowsProcessed(0)
{ }

WSMSGridder::MSData::~MSData()
{ }

WSMSGridder::WSMSGridder(ImageBufferAllocator* imageAllocator, size_t threadCount, double memFraction, double absMemLimit) :
	InversionAlgorithm(),
	_phaseCentreRA(0.0), _phaseCentreDec(0.0),
	_phaseCentreDL(0.0), _phaseCentreDM(0.0),
	_denormalPhaseCentre(false), _hasFrequencies(false),
	_freqHigh(0.0), _freqLow(0.0),
	_bandStart(0.0), _bandEnd(0.0),
	_beamSize(0.0),
	_totalWeight(0.0),
	_startTime(0.0),
	_gridMode(WStackingGridder::NearestNeighbour),
	_cpuCount(threadCount),
	_laneBufferSize(std::max<size_t>(_cpuCount*2,1024)),
	_imageBufferAllocator(imageAllocator),
	_trimWidth(0), _trimHeight(0),
	_actualInversionWidth(0), _actualInversionHeight(0),
	_actualPixelSizeX(0), _actualPixelSizeY(0)
{
	long int pageCount = sysconf(_SC_PHYS_PAGES), pageSize = sysconf(_SC_PAGE_SIZE);
	_memSize = (int64_t) pageCount * (int64_t) pageSize;
	double memSizeInGB = (double) _memSize / (1024.0*1024.0*1024.0);
	if(memFraction == 1.0 && absMemLimit == 0.0) {
		std::cout << "Detected " << round(memSizeInGB*10.0)/10.0 << " GB of system memory, usage not limited.\n";
	}
	else {
		double limitInGB = memSizeInGB*memFraction;
		if(absMemLimit!=0.0 && limitInGB > absMemLimit)
			limitInGB = absMemLimit;
		std::cout << "Detected " << round(memSizeInGB*10.0)/10.0 << " GB of system memory, usage limited to " << round(limitInGB*10.0)/10.0 << " GB (frac=" << round(memFraction*1000.0)/10.0 << "%, ";
		if(absMemLimit == 0.0)
			std::cout << "no limit)\n";
		else
			std::cout << "limit=" << round(absMemLimit*10.0)/10.0 << "GB)\n";
		
		_memSize = int64_t((double) pageCount * (double) pageSize * memFraction);
		if(absMemLimit!=0.0 && double(_memSize) > double(1024.0*1024.0*1024.0) * absMemLimit)
			_memSize = int64_t(double(absMemLimit) * double(1024.0*1024.0*1024.0));
	}
}
		
void WSMSGridder::initializeMeasurementSet(size_t msIndex, WSMSGridder::MSData& msData)
{
	MSProvider& msProvider = MeasurementSet(msIndex);
	msData.msProvider = &msProvider;
	casacore::MeasurementSet& ms(msProvider.MS());
	if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	/**
		* Read some meta data from the measurement set
		*/
	casacore::MSAntenna aTable = ms.antenna();
	size_t antennaCount = aTable.nrow();
	if(antennaCount == 0) throw std::runtime_error("No antennae in set");
	casacore::MPosition::ROScalarColumn antPosColumn(aTable, aTable.columnName(casacore::MSAntennaEnums::POSITION));
	casacore::MPosition ant1Pos = antPosColumn(0);
	
	msData.bandData = MultiBandData(ms.spectralWindow(), ms.dataDescription());
	if(Selection(msIndex).HasChannelRange())
	{
		msData.startChannel = Selection(msIndex).ChannelRangeStart();
		msData.endChannel = Selection(msIndex).ChannelRangeEnd();
		std::cout << "Selected channels: " << msData.startChannel << '-' << msData.endChannel << '\n';
		const BandData& firstBand = msData.bandData.FirstBand();
		if(msData.startChannel >= firstBand.ChannelCount() || msData.endChannel > firstBand.ChannelCount()
			|| msData.startChannel == msData.endChannel)
		{
			std::ostringstream str;
			str << "An invalid channel range was specified! Measurement set only has " << firstBand.ChannelCount() << " channels, requested imaging range is " << msData.startChannel << " -- " << msData.endChannel << '.';
			throw std::runtime_error(str.str());
		}
	}
	else {
		msData.startChannel = 0;
		msData.endChannel = msData.bandData.FirstBand().ChannelCount();
	}
	casacore::MEpoch::ROScalarColumn timeColumn(ms, ms.columnName(casacore::MSMainEnums::TIME));
	const MultiBandData selectedBand = msData.SelectedBand();
	if(_hasFrequencies)
	{
		_freqLow = std::min(_freqLow, selectedBand.LowestFrequency());
		_freqHigh = std::max(_freqHigh, selectedBand.HighestFrequency());
		_bandStart = std::min(_bandStart, selectedBand.BandStart());
		_bandEnd = std::max(_bandEnd, selectedBand.BandEnd());
		_startTime = std::min(_startTime, msProvider.StartTime());
	} else {
		_freqLow = selectedBand.LowestFrequency();
		_freqHigh = selectedBand.HighestFrequency();
		_bandStart = selectedBand.BandStart();
		_bandEnd = selectedBand.BandEnd();
		_startTime = msProvider.StartTime();
		_hasFrequencies = true;
	}
	
	casacore::MSField fTable(ms.field());
	casacore::MDirection::ROScalarColumn phaseDirColumn(fTable, fTable.columnName(casacore::MSFieldEnums::PHASE_DIR));
	casacore::MDirection phaseDir = phaseDirColumn(Selection(msIndex).FieldId());
	casacore::MEpoch curtime = timeColumn(0);
	casacore::MeasFrame frame(ant1Pos, curtime);
	casacore::MDirection::Ref j2000Ref(casacore::MDirection::J2000, frame);
	casacore::MDirection j2000 = casacore::MDirection::Convert(phaseDir, j2000Ref)();
	casacore::Vector<casacore::Double> j2000Val = j2000.getValue().get();
	_phaseCentreRA = j2000Val[0];
	_phaseCentreDec = j2000Val[1];
	if(fTable.keywordSet().isDefined("WSCLEAN_DL"))
		_phaseCentreDL = fTable.keywordSet().asDouble(casacore::RecordFieldId("WSCLEAN_DL"));
	else _phaseCentreDL = 0.0;
	if(fTable.keywordSet().isDefined("WSCLEAN_DM"))
		_phaseCentreDM = fTable.keywordSet().asDouble(casacore::RecordFieldId("WSCLEAN_DM"));
	else _phaseCentreDM = 0.0;

	_denormalPhaseCentre = _phaseCentreDL != 0.0 || _phaseCentreDM != 0.0;
	if(_denormalPhaseCentre)
		std::cout << "Set has denormal phase centre: dl=" << _phaseCentreDL << ", dm=" << _phaseCentreDM << '\n';
	
	std::cout << "Determining min and max w & theoretical beam size... " << std::flush;
	msData.maxW = 0.0;
	msData.minW = 1e100;
	double maxBaseline = 0.0;
	std::vector<float> weightArray(selectedBand.MaxChannels());
	msProvider.Reset();
	while(msProvider.CurrentRowAvailable())
	{
		size_t dataDescId;
		double uInM, vInM, wInM;
		msProvider.ReadMeta(uInM, vInM, wInM, dataDescId);
		const BandData& curBand = selectedBand[dataDescId];
		double wHi = fabs(wInM / curBand.SmallestWavelength());
		double wLo = fabs(wInM / curBand.LongestWavelength());
		double baselineInM = sqrt(uInM*uInM + vInM*vInM + wInM*wInM);
		double halfWidth = 0.5*ImageWidth(), halfHeight = 0.5*ImageHeight();
		if(wHi > msData.maxW || wLo < msData.minW || baselineInM / curBand.SmallestWavelength() > maxBaseline)
		{
			msProvider.ReadWeights(weightArray.data());
			const float* weightPtr = weightArray.data();
			for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
			{
				if(*weightPtr != 0.0)
				{
					const double wavelength = curBand.ChannelWavelength(ch);
					double
						uInL = uInM/wavelength, vInL = vInM/wavelength,
						wInL = wInM/wavelength,
						x = uInL * PixelSizeX() * ImageWidth(),
						y = vInL * PixelSizeY() * ImageHeight(),
						imagingWeight = this->PrecalculatedWeightInfo()->GetWeight(uInL, vInL);
					if(imagingWeight != 0.0)
					{
						if(floor(x) > -halfWidth  && ceil(x) < halfWidth &&
							floor(y) > -halfHeight && ceil(y) < halfHeight)
						{
							msData.maxW = std::max(msData.maxW, fabs(wInL));
							msData.minW = std::min(msData.minW, fabs(wInL));
							maxBaseline = std::max(maxBaseline, baselineInM / wavelength);
						}
					}
				}
				++weightPtr;
			}
		}
		
		msProvider.NextRow();
	}
	if(msData.minW == 1e100)
	{
		msData.minW = 0.0;
		msData.maxW = 0.0;
	}
	_beamSize = 1.0 / maxBaseline;
	std::cout << "DONE (w=[" << msData.minW << ":" << msData.maxW << "] lambdas, maxuvw=" << maxBaseline << " lambda, beam=" << Angle::ToNiceString(_beamSize) << ")\n";
	if(HasWLimit()) {
		msData.maxW *= (1.0 - WLimit());
		if(msData.maxW < msData.minW) msData.maxW = msData.minW;
	}

	if(_trimWidth == 0 || _trimHeight ==  0)
	{
		_trimWidth = ImageWidth();
		_trimHeight = ImageHeight();
	}
	
	_actualInversionWidth = ImageWidth();
	_actualInversionHeight = ImageHeight();
	_actualPixelSizeX = PixelSizeX();
	_actualPixelSizeY = PixelSizeY();
	
	if(SmallInversion())
	{
		double totalWidth = _actualInversionWidth * _actualPixelSizeX, totalHeight = _actualInversionHeight * _actualPixelSizeY;
		// Calc min res based on Nyquist sampling rate
		size_t minResX = size_t(ceil(totalWidth*2 / _beamSize));
		if(minResX%4 != 0) minResX += 4 - (minResX%4);
		size_t minResY = size_t(ceil(totalHeight*2 / _beamSize));
		if(minResY%4 != 0) minResY += 4 - (minResY%4);
		if(minResX < _actualInversionWidth || minResY < _actualInversionHeight)
		{
			_actualInversionWidth = std::max(std::min(minResX, _actualInversionWidth), size_t(32));
			_actualInversionHeight = std::max(std::min(minResY, _actualInversionHeight), size_t(32));
			std::cout << "Setting small inversion image size of " << _actualInversionWidth << " x " << _actualInversionHeight << "\n";
			_actualPixelSizeX = totalWidth / _actualInversionWidth;
			_actualPixelSizeY = totalHeight / _actualInversionHeight;
		}
		else {
			std::cout << "Small inversion enabled, but inversion resolution already smaller than beam size: not using optimization.\n";
		}
	}
	
	if(Verbose() || !HasWGridSize())
	{
		double
			maxL = _trimWidth * PixelSizeX() * 0.5 + fabs(_phaseCentreDL),
			maxM = _trimHeight * PixelSizeY() * 0.5 + fabs(_phaseCentreDM),
			lmSq = maxL * maxL + maxM * maxM;
		double cMinW = IsComplex() ? -msData.maxW : msData.minW;
		double radiansForAllLayers;
		if(lmSq < 1.0)
			radiansForAllLayers = 2 * M_PI * (msData.maxW - cMinW) * (1.0 - sqrt(1.0 - lmSq));
		else
			radiansForAllLayers = 2 * M_PI * (msData.maxW - cMinW);
		size_t suggestedGridSize = size_t(ceil(radiansForAllLayers));
		if(suggestedGridSize == 0) suggestedGridSize = 1;
		if(suggestedGridSize < _cpuCount)
		{
			// When nwlayers is lower than the nr of cores, we cannot parallellize well. 
			// However, we don't want extra w-layers if we are low on mem, as that might slow down the process
			double memoryRequired = double(_cpuCount) * double(sizeof(double))*double(_actualInversionWidth*_actualInversionHeight);
			if(4.0 * memoryRequired < double(_memSize))
			{
				std::cout <<
					"The theoretically suggested number of w-layers (" << suggestedGridSize << ") is less than the number of availables\n"
					"cores (" << _cpuCount << "). Changing suggested number of w-layers to " << _cpuCount << ".\n";
				suggestedGridSize = _cpuCount;
			}
			else {
				std::cout <<
					"The theoretically suggested number of w-layers (" << suggestedGridSize << ") is less than the number of availables\n"
					"cores (" << _cpuCount << "), but there is not enough memory available to increase the number of w-layers.\n"
					"Not all cores can be used efficiently.\n";
			}
		}
		if(Verbose())
			std::cout << "Suggested number of w-layers: " << ceil(suggestedGridSize) << '\n';
		if(!HasWGridSize())
			SetWGridSize(suggestedGridSize);
	}
}

void WSMSGridder::countSamplesPerLayer(MSData& msData)
{
	std::vector<size_t> sampleCount(WGridSize());
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
				++sampleCount[wLayerIndex];
		}
		++msData.matchingRows;
		msData.msProvider->NextRow();
	}
	std::cout << "Visibility count per layer: ";
	for(std::vector<size_t>::const_iterator i=sampleCount.begin(); i!=sampleCount.end(); ++i)
	{
		std::cout << *i << ' ';
	}
	std::cout << '\n';
}

void WSMSGridder::gridMeasurementSet(MSData &msData)
{
	const MultiBandData selectedBand(msData.SelectedBand());
	_gridder->PrepareBand(selectedBand);
	std::vector<std::complex<float>> modelBuffer(selectedBand.MaxChannels());
	std::vector<float> weightBuffer(selectedBand.MaxChannels());
	
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
	
	InversionWorkItem newItem;
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
			newItem.u = uInMeters;
			newItem.v = vInMeters;
			newItem.w = wInMeters;
			newItem.dataDescId = dataDescId;
			
			if(DoImagePSF())
			{
				msData.msProvider->ReadWeights(newItem.data);
				if(_denormalPhaseCentre)
				{
					double lmsqrt = sqrt(1.0-_phaseCentreDL*_phaseCentreDL- _phaseCentreDM*_phaseCentreDM);
					double shiftFactor = 2.0*M_PI* (newItem.w * (lmsqrt-1.0));
					rotateVisibilities(curBand, shiftFactor, newItem.data);
				}
			}
			else {
				msData.msProvider->ReadData(newItem.data);
			}
			
			if(DoSubtractModel())
			{
				msData.msProvider->ReadModel(modelBuffer.data());
				std::complex<float>* modelIter = modelBuffer.data();
				for(std::complex<float>* iter = newItem.data; iter!=newItem.data+curBand.ChannelCount(); ++iter)
				{
					*iter -= *modelIter;
					modelIter++;
				}
			}
			msData.msProvider->ReadWeights(weightBuffer.data());
			switch(VisibilityWeightingMode())
			{
				case NormalVisibilityWeighting:
					// The MS provider has already preweighted the
					// visibilities for their weight, so we do not
					// have to do anything.
					break;
				case SquaredVisibilityWeighting:
					for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
						newItem.data[ch] *= weightBuffer[ch];
					break;
				case UnitVisibilityWeighting:
					for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
					{
						if(weightBuffer[ch] == 0.0)
							newItem.data[ch] = 0.0;
						else
							newItem.data[ch] /= weightBuffer[ch];
					}
					break;
			}
			switch(Weighting().Mode())
			{
				case WeightMode::UniformWeighted:
				case WeightMode::BriggsWeighted:
				case WeightMode::NaturalWeighted:
				{
					std::complex<float>* dataIter = newItem.data;
					float* weightIter = weightBuffer.data();
					for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
					{
						double
							u = newItem.u / curBand.ChannelWavelength(ch),
							v = newItem.v / curBand.ChannelWavelength(ch),
							weight = PrecalculatedWeightInfo()->GetWeight(u, v);
						*dataIter *= weight;
						_totalWeight += weight * *weightIter;
						++dataIter;
						++weightIter;
					}
				} break;
				case WeightMode::DistanceWeighted:
				{
					float* weightIter = weightBuffer.data();
					double mwaWeight = sqrt(newItem.u*newItem.u + newItem.v*newItem.v + newItem.w*newItem.w);
					for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
					{
						_totalWeight += *weightIter * mwaWeight;
						++weightIter;
					}
				} break;
			}
			
			InversionWorkSample sampleData;
			for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
			{
				double wavelength = curBand.ChannelWavelength(ch);
				sampleData.sample = newItem.data[ch];
				sampleData.uInLambda = newItem.u / wavelength;
				sampleData.vInLambda = newItem.v / wavelength;
				sampleData.wInLambda = newItem.w / wavelength;
				size_t cpu = _gridder->WToLayer(sampleData.wInLambda) % _cpuCount;
				//std::cout << cpu << ' ' << lanes[cpu].size() << '\n';
				bufferedLanes[cpu].write(sampleData);
			}
			
			++rowsRead;
		}
		
		msData.msProvider->NextRow();
	}
	
	for(size_t i=0; i!=_cpuCount; ++i)
		bufferedLanes[i].write_end();
	
	if(Verbose())
		std::cout << "Rows that were required: " << rowsRead << '/' << msData.matchingRows << '\n';
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
		newItem.data = new std::complex<float>[selectedBandData[dataIds[i]].ChannelCount()];
		newItem.rowId = rowIds[i];
				
		bufferedCalcLane.write(newItem);
	}
	if(Verbose())
		std::cout << "Rows that were required: " << rowsProcessed << '/' << msData.matchingRows << '\n';
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
		_gridder->SampleData(item.data, item.dataDescId, item.u, item.v, item.w);
		
		writeBuffer.write(item);
	}
}

void WSMSGridder::predictWriteThread(ao::lane<PredictionWorkItem>* predictionWorkLane, const MSData* msData)
{
	lane_read_buffer<PredictionWorkItem> buffer(predictionWorkLane, std::min(_laneBufferSize, predictionWorkLane->capacity()));
	PredictionWorkItem workItem;
	while(buffer.read(workItem))
	{
		msData->msProvider->WriteModel(workItem.rowId, workItem.data);
		delete[] workItem.data;
	}
}

void WSMSGridder::Invert()
{
	MSData* msDataVector = new MSData[MeasurementSetCount()];
	_hasFrequencies = false;
	for(size_t i=0; i!=MeasurementSetCount(); ++i)
		initializeMeasurementSet(i, msDataVector[i]);
	
	double minW = msDataVector[0].minW;
	double maxW = msDataVector[0].maxW;
	for(size_t i=1; i!=MeasurementSetCount(); ++i)
	{
		if(msDataVector[i].minW < minW) minW = msDataVector[i].minW;
		if(msDataVector[i].maxW > maxW) maxW = msDataVector[i].maxW;
	}
	
	_gridder = std::unique_ptr<WStackingGridder>(new WStackingGridder(_actualInversionWidth, _actualInversionHeight, _actualPixelSizeX, _actualPixelSizeY, _cpuCount, _imageBufferAllocator, AntialiasingKernelSize(), OverSamplingFactor()));
	_gridder->SetGridMode(_gridMode);
	if(_denormalPhaseCentre)
		_gridder->SetDenormalPhaseCentre(_phaseCentreDL, _phaseCentreDM);
	_gridder->SetIsComplex(IsComplex());
	//_imager->SetImageConjugatePart(Polarization() == Polarization::YX && IsComplex());
	_gridder->PrepareWLayers(WGridSize(), double(_memSize)*(7.0/10.0), minW, maxW);
	
	if(Verbose())
	{
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
			countSamplesPerLayer(msDataVector[i]);
	}
	
	_totalWeight = 0.0;
	for(size_t pass=0; pass!=_gridder->NPasses(); ++pass)
	{
		std::cout << "Gridding pass " << pass << "... ";
		if(Verbose()) std::cout << '\n';
		else std::cout << std::flush;
		
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
		
		std::cout << "Fourier transforms...\n";
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
		
		std::cout << "Total rows read: " << totalRowsRead;
		if(totalMatchingRows != 0)
			std::cout << " (overhead: " << std::max(0.0, round(totalRowsRead * 100.0 / totalMatchingRows - 100.0)) << "%)";
		std::cout << '\n';
	}
	
	if(NormalizeForWeighting())
		_gridder->FinalizeImage(1.0/_totalWeight, false);
	else {
		std::cout << "Not dividing by normalization factor of " << _totalWeight << ".\n";
		_gridder->FinalizeImage(1.0, true);
	}
	
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
	
	if(_trimWidth != ImageWidth() || _trimHeight != ImageHeight())
	{
		std::cout << "Trimming " << ImageWidth() << " x " << ImageHeight() << " -> " << _trimWidth << " x " << _trimHeight << '\n';
		// Perform trimming
		
		double *trimmed = _imageBufferAllocator->Allocate(_trimWidth * _trimHeight);
		ImageOperations::Trim(trimmed, _trimWidth, _trimHeight, _gridder->RealImage(), ImageWidth(), ImageHeight());
		_gridder->ReplaceRealImageBuffer(trimmed);
		
		if(IsComplex())
		{
			double *trimmedImag = _imageBufferAllocator->Allocate(_trimWidth * _trimHeight);
			ImageOperations::Trim(trimmedImag, _trimWidth, _trimHeight, _gridder->ImaginaryImage(), ImageWidth(), ImageHeight());
			_gridder->ReplaceImaginaryImageBuffer(trimmedImag);
		}
	}
	
	delete[] msDataVector;
}

void WSMSGridder::Predict(double* real, double* imaginary)
{
	if(imaginary==0 && IsComplex())
		throw std::runtime_error("Missing imaginary in complex prediction");
	if(imaginary!=0 && !IsComplex())
		throw std::runtime_error("Imaginary specified in non-complex prediction");
	
	MSData* msDataVector = new MSData[MeasurementSetCount()];
	_hasFrequencies = false;
	for(size_t i=0; i!=MeasurementSetCount(); ++i)
		initializeMeasurementSet(i, msDataVector[i]);
	
	double minW = msDataVector[0].minW;
	double maxW = msDataVector[0].maxW;
	for(size_t i=1; i!=MeasurementSetCount(); ++i)
	{
		if(msDataVector[i].minW < minW) minW = msDataVector[i].minW;
		if(msDataVector[i].maxW > maxW) maxW = msDataVector[i].maxW;
	}
	
	_gridder = std::unique_ptr<WStackingGridder>(new WStackingGridder(_actualInversionWidth, _actualInversionHeight, _actualPixelSizeX, _actualPixelSizeY, _cpuCount, _imageBufferAllocator, AntialiasingKernelSize(), OverSamplingFactor()));
	_gridder->SetGridMode(_gridMode);
	if(_denormalPhaseCentre)
		_gridder->SetDenormalPhaseCentre(_phaseCentreDL, _phaseCentreDM);
	_gridder->SetIsComplex(IsComplex());
	//_imager->SetImageConjugatePart(Polarization() == Polarization::YX && IsComplex());
	_gridder->PrepareWLayers(WGridSize(), double(_memSize)*(7.0/10.0), minW, maxW);
	
	if(Verbose())
	{
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
			countSamplesPerLayer(msDataVector[i]);
	}
	
	ImageBufferAllocator::Ptr untrimmedReal, untrimmedImag;
	if(_trimWidth != ImageWidth() || _trimHeight != ImageHeight())
	{
		std::cout << "Untrimming " << _trimWidth << " x " << _trimHeight << " -> " << ImageWidth() << " x " << ImageHeight() << '\n';
		// Undo trimming (i.e., extend with zeros)
		// The input is of size _trimWidth x _trimHeight
		// This will make the model image of size ImageWidth() x ImageHeight()
		_imageBufferAllocator->Allocate(ImageWidth() * ImageHeight(), untrimmedReal);
		ImageOperations::Untrim(untrimmedReal.data(), ImageWidth(), ImageHeight(), real, _trimWidth, _trimHeight);
		real = untrimmedReal.data();
		
		if(IsComplex())
		{
			_imageBufferAllocator->Allocate(ImageWidth() * ImageHeight(), untrimmedImag);
			ImageOperations::Untrim(untrimmedImag.data(), ImageWidth(), ImageHeight(), imaginary, _trimWidth, _trimHeight);
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
		std::cout << "Fourier transforms for pass " << pass << "... ";
		if(Verbose()) std::cout << '\n';
		else std::cout << std::flush;
		if(imaginary == 0)
			_gridder->InitializePrediction(real);
		else
			_gridder->InitializePrediction(real, imaginary);
		
		_gridder->StartPredictionPass(pass);
		
		std::cout << "Predicting...\n";
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
	
	std::cout << "Total rows written: " << totalRowsWritten;
	if(totalMatchingRows != 0)
		std::cout << " (overhead: " << std::max(0.0, round(totalRowsWritten * 100.0 / totalMatchingRows - 100.0)) << "%)";
	std::cout << '\n';
	delete[] msDataVector;
}

void WSMSGridder::rotateVisibilities(const BandData &bandData, double shiftFactor, std::complex<float>* dataIter)
{
	for(unsigned ch=0; ch!=bandData.ChannelCount(); ++ch)
	{
		const double wShiftRad = shiftFactor / bandData.ChannelWavelength(ch);
		double rotSinD, rotCosD;
		sincos(wShiftRad, &rotSinD, &rotCosD);
		float rotSin = rotSinD, rotCos = rotCosD;
		std::complex<float> v = *dataIter;
		*dataIter = std::complex<float>(
			v.real() * rotCos  -  v.imag() * rotSin,
			v.real() * rotSin  +  v.imag() * rotCos);
		++dataIter;
	}
}
