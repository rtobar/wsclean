#include "wsmsgridder.h"

#include "logger.h"
#include "../buffered_lane.h"

#include "../msproviders/msprovider.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/measures/Measures/MBaseline.h>
#include <casacore/measures/Measures/MCBaseline.h>
#include <casacore/measures/Measures/MeasFrame.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/Muvw.h>
#include <casacore/measures/Measures/MPosition.h>

#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>

#include <vector>

static casacore::Muvw calculateUVW(const casacore::MPosition &antennaPos, const casacore::MPosition &refPos, const casacore::MEpoch &time, const casacore::MDirection &direction)
{
	const casacore::Vector<double> posVec = antennaPos.getValue().getVector();
	const casacore::Vector<double> refVec = refPos.getValue().getVector();
	casacore::MVPosition relativePos(posVec[0]-refVec[0], posVec[1]-refVec[1], posVec[2]-refVec[2]);
	casacore::MeasFrame frame(time, refPos, direction);
	casacore::MBaseline baseline(casacore::MVBaseline(relativePos), casacore::MBaseline::Ref(casacore::MBaseline::ITRF, frame));
	casacore::MBaseline j2000Baseline = casacore::MBaseline::Convert(baseline, casacore::MBaseline::J2000)();
	casacore::MVuvw uvw(j2000Baseline.getValue(), direction.getValue());
	return casacore::Muvw(uvw, casacore::Muvw::J2000);
}

/**
 * This is a hacked version of WSMSGridder::gridMeasurementSet, that takes the
 * LOFAR_FULL_RES_FLAG column into account to determine the UVWs
 * position. It's separated from the rest because it is an experimental
 * test.
 */
void WSMSGridder::gridLOFARCentroidMeasurementSet(MSData& msData)
{
	const MultiBandData selectedBand(msData.SelectedBand());
	_gridder->PrepareBand(selectedBand);
	ao::uvector<std::complex<float>> modelBuffer(selectedBand.MaxChannels());
	ao::uvector<float> weightBuffer(selectedBand.MaxChannels());
	ao::uvector<bool> isSelected(selectedBand.MaxChannels());
	
	casacore::ArrayColumn<unsigned char>
		fullresFlagCol(msData.msProvider->MS(), "LOFAR_FULL_RES_FLAG");
	casacore::ScalarColumn<double>
		intervalCol(msData.msProvider->MS(), casacore::MeasurementSet::columnName(casacore::MeasurementSet::INTERVAL));
	casacore::IPosition fullresShape = fullresFlagCol.shape(0);
	size_t freqFlagFactor = fullresShape[0], timeFlagFactor = fullresShape[1];
	casacore::Array<unsigned char> fullResArr(fullresShape);
	Logger::Debug << "Applying centroids from full res of " << _fullResWidth << " channels (" << freqFlagFactor*8 << " bits) and " << timeFlagFactor << " timesteps.\n";
	
	// Samples of the same w-layer are collected in a buffer
	// before they are written into the lane. This is done because writing
	// to a lane is reasonably slow; it requires holding a mutex. Without
	// these buffers, writing the lane was a bottleneck and multithreading
	// did not help. I think.
	std::vector<lane_write_buffer<InversionWorkSample>> bufferedLanes(_cpuCount);
	size_t bufferSize = std::max<size_t>(8u, _inversionCPULanes[0].capacity()/8);
	bufferSize = std::min<size_t>(128, std::min(bufferSize, _inversionCPULanes[0].capacity()));
	for(size_t i=0; i!=_cpuCount; ++i)
	{
		bufferedLanes[i].reset(&_inversionCPULanes[i], bufferSize);
	}
	
	InversionRow newItem;
	ao::uvector<std::complex<float>> newItemData(selectedBand.MaxChannels());
	newItem.data = newItemData.data();
	
	// Read antenna table
	casacore::MSAntenna antTable = msData.msProvider->MS().antenna();
	size_t nAntenna = antTable.nrow();
	std::vector<casacore::MPosition> antennaPositions(nAntenna);
	casacore::MPosition::ScalarColumn
		posCol(antTable, casacore::MSAntenna::columnName(casacore::MSAntenna::POSITION));
	for(size_t row=0; row!=nAntenna; ++row)
		antennaPositions[row] = posCol(row);
	
	// Read phase center
	casacore::MSField fieldTable = msData.msProvider->MS().field();
	casacore::MDirection::ScalarColumn
		dirCol(fieldTable, casacore::MSField::columnName(casacore::MSField::PHASE_DIR));
	casacore::MDirection direction = dirCol(0);
	
	// Array of uvws[antenna][timestep]
	std::vector<casacore::MVuvw> medianUVWs(nAntenna);
	std::vector<std::vector<casacore::MVuvw>> uvws(nAntenna);
	std::vector<std::array<double, 3>> channelCentroids;
	std::vector<size_t> channelSumCounts;
		
	std::vector<size_t> idToMSRow;
	msData.msProvider->MakeIdToMSRowMapping(idToMSRow);
	
	double uDeltaRMS = 0.0, vDeltaRMS = 0.0, wDeltaRMS = 0.0;
	size_t deltaCount = 0;
	
	size_t rowsRead = 0;
	msData.msProvider->Reset();
	double curTime = -1.0;
	while(msData.msProvider->CurrentRowAvailable())
	{
		MSProvider::MetaData metaData;
		msData.msProvider->ReadMeta(metaData);
		const size_t
			dataDescId = metaData.dataDescId;
		const double
			uInMeters = metaData.uInM,
			vInMeters = metaData.vInM,
			wInMeters = metaData.wInM;
		const BandData& curBand(selectedBand[dataDescId]);
		
		size_t row = idToMSRow[msData.msProvider->RowId()];
		fullresFlagCol.get(row, fullResArr, false);
		if(curTime != metaData.time)
		{
			curTime = metaData.time;
			std::vector<casacore::MEpoch> timesteps;
			timesteps.reserve(timeFlagFactor);
			double interval = intervalCol(row);
			double start = metaData.time - interval*0.5;
			for(size_t t=0; t!=timeFlagFactor; ++t)
			{
				double timestep = start + interval*(double(t)+0.5)/timeFlagFactor;
				timesteps.emplace_back(casacore::MVEpoch(timestep/86400.0), casacore::MEpoch::UTC);
				// Logger::Debug << timesteps.back() << '\n';
			}
			
			// Calculate high-res UVWs
			for(size_t antIndex=0; antIndex!=nAntenna; ++antIndex)
			{
				uvws[antIndex].resize(timeFlagFactor);
				for(size_t timeIndex=0; timeIndex!=timeFlagFactor; ++timeIndex)
				{
					uvws[antIndex][timeIndex] =
						calculateUVW(antennaPositions[antIndex], antennaPositions[0], timesteps[timeIndex], direction).getValue();
				}
				casacore::MEpoch curEpoch(casacore::MVEpoch(curTime/86400.0), casacore::MEpoch::UTC);
				medianUVWs[antIndex] = calculateUVW(antennaPositions[antIndex], antennaPositions[0], curEpoch, direction).getValue();
			}
		}
		
		size_t nChan = curBand.ChannelCount();
		channelCentroids.assign(nChan, std::array<double, 3>{{ 0.0, 0.0, 0.0 }});
		channelSumCounts.assign(nChan, 0);
		const unsigned char* data = fullResArr.cbegin();
		size_t bit = 0;
		for(size_t t=0; t!=timeFlagFactor; ++t)
		{
			for(size_t fc=0; fc!=nChan; ++fc)
			{
				double
					width = curBand.ChannelWidth(fc),
					freq0 = curBand.ChannelFrequency(fc) - width*0.5;
				size_t fco = fc*(_fullResWidth/nChan);
				for(size_t f=0; f!=_fullResWidth/nChan; ++f)
				{
					// freq0 is at the start of the averaged channel. In case channels were skipped during averaging, freq0 is
					// past the skipped channels. 
					double wavelength = BandData::FrequencyToLambda(freq0 + width * (double(fco+f)+0.5) / _fullResWidth);
					//Logger::Debug << "Freq=" << freq + width * (double(fco+f)+0.5) / _fullResWidth << '\n';
					bool value = *data & (1<<bit);
					//Logger::Debug << (value ? 'X' : ' ');
					if(value)
					{
						channelCentroids[fc][0] += (uvws[metaData.antenna2][t].getVector()[0]
							- uvws[metaData.antenna1][t].getVector()[0]) / wavelength;
						channelCentroids[fc][1] += (uvws[metaData.antenna2][t].getVector()[1]
							- uvws[metaData.antenna1][t].getVector()[1]) / wavelength;
						channelCentroids[fc][2] += (uvws[metaData.antenna2][t].getVector()[2]
							- uvws[metaData.antenna1][t].getVector()[2]) / wavelength;
						channelSumCounts[fc]++;
					}
					bit++;
					if(bit == 8) {
						++data;
						bit = 0;
					}
				}
				//Logger::Debug << '\n';
			}
			if(bit!=0)
			{
				++data;
				bit = 0;
			}
		}
		for(size_t fc=0; fc!=nChan; ++fc)
		{
			double wavelength = curBand.ChannelWavelength(fc);
			double medU = (medianUVWs[metaData.antenna2].getVector()[0] - medianUVWs[metaData.antenna1].getVector()[0]) / wavelength;
			double medV = (medianUVWs[metaData.antenna2].getVector()[1] - medianUVWs[metaData.antenna1].getVector()[1]) / wavelength;
			double medW = (medianUVWs[metaData.antenna2].getVector()[2] - medianUVWs[metaData.antenna1].getVector()[2]) / wavelength;
			
			// We calculate the delta u,v,w and add it to the measurement uvw. This is because
			// measurement sets sometimes not have accurate uvw values, but if they were already
			// used for calibration, we shouldn't change them "absolutely", but rather shift them
			// relatively to match the centroid.
			if(channelSumCounts[fc] == 0)
			{
				for(double& v : channelCentroids[fc])
					v = 0.0;
			}
			else {
				channelCentroids[fc][0] = channelCentroids[fc][0] / double(channelSumCounts[fc]) - medU;
				channelCentroids[fc][1] = channelCentroids[fc][1] / double(channelSumCounts[fc]) - medV;
				channelCentroids[fc][2] = channelCentroids[fc][2] / double(channelSumCounts[fc]) - medW;
			}
			uDeltaRMS += channelCentroids[fc][0]*channelCentroids[fc][0];
			vDeltaRMS += channelCentroids[fc][1]*channelCentroids[fc][1];
			wDeltaRMS += channelCentroids[fc][2]*channelCentroids[fc][2];
			++deltaCount;
			/*Logger::Debug << "["
				<< uInMeters / wavelength << " , "
				<< vInMeters / wavelength << " , "
				<< wInMeters / wavelength << "] -> ["
				<< medU << " , "
				<< medV << " , "
				<< medW << "] -> ["
				<< channelCentroids[fc][0] << " , "
				<< channelCentroids[fc][1] << " , "
				<< channelCentroids[fc][2] << "]\n";*/
		}
		
		const double
			w1 = wInMeters / curBand.LongestWavelength() + channelCentroids.front()[2],
			w2 = wInMeters / curBand.LongestWavelength() + channelCentroids.back()[2];
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
				double w = newItem.uvw[2] / curBand.ChannelWavelength(ch) + channelCentroids[ch][2];
				isSelected[ch] = _gridder->IsInLayerRange(w);
			}

			readAndWeightVisibilities<1>(*msData.msProvider, newItem, curBand, weightBuffer.data(), modelBuffer.data(), isSelected.data());
			
			InversionWorkSample sampleData;
			for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
			{
				double wavelength = curBand.ChannelWavelength(ch);
				sampleData.sample = newItem.data[ch];
				sampleData.uInLambda = newItem.uvw[0] / wavelength + channelCentroids[ch][0];
				sampleData.vInLambda = newItem.uvw[1] / wavelength + channelCentroids[ch][1];
				sampleData.wInLambda = newItem.uvw[2] / wavelength + channelCentroids[ch][2];
				size_t cpu = _gridder->WToLayer(sampleData.wInLambda) % _cpuCount;
				bufferedLanes[cpu].write(sampleData);
			}
			
			++rowsRead;
		}
		
		msData.msProvider->NextRow();
	}
	
	for(lane_write_buffer<InversionWorkSample>& buflane : bufferedLanes)
		buflane.write_end();
	
	if(Verbose())
		Logger::Info << "Rows that were required: " << rowsRead << '/' << msData.matchingRows << '\n';
	msData.totalRowsProcessed += rowsRead;
	
	uDeltaRMS = sqrt(uDeltaRMS/deltaCount);
	vDeltaRMS = sqrt(vDeltaRMS/deltaCount);
	wDeltaRMS = sqrt(wDeltaRMS/deltaCount);
	Logger::Info << "RMS of centroid deltas: " << uDeltaRMS << ", " << vDeltaRMS << ", " << wDeltaRMS << " wavelengths.\n";
}

