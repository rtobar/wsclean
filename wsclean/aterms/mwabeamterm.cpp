#include "mwabeamterm.h"

#include "../units/imagecoordinates.h"

#include "../wsclean/logger.h"

#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/Measures/MEpoch.h>

MWABeamTerm::MWABeamTerm(casacore::MeasurementSet& ms, size_t width, size_t height, double ra, double dec, double dl, double dm, double phaseCentreDL, double phaseCentreDM) :
	_width(width),
	_height(height),
	_phaseCentreRA(ra), _phaseCentreDec(dec),
	_dl(dl), _dm(dm),
	_phaseCentreDL(phaseCentreDL),
	_phaseCentreDM(phaseCentreDM),
	_frequencyInterpolation(true)
{
	casacore::MSAntenna aTable = ms.antenna();
	if(aTable.nrow() == 0) throw std::runtime_error("No antennae in set");
	_nStations = aTable.nrow();
	
	casacore::MPosition::ROScalarColumn antPosColumn(aTable, aTable.columnName(casacore::MSAntennaEnums::POSITION));
	_arrayPos = antPosColumn(0);
		
	casacore::Table mwaTilePointing = ms.keywordSet().asTable("MWA_TILE_POINTING");
	casacore::ROArrayColumn<int> delaysCol(mwaTilePointing, "DELAYS");
	casacore::Array<int> delaysArr = delaysCol(0);
	casacore::Array<int>::contiter delaysArrPtr = delaysArr.cbegin();
	for(int i=0; i!=16; ++i)
		_delays[i] = delaysArrPtr[i];
	
	Logger::Debug << "MWA beam delays: [";
	for(int i=0; i!=16; ++i)
	{
		Logger::Debug << _delays[i];
		if(i != 15) Logger::Debug << ',';
	}
	Logger::Debug << "]\n";
}

bool MWABeamTerm::calculateBeam(std::complex<float>* buffer, double time, double frequency)
{
	casacore::MEpoch timeEpoch(casacore::Quantity(time, "s"));
	casacore::MeasFrame frame(_arrayPos, timeEpoch);
	
	const casacore::MDirection::Ref hadecRef(casacore::MDirection::HADEC, frame);
	const casacore::MDirection::Ref azelgeoRef(casacore::MDirection::AZELGEO, frame);
	const casacore::MDirection::Ref j2000Ref(casacore::MDirection::J2000, frame);
	casacore::MDirection::Convert
		j2000ToHaDecRef(j2000Ref, hadecRef),
		j2000ToAzelGeoRef(j2000Ref, azelgeoRef);
	casacore::MPosition wgs = casacore::MPosition::Convert(_arrayPos, casacore::MPosition::WGS84)();
	double arrLatitude = wgs.getValue().getLat();
	
	casacore::MDirection zenith(casacore::MVDirection(0.0, 0.0, 1.0), azelgeoRef);
	casacore::MDirection zenithHaDec = casacore::MDirection::Convert(zenith, hadecRef)();
	double zenithHa = zenithHaDec.getAngle().getValue()[0];
	double zenithDec = zenithHaDec.getAngle().getValue()[1];
	
	if(!_tileBeam)
	{
		_tileBeam.reset(new TileBeamBase<TileBeam2016>(
			_delays, _frequencyInterpolation, _searchPath));
	}
	std::complex<float>* bufferPtr = buffer;
	for(size_t y=0; y!=_height; ++y)
	{
		for(size_t x=0; x!=_width; ++x)
		{
			double l, m, ra, dec;
			ImageCoordinates::XYToLM(x, y, _dl, _dm, _width, _height, l, m);
			l += _phaseCentreDL; m += _phaseCentreDM;
			ImageCoordinates::LMToRaDec(l, m, _phaseCentreRA, _phaseCentreDec, ra, dec);
			
			std::complex<double> gain[4];
			_tileBeam->ArrayResponse(ra, dec, j2000Ref, j2000ToHaDecRef, j2000ToAzelGeoRef, arrLatitude, zenithHa, zenithDec, frequency, gain);
			
			for(size_t i=0; i!=4; ++i)
			{
				*bufferPtr = gain[i];
				++bufferPtr;
			}
		}
	}
	
	// Now, copy the beam response to all antennas
	size_t nValues = _width*_height*4;
	for(size_t station=1; station!=_nStations; ++station)
		std::copy_n(buffer, nValues, buffer + nValues*station);
	
	saveATermsIfNecessary(buffer, _nStations, _width, _height);
	
	return true;
}
