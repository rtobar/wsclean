#ifndef LOFAR_BEAM_TERM_H
#define LOFAR_BEAM_TERM_H

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <complex>

#include "atermstub.h"

#ifdef HAVE_LOFAR_BEAM

#include <StationResponse/LofarMetaDataUtil.h>

using namespace LOFAR::StationResponse;

class LofarBeamTerm
{
public:
	LofarBeamTerm(casacore::MeasurementSet& ms, size_t width, size_t height, double dl, double dm, double phaseCentreDL, double phaseCentreDM, bool useDifferentialBeam);
	
	void Calculate(std::complex<float>* buffer, double time, double frequency);
	
private:
	void calcThread(struct LofarBeamTermThreadData* data);
	
	std::vector<LOFAR::StationResponse::Station::Ptr> _stations;
	size_t _width, _height;
	double _subbandFrequency, _phaseCentreRA, _phaseCentreDec, _dl, _dm, _phaseCentreDL, _phaseCentreDM;
	casacore::MDirection _delayDir, _referenceDir, _tileBeamDir;
	casacore::MPosition _arrayPos;
	bool _useDifferentialBeam;
};

#else
using LofarBeamTerm = ATermStub;
#endif // HAVE_LOFAR_BEAM

#endif 
