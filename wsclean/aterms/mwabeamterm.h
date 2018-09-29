#ifndef MWA_BEAM_TERM_H
#define MWA_BEAM_TERM_H

#include <thread>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MPosition.h>

#include <complex>

#include "../mwa/tilebeambase.h"
#include "../mwa/tilebeam2016.h"

#include "../lane.h"
#include "../matrix2x2.h"

#include "atermbeam.h"
#include "atermstub.h"

class MWABeamTerm : public ATermBeam
{
public:
	MWABeamTerm(casacore::MeasurementSet& ms, size_t width, size_t height, double ra, double dec, double dl, double dm, double phaseCentreDL, double phaseCentreDM);
	
private:
	virtual bool calculateBeam(std::complex<float>* buffer, double time, double frequency) final override;

	size_t _width, _height, _nStations;
	double _subbandFrequency, _phaseCentreRA, _phaseCentreDec, _dl, _dm, _phaseCentreDL, _phaseCentreDM;
	casacore::MPosition _arrayPos;
	double _delays[16];
	bool _frequencyInterpolation;
	std::unique_ptr<TileBeamBase<TileBeam2016>> _tileBeam;
};

#endif 

