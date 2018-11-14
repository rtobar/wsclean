#ifndef LOFAR_BEAM_TERM_H
#define LOFAR_BEAM_TERM_H

#include <thread>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <complex>

#include "atermstub.h"
#include "atermbeam.h"

#ifdef HAVE_LOFAR_BEAM

#include <StationResponse/LofarMetaDataUtil.h>

#include "../lane.h"
#include "../matrix2x2.h"

using namespace LOFAR::StationResponse;

class LofarBeamTerm : public ATermBeam
{
public:
	LofarBeamTerm(casacore::MeasurementSet& ms, size_t width, size_t height, double dl, double dm, double phaseCentreDL, double phaseCentreDM, const std::string& dataColumnName);
	
	void SetUseDifferentialBeam(bool useDiffBeam)
	{
		_useDifferentialBeam = useDiffBeam;
	}
	
	void SetUseChannelFrequency(bool useChannelFrequency)
	{
		_useChannelFrequency = useChannelFrequency;
	}
	
private:
	bool calculateBeam(std::complex<float>* buffer, double time, double frequency) final override;

	void calcThread(std::complex<float>* buffer, double time, double frequency);
	
	std::vector<LOFAR::StationResponse::Station::Ptr> _stations;
	size_t _width, _height;
	double _subbandFrequency, _phaseCentreRA, _phaseCentreDec, _dl, _dm, _phaseCentreDL, _phaseCentreDM;
	casacore::MDirection _delayDir, _preappliedBeamDir, _tileBeamDir;
	casacore::MPosition _arrayPos;
	bool _useDifferentialBeam, _useChannelFrequency;
	vector3r_t _l_vector_itrf;
	vector3r_t _m_vector_itrf;
	vector3r_t _n_vector_itrf;
	std::vector<MC2x2F> _inverseCentralGain;
	vector3r_t _station0, _tile0;
	
	ao::lane<size_t> *_lane;
	size_t _nThreads;
	std::vector<std::thread> _threads;
};

#else
using LofarBeamTerm = ATermStub;
#endif // HAVE_LOFAR_BEAM

#endif 
