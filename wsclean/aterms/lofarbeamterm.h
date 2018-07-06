#ifndef LOFAR_BEAM_TERM_H
#define LOFAR_BEAM_TERM_H

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <complex>

#include "atermstub.h"
#include "atermbase.h"

#ifdef HAVE_LOFAR_BEAM

#include <StationResponse/LofarMetaDataUtil.h>

using namespace LOFAR::StationResponse;

class LofarBeamTerm : public ATermBase
{
public:
	LofarBeamTerm(casacore::MeasurementSet& ms, size_t width, size_t height, double dl, double dm, double phaseCentreDL, double phaseCentreDM, bool useDifferentialBeam);
	
	virtual void Calculate(std::complex<float>* buffer, double time, double frequency) final override;
		
	void StoreATerms(const std::string& filename, std::complex<float>* buffer);
	
	/**
	 * Set whether a fits image with the a-terms should be written to disk
	 * every time they are calculated.
	 * @param saveATerms Fits images are saved when set to true.
	 */
	void SetSaveATerms(bool saveATerms);

private:
	void calcThread(struct LofarBeamTermThreadData* data);
	
	std::vector<LOFAR::StationResponse::Station::Ptr> _stations;
	size_t _width, _height;
	double _subbandFrequency, _phaseCentreRA, _phaseCentreDec, _dl, _dm, _phaseCentreDL, _phaseCentreDM;
	casacore::MDirection _delayDir, _referenceDir, _tileBeamDir;
	casacore::MPosition _arrayPos;
	bool _useDifferentialBeam, _saveATerms;
};

#else
using LofarBeamTerm = ATermStub;
#endif // HAVE_LOFAR_BEAM

#endif 
