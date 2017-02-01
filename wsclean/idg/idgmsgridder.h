#ifndef IDG_MS_GRIDDER_H
#define IDG_MS_GRIDDER_H

#ifdef HAVE_IDG

#include "../wsclean/msgridderbase.h"

//#include "interface.h"
#include <idg.h>

#include "../lane.h"
#include "../uvector.h"

#include <boost/thread/mutex.hpp>

class IdgMsGridder : public MSGridderBase
{
public:
	IdgMsGridder();
	
	virtual ~IdgMsGridder();
	
	virtual void Invert();
	
	virtual void Predict(double* real);
	
	virtual void Predict(double* real, double* imaginary);
	
	virtual double* ImageRealResult();
	
	virtual double* ImageImaginaryResult();
	
	virtual double BeamSize() const { return 0.0; }
	
	virtual void GetGriddingCorrectionImage(double* image) const;
	
	virtual bool HasGriddingCorrectionImage() const;
	
private:
	virtual size_t getSuggestedWGridSize() const   {
		return 1; // TODO
	}
		
	void constructGridders(const MultiBandData& selectedBands, size_t nStations, bool constructDegridders);
	
	void gridMeasurementSet(MSGridderBase::MSData& msData);
	void gridThreadFunction();
	
	void predictMeasurementSet(MSGridderBase::MSData& msData);
	void predictCalcThreadFunction();
	void predictWriteThreadFunction(boost::mutex* mutex);
	
	struct IDGInversionRow : public MSGridderBase::InversionRow {
		size_t antenna1, antenna2, timeIndex;
	};
	struct IDGPredictionRow {
		double uvw[3];
		size_t dataDescId, antenna1, antenna2, timeIndex, rowId;
	};
	struct IDGRowForWriting {
		std::complex<float>* data;
		size_t rowId;
	};
	
	std::vector<idg::GridderPlan*> _gridderPlans;
	std::vector<idg::DegridderPlan*> _degridderPlans;
	size_t _subgridSize;
	ao::uvector<std::complex<double>> _grid;
	ao::uvector<double> _image;
	ao::uvector<float> _taper_subgrid;
	ao::uvector<float> _taper_grid;
	ao::lane<IDGInversionRow> _inversionLane;
	ao::lane<IDGPredictionRow> _predictionCalcLane;
	ao::lane<IDGRowForWriting> _predictionWriteLane;
	MSProvider* _outputProvider;
	MultiBandData _selectedBands;
};

#else

#include "../wsclean/unavailablegridder.h"

#define IdgMsGridder UnavailableGridder

#endif // HAVE IDG

#endif
