#ifndef IDG_MS_GRIDDER_H
#define IDG_MS_GRIDDER_H

#ifdef HAVE_IDG

#include "../wsclean/msgridderbase.h"

#include "interface.h"

#include "../lane.h"
#include "../uvector.h"

#include <boost/thread/mutex.hpp>

class IdgMsGridder : public MSGridderBase
{
public:
	IdgMsGridder();
	
	virtual ~IdgMsGridder();
	
	virtual void Invert() final override;
	
	virtual void Predict(double* real) final override;
	
	virtual void Predict(double* real, double* imaginary) final override;
	
	virtual double* ImageRealResult() final override;
	
	virtual double* ImageImaginaryResult() final override;
	
	virtual double BeamSize() const final override { return 0.0; }
	
	virtual void GetGriddingCorrectionImage(double* image) const final override;
	
	virtual bool HasGriddingCorrectionImage() const final override;
	
private:
	virtual size_t getSuggestedWGridSize() const override final {
		return 1; // TODO
	}
		
	void constructGridders(const MultiBandData& selectedBands, size_t nStations);
	
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
	
	std::vector<HighLevelGridderInterface*> _interfaces;
	size_t _kernelSize;
	ao::uvector<std::complex<double>> _grid;
	ao::uvector<double> _image;
	ao::uvector<double> _kernel;
	ao::lane<IDGInversionRow> _inversionLane;
	ao::lane<IDGPredictionRow> _predictionCalcLane;
	ao::lane<IDGRowForWriting> _predictionWriteLane;
	MSProvider* _outputProvider;
	MultiBandData _selectedBands;
};

#else

#include "../wsclean/unavailablegridder.h"

using IdgMsGridder=UnavailableGridder;

#endif // HAVE IDG

#endif
