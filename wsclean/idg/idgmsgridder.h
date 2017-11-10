#ifndef IDG_MS_GRIDDER_H
#define IDG_MS_GRIDDER_H

#ifdef HAVE_IDG

#include "../wsclean/msgridderbase.h"

//#include "interface.h"
#include <idg-api.h>

#include "../lane.h"
#include "../uvector.h"

#include <boost/thread/mutex.hpp>

class IdgMsGridder : public MSGridderBase
{
public:
	IdgMsGridder(const class WSCleanSettings& settings);
	
	virtual ~IdgMsGridder() final override;
	
	virtual void Invert();
	
	virtual void Predict(double* real);
	
	virtual void Predict(double* real, double* imaginary);
	
	virtual double* ImageRealResult();
	
	virtual double* ImageImaginaryResult();
	
	virtual void GetGriddingCorrectionImage(double* image) const;
	
	virtual bool HasGriddingCorrectionImage() const;

	/**
	 * Stub functie die te implementeren methode spiegelt.
	 * 
	 * De aterms moeten dimensies n_stations * subgridsize * subgridsize * 2 * 2 hebben (2x2 is Jones matrix). 
	 * 
	 * Werking van set_aterm:
	 * - Voegt aan m_aterm een slice toe met een kopie van aterms
	 * - Voegt aan m_atermoffsets de timeIndex toe, zodat de gridder weet dat op dat timeslot de aterm verandert
	 * De functie set_aterm roep je aan het gridden begint. Dus bijvoorbeeld op het moment dat je de visibility
	 * voert waarom de a-term verandert. Maar eerder mag ook. Iets later eventueel ook, maar dan voordat je een
	 * visibility van een volgend tijdslot pusht (want dan kan hij besluiten te gaan flushen).
	 */
	void set_aterm(size_t timeIndex, const std::complex<float>* aTerms) { }
	
private:
	virtual size_t getSuggestedWGridSize() const   {
		return 1; // TODO
	}
		
	void gridMeasurementSet(MSGridderBase::MSData& msData);
	void gridThreadFunction();
	
	void predictMeasurementSet(MSGridderBase::MSData& msData);
	void readConfiguration();
	
	void setIdgType();
	
	struct IDGInversionRow : public MSGridderBase::InversionRow {
		size_t antenna1, antenna2, timeIndex;
	};
	struct IDGPredictionRow {
		double uvw[3];
		size_t dataDescId, antenna1, antenna2, timeIndex, rowId;
	};
	void predictRow(IDGPredictionRow& row);
	void computePredictionBuffer(size_t dataDescId);
	
	std::unique_ptr<idg::api::BufferSet> _bufferset;
	size_t _subgridSize;
	ao::uvector<double> _image;
	ao::uvector<float> _taper_subgrid;
	ao::uvector<float> _taper_grid;
	MSProvider* _outputProvider;
	MultiBandData _selectedBands;
	const WSCleanSettings& _settings;
	idg::api::Type _proxyType;
	int _buffersize;
	idg::api::options_type _options;
};

void init_optimal_taper_1D(int subgridsize, int gridsize, float kernelsize, float padding, float* taper_subgrid, float* taper_grid);
void init_optimal_gridding_taper_1D(int subgridsize, int gridsize, float kernelsize, float* taper_subgrid, float* taper_grid);

#else

#include "../wsclean/unavailablegridder.h"

#define IdgMsGridder UnavailableGridder

#endif // HAVE IDG

#endif
