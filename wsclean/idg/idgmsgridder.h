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

private:

    class AverageBeam : public AverageBeamBase
    {
    public:
        AverageBeam() {std::cout << "Constructing AverageBeam" << std::endl;}
        bool is_empty() {return (!m_scalar_beam || !m_matrix_inverse_beam);}
        void set_scalar_beam(std::shared_ptr<std::vector<float>> scalar_beam) {m_scalar_beam = scalar_beam;}
        void set_matrix_inverse_beam(std::shared_ptr<std::vector<std::complex<float>>> matrix_inverse_beam) {m_matrix_inverse_beam = matrix_inverse_beam;}
        std::shared_ptr<std::vector<float>> get_scalar_beam() {return m_scalar_beam;}
        std::shared_ptr<std::vector<std::complex<float>>> get_matrix_inverse_beam() {return m_matrix_inverse_beam;}
    private:
        std::shared_ptr<std::vector<float>> m_scalar_beam;
        std::shared_ptr<std::vector<std::complex<float>>> m_matrix_inverse_beam;
    };

    std::shared_ptr<AverageBeam> _average_beam;

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
