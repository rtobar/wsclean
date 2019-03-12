#ifndef IDG_MS_GRIDDER_H
#define IDG_MS_GRIDDER_H

#ifdef HAVE_IDG

#include "../wsclean/imagebufferallocator.h"
#include "../wsclean/msgridderbase.h"

//#include "interface.h"
#include <idg-api.h>

#include "../lane.h"
#include "../uvector.h"

#include <boost/thread/mutex.hpp>

class IdgMsGridder : public MSGridderBase
{
public:
	IdgMsGridder(const class WSCleanSettings& settings, ImageBufferAllocator& allocator);
	
	virtual ~IdgMsGridder() final override { }
	
	virtual void Invert() final override;
	
	virtual void Predict(ImageBufferAllocator::Ptr real) final override;
	
	virtual void Predict(ImageBufferAllocator::Ptr real, ImageBufferAllocator::Ptr imaginary) final override;
	
	virtual ImageBufferAllocator::Ptr ImageRealResult() final override;
	
	virtual ImageBufferAllocator::Ptr ImageImaginaryResult() final override;
	
	virtual void GetGriddingCorrectionImage(double* image) const final override;
	
	virtual bool HasGriddingCorrectionImage() const final override;
	
	void SavePBCorrectedImages(class FitsWriter& writer, const class ImageFilename& filename, const std::string& filenameKind, class ImageBufferAllocator& allocator) const;
	
	void SaveBeamImage(const class ImagingTableEntry& entry, class ImageFilename& filename) const;

private:
	class AverageBeam : public AverageBeamBase
	{
	public:
		AverageBeam() { }
		bool Empty() { return (!_scalarBeam || !_matrixInverseBeam); }
		void SetScalarBeam(const std::shared_ptr<std::vector<float>>& scalarBeam) { _scalarBeam = scalarBeam; }
		void SetMatrixInverseBeam(const std::shared_ptr<std::vector<std::complex<float>>>& matrixInverseBeam) { _matrixInverseBeam = matrixInverseBeam; }
		
		/**
		 * The image resulting from IDG gridding is multiplied by the scalar beam. It is the result of this multiplication that is
		 * returned to WSClean. The scalar beam is chosen such that the Stokes I image is flat noise. There is no guarantee that
		 * Stokes Q, U, V will be flat noise, nor that there is no correlation between the noise in I,Q,U,V, but for all practical
		 * purposes they can be treated as such.
		 */
		std::shared_ptr<std::vector<float>>& ScalarBeam() { return _scalarBeam; }
		
		/**
		 * The matrix inverse beam is applied while gridding. It is the inverse of the mean square matrix beam.
		 */
		std::shared_ptr<std::vector<std::complex<float>>>& MatrixInverseBeam() { return _matrixInverseBeam; }
		
	private:
		std::shared_ptr<std::vector<float>> _scalarBeam;
		std::shared_ptr<std::vector<std::complex<float>>> _matrixInverseBeam;
	};

	AverageBeam* _averageBeam;

	virtual size_t getSuggestedWGridSize() const override final {
		return 1; // TODO
	}
		
	void gridMeasurementSet(MSGridderBase::MSData& msData);
	void gridThreadFunction();
	
	void predictMeasurementSet(MSGridderBase::MSData& msData);
	void readConfiguration();
	
	void setIdgType();
	
	std::unique_ptr<class ATermBase> getATermMaker(MSGridderBase::MSData& msData);
	
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
	ImageBufferAllocator& _allocator;
};

void init_optimal_taper_1D(int subgridsize, int gridsize, float kernelsize, float padding, float* taper_subgrid, float* taper_grid);
void init_optimal_gridding_taper_1D(int subgridsize, int gridsize, float kernelsize, float* taper_subgrid, float* taper_grid);

#else

#include "../wsclean/unavailablegridder.h"

#define IdgMsGridder UnavailableGridder

#endif // HAVE IDG

#endif
