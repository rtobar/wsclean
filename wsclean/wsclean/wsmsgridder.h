#ifndef WS_MS_GRIDDER_H
#define WS_MS_GRIDDER_H

#include "msgridderbase.h"
#include "wstackinggridder.h"

#include "../lane.h"
#include "../multibanddata.h"

#include <complex>
#include <memory>

#include <casacore/casa/Arrays/Array.h>
#include <casacore/tables/Tables/ArrayColumn.h>

#include <thread>

namespace casacore {
	class MeasurementSet;
}
class ImageBufferAllocator;

class WSMSGridder : public MSGridderBase
{
	public:
		WSMSGridder(class ImageBufferAllocator* imageAllocator, size_t threadCount, double memFraction, double absMemLimit);
	
		virtual void Invert() final override;
		
		virtual void Predict(ImageBufferAllocator::Ptr image) final override { Predict(std::move(image), nullptr); }
		virtual void Predict(ImageBufferAllocator::Ptr real, ImageBufferAllocator::Ptr imaginary) final override;
		
		virtual ImageBufferAllocator::Ptr ImageRealResult() final override { return _gridder->RealImage(); }
		virtual ImageBufferAllocator::Ptr ImageImaginaryResult() final override {
			if(!IsComplex())
				throw std::runtime_error("No imaginary result available for non-complex inversion");
			return _gridder->ImaginaryImage();
		}
		virtual bool HasGriddingCorrectionImage() const final override { return GridMode() != NearestNeighbourGridding; }
		virtual void GetGriddingCorrectionImage(double *image) const final override { _gridder->GetGriddingCorrectionImage(image); }
		
		virtual size_t ActualInversionWidth() const final override { return _actualInversionWidth; }
		virtual size_t ActualInversionHeight() const final override { return _actualInversionHeight; }
		
		virtual void FreeImagingData() final override
		{
			_gridder.reset();
		}
		
	private:
		struct InversionWorkSample
		{
			double uInLambda, vInLambda, wInLambda;
			std::complex<float> sample;
		};
		struct PredictionWorkItem
		{
			double u, v, w;
			std::unique_ptr<std::complex<float>[]> data;
			size_t rowId, dataDescId;
		};
		
		void gridMeasurementSet(MSData& msData);
		void countSamplesPerLayer(MSData& msData);
		virtual size_t getSuggestedWGridSize() const final override;

		void predictMeasurementSet(MSData& msData);

		void workThread(ao::lane<InversionRow>* workLane)
		{
			InversionRow workItem;
			while(workLane->read(workItem))
			{
				_gridder->AddData(workItem.data, workItem.dataDescId, workItem.uvw[0], workItem.uvw[1], workItem.uvw[2]);
				delete[] workItem.data;
			}
		}
		
		void startInversionWorkThreads(size_t maxChannelCount);
		void finishInversionWorkThreads();
		void workThreadPerSample(ao::lane<InversionWorkSample>* workLane);
		
		void predictCalcThread(ao::lane<PredictionWorkItem>* inputLane, ao::lane<PredictionWorkItem>* outputLane);
		void predictWriteThread(ao::lane<PredictionWorkItem>* samplingWorkLane, const MSData* msData);

		std::unique_ptr<WStackingGridder> _gridder;
		std::vector<ao::lane<InversionWorkSample>> _inversionCPULanes;
		std::vector<std::thread> _threadGroup;
		size_t _cpuCount, _laneBufferSize;
		int64_t _memSize;
		ImageBufferAllocator* _imageBufferAllocator;
};

#endif
