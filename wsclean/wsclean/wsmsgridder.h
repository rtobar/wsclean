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
	
		virtual void Invert();
		
		virtual void Predict(double* image) { Predict(image, 0); }
		virtual void Predict(double* real, double* imaginary);
		
		virtual double *ImageRealResult() { return _gridder->RealImage(); }
		virtual double *ImageImaginaryResult() {
			if(!IsComplex())
				throw std::runtime_error("No imaginary result available for non-complex inversion");
			return _gridder->ImaginaryImage();
		}
		virtual bool HasGriddingCorrectionImage() const { return GridMode() != NearestNeighbourGridding; }
		virtual void GetGriddingCorrectionImage(double *image) const { _gridder->GetGriddingCorrectionImage(image); }
		
		virtual size_t ActualInversionWidth() const { return _actualInversionWidth; }
		virtual size_t ActualInversionHeight() const { return _actualInversionHeight; }
		
		virtual void FreeImagingData()
		{
			_gridder.reset();
		}
		
		void SetGridAtLOFARCentroid(bool gridAtLOFARCentroid, size_t fullResWidth)
		{
			_gridAtLOFARCentroid = gridAtLOFARCentroid;
			_fullResWidth = fullResWidth;
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
		void gridLOFARCentroidMeasurementSet(MSData& msData);
		void countSamplesPerLayer(MSData& msData);
		virtual size_t getSuggestedWGridSize() const  ;

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
		bool _gridAtLOFARCentroid;
		size_t _fullResWidth;
		ImageBufferAllocator* _imageBufferAllocator;
};

#endif
