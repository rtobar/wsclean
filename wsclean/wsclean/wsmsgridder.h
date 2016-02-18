#ifndef WS_MS_GRIDDER_H
#define WS_MS_GRIDDER_H

#include "inversionalgorithm.h"
#include "wstackinggridder.h"

#include "../lane.h"
#include "../multibanddata.h"

#include <complex>
#include <memory>

#include <casacore/casa/Arrays/Array.h>
#include <casacore/tables/Tables/ArrayColumn.h>

#include <boost/thread/thread.hpp>

namespace casacore {
	class MeasurementSet;
}
class ImageBufferAllocator;

class WSMSGridder : public InversionAlgorithm
{
	public:
		WSMSGridder(class ImageBufferAllocator* imageAllocator, size_t threadCount, double memFraction, double absMemLimit);
	
		virtual void Invert();
		
		virtual void Predict(double* image) { Predict(image, 0); }
		virtual void Predict(double* real, double* imaginary);
		
		virtual double *ImageRealResult() const { return _gridder->RealImage(); }
		virtual double *ImageImaginaryResult() const {
			if(!IsComplex())
				throw std::runtime_error("No imaginary result available for non-complex inversion");
			return _gridder->ImaginaryImage();
		}
		virtual double PhaseCentreRA() const { return _phaseCentreRA; }
		virtual double PhaseCentreDec() const { return _phaseCentreDec; }
		virtual double HighestFrequencyChannel() const { return _freqHigh; }
		virtual double LowestFrequencyChannel() const { return _freqLow; }
		virtual double BandStart() const { return _bandStart; }
		virtual double BandEnd() const { return _bandEnd; }
		virtual double BeamSize() const { return _beamSize; }
		virtual double StartTime() const { return _startTime; }
		virtual bool HasDenormalPhaseCentre() const { return _denormalPhaseCentre; }
		virtual double PhaseCentreDL() const { return _phaseCentreDL; }
		virtual double PhaseCentreDM() const { return _phaseCentreDM; }
		virtual double ImageWeight() const { return _totalWeight; }
		
		enum WStackingGridder::GridModeEnum GridMode() const { return _gridMode; }
		void SetGridMode(WStackingGridder::GridModeEnum gridMode) { _gridMode = gridMode; }
		
		virtual bool HasGriddingCorrectionImage() const { return _gridMode != WStackingGridder::NearestNeighbourGridding; }
		virtual void GetGriddingCorrectionImage(double *image) const { _gridder->GetGriddingCorrectionImage(image); }
		
		size_t ActualInversionWidth() const { return _actualInversionWidth; }
		size_t ActualInversionHeight() const { return _actualInversionHeight; }
		
		void SetTrimSize(size_t trimWidth, size_t trimHeight) {
			_trimWidth = trimWidth;
			_trimHeight = trimHeight;
		}
		void SetNWSize(size_t nwWidth, size_t nwHeight) {
			_nwWidth = nwWidth;
			_nwHeight = nwHeight;
		}
		
		virtual void FreeImagingData()
		{
			_gridder.reset();
		}
	private:
		struct InversionWorkItem
		{
			double u, v, w;
			size_t dataDescId;
			std::complex<float> *data;
		};
		struct InversionWorkSample
		{
			double uInLambda, vInLambda, wInLambda;
			std::complex<float> sample;
		};
		struct PredictionWorkItem
		{
			double u, v, w;
			std::complex<float> *data;
			size_t rowId, dataDescId;
		};
		
		struct MSData
		{
			public:
				MSData();
				~MSData();
				class MSProvider *msProvider;
				MultiBandData bandData;
				size_t startChannel, endChannel;
				size_t matchingRows, totalRowsProcessed;
				double minW, maxW, maxBaseline;
				size_t rowStart, rowEnd;
			
				MultiBandData SelectedBand() const { return MultiBandData(bandData, startChannel, endChannel); }
			private:
				MSData(const MSData &source);
				
				void operator=(const MSData &source);
		};
		
		void initializeMeasurementSet(size_t msIndex, MSData &msData);
		void calculateMetaData(const MSData* msDataVector);
		void gridMeasurementSet(MSData &msData);
		void countSamplesPerLayer(MSData &msData);

		void predictMeasurementSet(MSData &msData);

		void workThread(ao::lane<InversionWorkItem>* workLane)
		{
			InversionWorkItem workItem;
			while(workLane->read(workItem))
			{
				_gridder->AddData(workItem.data, workItem.dataDescId, workItem.u, workItem.v, workItem.w);
				delete[] workItem.data;
			}
		}
		
		void startInversionWorkThreads(size_t maxChannelCount);
		void finishInversionWorkThreads();
		void workThreadPerSample(ao::lane<InversionWorkSample>* workLane);
		
		void predictCalcThread(ao::lane<PredictionWorkItem>* inputLane, ao::lane<PredictionWorkItem>* outputLane);
		void predictWriteThread(ao::lane<PredictionWorkItem>* samplingWorkLane, const MSData* msData);
		static void rotateVisibilities(const BandData &bandData, double shiftFactor, std::complex<float>* dataIter);

		std::unique_ptr<WStackingGridder> _gridder;
		std::unique_ptr<ao::lane<InversionWorkItem>> _inversionWorkLane;
		std::unique_ptr<ao::lane<InversionWorkSample>[]> _inversionCPULanes;
		std::unique_ptr<boost::thread_group> _threadGroup;
		double _phaseCentreRA, _phaseCentreDec, _phaseCentreDL, _phaseCentreDM;
		bool _denormalPhaseCentre, _hasFrequencies;
		double _maxW, _minW;
		double _freqHigh, _freqLow;
		double _bandStart, _bandEnd;
		double _beamSize;
		double _totalWeight;
		double _startTime;
		WStackingGridder::GridModeEnum _gridMode;
		size_t _cpuCount, _laneBufferSize;
		int64_t _memSize;
		ImageBufferAllocator* _imageBufferAllocator;
		size_t _trimWidth, _trimHeight;
		size_t _nwWidth, _nwHeight;
		size_t _actualInversionWidth, _actualInversionHeight;
		double _actualPixelSizeX, _actualPixelSizeY;
};

#endif
