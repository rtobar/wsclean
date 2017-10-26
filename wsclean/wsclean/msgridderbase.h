#ifndef MS_GRIDDER_BASE_H
#define MS_GRIDDER_BASE_H

#include "measurementsetgridder.h"

#include "../multibanddata.h"
#include "../uvector.h"

class MSGridderBase : public MeasurementSetGridder
{
public:
	MSGridderBase();
	~MSGridderBase();
	
	virtual double StartTime() const final override { return _startTime; }
	virtual double PhaseCentreRA() const final override { return _phaseCentreRA; }
	virtual double PhaseCentreDec() const final override { return _phaseCentreDec; }
	virtual double PhaseCentreDL() const final override { return _phaseCentreDL; }
	virtual double PhaseCentreDM() const final override { return _phaseCentreDM; }
	virtual bool HasDenormalPhaseCentre() const final override { return _denormalPhaseCentre; }
	virtual double ImageWeight() const final override { return _totalWeight*2.0; }
	virtual double NormalizationFactor() const final override {
		return NormalizeForWeighting() ? _totalWeight*2.0 : 1.0;
	}
	virtual double BeamSize() const final override { return _theoreticalBeamSize; }
	
	/**
	 * This is the sum of the weights as given by the measurement set, before the
	 * image weighting is applied.
	 */
	double VisibilityWeightSum() const { return _visibilityWeightSum; }
	/**
	 * The number of visibilities that were gridded.
	 */
	size_t GriddedVisibilityCount() const { return _griddedVisibilityCount; }
	/**
	 * The maximum weight, after having applied the imaging weighting.
	 */
	double MaxGriddedWeight() const { return _maxGriddedWeight; }
	/**
	 * The effective number of visibilities, taking into account imaging weighting
	 * and visibility weighting. This number is relative to the "best" visibility:
	 * if one visibility with a weight of 10 and 5 visibilities with
	 * a weight of 4 were gridded, the effective number of visibilities is
	 * (10 + 5 x 4) / 10 = 3
	 */
	double EffectiveGriddedVisibilityCount() const { return totalWeight()/MaxGriddedWeight(); }
	
	static void GetPhaseCentreInfo(casacore::MeasurementSet& ms, size_t fieldId, double& ra, double& dec, double& dl, double& dm);
	
	const std::string& TelescopeName() const { return _telescopeName; }
	
	const std::string& Observer() const { return _observer; }
	
	const std::string& FieldName() const { return _fieldName; }
	
	struct MetaDataCache
	{
		struct Entry {
			double minW, maxW, maxBaselineUVW, maxBaselineInM;
		};
		std::vector<Entry> msDataVector;
	};
	
	void SetMetaDataCache(MetaDataCache* cache) { _metaDataCache = cache; }
	
protected:
	struct MSData
	{
		public:
			MSData();
			~MSData();
			class MSProvider *msProvider;
			size_t msIndex;
			MultiBandData bandData;
			size_t startChannel, endChannel;
			size_t matchingRows, totalRowsProcessed;
			double minW, maxW, maxBaselineUVW, maxBaselineInM;
			size_t rowStart, rowEnd;
		
			MultiBandData SelectedBand() const { return MultiBandData(bandData, startChannel, endChannel); }
		private:
			MSData(const MSData &source);
			
			void operator=(const MSData &source);
	};

	struct InversionRow
	{
		double uvw[3];
		size_t dataDescId, rowId;
		std::complex<float>* data;
	};
		
	void resetMetaData()
	{
		_hasFrequencies = false;
	}
	
	void calculateMSLimits(const MultiBandData& selectedBand, double startTime)
	{
		if(_hasFrequencies)
		{
			_freqLow = std::min(_freqLow, selectedBand.LowestFrequency());
			_freqHigh = std::max(_freqHigh, selectedBand.HighestFrequency());
			_bandStart = std::min(_bandStart, selectedBand.BandStart());
			_bandEnd = std::max(_bandEnd, selectedBand.BandEnd());
			_startTime = std::min(_startTime, startTime);
		} else {
			_freqLow = selectedBand.LowestFrequency();
			_freqHigh = selectedBand.HighestFrequency();
			_bandStart = selectedBand.BandStart();
			_bandEnd = selectedBand.BandEnd();
			_startTime = startTime;
			_hasFrequencies = true;
		}
	}
	
	template<size_t NPolInMSProvider>
	void calculateWLimits(MSGridderBase::MSData& msData);
	
	void initializeMeasurementSet(MSGridderBase::MSData& msData, MetaDataCache::Entry& cacheEntry, bool isCacheInitialized);
	
	void calculateOverallMetaData(const MSData* msDataVector);
	
	template<size_t PolarizationCount>
	void readAndWeightVisibilities(MSProvider& msProvider, InversionRow& rowData, const BandData& curBand, float* weightBuffer, std::complex<float>* modelBuffer, const bool* isSelected);

	double _maxW, _minW;
	double _theoreticalBeamSize;
	size_t _actualInversionWidth, _actualInversionHeight;
	double _actualPixelSizeX, _actualPixelSizeY;
	
	virtual size_t getSuggestedWGridSize() const = 0;
	
	void resetVisibilityCounters() {
		_griddedVisibilityCount = 0;
		_totalWeight = 0.0;
		_maxGriddedWeight = 0.0;
	}
	
	double totalWeight() const { return _totalWeight; }
	
	void initializeMSDataVector(std::vector<MSData>& msDataVector, size_t nPolInMSProvider);
	
private:
	template<size_t PolarizationCount>
	static void rotateVisibilities(const BandData &bandData, double shiftFactor, std::complex<float>* dataIter);
	
	void initializePhaseCentre(casacore::MeasurementSet& ms, size_t fieldId);
	
	void initializeBandData(casacore::MeasurementSet& ms, MSGridderBase::MSData& msData);
	
	void initializeMetaData(casacore::MeasurementSet& ms, size_t fieldId);
		
	bool _hasFrequencies;
	double _freqHigh, _freqLow;
	double _bandStart, _bandEnd;
	double _startTime;
	struct MetaDataCache* _metaDataCache;
	
	double _phaseCentreRA, _phaseCentreDec, _phaseCentreDL, _phaseCentreDM;
	bool _denormalPhaseCentre;
	std::string _telescopeName, _observer, _fieldName;
	
	size_t _griddedVisibilityCount;
	double _totalWeight;
	double _maxGriddedWeight;
	double _visibilityWeightSum;
	
	ao::uvector<float> _scratchWeights;
};

#endif
