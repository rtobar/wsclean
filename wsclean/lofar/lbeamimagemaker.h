#ifndef BEAM_IMAGE_MAKER
#define BEAM_IMAGE_MAKER

#include <set>

#include "../polarization.h"
#include "../uvector.h"
#include "../wsclean/imagingtable.h"
#include "../wsclean/imagebufferallocator.h"
#include "../wsclean/primarybeamimageset.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/measures/Measures/MDirection.h>

#ifdef HAVE_LOFAR_BEAM
#include <StationResponse/Station.h>
#endif

class LBeamImageMaker
{
public:
	LBeamImageMaker(const ImagingTableEntry* tableEntry, ImageBufferAllocator* allocator) :
	_tableEntry(tableEntry), _allocator(allocator),
	_undersample(8), _secondsBeforeBeamUpdate(1800),
	_useDifferentialBeam(false)
	{
	}
	
	void AddMS(class MSProvider* msProvider, const MSSelection* selection, size_t msIndex)
	{
		_msProviders.push_back(MSProviderInfo(msProvider, selection, msIndex));
	}
	
	void SetImageWeight(const class ImageWeightCache* imageWeightCache)
	{
		_imageWeightCache = imageWeightCache;
	}
	
	void SetImageDetails(size_t width, size_t height, double pixelSizeX, double pixelSizeY, double phaseCentreRA, double phaseCentreDec, double phaseCentreDL, double phaseCentreDM)
	{
		_width = width;
		_height = height;
		_pixelSizeX = pixelSizeX;
		_pixelSizeY = pixelSizeY;
		_phaseCentreRA = phaseCentreRA;
		_phaseCentreDec = phaseCentreDec;
		_phaseCentreDL = phaseCentreDL;
		_phaseCentreDM = phaseCentreDM;
	}
	
	void Make(PrimaryBeamImageSet& beamImages);
	
	void SetUseDifferentialBeam(bool useDifferentialBeam) {
		_useDifferentialBeam = useDifferentialBeam;
	}
	
	void SetUndersampling(size_t undersamplingFactor)
	{
		_undersample = undersamplingFactor;
	}
	
private:
#ifdef HAVE_LOFAR_BEAM
	class WeightMatrix
	{
	public:
		explicit WeightMatrix(size_t nAntenna) : _nAntenna(nAntenna), _weights(nAntenna*nAntenna, 0)
		{ }
		double& Value(size_t a1, size_t a2)
		{
			if(a1 < a2)
				return _weights[a1*_nAntenna + a2];
			else
				return _weights[a1 + a2*_nAntenna];
		}
		const double& Value(size_t a1, size_t a2) const
		{
			if(a1 < a2)
				return _weights[a1*_nAntenna + a2];
			else
				return _weights[a1 + a2*_nAntenna];
		}
	private:
		size_t _nAntenna;
		ao::uvector<double> _weights;
	};
	
	void makeBeamForMS(PrimaryBeamImageSet& beamImages, MSProvider& msProvider, const ImagingTableEntry::MSInfo& msInfo, const MSSelection& selection, double centralFrequency);

	void makeBeamSnapshot(const std::vector<LOFAR::StationResponse::Station::Ptr>& stations, const ao::uvector<double>& weights, const WeightMatrix& baselineWeights, double** imgPtr, double time, double frequency, double subbandFrequency, const casacore::MeasFrame& frame);
	
	void calculateStationWeights(const class ImageWeights& imageWeights, double& totalWeight, ao::uvector<double>& weights, WeightMatrix& baselineWeights, MSProvider& msProvider, const MSSelection& selection, double endTime);
	
	void logWeights(casacore::MeasurementSet& ms, const ao::uvector<double>& weights);
#endif
	
	struct MSProviderInfo
	{
		MSProviderInfo(MSProvider* _provider, const MSSelection* _selection, size_t _msIndex) :
			provider(_provider), selection(_selection), msIndex(_msIndex)
		{ }
		MSProvider* provider;
		const MSSelection* selection;
		size_t msIndex;
	};
	
	const ImagingTableEntry* _tableEntry;
	std::vector<MSProviderInfo> _msProviders;
	
	const class ImageWeightCache* _imageWeightCache;
	class ImageBufferAllocator* _allocator;
	
	size_t _width, _height, _sampledWidth, _sampledHeight;
	size_t _undersample, _secondsBeforeBeamUpdate;
	double _pixelSizeX, _pixelSizeY, _phaseCentreRA, _phaseCentreDec, _phaseCentreDL, _phaseCentreDM;
	double _sPixelSizeX, _sPixelSizeY, _totalWeightSum;
 	bool _useDifferentialBeam;
	casacore::MDirection _delayDir, _referenceDir, _tileBeamDir;
};

#endif
