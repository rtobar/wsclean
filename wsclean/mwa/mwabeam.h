#ifndef MWA_BEAM
#define MWA_BEAM

#include <set>

#include "../polarization.h"
#include "../uvector.h"
#include "../wsclean/imagingtable.h"
#include "../wsclean/imagebufferallocator.h"
#include "../wsclean/primarybeamimageset.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/measures/Measures/MDirection.h>

class MWABeam
{
public:
	MWABeam(const ImagingTableEntry* tableEntry, ImageBufferAllocator* allocator) :
	_tableEntry(tableEntry), _allocator(allocator),
	_undersample(8), _secondsBeforeBeamUpdate(1800),
	_frequencyInterpolation(true)
	{
	}
	
	void AddMS(class MSProvider* msProvider, const MSSelection* selection, size_t msIndex)
	{
		_msProviders.push_back(MSProviderInfo(msProvider, selection, msIndex));
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
	
	void SetUndersampling(size_t undersamplingFactor)
	{
		_undersample = undersamplingFactor;
	}
	
private:
	void makeBeamForMS(PrimaryBeamImageSet& beamImages, MSProvider& msProvider, double centralFrequency);

	void makeBeamSnapshot(double** imgPtr, double frequency, casacore::MeasFrame frame, casacore::MPosition arrayPos);
	
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
	
	class ImageBufferAllocator* _allocator;
	
	size_t _width, _height, _sampledWidth, _sampledHeight;
	size_t _undersample, _secondsBeforeBeamUpdate;
	double _pixelSizeX, _pixelSizeY, _phaseCentreRA, _phaseCentreDec, _phaseCentreDL, _phaseCentreDM;
	double _sPixelSizeX, _sPixelSizeY, _totalWeightSum;
	casacore::MDirection _delayDir, _referenceDir, _tileBeamDir;
	
	double _delays[16];
	bool _frequencyInterpolation;
};

#endif
