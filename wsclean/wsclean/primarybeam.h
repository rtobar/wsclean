#ifndef PRIMARY_BEAM_H
#define PRIMARY_BEAM_H

#include <string>

#include "imagefilename.h"
#include "imagingtable.h"
#include "primarybeamimageset.h"
#include "wscleansettings.h"

#include "../polarization.h"

#include "../aterms/telescope.h"

class PrimaryBeam
{
public:
	PrimaryBeam(const WSCleanSettings& settings) :
		_settings(settings),
		_phaseCentreRA(0.0), _phaseCentreDec(0.0), _phaseCentreDL(0.0), _phaseCentreDM(0.0)
	{ }
	
	void SetPhaseCentre(double ra, double dec, double dl, double dm)
	{
		_phaseCentreRA = ra;
		_phaseCentreDec = dec;
		_phaseCentreDL = dl;
		_phaseCentreDM = dm;
	}
	
	void CorrectImages(const ImageFilename& imageName, std::vector<double*>& images, ImageBufferAllocator& allocator)
	{
		PrimaryBeamImageSet beamImages(_settings.trimmedImageWidth, _settings.trimmedImageHeight, allocator);
		load(beamImages, imageName, _settings);
		if(_settings.polarizations.size() == 1 && *_settings.polarizations.begin() == Polarization::StokesI)
		{
			beamImages.ApplyStokesI(images[0]);
		}
		else if(_settings.polarizations.size() == 4 && Polarization::HasFullStokesPolarization(_settings.polarizations))
		{
			beamImages.ApplyFullStokes(images.data());
		}
	}
	
	void Load(PrimaryBeamImageSet& beamImages, const ImageFilename& imageName)
	{
		load(beamImages, imageName, _settings);
	}
	
	void AddMS(class MSProvider* msProvider, const MSSelection& selection)
	{
		_msProviders.push_back(std::make_pair(msProvider, selection));
	}
	
	void MakeBeamImages(const ImageFilename& imageName, const ImagingTableEntry& entry, std::shared_ptr<class ImageWeights> imageWeights, ImageBufferAllocator& allocator);
	
	void CorrectImages(class FitsWriter& writer, const ImageFilename& imageName, const std::string& filenameKind, ImageBufferAllocator& allocator);
	
private:
	const WSCleanSettings& _settings;
	std::vector<std::pair<MSProvider*, MSSelection>> _msProviders;
	double _phaseCentreRA, _phaseCentreDec, _phaseCentreDL, _phaseCentreDM;
	
	static void load(PrimaryBeamImageSet& beamImages, const ImageFilename& imageName, const WSCleanSettings& settings);
	
	void makeLOFARImage(PrimaryBeamImageSet& beamImages, const ImagingTableEntry& entry, std::shared_ptr<class ImageWeights> imageWeights, ImageBufferAllocator& allocator);
	
	void makeMWAImage(PrimaryBeamImageSet& beamImages, const ImagingTableEntry& entry, ImageBufferAllocator& allocator);
	
	void makeATCAImage(PrimaryBeamImageSet& beamImages, const ImagingTableEntry& entry);
};

#endif
