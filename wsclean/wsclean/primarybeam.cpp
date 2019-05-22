#include "primarybeam.h"

#include "../fitswriter.h"
#include "../matrix2x2.h"
#include "../fitsreader.h"

#include "../lofar/lbeamimagemaker.h"

#include "../mwa/mwabeam.h"

#include "../atca/atcabeam.h"

#include <boost/filesystem/operations.hpp>

#include <boost/algorithm/string/case_conv.hpp>

void PrimaryBeam::MakeBeamImages(const ImageFilename& imageName, const ImagingTableEntry& entry, std::shared_ptr<ImageWeights> imageWeights, ImageBufferAllocator& allocator)
{
	bool useExistingBeam = false;
	if(_settings.reusePrimaryBeam)
	{
		ImageFilename firstPolName(imageName);
		firstPolName.SetPolarization(Polarization::XX);
		firstPolName.SetIsImaginary(false);
		std::string f(firstPolName.GetBeamPrefix(_settings) + ".fits");
		if(boost::filesystem::exists(f))
		{
			FitsReader reader(f);
			if(reader.ImageWidth()==_settings.trimmedImageWidth && reader.ImageHeight()==_settings.trimmedImageHeight)
			{
				useExistingBeam = true;
				Logger::Info << "File '" << f << "' exists on disk -- reusing files for primary beam.\n";
			}
			else {
				Logger::Info << "File '" << f << "' exists on disk but has different dimensions. Beam will be recreated.\n";
			}
		}
		else {
			Logger::Info << "Primary beam not yet available (file '" << f << "' does not exist). Beam will be created.\n";
		}
	}
	if(!useExistingBeam)
	{
		Logger::Info << " == Constructing primary beam ==\n";
		
		PrimaryBeamImageSet beamImages(_settings.trimmedImageWidth, _settings.trimmedImageHeight, allocator);
		beamImages.SetToZero();
		
		switch(Telescope::GetType(_msProviders.front().first->MS()))
		{
			case Telescope::LOFAR:
			case Telescope::AARTFAAC:
				makeLOFARImage(beamImages, entry, imageWeights, allocator);
				break;
			case Telescope::MWA:
				makeMWAImage(beamImages, entry, allocator);
				break;
			case Telescope::ATCA:
				makeATCAImage(beamImages, entry);
				break;
			default:
				throw std::runtime_error("Can't make beam for this telescope");
		}
		
		// Save the beam images as fits files
		PolarizationEnum
			linPols[4] = { Polarization::XX, Polarization::XY, Polarization::YX, Polarization::YY };
		FitsWriter writer;
		writer.SetImageDimensions(_settings.trimmedImageWidth, _settings.trimmedImageHeight, _phaseCentreRA, _phaseCentreDec, _settings.pixelScaleX, _settings.pixelScaleY);
		writer.SetPhaseCentreShift(_phaseCentreDL, _phaseCentreDM);
		for(size_t i=0; i!=8; ++i)
		{
			PolarizationEnum p = linPols[i/2];
			ImageFilename polName(imageName);
			polName.SetPolarization(p);
			polName.SetIsImaginary(i%2 != 0);
			writer.SetPolarization(p);
			writer.SetFrequency(entry.CentralFrequency(), entry.bandEndFrequency - entry.bandStartFrequency);
			writer.Write<double>(polName.GetBeamPrefix(_settings) + ".fits", beamImages[i].data());
		}
	}
}

void PrimaryBeam::CorrectImages(FitsWriter& writer, const ImageFilename& imageName, const std::string& filenameKind, ImageBufferAllocator& allocator)
{
	PrimaryBeamImageSet beamImages(_settings.trimmedImageWidth, _settings.trimmedImageHeight, allocator);
	load(beamImages, imageName, _settings);
	if(_settings.polarizations.size() == 1 || filenameKind == "psf")
	{
		PolarizationEnum pol = *_settings.polarizations.begin();
		
		if(pol == Polarization::StokesI)
		{
			ImageFilename stokesIName(imageName);
			stokesIName.SetPolarization(pol);
			std::string prefix;
			if(filenameKind == "psf")
				prefix = stokesIName.GetPSFPrefix(_settings);
			else
				prefix = stokesIName.GetPrefix(_settings);
			FitsReader reader(prefix + "-" + filenameKind + ".fits");
			ImageBufferAllocator::Ptr image;
			allocator.Allocate(reader.ImageWidth() * reader.ImageHeight(), image);
			reader.Read(image.data());
			
			beamImages.ApplyStokesI(image.data());
			writer.Write(prefix + "-" + filenameKind + "-pb.fits", image.data());
		}
		else {
			throw std::runtime_error("Primary beam correction is requested, but this is not supported when imaging a single polarization that is not Stokes I. Either image all four polarizations or turn off beam correction.");
		}
	}
	else if(Polarization::HasFullStokesPolarization(_settings.polarizations))
	{
		ImageBufferAllocator::Ptr images[4];
		std::unique_ptr<FitsReader> reader;
		for(size_t polIndex = 0; polIndex != 4; ++polIndex)
		{
			PolarizationEnum pol = Polarization::IndexToStokes(polIndex);
			ImageFilename name(imageName);
			name.SetPolarization(pol);
			reader.reset(new FitsReader(name.GetPrefix(_settings) + "-" + filenameKind + ".fits"));
			allocator.Allocate(reader->ImageWidth() * reader->ImageHeight(), images[polIndex]);
			reader->Read(images[polIndex].data());
		}
		
		double* imagePtrs[4] = { images[0].data(), images[1].data(), images[2].data(), images[3].data() };
		beamImages.ApplyFullStokes(imagePtrs);
		for(size_t polIndex = 0; polIndex != 4; ++polIndex)
		{
			PolarizationEnum pol = Polarization::IndexToStokes(polIndex);
			ImageFilename name(imageName);
			name.SetPolarization(pol);
			writer.SetPolarization(pol);
			writer.Write(name.GetPrefix(_settings) + "-" + filenameKind + "-pb.fits", images[polIndex].data());
		}
	}
	else {
		throw std::runtime_error("Primary beam correction can only be performed on Stokes I or when imaging all four polarizations.");
	}
}

void PrimaryBeam::load(PrimaryBeamImageSet& beamImages, const ImageFilename& imageName, const WSCleanSettings& settings)
{
	PolarizationEnum
		linPols[4] = { Polarization::XX, Polarization::XY, Polarization::YX, Polarization::YY };
	for(size_t i=0; i!=8; ++i)
	{
		PolarizationEnum p = linPols[i/2];
		ImageFilename polName(imageName);
		polName.SetPolarization(p);
		polName.SetIsImaginary(i%2 != 0);
		FitsReader reader(polName.GetBeamPrefix(settings) + ".fits");
		reader.Read(beamImages[i].data());
	}
}

void PrimaryBeam::makeLOFARImage(PrimaryBeamImageSet& beamImages, const ImagingTableEntry& entry, std::shared_ptr<ImageWeights> imageWeights, ImageBufferAllocator& allocator)
{
	LBeamImageMaker lbeam(&entry, &allocator);
	for(size_t i=0; i!=_msProviders.size(); ++i)
		lbeam.AddMS(_msProviders[i].first, &_msProviders[i].second, i);
	lbeam.SetUseDifferentialBeam(_settings.useDifferentialLofarBeam);
	lbeam.SetImageDetails(_settings.trimmedImageWidth, _settings.trimmedImageHeight, _settings.pixelScaleX, _settings.pixelScaleY, _phaseCentreRA, _phaseCentreDec, _phaseCentreDL, _phaseCentreDM);
	lbeam.SetImageWeight(std::move(imageWeights));
	lbeam.SetUndersampling(_settings.primaryBeamUndersampling);
	lbeam.Make(beamImages);
}

void PrimaryBeam::makeMWAImage(PrimaryBeamImageSet& beamImages, const ImagingTableEntry& entry, ImageBufferAllocator& allocator)
{
	MWABeam mwaBeam(&entry, &allocator);
	for(size_t i=0; i!=_msProviders.size(); ++i)
		mwaBeam.AddMS(_msProviders[i].first, &_msProviders[i].second, i);
	mwaBeam.SetImageDetails(_settings.trimmedImageWidth, _settings.trimmedImageHeight, _settings.pixelScaleX, _settings.pixelScaleY, _phaseCentreRA, _phaseCentreDec, _phaseCentreDL, _phaseCentreDM);
	mwaBeam.SetUndersampling(_settings.primaryBeamUndersampling);
	if(!_settings.mwaPath.empty())
		mwaBeam.SetSearchPath(_settings.mwaPath);
	mwaBeam.Make(beamImages);
}

void PrimaryBeam::makeATCAImage(PrimaryBeamImageSet& beamImages, const ImagingTableEntry& entry)
{
	Logger::Info << "Calculating ATCA primary beam...\n";
	ATCABeam::Band band = ATCABeam::GetBand(entry.CentralFrequency() * 1e-9);
	VoltagePattern vp = ATCABeam::CalculateVoltagePattern(band);
	ATCABeam::Calculate(beamImages, _settings.trimmedImageWidth, _settings.trimmedImageHeight, _settings.pixelScaleX, _settings.pixelScaleY, _phaseCentreRA, _phaseCentreDec, _phaseCentreDL, _phaseCentreDM, entry.CentralFrequency(), vp);
}
