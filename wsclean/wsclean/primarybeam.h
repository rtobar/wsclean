#ifndef PRIMARY_BEAM_H
#define PRIMARY_BEAM_H

#include <string>
#include <boost/filesystem/operations.hpp>

#include "imagefilename.h"
#include "imagingtable.h"
#include "primarybeamimageset.h"
#include "wscleansettings.h"

#include "../polarization.h"

#include "../lofar/lbeamimagemaker.h"
#include "../fitswriter.h"
#include "../matrix2x2.h"
#include "../fitsreader.h"

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
	
	void MakeBeamImages(const ImageFilename& imageName, const ImagingTableEntry& entry, const ImageWeightCache* imageWeightCache, ImageBufferAllocator& allocator)
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
			
			// TODO Find out what telescope array is used in the (first) measurement set
			
			PrimaryBeamImageSet beamImages(_settings.trimmedImageWidth, _settings.trimmedImageHeight, allocator);
			beamImages.SetToZero();
			
			makeLOFARImage(beamImages, entry, imageWeightCache, allocator);
		
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
	
	void CorrectImages(const ImageFilename& imageName, std::vector<double*>& images, ImageBufferAllocator& allocator)
	{
		PrimaryBeamImageSet beamImages(_settings.trimmedImageWidth, _settings.trimmedImageHeight, allocator);
		load(beamImages, imageName);
		if(_settings.polarizations.size() == 1 && *_settings.polarizations.begin() == Polarization::StokesI)
		{
			applyStokesI(images[0], beamImages);
		}
		else if(_settings.polarizations.size() == 4 && Polarization::HasFullStokesPolarization(_settings.polarizations))
		{
			applyFullStokes(images.data(), beamImages);
		}
	}
	
	void Load(PrimaryBeamImageSet& beamImages, const ImageFilename& imageName)
	{
		load(beamImages, imageName);
	}
	
	void CorrectImages(FitsWriter& writer, const ImageFilename& imageName, const std::string& filenameKind, ImageBufferAllocator& allocator)
	{
		PrimaryBeamImageSet beamImages(_settings.trimmedImageWidth, _settings.trimmedImageHeight, allocator);
		load(beamImages, imageName);
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
				
				applyStokesI(image.data(), beamImages);
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
			applyFullStokes(imagePtrs, beamImages);
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
	
	void AddMS(class MSProvider* msProvider, const MSSelection& selection)
	{
		_msProviders.push_back(std::make_pair(msProvider, selection));
	}
private:
	const WSCleanSettings& _settings;
	std::vector<std::pair<MSProvider*, MSSelection>> _msProviders;
	double _phaseCentreRA, _phaseCentreDec, _phaseCentreDL, _phaseCentreDM;
	
	void load(PrimaryBeamImageSet& beamImages, const ImageFilename& imageName) {
		PolarizationEnum
			linPols[4] = { Polarization::XX, Polarization::XY, Polarization::YX, Polarization::YY };
		for(size_t i=0; i!=8; ++i)
		{
			PolarizationEnum p = linPols[i/2];
			ImageFilename polName(imageName);
			polName.SetPolarization(p);
			polName.SetIsImaginary(i%2 != 0);
			FitsReader reader(polName.GetBeamPrefix(_settings) + ".fits");
			reader.Read(beamImages[i].data());
		}
	}
	
	void makeLOFARImage(PrimaryBeamImageSet& beamImages, const ImagingTableEntry& entry, const ImageWeightCache* imageWeightCache, ImageBufferAllocator& allocator)
	{
		LBeamImageMaker lbeam(&entry, &allocator);
		for(size_t i=0; i!=_msProviders.size(); ++i)
			lbeam.AddMS(_msProviders[i].first, &_msProviders[i].second, i);
		lbeam.SetUseDifferentialBeam(_settings.useDifferentialLofarBeam);
		lbeam.SetImageDetails(_settings.trimmedImageWidth, _settings.trimmedImageHeight, _settings.pixelScaleX, _settings.pixelScaleY, _phaseCentreRA, _phaseCentreDec, _phaseCentreDL, _phaseCentreDM);
		lbeam.SetImageWeight(imageWeightCache);
		lbeam.Make(beamImages);
	}
	
	void applyStokesI(double* stokesI, PrimaryBeamImageSet& beamImages) const
	{
		// If Iu is uncorrected and Ic is corrected:
		// Iu = B Ic B^*
		// when I is unpolarized (diagonal, constant)
		// Iu = Ic B B^*
		// Ic = Iu (B B^*)^-1
		// Since we have measured Iu_xx + Iu_yy, and want to know Ic_xx + Ic_yy, let B2 = (B B^*)^-1 and Iu_xx = Iu_yy :
		// Ic_xx + Ic_yy = Iu_xx B2_xx + Iu_yy B2_yy = (Iu_xx + Iu_yy) (B2_xx + B2_yy)
		
		size_t size = _settings.trimmedImageWidth * _settings.trimmedImageHeight;
		for(size_t j=0; j!=size; ++j)
		{
			MC2x2 val, squared;
			val[0] = std::complex<double>(beamImages[0][j], beamImages[1][j]);
			val[1] = std::complex<double>(beamImages[2][j], beamImages[3][j]);
			val[2] = std::complex<double>(beamImages[4][j], beamImages[5][j]);
			val[3] = std::complex<double>(beamImages[6][j], beamImages[7][j]);
			MC2x2::ATimesHermB(squared, val, val);
			if(squared.Invert())
				stokesI[j] = stokesI[j] * 0.5 * (squared[0].real() + squared[3].real());
			else
				stokesI[j] = std::numeric_limits<double>::quiet_NaN();
		}
	}
	
	void applyFullStokes(double* images[4], PrimaryBeamImageSet& beamImages) const
	{
		size_t size = _settings.trimmedImageWidth * _settings.trimmedImageHeight;
		for(size_t j=0; j!=size; ++j)
		{
			MC2x2 beamVal;
			beamVal[0] = std::complex<double>(beamImages[0][j], beamImages[1][j]);
			beamVal[1] = std::complex<double>(beamImages[2][j], beamImages[3][j]);
			beamVal[2] = std::complex<double>(beamImages[4][j], beamImages[5][j]);
			beamVal[3] = std::complex<double>(beamImages[6][j], beamImages[7][j]);
			if(beamVal.Invert())
			{
				double stokesVal[4] = { images[0][j], images[1][j], images[2][j], images[3][j] };
				MC2x2 linearVal, scratch;
				Polarization::StokesToLinear(stokesVal, linearVal.Data());
				MC2x2::ATimesB(scratch, beamVal, linearVal);
				MC2x2::ATimesHermB(linearVal, scratch, beamVal);
				Polarization::LinearToStokes(linearVal.Data(), stokesVal);
				for(size_t p=0; p!=4; ++p)
					images[p][j] = stokesVal[p];
			}
			else {
				for(size_t p=0; p!=4; ++p)
					images[p][j] = std::numeric_limits<double>::quiet_NaN();
			}
		}
	}
};

#endif
