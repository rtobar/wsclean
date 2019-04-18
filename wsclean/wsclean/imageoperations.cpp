#include "imageoperations.h"

#include "logger.h"
#include "wscleansettings.h"
#include "wscfitswriter.h"

#include "../fitsreader.h"
#include "../gaussianfitter.h"
#include "../modelrenderer.h"

#include "../units/angle.h"

void ImageOperations::FitBeamSize(const WSCleanSettings& settings, double& bMaj, double& bMin, double& bPA, const double* image, double beamEstimate)
{
	GaussianFitter beamFitter;
	Logger::Info << "Fitting beam... ";
	Logger::Info.Flush();
	if(settings.circularBeam)
	{
		bMaj = beamEstimate;
		beamFitter.Fit2DCircularGaussianCentred(
			image,
			settings.trimmedImageWidth, settings.trimmedImageHeight,
			bMaj,
			settings.beamFittingBoxSize);
		bMin = bMaj;
		bPA = 0.0;
	}
	else {
		beamFitter.Fit2DGaussianCentred(
			image,
			settings.trimmedImageWidth, settings.trimmedImageHeight,
			beamEstimate,
			bMaj, bMin, bPA,
			settings.beamFittingBoxSize);
	}
	bMaj = bMaj*0.5*(settings.pixelScaleX + settings.pixelScaleY);
	bMin = bMin*0.5*(settings.pixelScaleX + settings.pixelScaleY);
}

void ImageOperations::DetermineBeamSize(const WSCleanSettings& settings, double& bMaj, double& bMin, double& bPA, const double* image, double theoreticBeam)
{
	double theoreticBeamWithTaper = theoreticBeam;
	if(settings.gaussianTaperBeamSize != 0.0)
	{
		if(settings.gaussianTaperBeamSize > theoreticBeamWithTaper)
		{
			theoreticBeamWithTaper = settings.gaussianTaperBeamSize;
			Logger::Debug << "Beam is tapered; using " << Angle::ToNiceString(theoreticBeamWithTaper) << " as initial value in PSF fitting.\n";
		}
	}
	if(settings.manualBeamMajorSize != 0.0)
	{
		bMaj = settings.manualBeamMajorSize;
		bMin = settings.manualBeamMinorSize;
		bPA = settings.manualBeamPA;
	} else if(settings.fittedBeam)
	{
		FitBeamSize(settings, bMaj, bMin, bPA, image, theoreticBeamWithTaper*2.0/(settings.pixelScaleX + settings.pixelScaleY));
		Logger::Info << "major=" << Angle::ToNiceString(bMaj) << ", minor=" <<
		Angle::ToNiceString(bMin) << ", PA=" << Angle::ToNiceString(bPA) << ", theoretical=" <<
		Angle::ToNiceString(theoreticBeamWithTaper)<< ".\n";
	}
	else if(settings.theoreticBeam) {
		bMaj = theoreticBeamWithTaper;
		bMin = theoreticBeamWithTaper;
		bPA = 0.0;
		Logger::Info << "Beam size is " << Angle::ToNiceString(theoreticBeamWithTaper) << '\n';
	} else {
		bMaj = std::numeric_limits<double>::quiet_NaN();
		bMin = std::numeric_limits<double>::quiet_NaN();
		bPA = std::numeric_limits<double>::quiet_NaN();
	}
}

void ImageOperations::MakeMFSImage(const WSCleanSettings& settings, const std::vector<OutputChannelInfo>& infoPerChannel, OutputChannelInfo& mfsInfo, const string& suffix, size_t intervalIndex, PolarizationEnum pol, bool isImaginary, bool isPSF)
{
	double lowestFreq = 0.0, highestFreq = 0.0;
	const size_t size = settings.trimmedImageWidth * settings.trimmedImageHeight;
	ao::uvector<double> mfsImage(size, 0.0), addedImage(size), weightImage(size, 0.0);
	double weightSum = 0.0;
	FitsWriter writer;
	for(size_t ch=0; ch!=settings.channelsOut; ++ch)
	{
		std::string prefixStr = isPSF ?
			ImageFilename::GetPSFPrefix(settings, ch, intervalIndex) :
			ImageFilename::GetPrefix(settings, pol, ch, intervalIndex, isImaginary);
		const std::string name(prefixStr + '-' + suffix);
		FitsReader reader(name);
		if(ch == 0)
		{
			WSCFitsWriter wscWriter(reader);
			writer = wscWriter.Writer();
			lowestFreq = reader.Frequency() - reader.Bandwidth()*0.5;
			highestFreq = reader.Frequency() + reader.Bandwidth()*0.5;
		}
		else {
			lowestFreq = std::min(lowestFreq, reader.Frequency() - reader.Bandwidth()*0.5);
			highestFreq = std::max(highestFreq, reader.Frequency() + reader.Bandwidth()*0.5);
		}
		double weight;
		if(!reader.ReadDoubleKeyIfExists("WSCIMGWG", weight))
		{
			Logger::Error << "Error: image " << name << " did not have the WSCIMGWG keyword.\n";
			weight = 0.0;
		}
		weightSum += weight;
		reader.Read(addedImage.data());
		for(size_t i=0; i!=size; ++i)
		{
			if(std::isfinite(addedImage[i])) {
				mfsImage[i] += addedImage[i] * weight;
				weightImage[i] += weight;
			}
		}
	}
	for(size_t i=0; i!=size; ++i)
		mfsImage[i] /= weightImage[i];
	
	if(isPSF)
	{
		double smallestTheoreticBeamSize = 0.0;
		for(std::vector<OutputChannelInfo>::const_reverse_iterator r = infoPerChannel.rbegin();
				r != infoPerChannel.rend(); ++r)
		{
			if(std::isfinite(r->theoreticBeamSize)) {
				smallestTheoreticBeamSize = r->theoreticBeamSize;
				break;
			}
		}
		
		ImageOperations::DetermineBeamSize(settings, mfsInfo.beamMaj, mfsInfo.beamMin, mfsInfo.beamPA, mfsImage.data(), smallestTheoreticBeamSize);
	}
	if(std::isfinite(mfsInfo.beamMaj))
		writer.SetBeamInfo(mfsInfo.beamMaj, mfsInfo.beamMin, mfsInfo.beamPA);
	else
		writer.SetNoBeamInfo();
	
	std::string mfsName(ImageFilename::GetMFSPrefix(settings, pol, intervalIndex, isImaginary, isPSF) + '-' + suffix);
	Logger::Info << "Writing " << mfsName << "...\n";
	writer.SetFrequency((lowestFreq+highestFreq)*0.5, highestFreq-lowestFreq);
	writer.SetExtraKeyword("WSCIMGWG", weightSum);
	writer.RemoveExtraKeyword("WSCCHANS");
	writer.RemoveExtraKeyword("WSCCHANE");
	writer.Write(mfsName, mfsImage.data());
}

void ImageOperations::RenderMFSImage(const WSCleanSettings& settings, const OutputChannelInfo& mfsInfo, size_t intervalIndex, PolarizationEnum pol, bool isImaginary, bool isPBCorrected)
{
	const size_t size = settings.trimmedImageWidth * settings.trimmedImageHeight;
	
	std::string mfsPrefix(ImageFilename::GetMFSPrefix(settings, pol, intervalIndex, isImaginary, false));
	std::string postfix = isPBCorrected ? "-pb.fits" : ".fits";
	FitsReader residualReader(mfsPrefix + "-residual" + postfix);
	FitsReader modelReader(mfsPrefix + "-model" + postfix);
	ao::uvector<double> image(size), modelImage(size);
	residualReader.Read(image.data());
	modelReader.Read(modelImage.data());
	
	double beamMaj = mfsInfo.beamMaj;
	double beamMin, beamPA;
	std::string beamStr;
	if(std::isfinite(beamMaj))
	{
		beamMin = mfsInfo.beamMin;
		beamPA = mfsInfo.beamPA;
		beamStr = "(beam=" + Angle::ToNiceString(beamMin) + "-" +
		Angle::ToNiceString(beamMaj) + ", PA=" +
		Angle::ToNiceString(beamPA) + ")";
	}
	else {
		beamStr = "(beam is neither fitted nor estimated -- using delta scales!)";
		beamMaj = 0.0; beamMin = 0.0; beamPA = 0.0;
	}
	Logger::Info << "Rendering sources to restored image " + beamStr + "... ";
	Logger::Info.Flush();
	bool hasWarned = false;
	for(double& v : modelImage)
	{
		if(!std::isfinite(v))
		{
			if(!hasWarned)
			{
				Logger::Warn <<
					"\nWarning: Some beam corrected model components are NaN. These won't be restored on the pb mf image.\n"
					"This can be caused by an undefined or zero beam inside the FOV. Be sure to check the individual model\n"
					"images to check if the mf image is as expected.\n";
				hasWarned = true;
			}
			v = 0.0;
		}
	}
	ModelRenderer::Restore(image.data(), modelImage.data(), settings.trimmedImageWidth, settings.trimmedImageHeight, beamMaj, beamMin, beamPA, settings.pixelScaleX, settings.pixelScaleY);
	Logger::Info << "DONE\n";
	
	Logger::Info << "Writing " << mfsPrefix << "-image" << postfix << "...\n";
	FitsWriter imageWriter(residualReader);
	imageWriter.Write(mfsPrefix + "-image" + postfix, image.data());
}

