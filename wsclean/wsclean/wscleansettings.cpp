#include "wscleansettings.h"

#include <sstream>

void WSCleanSettings::Validate() const
{
	if(mode == ImagingMode)
	{
		if(trimmedImageWidth == 0 && trimmedImageHeight == 0)
			throw std::runtime_error("Image size has not been set");
		
		if(trimmedImageWidth == 0 || trimmedImageHeight == 0)
			throw std::runtime_error("Invalid image size given: one of the dimensions was zero.");
		
		if(pixelScaleX == 0.0 && pixelScaleY == 0.0)
			throw std::runtime_error("Pixel scale has not been set");
		
		if(pixelScaleX == 0.0 || pixelScaleY == 0.0)
			throw std::runtime_error("Invalid pixel scale given: one direction was set to zero");
	}
	
	// antialiasingKernelSize should be odd
	if(antialiasingKernelSize%2 == 0)
	{
		std::stringstream s;
		s << "Bad anti-aliasing kernel size given of " << antialiasingKernelSize << ". The kernel size has to be odd.";
		throw std::runtime_error(s.str());
	}
	
	if(useIDG)
	{
		bool stokesIOnly = polarizations.size()==1 && *polarizations.begin() == Polarization::StokesI;
		bool allStokes = Polarization::HasFullStokesPolarization(polarizations) &&
			polarizations.size() == 4;
		if(!allStokes && !stokesIOnly)
		{
			throw std::runtime_error("When using IDG, it is only possible to either image Stokes I or to image all 4 Stokes polarizations: use -pol i or -pol iquv");
		}
		if(allStokes && !joinedPolarizationCleaning && deconvolutionIterationCount!=0)
			throw std::runtime_error("Cleaning IDG images with multiple polarizations is only possible in joined polarization mode");
	}
	
	if(baselineDependentAveragingInWavelengths != 0.0)
	{
		if(forceNoReorder)
			throw std::runtime_error("Baseline dependent averaging can not be performed without reordering");
		if(modelUpdateRequired)
			throw std::runtime_error("Baseline dependent averaging can not update the model column (yet) -- you have to add -no-update-model-required.");
	}
	
	if(simulateNoise)
	{
		if(forceNoReorder)
			throw std::runtime_error("Noise simulation can not be performed without reordering");
	}
	
	if(channelsOut == 0)
		throw std::runtime_error("You have specified 0 output channels -- at least one output channel is required.");
	
	if(joinedFrequencyCleaning && channelsOut == 1)
		throw std::runtime_error("Joined frequency cleaning was requested, but only one output channel is being requested. Did you forget -channels-out?");
	
	if(forceReorder && forceNoReorder)
		throw std::runtime_error("Can not both force reordering and force not reordering!");
	
	if(deconvolutionChannelCount != 0 && deconvolutionChannelCount != channelsOut && spectralFittingMode == NoSpectralFitting)
		throw std::runtime_error("You have requested to deconvolve with a decreased number of channels (-deconvolution-channels), but you have not enabled spectral fitting. You should specify an interpolation function by enabling spectral fitting in order to interpolate the deconvolved channels back to the full number of channels. The most useful and common spectral fitting function is -fit-spectral-pol.");
	
	if(savePsfPb && !applyPrimaryBeam)
		throw std::runtime_error("You can not save the primary-beam corrected PSF without enabling primary beam correction: add -apply-primary-beam to your commandline.");
	
	if(saveSourceList && (polarizations.size()!=1 || (*polarizations.begin())!=Polarization::StokesI))
		throw std::runtime_error("Saving a source list currently only works for Stokes I imaging");
	
	if(saveSourceList && deconvolutionIterationCount==0)
		throw std::runtime_error("A source list cannot be saved without cleaning");
	
	checkPolarizations();
}

void WSCleanSettings::checkPolarizations() const
{
	bool hasXY = polarizations.count(Polarization::XY)!=0;
	bool hasYX = polarizations.count(Polarization::YX)!=0;
	if(joinedPolarizationCleaning)
	{
		if(polarizations.size() == 1)
			throw std::runtime_error("Joined/linked polarization cleaning requested, but only one polarization is being imaged. Specify multiple polarizations, or do not request to join the polarizations");
	}
	else {
		if((hasXY || hasYX) && deconvolutionIterationCount !=0)
			throw std::runtime_error("You are imaging XY and/or YX polarizations and have enabled cleaning (niter!=0). This is not possible -- you have to specify '-join-polarizations' or disable cleaning.");
	}
	
	for(PolarizationEnum p : linkedPolarizations)
	{
		if(polarizations.count(p) == 0)
		{
			std::ostringstream str;
			str << "Linked polarization cleaning was requested for polarization "
				<< Polarization::TypeToFullString(p)
				<< ", but this polarization is not imaged";
			throw std::runtime_error(str.str());
		}
	}
	
	if((hasXY && !hasYX) || (!hasXY && hasYX))
		throw std::runtime_error("You are imaging only one of the XY or YX polarizations. This is not possible -- you have to specify both XY and YX polarizations (the output of imaging both polarizations will be the XY and imaginary XY images).");
	if(IsSpectralFittingEnabled())
	{
		if(joinedPolarizationCleaning)
			throw std::runtime_error("You have requested spectral fitting, but you are joining multiple polarizations. This is not supported. You probably want to turn off the joining of polarizations (leave out -join-polarizations).");
		if(!joinedFrequencyCleaning)
			throw std::runtime_error("You have requested spectral fitting, but you are not joining channels. This is not possible: you probably want to turn channel joining on (add -join-channels).");
	}
	
	if(autoDeconvolutionThreshold && autoMask)
	{
		if(autoDeconvolutionThresholdSigma >= autoMaskSigma)
			throw std::runtime_error("The auto-masking threshold was smaller or equal to the auto-threshold. This does not make sense. Did you accidentally reverse the auto-mask and auto-threshold values?");
	}
}

void WSCleanSettings::RecalculatePaddedDimensions()
{
	paddedImageWidth = (size_t) ceil(trimmedImageWidth * imagePadding);
	paddedImageHeight = (size_t) ceil(trimmedImageHeight * imagePadding);
	// Make the width and height divisable by four.
	paddedImageWidth += (4-(paddedImageWidth%4))%4;
	paddedImageHeight += (4-(paddedImageHeight%4))%4;
	if(trimmedImageWidth!=0 && trimmedImageHeight!=0)
	{
		Logger::Debug << "Using image size of " << trimmedImageWidth << " x " << trimmedImageHeight << ", padded to " << paddedImageWidth << " x " << paddedImageHeight << ".\n";
	}
}
