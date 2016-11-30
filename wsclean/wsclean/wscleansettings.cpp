#include "wscleansettings.h"

#include <sstream>

void WSCleanSettings::Validate() const
{
	if(mode == ImagingMode)
	{
		if(untrimmedImageWidth == 0 && untrimmedImageHeight == 0)
			throw std::runtime_error("Image size has not been set");
		
		if(untrimmedImageWidth == 0 || untrimmedImageHeight == 0)
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
		throw std::runtime_error("Joined frequency cleaning was requested, but only one output channel is being requested. Did you forget -channelsout?");
	
	if(forceReorder && forceNoReorder)
		throw std::runtime_error("Can not both force reordering and force not reordering!");
	
	checkPolarizations();
}

void WSCleanSettings::checkPolarizations() const
{
	bool hasXY = polarizations.count(Polarization::XY)!=0;
	bool hasYX = polarizations.count(Polarization::YX)!=0;
	if(joinedPolarizationCleaning)
	{
		if(polarizations.size() == 1)
			throw std::runtime_error("Joined polarization cleaning requested, but only one polarization is being imaged. Specify multiple polarizatons, or do not request to join the polarizations");
		else if(polarizations.size() != 2 && polarizations.size() != 4)
			throw std::runtime_error("Joined polarization cleaning requested, but neither 2 or 4 polarizations are imaged that are suitable for this");
	}
	else {
		if((hasXY || hasYX) && deconvolutionIterationCount !=0)
			throw std::runtime_error("You are imaging XY and/or YX polarizations and have enabled cleaning (niter!=0). This is not possible -- you have to specify '-joinpolarizations' or disable cleaning.");
	}
	if((hasXY && !hasYX) || (!hasXY && hasYX))
		throw std::runtime_error("You are imaging only one of XY or YX polarizations. This is not possible -- you have to specify both XY and YX polarizations (the output of imaging both polarizations will be the XY and imaginary XY images).");
	if(IsSpectralFittingEnabled())
	{
		if(joinedPolarizationCleaning)
			throw std::runtime_error("You have requested spectral fitting, but you are joining multiple polarizations. This is not supported. You probably want to turn off the joining of polarizations (leave out -joinpolarizations).");
		if(!joinedFrequencyCleaning)
			throw std::runtime_error("You have requested spectral fitting, but you are not joining channels. This is not possible: you probably want to turn channel joining on (add -joinchannels).");
	}
}

void WSCleanSettings::setDimensions()
{
	if(trimmedImageWidth==0 && trimmedImageHeight==0)
	{
		trimmedImageWidth = untrimmedImageWidth;
		trimmedImageHeight = untrimmedImageHeight;
	}
	else if(trimmedImageWidth > untrimmedImageWidth || trimmedImageHeight > untrimmedImageHeight)
	{
		throw std::runtime_error("Error in specified trim dimensions: at least one dimension of the trimmed image is larger than in the untrimmed image");
	}
}

