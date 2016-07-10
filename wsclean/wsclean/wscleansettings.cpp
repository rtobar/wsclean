#include "wscleansettings.h"

#include <sstream>

void WSCleanSettings::Validate()
{
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
}
