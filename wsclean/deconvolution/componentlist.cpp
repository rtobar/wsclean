#include "componentlist.h"

#include "../ndppp.h"

#include "../wsclean/imagefilename.h"
#include "../wsclean/imagingtable.h"
#include "../wsclean/wscleansettings.h"

#include "../multiscale/multiscalealgorithm.h"

void ComponentList::Write(const MultiScaleAlgorithm& multiscale, const WSCleanSettings& settings, long double pixelScaleX, long double pixelScaleY, long double phaseCentreRA, long double phaseCentreDec)
{
	if(_componentsAddedSinceLastMerge != 0)
		MergeDuplicates();
	
	const SpectralFitter& fitter = multiscale.Fitter();
	if(fitter.Mode() == NoSpectralFitting && _nFrequencies>1)
		throw std::runtime_error("Can't write component list, because you have not specified a spectral fitting method. You probably want to add '-fit-spectral-pol'.");
	
	std::string filename = settings.prefixName + "-components.txt";
	std::ofstream file(filename);
	std::string spectralFittingMode;
	switch(fitter.Mode())
	{
		case NoSpectralFitting:
			spectralFittingMode = "None";
			break;
		case PolynomialSpectralFitting:
			spectralFittingMode = "Polynomial";
			break;
		case LogPolynomialSpectralFitting:
			spectralFittingMode = "LogPolynomial";
			break;
	}
	NDPPP::WriteHeaderForSpectralTerms(file, fitter.ReferenceFrequency(), spectralFittingMode);
	ao::uvector<double> spectrum(_nFrequencies);
	ao::uvector<double> terms;
	for(size_t scaleIndex=0; scaleIndex!=_nScales; ++scaleIndex)
	{
		ScaleList& list = _listPerScale[scaleIndex];
		size_t componentIndex = 0;
		const double
			scale = multiscale.ScaleSize(scaleIndex),
			// Using the FWHM formula for a Gaussian
			fwhm = 2.0L * sqrtl(2.0L * logl(2.0L)) * MultiScaleTransforms::GaussianSigma(scale),
			scaleFWHML = fwhm * pixelScaleX * (180.0*60.0*60.0/ M_PI),
			scaleFWHMM = fwhm * pixelScaleY * (180.0*60.0*60.0/ M_PI);
		size_t valueIndex = 0;
		for(size_t index=0; index!=list.positions.size(); ++index)
		{
			const size_t x = list.positions[index].x;
			const size_t y = list.positions[index].y;
			for(size_t frequency = 0; frequency != _nFrequencies; ++frequency)
			{
				spectrum[frequency] = list.values[valueIndex];
				++valueIndex;
			}
			if(_nFrequencies == 1)
				terms.assign(1, spectrum[0]);
			else
				fitter.Fit(terms, spectrum.data());
			long double l, m;
			ImageCoordinates::XYToLM<long double>(x, y, pixelScaleX, pixelScaleY, _width, _height, l, m);
			long double ra, dec;
			ImageCoordinates::LMToRaDec(l, m, phaseCentreRA, phaseCentreDec, ra, dec);
			std::ostringstream name;
			name << 's' << scaleIndex << 'c' << componentIndex;
			if(scale == 0.0)
					NDPPP::WritePolynomialPointComponent(file, name.str(), ra, dec, terms);
			else {
					NDPPP::WritePolynomialGausssianComponent(file, name.str(), ra, dec, terms, scaleFWHML, scaleFWHMM, 0.0);
				++componentIndex;
			}
		}
	}	
}

