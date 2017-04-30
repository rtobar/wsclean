#include "componentlist.h"

#include "../ndppp.h"

#include "../wsclean/imagefilename.h"
#include "../wsclean/imagingtable.h"
#include "../wsclean/wscleansettings.h"

#include "../multiscale/multiscalealgorithm.h"

void ComponentList::Write(const MultiScaleAlgorithm& multiscale, const WSCleanSettings& settings, long double pixelScaleX, long double pixelScaleY, long double phaseCentreRA, long double phaseCentreDec)
{
	ao::uvector<double> scaleSizes(_nScales);
	for(size_t scaleIndex=0; scaleIndex!=_nScales; ++scaleIndex)
		scaleSizes[scaleIndex] = multiscale.ScaleSize(scaleIndex);
	write(multiscale, scaleSizes, settings, pixelScaleX, pixelScaleY, phaseCentreRA, phaseCentreDec);
}

void ComponentList::WriteSingleScale(const class DeconvolutionAlgorithm& algorithm, const class WSCleanSettings& settings, long double pixelScaleX, long double pixelScaleY, long double phaseCentreRA, long double phaseCentreDec)
{
	ao::uvector<double> scaleSizes(1, 0);
	write(algorithm, scaleSizes, settings, pixelScaleX, pixelScaleY, phaseCentreRA, phaseCentreDec);
}

void ComponentList::write(const class DeconvolutionAlgorithm& algorithm, const ao::uvector<double>& scaleSizes, const class WSCleanSettings& settings, long double pixelScaleX, long double pixelScaleY, long double phaseCentreRA, long double phaseCentreDec)
{
	if(_componentsAddedSinceLastMerge != 0)
		MergeDuplicates();
	
	const SpectralFitter& fitter = algorithm.Fitter();
	if(fitter.Mode() == NoSpectralFitting && _nFrequencies>1)
		throw std::runtime_error("Can't write component list, because you have not specified a spectral fitting method. You probably want to add '-fit-spectral-pol'.");
	
	std::string filename = settings.prefixName + "-sources.txt";
	std::ofstream file(filename);
  bool useLogSI = false;
	switch(fitter.Mode())
	{
		case NoSpectralFitting:
		case PolynomialSpectralFitting:
			useLogSI = false;
			break;
		case LogPolynomialSpectralFitting:
			useLogSI = true;
			break;
	}
	NDPPP::WriteHeaderForSpectralTerms(file, fitter.ReferenceFrequency());
	ao::uvector<double> terms;
	for(size_t scaleIndex=0; scaleIndex!=_nScales; ++scaleIndex)
	{
		ScaleList& list = _listPerScale[scaleIndex];
		size_t componentIndex = 0;
		const double
			scale = scaleSizes[scaleIndex],
			// Using the FWHM formula for a Gaussian
			fwhm = 2.0L * sqrtl(2.0L * logl(2.0L)) * MultiScaleTransforms::GaussianSigma(scale),
			scaleFWHML = fwhm * pixelScaleX * (180.0*60.0*60.0/ M_PI),
			scaleFWHMM = fwhm * pixelScaleY * (180.0*60.0*60.0/ M_PI);
		size_t valueIndex = 0;
		for(size_t index=0; index!=list.positions.size(); ++index)
		{
			const size_t x = list.positions[index].x;
			const size_t y = list.positions[index].y;
      ao::uvector<double> spectrum(_nFrequencies);
			for(size_t frequency = 0; frequency != _nFrequencies; ++frequency)
			{
				spectrum[frequency] = list.values[valueIndex];
				++valueIndex;
			}
			if(_nFrequencies == 1)
				terms.assign(1, spectrum[0]);
			else
				fitter.Fit(terms, spectrum.data());
      double stokesI = terms[0];
      terms.erase(terms.begin());
			long double l, m;
			ImageCoordinates::XYToLM<long double>(x, y, pixelScaleX, pixelScaleY, _width, _height, l, m);
			long double ra, dec;
			ImageCoordinates::LMToRaDec(l, m, phaseCentreRA, phaseCentreDec, ra, dec);
			std::ostringstream name;
			name << 's' << scaleIndex << 'c' << componentIndex;
			if(scale == 0.0)
					NDPPP::WritePolynomialPointComponent(file, name.str(), ra, dec, stokesI, useLogSI, terms, fitter.ReferenceFrequency());
			else {
					NDPPP::WritePolynomialGaussianComponent(file, name.str(), ra, dec, stokesI, useLogSI, terms, fitter.ReferenceFrequency(), scaleFWHML, scaleFWHMM, 0.0);
			}
			++componentIndex;
		}
	}	
}

void ComponentList::loadFromImageSet(ImageSet& imageSet, size_t scaleIndex)
{
	_listPerScale.clear();
	_componentsAddedSinceLastMerge = 0;
	for(size_t y=0; y!=_height; ++y)
	{
		const size_t rowIndex = y*_width;
		for(size_t x=0; x!=_width; ++x)
		{
			const size_t posIndex = rowIndex + x;
			bool isNonZero = false;
			for(size_t imageIndex=0; imageIndex!=imageSet.size(); ++imageIndex)
			{
				if(imageSet[imageIndex][posIndex] != 0.0)
				{
					isNonZero = true;
					break;
				}
			}
			if(isNonZero)
			{
				_listPerScale[scaleIndex].positions.emplace_back(x, y);
				for(size_t imageIndex=0; imageIndex!=imageSet.size(); ++imageIndex)
					_listPerScale[scaleIndex].values.push_back(imageSet[imageIndex][posIndex]);
			}
		}
	}
}

void ComponentList::CorrectForBeam()
{
  // TODO
}
