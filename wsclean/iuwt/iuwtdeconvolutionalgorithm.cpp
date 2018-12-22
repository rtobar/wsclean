#include "iuwtdeconvolutionalgorithm.h"

#include "imageanalysis.h"

#include "../fftconvolver.h"
#include "../fftwmanager.h"
#include "../gaussianfitter.h"
#include "../modelrenderer.h"
#include "../threadpool.h"

#include "../deconvolution/imageset.h"
#include "../image.h"

#include <algorithm>
#include <iostream>

#include <boost/numeric/conversion/bounds.hpp>

IUWTDeconvolutionAlgorithm::IUWTDeconvolutionAlgorithm(FFTWManager& fftwManager, size_t width, size_t height, double gain, double mGain, double cleanBorder, bool allowNegativeComponents, const bool* mask, double absoluteThreshold, double thresholdSigmaLevel, double tolerance, bool useSNRTest) :
	_fftwManager(fftwManager),
	_width(width), _height(height),
	_gain(gain), _mGain(mGain), _cleanBorder(cleanBorder),
	_mask(mask),
	_absoluteThreshold(absoluteThreshold),
	_thresholdSigmaLevel(thresholdSigmaLevel),
	_tolerance(tolerance), _allowNegativeComponents(allowNegativeComponents),
	_useSNRTest(useSNRTest)
{ }

void IUWTDeconvolutionAlgorithm::measureRMSPerScale(const double* image, const double* convolvedImage, double* scratch, size_t endScale, std::vector<ScaleResponse>& psfResponse)
{
	IUWTDecomposition imageIUWT(endScale, _width, _height);
	imageIUWT.Decompose(*_threadPool, image, scratch, false);

	_psfMaj = 2.0; _psfMin = 2.0; _psfPA = 0.0;
	double fl = 0.0;
	GaussianFitter fitter;
	fitter.Fit2DGaussianCentred(image, _width, _height, 2.0, _psfMaj, _psfMin, _psfPA);
	_psfVolume = (M_PI/4.0) * _psfMaj * _psfMin / M_LOG2E;
	
	double v=1.0, x=_width/2, y=_height/2;
	double bMaj = _psfMaj, bMin = _psfMin, bPA = _psfPA;
	fitter.Fit2DGaussianFull(image, _width, _height, v, x, y, bMaj, bMin, bPA, &fl);
	
	psfResponse.resize(endScale);
	for(size_t scale=0; scale!=endScale; ++scale)
	{
		psfResponse[scale].rms = rms(imageIUWT[scale].Coefficients());
		psfResponse[scale].peakResponse = centralPeak(imageIUWT[scale].Coefficients());
		bMaj = 2.0; bMin = 2.0; bPA = 0.0;
		v=1.0; x=_width/2; y=_height/2;
		fitter.Fit2DGaussianFull(imageIUWT[scale].Coefficients().data(), _width, _height, v, x, y, bMaj, bMin, bPA, &fl);
		psfResponse[scale].bMaj = bMaj;
		psfResponse[scale].bMin = bMin;
		psfResponse[scale].bPA = bPA;
		
	}
	
	imageIUWT.Decompose(*_threadPool, imageIUWT[1].Coefficients().data(), scratch, false);
	for(size_t scale=0; scale!=endScale; ++scale)
	{
		psfResponse[scale].peakResponseToNextScale = centralPeak(imageIUWT[scale].Coefficients());
	}
	
	imageIUWT.Decompose(*_threadPool, convolvedImage, scratch, false);
	
	for(size_t scale=0; scale!=endScale; ++scale)
	{
		psfResponse[scale].convolvedPeakResponse = centralPeak(imageIUWT[scale].Coefficients());
	}
	
	ao::uvector<double> thresholds(imageIUWT.NScales());
	for(size_t i=0; i!=imageIUWT.NScales(); ++i)
	{
		thresholds[i] = psfResponse[0].convolvedPeakResponse*_tolerance;
	}
	IUWTMask mask(imageIUWT.NScales(), _width, _height);
	ImageAnalysis::Component component(_width/2, _height/2, 0);
	size_t areaSize;
	ImageAnalysis::Floodfill(imageIUWT, mask, thresholds, 0, std::min<size_t>(endScale, 2), component, 0.0, areaSize);
	ao::uvector<bool> markedMask0(mask[0].size(), false);
	ImageAnalysis::Component2D c2D(_width/2, _height/2);
	double threshold = psfResponse[0].convolvedPeakResponse * _tolerance;
	ImageAnalysis::FloodFill2D(imageIUWT[0].Coefficients().data(), markedMask0.data(), threshold, c2D, _width, _height, psfResponse[0].convolvedArea);
}

double IUWTDeconvolutionAlgorithm::mad(const double* dest)
{
	ao::uvector<double> v(_width*_height);
	for(size_t i=0; i!=_width*_height; ++i)
		v[i] = std::fabs(dest[i]);
	size_t mid = (_width*_height)/2;
	std::nth_element(v.begin(), v.begin()+mid, v.end());
	return v[mid] / 0.674559;
}

double IUWTDeconvolutionAlgorithm::getMaxAbsWithoutMask(const ao::uvector<double>& data, size_t& x, size_t& y, size_t width)
{
	size_t height = data.size()/width;
	size_t xBorder = _cleanBorder*width;
	size_t yBorder = _cleanBorder*height;
	size_t minX = xBorder, maxX = width - xBorder;
	size_t minY = yBorder, maxY = height - yBorder;
	x = width;
	y = height;
	
	double maxVal = boost::numeric::bounds<double>::lowest();
	for(size_t yi=minY; yi!=maxY; ++yi) {
		const double* dataPtr = data.data() + yi*width;
		for(size_t xi=minX; xi!=maxX; ++xi) {
			double val = _allowNegativeComponents ? std::fabs(dataPtr[xi]) : dataPtr[xi];
			if(val > maxVal)
			{
				maxVal = val;
				x = xi;
				y = yi;
			}
		}
	}
	return maxVal;
}

double IUWTDeconvolutionAlgorithm::getMaxAbsWithMask(const ao::uvector<double>& data, size_t& x, size_t& y, size_t width)
{
	size_t height = data.size()/width;
	size_t xBorder = _cleanBorder*width;
	size_t yBorder = _cleanBorder*height;
	size_t minX = xBorder, maxX = width - xBorder;
	size_t minY = yBorder, maxY = height - yBorder;
	x = width;
	y = height;
	
	double maxVal = boost::numeric::bounds<double>::lowest();
	for(size_t yi=minY; yi!=maxY; ++yi) {
		const double* dataPtr = data.data() + yi*width;
		const bool* maskPtr = _mask + yi*width;
		for(size_t xi=minX; xi!=maxX; ++xi) {
			if(maskPtr[xi])
			{
				double val = _allowNegativeComponents ? std::fabs(dataPtr[xi]) : dataPtr[xi];
				if(val > maxVal)
				{
					maxVal = val;
					x = xi;
					y = yi;
				}
			}
		}
	}
	return maxVal;
}

double IUWTDeconvolutionAlgorithm::dotProduct(const ao::uvector<double>& lhs, const ao::uvector<double>& rhs)
{
	double sum = 0.0;
	for(size_t i=0; i!=lhs.size(); ++i)
		sum += lhs[i] * rhs[i];
	return sum;
}

void IUWTDeconvolutionAlgorithm::factorAdd(double* dest, const double* rhs, double factor, size_t width, size_t height)
{
	for(size_t i=0; i!=width*height; ++i)
		dest[i] += rhs[i] * factor;
}

void IUWTDeconvolutionAlgorithm::factorAdd(ao::uvector<double>& dest, const ao::uvector<double>& rhs, double factor)
{
	for(size_t i=0; i!=dest.size(); ++i)
		dest[i] += rhs[i] * factor;
}

void IUWTDeconvolutionAlgorithm::Subtract(double* dest, const ao::uvector<double>& rhs)
{
	for(size_t i=0; i!=rhs.size(); ++i)
		dest[i] -= rhs[i];
}

void IUWTDeconvolutionAlgorithm::boundingBox(size_t& x1, size_t& y1, size_t& x2, size_t& y2, const ao::uvector<double>& image, size_t width, size_t height)
{
	double mP = *std::max_element(image.begin(), image.end());
	double mN = *std::min_element(image.begin(), image.end());
	double m = std::max(mP, -mN);
	x1 = width; x2 = 0;
	y1 = height; y2 = 0;
	for(size_t y=0; y!=height; ++y)
	{
		const double* ptr = image.data() + y*width;
		for(size_t x=0; x!=x1; ++x)
		{
			if(std::fabs(ptr[x]) > m*0.01)
			{
				x1 = x;
				break;
			}
		}
		for(size_t x=width-1; x!=x2; --x)
		{
			if(std::fabs(ptr[x]) > m*0.01)
			{
				x2 = x;
				break;
			}
		}
	}
	x2++;
	for(size_t y=0; y!=height; ++y)
	{
		const double* ptr = image.data() + y*width;
		for(size_t x=0; x!=width; ++x)
		{
			if(std::fabs(ptr[x]) > m*0.01)
			{
				if(y1 > y)
					y1 = y;
				if(y2 < y)
					y2 = y+1;
			}
		}
	}
}

void IUWTDeconvolutionAlgorithm::adjustBox(size_t& x1, size_t& y1, size_t& x2, size_t& y2, size_t width, size_t height, int endScale)
{
	const int minBoxSize = std::max<int>(128, IUWTDecomposition::MinImageDimension(endScale)*3/2);
	
	int boxWidth = x2 - x1;
	int boxHeight = y2 - y1;
	int
		newX1 = x1 - 0.5*boxWidth,
		newX2 = x2 + 0.5*boxWidth,
		newY1 = y1 - 0.5*boxHeight,
		newY2 = y2 + 0.5*boxHeight;
	
	if(newX2 - newX1 < minBoxSize)
	{
		int mid = 0.5*(int(x1) + int(x2));
		newX1 = mid - minBoxSize/2;
		newX2 = mid + minBoxSize/2;
	}
	if(newY2 - newY1 < minBoxSize)
	{
		int mid = 0.5*(int(y1) + int(y2));
		newY1 = mid - minBoxSize/2;
		newY2 = mid + minBoxSize/2;
	}
	if(newX1 >= 0)
		x1 = newX1;
	else
		x1 = 0;
	if(newX2 < int(width))
		x2 = newX2;
	else
		x2 = width;
	if(newY1 >= 0)
		y1 = newY1;
	else
		y1 = 0;
	if(newY2 < int(height))
		y2 = newY2;
	else
		y2 = height;
	while((x2-x1)%8 != 0)
		x2--;
	while((y2-y1)%8 != 0)
		y2--;
}

void IUWTDeconvolutionAlgorithm::trim(ao::uvector<double>& dest, const double* source, size_t oldWidth, size_t /*oldHeight*/, size_t x1, size_t y1, size_t x2, size_t y2)
{
	// We do this so that dest and source can be the same image.
	if(dest.size() < (x2-x1) * (y2-y1))
		dest.resize((x2-x1) * (y2-y1));
	for(size_t y=y1; y!=y2; ++y)
	{
		const double* oldPtr = &source[y*oldWidth];
		double* newPtr = &dest[(y-y1)*(x2-x1)];
		for(size_t x=x1; x!=x2; ++x)
		{
			newPtr[x - x1] = oldPtr[x];
		}
	}
	dest.resize((x2-x1) * (y2-y1));
}

void IUWTDeconvolutionAlgorithm::untrim(ao::uvector<double>& image, size_t width, size_t height, size_t x1, size_t y1, size_t x2, size_t y2)
{
	image.resize(width*height, 0.0);
	size_t y=y2;
	while(y!=y1)
	{
		--y;
		double* newPtr = &image[y*width];
		double* oldPtr = &image[(y-y1)*(x2-x1)];
		size_t x=x2;
		while(x!=x1)
		{
			--x;
			newPtr[x] = oldPtr[x - x1];
		}
	}
	for(size_t y=0; y!=y1; ++y)
	{
		double* ptr = &image[y*width];
		for(size_t x=0; x!=width; ++x)
			ptr[x] = 0;
	}
	for(size_t y=y1; y!=y2; ++y)
	{
		double* ptr = &image[y*width];
		for(size_t x=0; x!=x1; ++x)
			ptr[x] = 0.0;
		for(size_t x=x2; x!=width; ++x)
			ptr[x] = 0.0;
	}
	for(size_t y=y2; y!=height; ++y)
	{
		double* ptr = &image[y*width];
		for(size_t x=0; x!=width; ++x)
			ptr[x] = 0;
	}
}

double IUWTDeconvolutionAlgorithm::sum(const ao::uvector<double>& img) const
{
	double s = 0.0;
	for(size_t i=0; i!=img.size(); ++i)
		s += img[i];
	return s;
}

double IUWTDeconvolutionAlgorithm::snr(const IUWTDecomposition& noisyImg, const IUWTDecomposition& model) const
{
	double mSum = 0.0, nSum = 0.0;
	for(size_t scale=0; scale!=noisyImg.NScales(); ++scale)
	{
		const ao::uvector<double>
			&n = noisyImg[scale].Coefficients(),
			&m = model[scale].Coefficients();
		for(size_t i=0; i!=n.size(); ++i)
		{
			mSum += m[i]*m[i];
			double d = m[i]-n[i];
			nSum += d*d;
		}
	}
	return mSum / nSum;
}

double IUWTDeconvolutionAlgorithm::rmsDiff(const ao::uvector<double>& a, const ao::uvector<double>& b)
{
	double sum = 0.0;
	for(size_t i=0; i!=a.size(); ++i)
	{
		double d = a[i]-b[i];
		sum += d*d;
	}
	return sqrt(sum/a.size());
}

double IUWTDeconvolutionAlgorithm::rms(const ao::uvector<double>& image)
{
	double sum = 0.0;
	for(size_t i=0; i!=image.size(); ++i)
	{
		double v = image[i];
		sum += v*v;
	}
	return sqrt(sum/image.size());
}

bool IUWTDeconvolutionAlgorithm::runConjugateGradient(IUWTDecomposition& iuwt, const IUWTMask& mask, ao::uvector<double>& maskedDirty, ao::uvector<double>& structureModel, ao::uvector<double>& scratch, const ao::uvector<double>& psfKernel, size_t width, size_t height)
{
	ao::uvector<double> gradient = maskedDirty;
	double modelSNR = 0.0;
	
	IUWTDecomposition initialDirtyIUWT(iuwt);
	
	for(size_t minorIter=0; minorIter!=20; ++minorIter)
	{
		// scratch = gradient (x) psf
		scratch = gradient;
		FFTConvolver::ConvolveSameSize(_fftwManager, scratch.data(), psfKernel.data(), width, height);
		
		// calc: IUWT gradient (x) psf
		iuwt.Decompose(*_threadPool, scratch.data(), scratch.data(), false);
		
		// calc: mask IUWT gradient (x) psf
		iuwt.ApplyMask(mask);
		
		// scratch = IUWT^-1 mask IUWT gradient (x) psf
		iuwt.Recompose(scratch, false);
		
		// stepsize = <residual, residual> / <gradient, scratch>
		double gradientDotScratch = dotProduct(gradient, scratch);
		if(gradientDotScratch == 0.0)
			return false;
		double stepSize = dotProduct(maskedDirty, maskedDirty) / gradientDotScratch;
		
		// model_i+1 = model_i + stepsize * gradient
		factorAdd(structureModel.data(), gradient.data(), stepSize, width, height);
		
		//For Dabbech's approach (see below) :
		//  ao::uvector<double> scratch2 = maskedDirty;
		
		double gradStepDen = dotProduct(maskedDirty, maskedDirty);
		if(gradStepDen == 0.0)
			return false;
		// residual_i+1 = residual_i - stepsize * scratch
		factorAdd(maskedDirty.data(), scratch.data(), -stepSize, width, height);
		
		// PyMORESANE uses this:
		// gradstep = <residual_i+1, residual_i+1> / <residual_i, residual_i>
		// double gradStep = dotProduct(maskedDirty, maskedDirty) / gradStepDen;
		// But in MORESANE's paper A. Dabbech says this:
		// gradstep = <residual_i+1 - residual_i, residual_i+1> / <residual_i, residual_i>
		// scratch = maskedDirty;
		// subtract(scratch, scratch2);
		// double gradStep = dotProduct(scratch, maskedDirty) / gradStepDen;
		double gradStep = dotProduct(maskedDirty, maskedDirty) / gradStepDen;
		
		// gradient_i+1 = residual_i+1 + gradstep * gradient_i
		scratch = gradient;
		gradient = maskedDirty;
		factorAdd(gradient.data(), scratch.data(), gradStep, width, height);
		
		// scratch = mask IUWT PSF (x) model
		scratch = structureModel;
		FFTConvolver::ConvolveSameSize(_fftwManager, scratch.data(), psfKernel.data(), width, height);
		iuwt.Decompose(*_threadPool, scratch.data(), scratch.data(), false);
		iuwt.ApplyMask(mask);
		
		double previousSNR = modelSNR;
		modelSNR = snr(iuwt, initialDirtyIUWT);
		if(modelSNR>100 && minorIter>2)
		{
			std::cout << "Converged after " << minorIter << " iterations.\n";
			return true;
		}
		else if(modelSNR < previousSNR && minorIter>5)
		{
			if(modelSNR > 3)
			{
				std::cout << "SNR decreased after " << minorIter << " iterations (SNR=" << modelSNR << ").\n";
				return true;
			}
		}
	}
	if(modelSNR <= 3.0)
	{
		std::cout << "Failed to converge (SNR=" << modelSNR << ").\n";
		structureModel.assign(width*height, 0.0);
		return false;
	}
	return true;
}

struct PointSource
{
	double x, y, flux;
	bool operator<(const PointSource& rhs) const
	{ return flux < rhs.flux; }
};

bool IUWTDeconvolutionAlgorithm::extractPointSources(const IUWTDecomposition& iuwt, const IUWTMask& mask, const double* dirty, double* model)
{
	size_t width = iuwt.Width(), height = iuwt.Height();
	IUWTMask markedMask(mask);
	bool pointSourcesWereFound = false;
	GaussianFitter posFitter;
	std::set<PointSource> sources;
	for(size_t y=0; y!=height; ++y)
	{
		for(size_t x=0; x!=width; ++x)
		{
			double flux;
			size_t maxX = x, maxY = y;
			if(markedMask[0][y*width + x] && ImageAnalysis::IsHighestOnScale0(iuwt, markedMask, maxX, maxY, iuwt.NScales(), flux))
			{
				size_t sourceArea;
				ao::uvector<bool> tmpMask(markedMask[0].size(), false);
				ImageAnalysis::Component2D c(maxX, maxY);
				ImageAnalysis::FloodFill2D(iuwt[0].Coefficients().data(), tmpMask.data(), flux*_tolerance, c, width, height, sourceArea);
				flux = flux / _psfResponse[0].peakResponse;
				if(sourceArea <= _psfResponse[0].convolvedArea*3/2 && flux>0)
				{
					// Fit the source
					double
						v = flux, xd = maxX, yd = maxY, fl = 0.0,
						bMaj = _psfResponse[0].bMaj, bMin = _psfResponse[0].bMin, bPA = _psfResponse[0].bPA;
					const double *d = dirty;
							//double *d = iuwt[0].Coefficients().data();
					posFitter.SetPosConstrained(2.0);
					posFitter.Fit2DGaussianFull(d, width, height, v, xd, yd, bMaj, bMin, bPA, &fl);
					bool badFit = std::fabs(maxX-xd) > 2.0 || std::fabs(maxY-yd) > 2.0 || v <= 0.0;
					if(!badFit)
					{
						double iFlux = v * bMaj * bMin * (M_PI/4.0) / (M_LOG2E * _psfVolume);
						std::cout << "Point source: " << flux << "/" << iFlux << " Jy (" << maxX << ',' << maxY << "), area=" << sourceArea << ", fit: v=" << v << ", x=" << xd << ", y=" << yd << ", ma=" << bMaj << ", mi=" << bMin << ", pa=" << bPA << "\n";
						//if(bMin * bMaj >= _psfMin * _psfMaj && bMin * bMaj < _psfMin * _psfMaj * 4.0)
						//	flux = 0.5*(iFlux + flux); // Compromise :-/
						if(std::fabs(flux) <= std::fabs(_rmses[0])*_thresholdSigmaLevel)
							std::cout << "Insignificant.\n";
						else {
							PointSource s;
							s.x = xd;
							s.y = yd;
							s.flux = flux;
							sources.insert(s);
						}
					}
				}
			}
		}
	}
	
	ao::uvector<bool>& m = markedMask[0];
	m.assign(m.size(), false);
	size_t nAcceptedSources = 0, nRejectedSources = 0;
	for(std::set<PointSource>::const_reverse_iterator sIter=sources.rbegin();
			sIter!=sources.rend(); ++sIter)
	{
		int halfBoxSize = ceil(_psfMin*2.0);
		size_t xi = round(sIter->x), yi = round(sIter->y);
		if(!m[xi + yi*width])
		{
			ModelRenderer::RenderInterpolatedSource(model, width, height, sIter->flux, sIter->x, sIter->y);
			pointSourcesWereFound = true;
			size_t xl = std::max<int>(int(xi)-halfBoxSize, 0);
			size_t yt = std::max<int>(int(yi)-halfBoxSize, 0);
			size_t xr = std::min<size_t>(xi+halfBoxSize, width);
			size_t yb = std::min<size_t>(yi+halfBoxSize, height);
			for(size_t y=yt; y!=yb; ++y)
			{
				bool* row = &m[y*width];
				for(size_t x=xl; x!=xr; ++x)
				{
					row[x] = true;
				}
			}
			++nAcceptedSources;
		}
		else {
			std::cout << "Source of " << sIter->flux << " Jy (" << sIter->x << ',' << sIter->y << ") too close (<" << halfBoxSize << " px) to other source.\n";
			++nRejectedSources;
		}
	}
	std::cout << "Subtracted sources: " << nAcceptedSources << ", rejected because of overlap: " << nRejectedSources << '\n';
	
	return pointSourcesWereFound;
}

void IUWTDeconvolutionAlgorithm::constrainedPSFConvolve(double* image, const double* psf, size_t width, size_t height)
{
	ao::uvector<double> smallerPsf(width*height, 0.0), kernel(width*height);
	size_t s = round(sqrt(_psfResponse[0].convolvedArea*25.0));
	size_t smallWidth = std::min(s, width);
	size_t smallHeight = std::min(s, height);
	std::cout << "Constrained PSF=" << smallWidth << " x " << smallHeight << '\n';
	size_t xMin = width/2 - smallWidth/2, xMax = width/2 + smallWidth/2;
	size_t yMin = height/2 - smallHeight/2, yMax = height/2 + smallHeight/2;
	for(size_t y=yMin; y!=yMax; ++y)
	{
		for(size_t x=xMin; x!=xMax; ++x)
		{
			smallerPsf[y*width + x] = psf[y*width + x];
		}
	}
	FFTConvolver::PrepareKernel(kernel.data(), smallerPsf.data(), width, height);
	FFTConvolver::ConvolveSameSize(_fftwManager, image, kernel.data(), width, height);
}

bool IUWTDeconvolutionAlgorithm::findAndDeconvolveStructure(IUWTDecomposition& iuwt, ao::uvector<double>& dirty, const ao::uvector<double>& psf, const ao::uvector<double>& psfKernel, ao::uvector<double>& scratch, ImageSet& structureModel, size_t curEndScale, size_t curMinScale, std::vector<IUWTDeconvolutionAlgorithm::ValComponent>& maxComponents)
{
	iuwt.Decompose(*_threadPool, dirty.data(), scratch.data(), false);
	ao::uvector<double> thresholds(curEndScale);
	_rmses.resize(curEndScale);
	for(size_t scale=0; scale!=curEndScale; ++scale)
	{
		double r = mad(iuwt[scale].Coefficients().data());
		_rmses[scale] = r;
		thresholds[scale] = r*(_thresholdSigmaLevel*4.0/5.0);
	}
	
	scratch = dirty;
	maxComponents.resize(curEndScale);
	for(size_t scale=0; scale!=curEndScale; ++scale)
	{
		size_t x, y;
		double maxAbsCoef = getMaxAbs(iuwt[scale].Coefficients(), x, y, _width);
		maxComponents[scale].x = x;
		maxComponents[scale].y = y;
		maxComponents[scale].scale = scale;
		maxComponents[scale].val = maxAbsCoef;
	}
	
	double maxVal = -1.0;
	size_t maxX = 0, maxY = 0;
	int maxValScale = -1;
	for(size_t scale=0; scale!=curEndScale; ++scale)
	{
		// Considerations for this section:
		// - Scale 0 should be chosen if the input corresponds to the PSF.
		//   Therefore, a peak on scale 1 should be at least:
		//   (PSF peak on scale 1) * (peak on scale 0) / (PSF (x) scale 1 peak response)
		//   Such that anything smaller than scale 1 will be considered scale 0.
		
		const ValComponent& val = maxComponents[scale];
		double absCoef = val.val/_psfResponse[scale].rms;
		//std::cout << scale << ">=" << curMinScale << " && " << absCoef << " > " << maxVal << " && " << val.val << " > " << _rmses[scale]*_thresholdLevel << "\n";
		if(scale>=curMinScale && absCoef > maxVal && val.val > _rmses[scale]*_thresholdSigmaLevel && val.val > _rmses[scale]/_rmses[0]*_absoluteThreshold)
		{
			maxX = val.x;
			maxY = val.y;
			maxValScale = scale;
			if(scale == 0)
			{
				double lowestRMS = std::min(_psfResponse[0].rms, _psfResponse[1].rms);
				maxVal = val.val/lowestRMS*_psfResponse[1].peakResponse/_psfResponse[0].peakResponseToNextScale;
			}
			else
				maxVal = absCoef;
		}
	}
	if(maxValScale == -1)
	{
		std::cout << "No significant pixel found.\n";
		return false;
	}
	
	maxVal = iuwt[maxValScale][maxX + maxY*_width];
	std::cout << "Most significant pixel: " << maxX << ',' << maxY << "=" << maxVal << " (" << maxVal/_rmses[maxValScale] << "σ) on scale " << maxValScale << '\n';
	
	if(std::fabs(maxVal) < thresholds[maxValScale])
	{
		std::cout << "Most significant pixel is in the noise, stopping.\n";
		return false;
	}

	double scaleMaxAbsVal = std::fabs(maxVal);
	for(size_t scale=0; scale!=curEndScale; ++scale)
	{
		if(thresholds[scale] < _tolerance * scaleMaxAbsVal)
		{
			thresholds[scale] = _tolerance * scaleMaxAbsVal;
		}
		if(maxVal < 0.0)
			thresholds[scale] = -thresholds[scale];
	}
	
	ImageAnalysis::Component maxComp(maxX, maxY, maxValScale);
	return fillAndDeconvolveStructure(iuwt, dirty, structureModel, scratch, psf, psfKernel, curEndScale, curMinScale, _width, _height, thresholds, maxComp, true, _mask);
}

bool IUWTDeconvolutionAlgorithm::fillAndDeconvolveStructure(IUWTDecomposition& iuwt, ao::uvector<double>& dirty, ImageSet& structureModelFull, ao::uvector<double>& scratch, const ao::uvector<double>& psf, const ao::uvector<double>& psfKernel, size_t curEndScale, size_t curMinScale, size_t width, size_t height, const ao::uvector<double>& thresholds, const ImageAnalysis::Component& maxComp, bool allowTrimming, const bool* priorMask)
{
	IUWTMask mask(curEndScale, width, height);
	size_t areaSize;
	ImageAnalysis::SelectStructures(iuwt, mask, thresholds, curMinScale, curEndScale, _cleanBorder, priorMask, areaSize);
	std::cout << "Flood-filled area contains " << areaSize << " significant components.\n";

	iuwt.ApplyMask(mask);
	iuwt.Recompose(scratch, false);
	
	// Find bounding box
	size_t x1, y1, x2, y2;
	boundingBox(x1, y1, x2, y2, scratch, width, height);
	adjustBox(x1, y1, x2, y2, width, height, maxComp.scale+1);
	if(allowTrimming && ((x2-x1)<width || (y2-y1)<height))
	{
		_curBoxXStart = x1; _curBoxXEnd = x2;
		_curBoxYStart = y1; _curBoxYEnd = y2;
		std::cout << "Bounding box: (" << x1 << ',' << y1 << ")-(" << x2 << ',' << y2 << ")\n";
		size_t newWidth = x2-x1, newHeight = y2-y1;
		trim(dirty, dirty, width, height, x1, y1, x2, y2);
		ao::uvector<double> smallPSF;
		
		trimPsf(smallPSF, psf.data(), width, height, newWidth, newHeight);
		
		ao::uvector<double> smallPSFKernel(smallPSF.size());
		FFTConvolver::PrepareKernel(smallPSFKernel.data(), smallPSF.data(), newWidth, newHeight);
		
		scratch.resize(dirty.size());
		
		int maxScale = std::max(IUWTDecomposition::EndScale(std::min(x2-x1, y2-y1)), maxComp.scale+1);
		if(maxScale < int(curEndScale))
		{
			std::cout << "Bounding box too small for largest scale of " << curEndScale << " -- ignoring scales>=" << maxScale << " in this iteration.\n";
			curEndScale = maxScale;
		}
		std::unique_ptr<IUWTDecomposition> trimmedIUWT(iuwt.CreateTrimmed(curEndScale, x1, y1, x2, y2));
		
		std::unique_ptr<ImageSet> trimmedStructureModel(structureModelFull.Trim(x1, y1, x2, y2, width));
		
		ao::uvector<bool> trimmedPriorMask;
		bool* trimmedPriorMaskPtr;
		if(priorMask == 0)
			trimmedPriorMaskPtr = 0;
		else
		{
			trimmedPriorMask.resize(newWidth * newHeight);
			trimmedPriorMaskPtr = trimmedPriorMask.data();
			Image::TrimBox(trimmedPriorMaskPtr, x1, y1, newWidth, newHeight, priorMask, width, height);
		}

		ImageAnalysis::Component newMaxComp(maxComp.x-x1, maxComp.y-y1, maxComp.scale);
		bool result = fillAndDeconvolveStructure(*trimmedIUWT, dirty, *trimmedStructureModel, scratch, smallPSF, smallPSFKernel, curEndScale, curMinScale, x2-x1, y2-y1, thresholds, newMaxComp, false, trimmedPriorMaskPtr);
		for(size_t i=0; i!=structureModelFull.size(); ++i)
		{
			memcpy(scratch.data(), (*trimmedStructureModel)[i], (y2-y1)*(x2-x1)*sizeof(double));
			untrim(scratch, width, height, x1, y1, x2, y2);
			memcpy(structureModelFull[i], scratch.data(), width*height*sizeof(double));
		}
		
		dirty.resize(scratch.size());
		_curBoxXStart = 0; _curBoxXEnd = width;
		_curBoxYStart = 0; _curBoxYEnd = height;
		return result;
	}
	else {
		if(curEndScale <= 3)
		{
			//bool pointSourcesWereFound = extractPointSources(iuwt, mask, dirty.data(), structureModel.data());
			//if(pointSourcesWereFound)
			//	return true;
		}
		
		// get undeconvolved dirty
		iuwt.Decompose(*_threadPool, dirty.data(), scratch.data(), false);
		
		iuwt.ApplyMask(mask);
		iuwt.Recompose(scratch, false);
		
		ao::uvector<double> maskedDirty = scratch;
		
		ao::uvector<double> structureModel(width*height, 0.0);
		bool success = runConjugateGradient(iuwt, mask, maskedDirty, structureModel, scratch, psfKernel, width, height);
		if(!success) return false;
		
		double rmsBefore = rms(dirty);
		scratch = structureModel;
		FFTConvolver::ConvolveSameSize(_fftwManager, scratch.data(), psfKernel.data(), width, height);
		maskedDirty = dirty; // we use maskedDirty as temporary
		factorAdd(maskedDirty.data(), scratch.data(), -_gain, width, height);
		double rmsAfter = rms(maskedDirty);
		if(rmsAfter > rmsBefore)
		{
			std::cout << "RMS got worse: " << rmsBefore << " -> " << rmsAfter << '\n';
			return false;
		}
		
		// TODO when only one image is available, this is not necessary
		performSubImageFitAll(iuwt, mask, structureModel, scratch, maskedDirty, maxComp, structureModelFull, psf.data(), dirty);
		
		return true;
	}
}

void IUWTDeconvolutionAlgorithm::performSubImageFitAll(IUWTDecomposition& iuwt, const IUWTMask& mask, const ao::uvector<double>& structureModel, ao::uvector<double>& scratchA, ao::uvector<double>& scratchB, const ImageAnalysis::Component& maxComp, ImageSet& fittedModel, const double* psf, const ao::uvector<double>& dirty)
{
	size_t width = iuwt.Width(), height = iuwt.Height();
	
	if(_dirtySet->size() == 1)
	{
		// With only one image, we don't have to refit
		fittedModel.Set(0, structureModel.data());
	}
	else {
		std::cout << "Fitting structure in images: " << std::flush;
		ao::uvector<double> correctionFactors;
		scratchA = dirty;
		performSubImageFitSingle(iuwt, mask, structureModel, scratchB, maxComp, psf, scratchA, 0, correctionFactors);
			
		fittedModel = 0.0;
		
		for(size_t imgIndex=0; imgIndex!=_dirtySet->size(); ++imgIndex)
		{
			std::cout << '.' << std::flush;
			const double* subPsf = _psfs[_dirtySet->PSFIndex(imgIndex)];
			
			trim(scratchA, (*_dirtySet)[imgIndex], _width, _height, _curBoxXStart, _curBoxYStart, _curBoxXEnd, _curBoxYEnd);
			
			ao::uvector<double> smallSubPsf;
			const double *subPsfData;
			if(_width != width || _height != height)
			{
				trimPsf(smallSubPsf, subPsf, _width, _height, width, height);
				subPsfData = smallSubPsf.data();
			}
			else {
				subPsfData = subPsf;
			}
		
			performSubImageFitSingle(iuwt, mask, structureModel, scratchB, maxComp, subPsfData, scratchA, fittedModel[imgIndex], correctionFactors);
		}
		std::cout << '\n';
	}
}

void IUWTDeconvolutionAlgorithm::performSubImageFitSingle(IUWTDecomposition& iuwt, const IUWTMask& mask, const ao::uvector<double>& structureModel, ao::uvector<double>& scratchB, const ImageAnalysis::Component& maxComp, const double* psf, ao::uvector<double>& subDirty, double* fittedSubModel, ao::uvector<double>& correctionFactors)
{
	size_t width = iuwt.Width(), height = iuwt.Height();
	
	ao::uvector<double> psfKernel(width*height);
	FFTConvolver::PrepareKernel(psfKernel.data(), psf, width, height);
		
	ao::uvector<double>& maskedDirty = scratchB;
	
	iuwt.Decompose(*_threadPool, subDirty.data(), subDirty.data(), false);
	iuwt.ApplyMask(mask);
	iuwt.Recompose(maskedDirty, false);
	ao::uvector<bool> mask2d(structureModel.size(), false);
	double peakLevel = std::fabs(structureModel[maxComp.y*width + maxComp.x]);
	size_t componentIndex = 0;
	for(size_t y=0; y!=height; ++y)
	{
		bool* maskRow = &mask2d[y*width];
		const double* modelRow = &structureModel[y*width];
		for(size_t x=0; x!=width; ++x)
		{
			if(!maskRow[x] && std::fabs(modelRow[x]) > peakLevel*1e-4)
			{
				std::vector<ImageAnalysis::Component2D> area;
				ImageAnalysis::Component2D comp(x, y);
				ImageAnalysis::FloodFill2D(structureModel.data(), mask2d.data(), peakLevel*1e-4, comp, width, height, area);
				// Find bounding box and copy active pixels to subDirty
				subDirty.assign(width*height, 0.0);
				size_t boxX1=width, boxX2=0, boxY1=height, boxY2=0;
				for(std::vector<ImageAnalysis::Component2D>::const_iterator a=area.begin(); a!=area.end(); ++a)
				{
					size_t index = a->x + a->y*width;
					boxX1 = std::min(a->x, boxX1);
					boxX2 = std::max(a->x ,boxX2);
					boxY1 = std::min(a->y, boxY1);
					boxY2 = std::max(a->y, boxY2);
					subDirty[index] = structureModel[index];
				}
				adjustBox(boxX1, boxY1, boxX2, boxY2, width, height, iuwt.NScales());
				
				double factor = performSubImageComponentFitBoxed(iuwt, mask, area, subDirty, maskedDirty, psf, psfKernel, boxX1, boxY1, boxX2, boxY2);
				
				// if no fittedSubModel was given, we just need to store the factors. Otherwise,
				// scale the deconvolved model and add it to the contaminated model.
				if(fittedSubModel != 0)
				{
					const double integratedFactor = correctionFactors[componentIndex];
					if(std::isfinite(factor) && std::isfinite(integratedFactor) && integratedFactor != 0.0)
					{
						for(std::vector<ImageAnalysis::Component2D>::const_iterator a=area.begin(); a!=area.end(); ++a)
						{
							size_t index = a->x + a->y*width;
							fittedSubModel[index] += structureModel[index]*factor/integratedFactor;
						}
					}
					++componentIndex;
				}
				else {
					correctionFactors.push_back(factor);
				}
			}
		}
	}
}

double IUWTDeconvolutionAlgorithm::performSubImageComponentFitBoxed(IUWTDecomposition& iuwt, const IUWTMask& mask, const std::vector<ImageAnalysis::Component2D>& area, ao::uvector<double>& model, ao::uvector<double>& maskedDirty, const double* psf, const ao::uvector<double>& psfKernel, size_t x1, size_t y1, size_t x2, size_t y2)
{
	const size_t width = iuwt.Width(), height = iuwt.Height();
	if(x1 > 0 || y1 > 0 || x2 < width || y2 < height)
	{
		size_t newWidth = x2-x1, newHeight = y2-y1;
		IUWTDecomposition smallIUWTW(iuwt.NScales(), newWidth, newHeight);
		std::unique_ptr<IUWTMask> smallMask(mask.CreateTrimmed(x1, y1, x2, y2));
		ao::uvector<double> smallModel;
		trim(smallModel, model, width, height, x1, y1, x2, y2);
		
		ao::uvector<double> smallPsf;
		trimPsf(smallPsf, psf, width, height, newWidth, newHeight);
		ao::uvector<double> smallPsfKernel(smallPsf.size());
		FFTConvolver::PrepareKernel(smallPsfKernel.data(), smallPsf.data(), newWidth, newHeight);
		
		ao::uvector<double> smallMaskedDirty;
		trim(smallMaskedDirty, maskedDirty, width, height, x1, y1, x2, y2);
		
		double factor = performSubImageComponentFit(smallIUWTW, *smallMask, area, smallModel, smallMaskedDirty, smallPsfKernel, x1, y1);
		return factor;
	}
	else {
		return performSubImageComponentFit(iuwt, mask, area, model, maskedDirty, psfKernel, 0, 0);
	}
}

double IUWTDeconvolutionAlgorithm::performSubImageComponentFit(IUWTDecomposition& iuwt, const IUWTMask& mask, const std::vector<ImageAnalysis::Component2D>& area, ao::uvector<double>& model, ao::uvector<double>& maskedDirty, const ao::uvector<double>& psfKernel, size_t xOffset, size_t yOffset)
{
	const size_t width = iuwt.Width(), height = iuwt.Height();
	// Calculate IUWT^-1 mask IUWT model (x) PSF
	FFTConvolver::ConvolveSameSize(_fftwManager, model.data(), psfKernel.data(), width, height);
	iuwt.Decompose(*_threadPool, model.data(), model.data(), false);
	iuwt.ApplyMask(mask);
	iuwt.Recompose(model, false);
	
	double modelSum = 0.0, dirtySum = 0.0;
	for(std::vector<ImageAnalysis::Component2D>::const_iterator a=area.begin(); a!=area.end(); ++a)
	{
		size_t index = (a->x-xOffset) + (a->y-yOffset)*width;
		modelSum += model[index];
		dirtySum += maskedDirty[index];
	}
	//std::cout << "factor=" << dirtySum << " / " << modelSum << " = " << dirtySum/modelSum << '\n';
	if(modelSum == 0.0 || !std::isfinite(dirtySum) || !std::isfinite(modelSum))
		return 0.0;
	else
		return dirtySum / modelSum;
}

double IUWTDeconvolutionAlgorithm::PerformMajorIteration(size_t& iterCounter, size_t nIter, ImageSet& modelSet, ImageSet& dirtySet, const ao::uvector<const double*>& psfs, bool& reachedMajorThreshold)
{
	FFTWManager::ThreadingScope fftwThreadsEnabled(_fftwManager);
	std::unique_ptr<ThreadPool> threadPool(new ThreadPool());
	_threadPool = &*threadPool;
	
	reachedMajorThreshold = false;
	if(iterCounter == nIter)
		return 0.0;
	
	_modelSet = &modelSet;
	_dirtySet = &dirtySet;
	_psfs = psfs;
	
	_curBoxXStart = 0; _curBoxXEnd = _width;
	_curBoxYStart = 0; _curBoxYEnd = _height;
		
	ao::uvector<double> dirty(_width * _height);
	dirtySet.GetLinearIntegrated(dirty.data());
	ao::uvector<double> psf(psfs[0], psfs[0] + _width*_height);
	dirtySet.GetIntegratedPSF(psf.data(), psfs);
	/*for(size_t i=1; i!=psfs.size(); ++i)
	{
		for(size_t j=0; j!=psf.size(); ++j)
			psf[j] += psfs[i][j];
	}
	for(size_t j=0; j!=psf.size(); ++j)
		psf[j] /= double(psfs.size());*/
	
	int maxScale = IUWTDecomposition::EndScale(std::min(_width, _height));
	int curEndScale = 2;
	
	// Prepare the PSF for convolutions later on
	ao::uvector<double> psfKernel(_width * _height);
	FFTConvolver::PrepareKernel(psfKernel.data(), psf.data(), _width, _height);
	
	std::cout << "Measuring PSF...\n";
	{
		ao::uvector<double> convolvedPSF(psf);
		ao::uvector<double> scratch(_width * _height);
		
		FFTConvolver::ConvolveSameSize(_fftwManager, convolvedPSF.data(), psfKernel.data(), _width, _height);
		measureRMSPerScale(psf.data(), convolvedPSF.data(), scratch.data(), maxScale, _psfResponse);
	}
	
	ImageSet structureModel(
		&modelSet.Table(), dirtySet.Allocator(),
		modelSet.ChannelsInDeconvolution(), modelSet.SquareJoinedChannels(),
		modelSet.LinkedPolarizations(),
		_width, _height);
	
	std::unique_ptr<IUWTDecomposition> iuwt(new IUWTDecomposition(curEndScale, _width, _height));
	
	ao::uvector<double> dirtyBeforeIteration;
	
	double maxValue = 0.0;
	size_t curMinScale = 0;
	reachedMajorThreshold = false;
	bool doContinue = true;
	std::vector<ValComponent> initialComponents;
	do
	{
		std::cout << "*** Deconvolution iteration " << iterCounter << " ***\n";
		dirtyBeforeIteration = dirty;
		FFTConvolver::PrepareKernel(psfKernel.data(), psf.data(), _width, _height);
		std::vector<ValComponent> maxComponents;
		ao::uvector<double> scratch(_width * _height);
		bool succeeded = findAndDeconvolveStructure(*iuwt, dirty, psf, psfKernel, scratch, structureModel, curEndScale, curMinScale, maxComponents);
		
		if(succeeded)
		{
			structureModel *= _gain;
			modelSet += structureModel;
		
			// Calculate: dirty = dirty - structureModel (x) psf
			for(size_t i=0; i!=dirtySet.size(); ++i)
			{
				scratch.assign(structureModel[i], structureModel[i] + _width*_height);
				size_t psfIndex = dirtySet.PSFIndex(i);
				FFTConvolver::PrepareKernel(psfKernel.data(), psfs[psfIndex], _width, _height);
				FFTConvolver::ConvolveSameSize(_fftwManager, scratch.data(), psfKernel.data(), _width, _height);
				Subtract(dirtySet[i], scratch);
			}
			dirtySet.GetLinearIntegrated(dirty.data());
			
			while(maxComponents.size() > initialComponents.size())
			{
				initialComponents.push_back(maxComponents[initialComponents.size()]);
			}
			maxValue = 0.0;
			for(size_t c=0; c!=initialComponents.size(); ++c)
			{
				std::cout << initialComponents[c].val << " now " << maxComponents[c].val << '\n';
				maxValue = std::max(maxValue, maxComponents[c].val);
				if(std::fabs(maxComponents[c].val) < std::fabs(initialComponents[c].val) * (1.0 - _mGain))
				{
					std::cout << "Scale " << c << " reached mGain (starting level: " << initialComponents[c].val << ", now: " << maxComponents[c].val << ").\n";
					reachedMajorThreshold = true;
				}
			}
			if(reachedMajorThreshold) break;
		}
		else {
			if(int(curMinScale)+1 < curEndScale)
			{
				++curMinScale;
				std::cout << "=> Min scale now " << curMinScale << '\n';
			}
			else {
				curMinScale = 0;
				if(curEndScale != maxScale)
				{
					++curEndScale;
					std::cout << "=> Scale now " << curEndScale << ".\n";
					iuwt.reset(new IUWTDecomposition(curEndScale, _width, _height));
				}
				else {
					std::cout << "Max scale reached: finished all scales, quiting.\n";
					doContinue = false;
				}
			}
			dirty = dirtyBeforeIteration;
		}
		
		++iterCounter;
	} while(iterCounter!=nIter && doContinue);
	return maxValue;
}
