#include "lsdeconvolution.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_multifit.h>

#define OUTPUT_LSD_DEBUG_INFO 1

struct LSDeconvolutionData
{
	LSDeconvolution* parent;
	gsl_multifit_fdfsolver* solver;
	ao::uvector<std::pair<size_t, size_t>> maskPositions;
	size_t width, height;
	double regularization;
	const double* dirty;
	const double* psf;
	
	static int fitting_func(const gsl_vector *xvec, void *data, gsl_vector *f)
	{
		const LSDeconvolutionData &lsData = *reinterpret_cast<LSDeconvolutionData*>(data);

		size_t
			i = 0,
			midX = lsData.width + (lsData.width/2),
			midY = lsData.height + (lsData.height/2);
#ifdef OUTPUT_LSD_DEBUG_INFO
		double rmsSum = 0.0, cost = 0.0, peakFlux = 0.0;
#endif
		// e = (y - sum modelledflux)^2 + mu modelledflux
		for(size_t y=0; y!=lsData.height; ++y)
		{
			for(size_t x=0; x!=lsData.width; ++x)
			{
				double modelledFlux = 0.0;
				for(size_t p=0; p!=lsData.maskPositions.size(); ++p)
				{
					int
						pX = lsData.maskPositions[p].first,
						pY = lsData.maskPositions[p].second;
					double pVal = gsl_vector_get(xvec, p);
					int psfX = (x + midX - pX) % lsData.width;
					int psfY = (y + midY - pY) % lsData.height;
					modelledFlux += lsData.psf[psfX + psfY*lsData.width] * pVal;
				}
#ifdef OUTPUT_LSD_DEBUG_INFO
				double pixelCost = (lsData.dirty[i] - modelledFlux) * (lsData.dirty[i] - modelledFlux);
				rmsSum += pixelCost;
				cost += pixelCost;
#endif
				
				gsl_vector_set(f, i, lsData.dirty[i] - modelledFlux);
				++i;
			}
		}
		double totalModelledFlux = 0.0;
		for(size_t p=0; p!=lsData.maskPositions.size(); ++p)
		{
			totalModelledFlux += std::fabs(gsl_vector_get(xvec, p));
#ifdef OUTPUT_LSD_DEBUG_INFO
			cost += std::fabs(gsl_vector_get(xvec, p))*lsData.regularization;
			peakFlux = std::max(peakFlux, gsl_vector_get(xvec, p));
#endif
		}
		//gsl_vector_set(f, lsData.width*lsData.height, sqrt(lsData.regularization*totalModelledFlux));
		gsl_vector_set(f, lsData.width*lsData.height, lsData.regularization*totalModelledFlux);
		
#ifdef OUTPUT_LSD_DEBUG_INFO
		Logger::Debug << "Current RMS=" << sqrt(rmsSum / (lsData.height*lsData.width)) << ", mean flux in model=" << totalModelledFlux/lsData.maskPositions.size() << ", peak=" << peakFlux << ", total cost=" << cost << '\n';
#endif
		
		return GSL_SUCCESS;
	}
	
	static int fitting_func_deriv(const gsl_vector *xvec, void *data, gsl_matrix *J)
	{
#ifdef OUTPUT_LSD_DEBUG_INFO
		//Logger::Info << "Calculating Jacobian... " << std::flush;
#endif
		const LSDeconvolutionData &lsData = *reinterpret_cast<LSDeconvolutionData*>(data);

		size_t
			i = 0,
			midX = lsData.width + (lsData.width/2),
			midY = lsData.height + (lsData.height/2);
		for(size_t y=0; y!=lsData.height; ++y)
		{
			for(size_t x=0; x!=lsData.width; ++x)
			{
				for(size_t p=0; p!=lsData.maskPositions.size(); ++p)
				{
					int
						pX = lsData.maskPositions[p].first,
						pY = lsData.maskPositions[p].second;
					//double pVal = gsl_vector_get(xvec, p);
					int psfX = (x + midX - pX) % lsData.width;
					int psfY = (y + midY - pY) % lsData.height;
					
					gsl_matrix_set(J, i, p, -lsData.psf[psfX + psfY*lsData.width]);
				}
				++i;
			}
		}
		for(size_t p=0; p!=lsData.maskPositions.size(); ++p)
		{
			// f = sqrt | pval |
			//   = sqrt sqrt(pval^2)
			// f'= 2pval * 0.5/sqrt(pval^2) * 0.5/sqrt(sqrt(pval^2))
			// f'= pval/( 2.0 * |pval| * sqrt(|pval|) )
			double dpval = gsl_vector_get(xvec, p);
			/*if(dpval < 0.0)
				dpval = -0.5/sqrt(-dpval);
			else if(dpval > 0.0)
				dpval = 0.5/sqrt(dpval);
			else
				dpval = 0.0;*/
			if(dpval < 0.0)
				dpval = -dpval;
			//else if(dpval > 0.0)
			//	dpval = dpval;
			//else
			//	dpval = 0.0;
			gsl_matrix_set(J, lsData.width*lsData.height, p, dpval * lsData.regularization);
		}
#ifdef OUTPUT_LSD_DEBUG_INFO
		//Logger::Info << "DONE\n";
#endif
		return GSL_SUCCESS;
	}
	
	static int fitting_func_both(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J)
	{
		fitting_func(x, data, f);
		fitting_func_deriv(x, data, J);
		return GSL_SUCCESS;
	}
};

LSDeconvolution::LSDeconvolution() : _allocator(nullptr), _data(new LSDeconvolutionData())
{
}

LSDeconvolution::~LSDeconvolution()
{
}

void LSDeconvolution::getMaskPositions(ao::uvector<std::pair<size_t, size_t>>& maskPositions, const bool* mask, size_t width, size_t height)
{
	const bool* maskPtr = mask;
	for(size_t y=0; y!=height; ++y)
	{
		for(size_t x=0; x!=width; ++x)
		{
			if(*maskPtr)
			{
				maskPositions.push_back(std::make_pair(x, y));
			}
			++maskPtr;
		}
	}
}

void LSDeconvolution::linearFit(double* dataImage, double* modelImage, const double* psfImage, size_t width, size_t height, bool& reachedMajorThreshold)
{
	ao::uvector<std::pair<size_t, size_t>> maskPositions;
	getMaskPositions(maskPositions, _cleanMask, width, height);
	Logger::Info << "Running LSDeconvolution with " << maskPositions.size() << " parameters.\n";
	
	// y = X c
	//   - y is vector of N,     N=number of data points (pixels in image)
	//     y_i = pixel value i
	//   - x is vector of N x M, M=number of parameters (pixels in mask)
	//     x_ij = (pixel value i) * (psf value j)
	
	size_t n = width*height;
	size_t pn = maskPositions.size();
	
	gsl_matrix* x = gsl_matrix_calloc(n, pn);
	gsl_vector* y = gsl_vector_alloc(n);
	gsl_vector* c = gsl_vector_calloc(pn);
	gsl_matrix* cov = gsl_matrix_alloc(pn, pn);
	double chisq;
	gsl_multifit_linear_workspace* work = gsl_multifit_linear_alloc(n, pn);
	
	for(size_t i=0; i!=n; ++i)
		gsl_vector_set(y, i, dataImage[i]);
	
	size_t
		i = 0,
		midX = width + (width/2),
		midY = height + (height/2);
	for(size_t yi=0; yi!=height; ++yi)
	{
		for(size_t xi=0; xi!=width; ++xi)
		{
			for(size_t p=0; p!=pn; ++p)
			{
				int
					pX = maskPositions[p].first,
					pY = maskPositions[p].second;
				int psfX = (xi + midX - pX) % width;
				int psfY = (yi + midY - pY) % height;
				
				gsl_matrix_set(x, i, p, psfImage[psfX + psfY*width]);
			}
			++i;
		}
	}
	
	Logger::Info << "psf(0,0) = " << psfImage[midX%width + (midY%height)*width] << "\n";
	Logger::Info << "Fitting... ";
	Logger::Info.Flush();
	int result = gsl_multifit_linear(x, y, c, cov, &chisq, work);
	Logger::Info << "result=" << gsl_strerror(result) << "\n";
	gsl_multifit_linear_free(work);
	
	for(size_t i=0; i!=n; ++n)
		modelImage[i] = 0.0;
	
	for(size_t p=0; p!=pn; ++p)
	{
		size_t
			pX = maskPositions[p].first,
			pY = maskPositions[p].second;
		modelImage[pY*width + pX] = gsl_vector_get(c, p);
	}
	
	for(size_t y=0; y!=height; ++y)
	{
		size_t index = y*width;
		for(size_t x=0; x!=width; ++x)
		{
			double val = dataImage[index];
			for(size_t p=0; p!=pn; ++p)
			{
				int
					pX = maskPositions[p].first,
					pY = maskPositions[p].second;
				double pVal = gsl_vector_get(c, p);
				int psfX = (x + midX - pX) % width;
				int psfY = (y + midY - pY) % height;
				val -= psfImage[psfX + psfY*width] * pVal;
			}
			dataImage[index] = val;
			++index;
		}
	}
	
	gsl_matrix_free(x);
	gsl_vector_free(y);
	gsl_vector_free(c);
	gsl_matrix_free(cov);
}

void LSDeconvolution::nonLinearFit(double* dataImage, double* modelImage, const double* psfImage, size_t width, size_t height, bool& reachedMajorThreshold)
{
	if(this->_cleanMask == 0)
		throw std::runtime_error("No mask available");
		
	getMaskPositions(_data->maskPositions, _cleanMask, width, height);
	size_t parameterCount = _data->maskPositions.size(), dataCount = width * height + 1;
	Logger::Info << "Running LSDeconvolution with " << parameterCount << " parameters.\n";
	
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
	_data->solver = gsl_multifit_fdfsolver_alloc (T, dataCount, parameterCount);
	_data->dirty = dataImage;
	_data->psf = psfImage;
	_data->width = width;
	_data->height = height;
	_data->parent = this;
	_data->regularization = 0.1;
	
	gsl_multifit_function_fdf fdf;
	fdf.f = &LSDeconvolutionData::fitting_func;
	fdf.df = &LSDeconvolutionData::fitting_func_deriv;
	fdf.fdf = &LSDeconvolutionData::fitting_func_both;
	fdf.n = dataCount;
	fdf.p = parameterCount;
	fdf.params = &*_data;
	
	ao::uvector<double> parameterArray(parameterCount, 0.0);
	gsl_vector_view initialVals = gsl_vector_view_array(parameterArray.data(), parameterCount);
	gsl_multifit_fdfsolver_set(_data->solver, &fdf, &initialVals.vector);

	int status;
	size_t iter = 0;
	do {
		iter++;
		status = gsl_multifit_fdfsolver_iterate (_data->solver);
		
		if(status)
			break;
		
		status = gsl_multifit_test_delta(_data->solver->dx, _data->solver->x, 1e-4, 1e-4);
		
  } while (status == GSL_CONTINUE && iter < 100);
	Logger::Info << "niter=" << iter << ", status=" << gsl_strerror(status) << "\n";
	
	for(size_t p=0; p!=parameterCount; ++p)
	{
		size_t
			pX = _data->maskPositions[p].first,
			pY = _data->maskPositions[p].second;
		modelImage[pY*width + pX] = gsl_vector_get(_data->solver->x, p);
	}
	
	size_t
		midX = width + (width/2),
		midY = height + (height/2);
	for(size_t y=0; y!=height; ++y)
	{
		size_t index = y*width;
		for(size_t x=0; x!=width; ++x)
		{
			double val = dataImage[index];
			for(size_t p=0; p!=parameterCount; ++p)
			{
				int
					pX = _data->maskPositions[p].first,
					pY = _data->maskPositions[p].second;
				double pVal = gsl_vector_get(_data->solver->x, p);
				int psfX = (x + midX - pX) % width;
				int psfY = (y + midY - pY) % height;
				val -= psfImage[psfX + psfY*width] * pVal;
			}
			dataImage[index] = val;
			++index;
		}
	}
	
	gsl_multifit_fdfsolver_free(_data->solver);
}
