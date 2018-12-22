#include "gaussianfitter.h"

void GaussianFitter::ToAnglesAndFWHM(double sx, double sy, double beta, double& ellipseMaj, double& ellipseMin, double& ellipsePA)
{
	//std::cout << "conv sx=" << sx << ", sy=" << sy << ", beta=" << beta << '\n';
	const long double sigmaToBeam = 2.0L * sqrtl(2.0L * logl(2.0L));
	const double betaFact = 1.0 - beta*beta;
	double cov[4];
	cov[0] = sx*sx / betaFact;
	cov[1] = beta * sx*sy / betaFact;
	cov[2] = cov[1];
	cov[3] = sy*sy / betaFact;
	
	double e1, e2, vec1[2], vec2[2];
	Matrix2x2::EigenValuesAndVectors(cov, e1, e2, vec1, vec2);
	if(std::isfinite(e1))
	{
		ellipseMaj = sqrt(std::fabs(e1)) * sigmaToBeam;
		ellipseMin = sqrt(std::fabs(e2)) * sigmaToBeam;
		if(ellipseMaj < ellipseMin)
		{
			std::swap(ellipseMaj, ellipseMin);
			vec1[0] = vec2[0];
			vec1[1] = vec2[1];
		}
		ellipsePA = -atan2(vec1[0], vec1[1]);
	}
	else {
		ellipseMaj = sqrt(std::fabs(sx)) * sigmaToBeam;
		ellipseMin = sqrt(std::fabs(sx)) * sigmaToBeam;
		ellipsePA = 0.0;
	}
}

void GaussianFitter::ToFWHM(double s, double& beamSize)
{
	const long double sigmaToBeam = 2.0L * sqrtl(2.0L * logl(2.0L));
	beamSize = s * sigmaToBeam;
}

void GaussianFitter::ToCovariance(double fwhmMaj, double fwhmMin, double positionAngle, double& sxsx, double& sxsy, double& sysy)
{
	const long double sigmaToBeam = 2.0L * sqrtl(2.0L * logl(2.0L));
	fwhmMaj /= sigmaToBeam;
	fwhmMin /= sigmaToBeam;
	const double
		saSq = fwhmMaj * fwhmMaj,
		sbSq = fwhmMin * fwhmMin,
		cosAngle = std::cos(positionAngle), sinAngle = std::sin(positionAngle),
		cosAngleSq = cosAngle*cosAngle, sinAngleSq = sinAngle*sinAngle;
	// The covariance matrix is Rotation x Scaling x Rotation^T:
	// sxsx sxsy    cos a -sin a  maj^2 0    cos a  sin a
	// sxsy sysy = (sin a  cos a)  0 min^2 (-sin a  cos a)
	//
	//             maj^2 cos a  -min^2 sin a     cos a  sin a)
	//          =( maj^2 sin a   min^2 cos a ) (-sin a  cos a)
	sxsx = saSq * cosAngleSq + sbSq * sinAngleSq;
	sxsy = saSq * sinAngle * cosAngle - sbSq * cosAngle * sinAngle;
	sysy = saSq * sinAngleSq + sbSq * cosAngleSq;
}

void GaussianFitter::FromCovariance(double sxsx, double sxsy, double sysy, double& fwhmMaj, double& fwhmMin, double& positionAngle)
{
	const long double sigmaToBeam = 2.0L * sqrtl(2.0L * logl(2.0L));
	double cov[4];
	cov[0] = sxsx;
	cov[1] = sxsy;
	cov[2] = sxsy;
	cov[3] = sysy;
	
	double e1, e2, vec1[2], vec2[2];
	Matrix2x2::EigenValuesAndVectors(cov, e1, e2, vec1, vec2);
	if(std::isfinite(e1))
	{
		fwhmMaj = sqrt(std::fabs(e1)) * sigmaToBeam;
		fwhmMin = sqrt(std::fabs(e2)) * sigmaToBeam;
		if(fwhmMaj < fwhmMin)
		{
			std::swap(fwhmMaj, fwhmMin);
			vec1[0] = vec2[0];
			vec1[1] = vec2[1];
		}
		positionAngle = atan2(vec1[1], vec1[0]);
	}
	else {
		fwhmMaj = sqrt(std::fabs(sxsx)) * sigmaToBeam;
		fwhmMin = sqrt(std::fabs(sxsx)) * sigmaToBeam;
		positionAngle = 0.0;
	}
}

void GaussianFitter::Fit2DGaussianCentred(const double* image, size_t width, size_t height, double beamEst, double& beamMaj, double& beamMin, double& beamPA, double boxScaleFactor, bool verbose)
{
	size_t prefSize = std::max<size_t>(std::ceil(boxScaleFactor), std::ceil(beamEst*boxScaleFactor));
	if(prefSize%2 != 0) ++prefSize;
	if(prefSize < width || prefSize < height)
	{
		size_t nIter = 0;
		bool boxWasLargeEnough;
		do {
			size_t boxWidth  = std::min(prefSize, width);
			size_t boxHeight = std::min(prefSize, height);
			if(verbose)
				std::cout << "Fit initial value:" << beamEst << "\n";
			fit2DGaussianCentredInBox(image, width, height, beamEst, beamMaj, beamMin, beamPA, boxWidth, boxHeight, verbose);
			if(verbose)
				std::cout << "Fit result:" << beamMaj << " x " << beamMin << " px, " << beamPA << " (box was " << boxWidth << " x " << boxHeight << ")\n";
			
			boxWasLargeEnough =
				(beamMaj*boxScaleFactor*0.8 < boxWidth || boxWidth>=width) &&
				(beamMaj*boxScaleFactor*0.8 < boxHeight || boxHeight>=height);
			if(!boxWasLargeEnough)
			{
				prefSize = std::max<size_t>(std::ceil(boxScaleFactor), std::ceil(beamMaj*boxScaleFactor));
				if(prefSize%2 != 0) ++prefSize;
				beamEst = std::max(beamMaj, beamEst);
			}
			++nIter;
		} while(!boxWasLargeEnough && nIter < 5);
	}
	else {
		if(verbose)
			std::cout << "Image is as large as the fitting box.\n";
		fit2DGaussianCentred(image, width, height, beamEst, beamMaj, beamMin, beamPA, verbose);
	}
}

void GaussianFitter::Fit2DCircularGaussianCentred(const double* image, size_t width, size_t height, double& beamSize, double boxScaleFactor)
{
	double initialValue = beamSize;
	size_t prefSize = std::max<size_t>(std::ceil(boxScaleFactor), std::ceil(beamSize*boxScaleFactor));
	if(prefSize%2 != 0) ++prefSize;
	if(prefSize < width || prefSize < height)
	{
		size_t boxWidth  = std::min(prefSize, width);
		size_t boxHeight = std::min(prefSize, height);
		size_t nIter = 0;
		bool boxWasLargeEnough;
		do {
			fit2DCircularGaussianCentredInBox(image, width, height, beamSize, boxWidth, boxHeight);
			
			boxWasLargeEnough =
				(beamSize*boxScaleFactor*0.8 < boxWidth || width>=boxWidth) &&
				(beamSize*boxScaleFactor*0.8 < boxHeight || height>=boxHeight);
			if(!boxWasLargeEnough)
			{
				prefSize = std::max<size_t>(std::ceil(boxScaleFactor), std::ceil(beamSize*boxScaleFactor));
				if(prefSize%2 != 0) ++prefSize;
				beamSize = std::max(initialValue, beamSize);
			}
			++nIter;
		} while(!boxWasLargeEnough && nIter < 5);
	}
	else {
		fit2DCircularGaussianCentred(image, width, height, beamSize);
	}
}

void GaussianFitter::Fit2DGaussianFull(const double* image, size_t width, size_t height, double& val, double& posX, double& posY, double& beamMaj, double& beamMin, double& beamPA, double* floorLevel)
{
	size_t prefSize = std::max<size_t>(10, std::ceil(beamMaj*10.0));
	if(prefSize%2 != 0) ++prefSize;
	if(prefSize < width || prefSize < height)
	{
		size_t xStart  = std::max<int>(0, int(round(posX)) - int(prefSize)/2);
		size_t xEnd    = std::min(width, size_t(round(posX)) + prefSize/2);
		size_t yStart  = std::max<int>(0, int(round(posY)) - int(prefSize)/2);
		size_t yEnd    = std::min(height, size_t(round(posY)) + prefSize/2);
		size_t nIter = 0;
		bool boxWasLargeEnough;
		do {
			fit2DGaussianWithAmplitudeInBox(image, width, height, val, posX, posY, beamMaj, beamMin, beamPA, floorLevel, xStart, xEnd, yStart, yEnd);
			
			size_t boxWidth = xEnd - xStart;
			size_t boxHeight = yEnd - yStart;
			boxWasLargeEnough =
				(beamMaj*4.0 < boxWidth || width>=boxWidth) &&
				(beamMaj*4.0 < boxHeight || height>=boxHeight);
			if(!boxWasLargeEnough)
			{
				prefSize = std::max<size_t>(10, std::ceil(beamMaj*10.0));
				if(prefSize%2 != 0) ++prefSize;
			}
			++nIter;
		} while(!boxWasLargeEnough && nIter < 5);
	}
	else {
		fit2DGaussianWithAmplitude(image, width, height, val, posX, posY, beamMaj, beamMin, beamPA, floorLevel);
	}
}

void GaussianFitter::fit2DGaussianCentredInBox(const double* image, size_t width, size_t height, double beamEst, double& beamMaj, double& beamMin, double& beamPA, size_t boxWidth, size_t boxHeight, bool verbose)
{
	size_t startX = (width-boxWidth)/2;
	size_t startY = (height-boxHeight)/2;
	ao::uvector<double> smallImage(boxWidth*boxHeight);
	for(size_t y=startY; y!=(height+boxHeight)/2; ++y)
	{
		memcpy(&smallImage[(y-startY)*boxWidth], &image[y*width + startX], sizeof(double)*boxWidth);
	}
	
	fit2DGaussianCentred(&smallImage[0], boxWidth, boxHeight, beamEst, beamMaj, beamMin, beamPA, verbose);
}

void GaussianFitter::fit2DCircularGaussianCentredInBox(const double* image, size_t width, size_t height, double& beamSize, size_t boxWidth, size_t boxHeight)
{
	size_t startX = (width-boxWidth)/2;
	size_t startY = (height-boxHeight)/2;
	ao::uvector<double> smallImage(boxWidth*boxHeight);
	for(size_t y=startY; y!=(height+boxHeight)/2; ++y)
	{
		memcpy(&smallImage[(y-startY)*boxWidth], &image[y*width + startX], sizeof(double)*boxWidth);
	}
	
	fit2DCircularGaussianCentred(&smallImage[0], boxWidth, boxHeight, beamSize);
}

void GaussianFitter::fit2DGaussianCentred(const double* image, size_t width, size_t height, double beamEst, double& beamMaj, double& beamMin, double& beamPA, bool verbose)
{
	_width = width;
	_height = height;
	_image = image;
	_scaleFactor = (width + height)/2;
	
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
	gsl_multifit_fdfsolver *solver = gsl_multifit_fdfsolver_alloc (T, _width*_height, 3);
	
	gsl_multifit_function_fdf fdf;
	fdf.f = &fitting_func_centered;
	fdf.df = &fitting_deriv_centered;
	fdf.fdf = &fitting_both_centered;
	fdf.n = _width*_height;
	fdf.p = 3;
	fdf.params = this;
	
	// Using the FWHM formula for a Gaussian:
	const long double sigmaToBeam = 2.0L * sqrtl(2.0L * logl(2.0L));
	double initialValsArray[3] = {
		beamEst/(_scaleFactor*double(sigmaToBeam)),
		beamEst/(_scaleFactor*double(sigmaToBeam)),
		0.0
	};
	gsl_vector_view initialVals = gsl_vector_view_array (initialValsArray, 3);
	gsl_multifit_fdfsolver_set (solver, &fdf, &initialVals.vector);

	int status;
	size_t iter = 0;
	do {
		if(verbose)
			std::cout << "Iteration " << iter << ": ";
		iter++;
		status = gsl_multifit_fdfsolver_iterate (solver);
		
		if(status)
			break;
		
		status = gsl_multifit_test_delta(solver->dx, solver->x, 1e-7, 1e-7);
		
	} while (status == GSL_CONTINUE && iter < 500);
	
	double
		sx = gsl_vector_get (solver->x, 0),
		sy = gsl_vector_get (solver->x, 1),
		beta = gsl_vector_get (solver->x, 2);
		
	gsl_multifit_fdfsolver_free(solver);
	
	ToAnglesAndFWHM(sx, sy, beta, beamMaj, beamMin, beamPA);
	beamMaj *= _scaleFactor;
	beamMin *= _scaleFactor;
}

void GaussianFitter::fit2DCircularGaussianCentred(const double* image, size_t width, size_t height, double& beamSize)
{
	_width = width;
	_height = height;
	_image = image;
	_scaleFactor = (width + height)/2;
	
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
	gsl_multifit_fdfsolver *solver = gsl_multifit_fdfsolver_alloc (T, _width*_height, 1);
	
	gsl_multifit_function_fdf fdf;
	fdf.f = &fitting_func_circular_centered;
	fdf.df = &fitting_deriv_circular_centered;
	fdf.fdf = &fitting_both_circular_centered;
	fdf.n = _width*_height;
	fdf.p = 1;
	fdf.params = this;
	
	// Using the FWHM formula for a Gaussian:
	const long double sigmaToBeam = 2.0L * sqrtl(2.0L * logl(2.0L));
	double initialValsArray[1] = { beamSize/(_scaleFactor*double(sigmaToBeam)) };
	gsl_vector_view initialVals = gsl_vector_view_array (initialValsArray, 1);
	gsl_multifit_fdfsolver_set (solver, &fdf, &initialVals.vector);

	int status;
	size_t iter = 0;
	do {
		iter++;
		status = gsl_multifit_fdfsolver_iterate (solver);
		
		if(status)
			break;
		
		status = gsl_multifit_test_delta(solver->dx, solver->x, 1e-7, 1e-7);
		
	} while (status == GSL_CONTINUE && iter < 500);
	
	double
		s = gsl_vector_get (solver->x, 0);			
	gsl_multifit_fdfsolver_free(solver);
	
	ToFWHM(s, beamSize);
	beamSize *= _scaleFactor;
}

/**
	* Fitting function for fit2DGaussianCentred(). Calculates the sum of the
	* squared errors(/residuals).
	*/
int GaussianFitter::fitting_func_centered(const gsl_vector *xvec, void *data, gsl_vector *f)
{
	GaussianFitter& fitter=*static_cast<GaussianFitter*>(data);
	double sx = gsl_vector_get(xvec, 0);
	double sy = gsl_vector_get(xvec, 1);
	double beta = gsl_vector_get(xvec, 2);
	const size_t width = fitter._width, height = fitter._height;
	int xMid = width/2, yMid = height/2;
	double scale = 1.0/fitter._scaleFactor;
	
	size_t dataIndex = 0;
	double errSum = 0.0;
	for(size_t yi=0; yi!=height; ++yi)
	{
		double y = (yi - yMid) * scale;
		for(size_t xi=0; xi!=width; ++xi)
		{
			double x = (xi - xMid) * scale;
			double e = err_centered(fitter._image[dataIndex], x, y, sx, sy, beta);
			errSum += e*e;
			gsl_vector_set(f, dataIndex, e);
			++dataIndex;
		}
	}
	//std::cout << "sx=" << sx << ", sy=" << sy << ", beta=" << beta << ", err=" << errSum << '\n';
	return GSL_SUCCESS;
}

int GaussianFitter::fitting_func_circular_centered(const gsl_vector *xvec, void *data, gsl_vector *f)
{
	GaussianFitter& fitter=*static_cast<GaussianFitter*>(data);
	double s = gsl_vector_get(xvec, 0);
	const size_t width = fitter._width, height = fitter._height;
	int xMid = width/2, yMid = height/2;
	double scale = 1.0/fitter._scaleFactor;
	
	size_t dataIndex = 0;
	double errSum = 0.0;
	for(size_t yi=0; yi!=height; ++yi)
	{
		double y = (yi - yMid) * scale;
		for(size_t xi=0; xi!=width; ++xi)
		{
			double x = (xi - xMid) * scale;
			double e = err_circular_centered(fitter._image[dataIndex], x, y, s);
			errSum += e*e;
			gsl_vector_set(f, dataIndex, e);
			++dataIndex;
		}
	}
	return GSL_SUCCESS;
}

/**
	* Derivative function belong with fit2DGaussianCentred().
	*/
int GaussianFitter::fitting_deriv_centered(const gsl_vector *xvec, void *data, gsl_matrix *J)
{
	GaussianFitter& fitter=*static_cast<GaussianFitter*>(data);
	const double sx = gsl_vector_get(xvec, 0);
	const double sy = gsl_vector_get(xvec, 1);
	const double beta = gsl_vector_get(xvec, 2);
	const size_t width = fitter._width, height = fitter._height;
	const int xMid = width/2, yMid = height/2;
	const double scale = 1.0 / fitter._scaleFactor;
	
	size_t dataIndex = 0;
	for(size_t yi=0; yi!=height; ++yi)
	{
		double y = (yi - yMid)*scale;
		for(size_t xi=0; xi!=width; ++xi)
		{
			double x = (xi - xMid)*scale;
			double expTerm = exp(-x*x/(2.0*sx*sx) + beta*x*y/(sx*sy) - y*y/(2.0*sy*sy));
			double dsx = (beta*x*y/(sx*sx*sy)+x*x/(sx*sx*sx)) * expTerm;
			double dsy = (beta*x*y/(sy*sy*sx)+y*y/(sy*sy*sy)) * expTerm;
			double dbeta = x*y/(sx*sy) * expTerm;
			gsl_matrix_set(J, dataIndex, 0, dsx);
			gsl_matrix_set(J, dataIndex, 1, dsy);
			gsl_matrix_set(J, dataIndex, 2, dbeta);
			++dataIndex;
		}
	}
	return GSL_SUCCESS;
}

int GaussianFitter::fitting_deriv_circular_centered(const gsl_vector *xvec, void *data, gsl_matrix *J)
{
	GaussianFitter& fitter=*static_cast<GaussianFitter*>(data);
	const double s = gsl_vector_get(xvec, 0);
	const size_t width = fitter._width, height = fitter._height;
	const int xMid = width/2, yMid = height/2;
	const double scale = 1.0 / fitter._scaleFactor;
	
	size_t dataIndex = 0;
	for(size_t yi=0; yi!=height; ++yi)
	{
		double y = (yi - yMid)*scale;
		for(size_t xi=0; xi!=width; ++xi)
		{
			double x = (xi - xMid)*scale;
			double expTerm = exp((-x*x - y*y)/(2.0*s*s));
			// derivative of exp((-x*x - y*y)/(2.0*s*s)) to s
			// = (-x*x - y*y)/2.0*-2/(s*s*s)
			// = (-x*x - y*y)/(-s*s*s)
			// = (x*x + y*y)/(s*s*s)
			double ds = ((x*x + y*y)/(s*s*s)) * expTerm;
			gsl_matrix_set(J, dataIndex, 0, ds);
			++dataIndex;
		}
	}
	return GSL_SUCCESS;
}

void GaussianFitter::fit2DGaussianWithAmplitudeInBox(const double* image, size_t width, size_t /*height*/, double& val, double& posX, double& posY, double& beamMaj, double& beamMin, double& beamPA, double* floorLevel, size_t xStart, size_t xEnd, size_t yStart, size_t yEnd)
{
	size_t boxWidth = xEnd - xStart;
	size_t boxHeight = yEnd - yStart;
	ao::uvector<double> smallImage(boxWidth*boxHeight);
	for(size_t y=yStart; y!=yEnd; ++y)
	{
		memcpy(&smallImage[(y-yStart)*boxWidth], &image[y*width + xStart], sizeof(double)*boxWidth);
	}
	
	posX -= xStart;
	posY -= yStart;
	fit2DGaussianWithAmplitude(&smallImage[0], boxWidth, boxHeight, val, posX, posY, beamMaj, beamMin, beamPA, floorLevel);
	posX += xStart;
	posY += yStart;
}

/**
	* Fits the position, size and amplitude of a Gaussian. If floorLevel is not
	* a nullptr, the floor (background level, or zero level) is fitted too. 
	*/
void GaussianFitter::fit2DGaussianWithAmplitude(const double* image, size_t width, size_t height, double& val, double& posX, double& posY, double& beamMaj, double& beamMin, double& beamPA, double* floorLevel)
{
	_width = width;
	_height = height;
	_image = image;
	_scaleFactor = (width + height)/2;
	
	if(floorLevel == 0)
		fit2DGaussianWithAmplitude(val, posX, posY, beamMaj, beamMin, beamPA);
	else
		fit2DGaussianWithAmplitudeWithFloor(val, posX, posY, beamMaj, beamMin, beamPA, *floorLevel);
}

void GaussianFitter::fit2DGaussianWithAmplitude(double& val, double& posX, double& posY, double& beamMaj, double& beamMin, double& beamPA)
{
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
	gsl_multifit_fdfsolver *solver = gsl_multifit_fdfsolver_alloc (T, _width*_height, 6);
	
	gsl_multifit_function_fdf fdf;
	fdf.f = &fitting_func_with_amplitude;
	fdf.df = &fitting_deriv_with_amplitude;
	fdf.fdf = &fitting_both_with_amplitude;
	fdf.n = _width*_height;
	fdf.p = 6;
	fdf.params = this;
	
	// Using the FWHM formula for a Gaussian:
	const long double sigmaToBeam = 2.0L * sqrtl(2.0L * logl(2.0L));
	_xInit = -(posX-_width/2)/_scaleFactor;
	_yInit = -(posY-_height/2)/_scaleFactor;
	double initialValsArray[6] = {
		val,
		_xInit,
		_yInit,
		beamMaj/(_scaleFactor*double(sigmaToBeam)),
		beamMaj/(_scaleFactor*double(sigmaToBeam)),
		0.0
	};
	gsl_vector_view initialVals = gsl_vector_view_array (initialValsArray, 6);
	gsl_multifit_fdfsolver_set (solver, &fdf, &initialVals.vector);

	int status;
	size_t iter = 0;
	do {
		iter++;
		status = gsl_multifit_fdfsolver_iterate (solver);
		
		if(status)
			break;
		
		status = gsl_multifit_test_delta(solver->dx, solver->x, 1e-7, 1e-7);
		
	} while (status == GSL_CONTINUE && iter < 500);
	
	val = gsl_vector_get (solver->x, 0);
	posX = -1.0*gsl_vector_get (solver->x, 1)*_scaleFactor + _width/2;
	posY = -1.0*gsl_vector_get (solver->x, 2)*_scaleFactor + _height/2;
	double
		sx = gsl_vector_get (solver->x, 3),
		sy = gsl_vector_get (solver->x, 4),
		beta = gsl_vector_get (solver->x, 5);
		
	gsl_multifit_fdfsolver_free(solver);
	
	ToAnglesAndFWHM(sx, sy, beta, beamMaj, beamMin, beamPA);
	beamMaj *= _scaleFactor;
	beamMin *= _scaleFactor;
}

void GaussianFitter::fit2DGaussianWithAmplitudeWithFloor(double& val, double& posX, double& posY, double& beamMaj, double& beamMin, double& beamPA, double& floorLevel)
{
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
	gsl_multifit_fdfsolver *solver = gsl_multifit_fdfsolver_alloc (T, _width*_height, 7);
	
	gsl_multifit_function_fdf fdf;
	fdf.f = &fitting_func_with_amplitude_and_floor;
	fdf.df = &fitting_deriv_with_amplitude_and_floor;
	fdf.fdf = &fitting_both_with_amplitude_and_floor;
	fdf.n = _width*_height;
	fdf.p = 7;
	fdf.params = this;
	
	// Using the FWHM formula for a Gaussian:
	const long double sigmaToBeam = 2.0L * sqrtl(2.0L * logl(2.0L));
	_xInit = -(posX-_width/2)/_scaleFactor;
	_yInit = -(posY-_height/2)/_scaleFactor;
	double initialValsArray[7] = {
		val,
		_xInit,
		_yInit,
		beamMaj/(_scaleFactor*double(sigmaToBeam)),
		beamMaj/(_scaleFactor*double(sigmaToBeam)),
		0.0,
		0.0
	};
	gsl_vector_view initialVals = gsl_vector_view_array (initialValsArray, 7);
	gsl_multifit_fdfsolver_set (solver, &fdf, &initialVals.vector);

	int status;
	size_t iter = 0;
	do {
		iter++;
		status = gsl_multifit_fdfsolver_iterate (solver);
		
		if(status)
			break;
		
		status = gsl_multifit_test_delta(solver->dx, solver->x, 1e-7, 1e-7);
		
	} while (status == GSL_CONTINUE && iter < 500);
	
	val = gsl_vector_get (solver->x, 0);
	posX = -1.0*gsl_vector_get (solver->x, 1)*_scaleFactor + _width/2;
	posY = -1.0*gsl_vector_get (solver->x, 2)*_scaleFactor + _height/2;
	double
		sx = gsl_vector_get (solver->x, 3),
		sy = gsl_vector_get (solver->x, 4),
		beta = gsl_vector_get (solver->x, 5);
	floorLevel = gsl_vector_get(solver->x, 6);
		
	gsl_multifit_fdfsolver_free(solver);
	
	ToAnglesAndFWHM(sx, sy, beta, beamMaj, beamMin, beamPA);
	beamMaj *= _scaleFactor;
	beamMin *= _scaleFactor;
}

int GaussianFitter::fitting_func_with_amplitude(const gsl_vector *xvec, void *data, gsl_vector *f)
{
	GaussianFitter& fitter=*static_cast<GaussianFitter*>(data);
	double
		v = gsl_vector_get (xvec, 0),
		xc = gsl_vector_get (xvec, 1),
		yc = gsl_vector_get (xvec, 2),
		sx = gsl_vector_get (xvec, 3),
		sy = gsl_vector_get (xvec, 4),
		beta = gsl_vector_get (xvec, 5);
	const size_t width = fitter._width, height = fitter._height;
	int xMid = width/2, yMid = height/2;
	double scale = 1.0/fitter._scaleFactor;
	
	size_t dataIndex = 0;
	double errSum = 0.0;
	for(int yi=0; yi!=int(height); ++yi)
	{
		double yS = yc + (yi - yMid) * scale;
		for(int xi=0; xi!=int(width); ++xi)
		{
			double xS = xc + (xi - xMid) * scale;
			double e = err_full(fitter._image[dataIndex], v, xS, yS, sx, sy, beta);
			errSum += e*e;
			gsl_vector_set(f, dataIndex, e);
			++dataIndex;
		}
	}
	//std::cout << "v=" << v << ", x=" << xc << ", y=" << yc << ", sx=" << sx << ", sy=" << sy << ", beta=" << beta << ", err=" << errSum << '\n';
	return GSL_SUCCESS;
}

int GaussianFitter::fitting_deriv_with_amplitude(const gsl_vector *xvec, void *data, gsl_matrix *J)
{
	GaussianFitter& fitter=*static_cast<GaussianFitter*>(data);
	const double scale = 1.0 / fitter._scaleFactor;
	double
		v = gsl_vector_get (xvec, 0),
		xc = gsl_vector_get (xvec, 1),
		yc = gsl_vector_get (xvec, 2),
		sx = gsl_vector_get (xvec, 3),
		sy = gsl_vector_get (xvec, 4),
		beta = gsl_vector_get (xvec, 5);
	if(fitter._posConstrained!=0.0 && (std::fabs(xc-fitter._xInit)>fitter._posConstrained*scale || std::fabs(yc-fitter._yInit)>fitter._posConstrained*scale))
	{
		std::cout << "GSL_EDOM\n";
		return GSL_EDOM;
	}
	const size_t width = fitter._width, height = fitter._height;
	int xMid = width/2, yMid = height/2;
	
	size_t dataIndex = 0;
	for(int yi=0; yi!=int(height); ++yi)
	{
		double y = yc + (yi - yMid)*scale;
		for(int xi=0; xi!=int(width); ++xi)
		{
			// TODO I need to go over the signs -- ds, dy, dsx, dsy in particular
			double x = xc + (xi - xMid)*scale;
			double expTerm = exp(-x*x/(2.0*sx*sx) + beta*x*y/(sx*sy) - y*y/(2.0*sy*sy));
			double dv = expTerm;
			expTerm *= v;
			double dx = (-beta*y/(sx*sy) - x/(sx*sx)) * expTerm;
			double dy = (-beta*x/(sy*sx) - y/(sy*sy)) * expTerm;
			double dsx = (beta*x*y/(sx*sx*sy)+x*x/(sx*sx*sx)) * expTerm;
			double dsy = (beta*x*y/(sy*sy*sx)+y*y/(sy*sy*sy)) * expTerm;
			double dbeta = x*y/(sx*sy) * expTerm;
			gsl_matrix_set(J, dataIndex, 0, dv);
			gsl_matrix_set(J, dataIndex, 1, dx);
			gsl_matrix_set(J, dataIndex, 2, dy);
			gsl_matrix_set(J, dataIndex, 3, dsx);
			gsl_matrix_set(J, dataIndex, 4, dsy);
			gsl_matrix_set(J, dataIndex, 5, dbeta);
			++dataIndex;
		}
	}
	//std::cout << "diff: v=" << v << ", x=" << xc << ", y=" << yc << ", sx=" << sx << ", sy=" << sy << ", beta=" << beta << '\n';
	return GSL_SUCCESS;
}

int GaussianFitter::fitting_func_with_amplitude_and_floor(const gsl_vector *xvec, void *data, gsl_vector *f)
{
	GaussianFitter& fitter=*static_cast<GaussianFitter*>(data);
	const double scale = 1.0/fitter._scaleFactor;
	double
		v = gsl_vector_get (xvec, 0),
		xc = gsl_vector_get (xvec, 1),
		yc = gsl_vector_get (xvec, 2),
		sx = gsl_vector_get (xvec, 3),
		sy = gsl_vector_get (xvec, 4),
		beta = gsl_vector_get (xvec, 5),
		fl = gsl_vector_get (xvec, 6);
	if(fitter._posConstrained!=0.0 && (std::fabs(xc-fitter._xInit)>fitter._posConstrained*scale || std::fabs(yc-fitter._yInit)>fitter._posConstrained*scale))
		return GSL_EDOM;
	const size_t width = fitter._width, height = fitter._height;
	int xMid = width/2, yMid = height/2;
	
	size_t dataIndex = 0;
	double errSum = 0.0;
	for(int yi=0; yi!=int(height); ++yi)
	{
		double yS = yc + (yi - yMid) * scale;
		for(int xi=0; xi!=int(width); ++xi)
		{
			double xS = xc + (xi - xMid) * scale;
			double e = err_full(fitter._image[dataIndex], v, xS, yS, sx, sy, beta) + fl;
			errSum += e*e;
			gsl_vector_set(f, dataIndex, e);
			++dataIndex;
		}
	}
	//std::cout << "v=" << v << ", x=" << xc << ", y=" << yc << ", sx=" << sx << ", sy=" << sy << ", beta=" << beta << ", err=" << errSum << '\n';
	return GSL_SUCCESS;
}

int GaussianFitter::fitting_deriv_with_amplitude_and_floor(const gsl_vector *xvec, void *data, gsl_matrix *J)
{
	GaussianFitter& fitter=*static_cast<GaussianFitter*>(data);
	double
		v = gsl_vector_get (xvec, 0),
		xc = gsl_vector_get (xvec, 1),
		yc = gsl_vector_get (xvec, 2),
		sx = gsl_vector_get (xvec, 3),
		sy = gsl_vector_get (xvec, 4),
		beta = gsl_vector_get (xvec, 5);
	const size_t width = fitter._width, height = fitter._height;
	int xMid = width/2, yMid = height/2;
	double scale = 1.0 / fitter._scaleFactor;
	
	size_t dataIndex = 0;
	for(int yi=0; yi!=int(height); ++yi)
	{
		double y = yc + (yi - yMid)*scale;
		for(int xi=0; xi!=int(width); ++xi)
		{
			double x = xc + (xi - xMid)*scale;
			double expTerm = exp(-x*x/(2.0*sx*sx) + beta*x*y/(sx*sy) - y*y/(2.0*sy*sy));
			double dv = expTerm;
			expTerm *= v;
			double dx = (-beta*y/(sx*sy) - x/(sx*sx)) * expTerm;
			double dy = (-beta*x/(sy*sx) - y/(sy*sy)) * expTerm;
			double dsx = (beta*x*y/(sx*sx*sy)+x*x/(sx*sx*sx)) * expTerm;
			double dsy = (beta*x*y/(sy*sy*sx)+y*y/(sy*sy*sy)) * expTerm;
			double dbeta = x*y/(sx*sy) * expTerm;
			double dfl = 1.0;
			gsl_matrix_set(J, dataIndex, 0, dv);
			gsl_matrix_set(J, dataIndex, 1, dx);
			gsl_matrix_set(J, dataIndex, 2, dy);
			gsl_matrix_set(J, dataIndex, 3, dsx);
			gsl_matrix_set(J, dataIndex, 4, dsy);
			gsl_matrix_set(J, dataIndex, 5, dbeta);
			gsl_matrix_set(J, dataIndex, 6, dfl);
			++dataIndex;
		}
	}
	//std::cout << "diff: v=" << v << ", x=" << xc << ", y=" << yc << ", sx=" << sx << ", sy=" << sy << ", beta=" << beta << '\n';
	return GSL_SUCCESS;
}
