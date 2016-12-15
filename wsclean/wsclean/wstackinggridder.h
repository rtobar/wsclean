//
// This file was written by André Offringa and is published
// under the GPL 3 license. 
//

#ifndef WSTACKING_GRIDDER_H
#define WSTACKING_GRIDDER_H

#ifndef AVOID_CASACORE
#include "../multibanddata.h"
#endif

#include "gridmodeenum.h"

// Forward declare boost::mutex to avoid header inclusion
namespace boost
{
	class mutex;
}

#include <cmath>
#include <cstring>
#include <complex>
#include <vector>
#include <stack>

class ImageBufferAllocator;

/**
 * This class grids and/or samples visibilities to/from UV space.
 * It also executes the FFT(s) required and performs the w-term correction.
 * It does this with the 'w-stacking' method (see
 * [Offringa et al., 2014](http://arxiv.org/abs/1407.1943)).
 * 
 * This class is quite generic in that it is agnostic to the underlying data format 
 * of the visibilities, and it is written to be of high performance. It is used
 * by the @ref WSMSGridder class to perform inversion and prediction. There, large
 * imaging operations are split up in multiple 'passes' over the data if necessary,
 * and in each pass a few w-layers are gridded/sampled. However, simpler applications
 * might chose not to support this when very large images/large w-terms are not a
 * concern, in which case things become easier.
 * 
 * This class can do both the calculation from visibilities to image, as well as going
 * from (model) image to visibilities. The first is referred to as "Inversion" in this
 * class, and the latter is called "Prediction".
 * 
 * The general call sequence for inversion (/imaging) is as follows:
 * 
 * - Construct an instance with @ref WStackingGridder()
 * - Set settings if necessary (grid mode, denormal phase centre, complex...)
 * - Call @ref PrepareWLayers();
 * - For each pass if multiple passes are necessary (or once otherwise) :
 *   - Call @ref StartInversionPass();
 *   - Add all samples with @ref AddDataSample();
 *   - Call @ref FinishInversionPass();
 * - Finally, call @ref FinalizeImage();
 * - Now, @ref RealImage() and optionally @ref ImaginaryImage() will return the
 *   image(s).
 * 
 * Alternatively, @ref AddData() can be used instead of @ref AddDataSample, to grid
 * several samples that only differ in frequency. To use @ref AddData(), it is necessary to
 * call @ref PrepareBand() first.
 * 
 * For prediction, the sequence is similar:
 * 
 * - Construct an instance with @ref WStackingGridder()
 * - Set settings if necessary
 * - Call @ref PrepareWLayers();
 * - For each pass if multiple passes are necessary (or once otherwise) :
 *   - Call @ref InitializePrediction();
 *   - Call @ref StartPredictionPass();
 *   - Get the predicted visibilities with @ref SampleDataSample();
 * 
 * Similar to Inversion, all values for a band can be predicted at once by calling
 * @ref SampleData() instead of @ref SampleDataSample() when the band has been
 * set beforehand with @ref PrepareBand().
 * 
 * Prediction does not require any finalisation calls.
 * 
 * @author André Offringa
 * @date 2013 (first version)
 * @sa [WSClean: an implementation of a fast, generic wide-field imager for radio astronomy](http://arxiv.org/abs/1407.1943)
 */
class WStackingGridder
{
	/**
	* @example wspredictionexample.cpp
	* This example demonstrates how the WStackingGridder can be used in a minimalistic
	* way to predict visibilities.
	*/
	public:
		/** Construct a new gridder with given settings.
		 * @param width The width of the image in pixels
		 * @param height The height of the image in pixels.
		 * @param pixelSizeX The angular width of a pixel in radians.
		 * @param pixelSizeY The angular height of a pixel in radians.
		 * @param fftThreadCount The number of threads used for FFTs.
		 *   See @ref NFFTThreads() for more info.
		 * @param allocator An ImageBufferAllocator that is used for allocating images.
		 *   See @ref ImageBufferAllocator for the rational for using this.
		 * @param kernelSize The full width and height in pixels of the kernel. Should be odd.
		 *   Large kernels cause smaller aliasing artefacts but are more expensive.
		 *   The default of 7 is a reasonable compromise between effectiveness and speed.
		 * @param overSamplingFactor The number of different horizontal and vertical kernels
		 *   that are precalculated at different rational positions. Larger is more accurate
		 *   but requires more memory and becomes slower, probably mainly due to cache misses.
		 */
		WStackingGridder(size_t width, size_t height, double pixelSizeX, double pixelSizeY, size_t fftThreadCount, ImageBufferAllocator* allocator, size_t kernelSize = 7, size_t overSamplingFactor = 63);
		
		WStackingGridder(const WStackingGridder&) = delete;
		
		WStackingGridder& operator=(const WStackingGridder&) = delete;
		
		/** De-allocate imagebuffers with the allocator and perform other clean-up. */
		~WStackingGridder();
		
		/** Initialize the inversion/prediction stage with the given number of w-layers.
		 * The number of w-layers sets the accuracy of w-term correction and can normally
		 * be calculated from the image dimensions and maximum w-term. 
		 * 
		 * After this call, @ref NPasses() can be called to acquire the number of required
		 * passes.
		 * 
		 * Note that for non-complex images (the default), @p minW and @p maxW can be
		 * set to the minimum absolute w-value, which results
		 * in finer w-gridding and thus better w-term correction.
		 * 
		 * The @p minW and @p maxW values need to be
		 * calculated before gridding, and this normally requires an extra iteration through
		 * the w-values of the data before the data can be gridded.
		 * 
		 * For imaging without w-term correction, the number of layers can be set to one
		 * (note that one w-layer disables w-term correcting -- do not set it to zero, as
		 * one uv-grid is always required). In that case, it is best to set @p minW to
		 * -@p maxW.
		 * 
		 * In cases where nwlayer > 1, the w-value @c w of all values should satisfy
		 * @p minW < abs(@c w) < @p maxW. When @p nWLayers == 1, this is not required.
		 * 
		 * @param nWLayers Number of uv grids at different w-values, should be >= 1.
		 * @param maxMem Allowed memory in bytes. The gridder will try to set the number
		 * of passes such that this value is not exceeded. Note that this is approximate.
		 * @param minW The smallest w-value to be inverted/predicted.
		 * @param maxW The largest w-value to be inverted/predicted.
		 */
		void PrepareWLayers(size_t nWLayers, double maxMem, double minW, double maxW);
		
#ifndef AVOID_CASACORE
		/**
		 * Initialize the inversion/prediction stage with a given band. This is
		 * required for calling methods that take a dataDescId, specifically @ref AddData() and
		 * @ref SampleData(). This call can be avoided
		 * by using the @ref AddDataSample() and @ref SampleDataSample methods instead.
		 * @param bandData The information about the bands. A MultiBandData links a
		 * dataDescId to a set of contiguous frequencies. This corresponds with the
		 * DATA_DESC_ID field in meaurement sets.
		 */
		void PrepareBand(const MultiBandData &bandData)
		{
			_bandData = bandData;
		}
#endif // AVOID_CASACORE
		
		/**
		 * Calculate on which layer this w-value belongs. Only valid once
		 * @ref PrepareWLayers() has been called.
		 * @param wInLambda The w-value, in units of number of wavelengths.
		 * @returns W-layer index.
		 * @see @ref LayerToW()
		 */
		size_t WToLayer(double wInLambda) const
		{
			if(_nWLayers == 1)
				return 0;
			else {
				if(_isComplex)
					return size_t(round((wInLambda + _maxW) * (_nWLayers-1) / (_maxW + _maxW)));
				else
					return size_t(round((fabs(wInLambda) - _minW) * (_nWLayers-1) / (_maxW - _minW)));
			}
		}
		
		/**
		 * Calculate the central w-value for the given w-layer index. This is the inverse
		 * of @ref WToLayer().
		 * @param layer W-layer index
		 * @returns Central w-value for specified w-layer in units of number of wavelengths.
		 */
		double LayerToW(size_t layer) const
		{
			if(_nWLayers == 1)
				return 0.0;
			else {
				if(_isComplex)
					// In this case, replace minW with -maxW :
					return layer * (_maxW + _maxW) / (_nWLayers-1) - _maxW;
				else
					return layer * (_maxW - _minW) / (_nWLayers-1) + _minW;
			}
		}
		
		/**
		 * Number of w-layers, as specified in @ref PrepareWLayers().
		 * @returns Number of w-layers.
		 */
		size_t NWLayers() const { return _nWLayers; }
		
		/**
		 * Determine whether specified w-value should be gridded in this pass.
		 * This method can only be called after calling @ref StartInversionPass()
		 * or @ref StartPredictionPass().
		 * @param wInLambda W-value in units of number of wavelengths.
		 * @returns true if the w-value is gridded on one of the w-layers that are processed in this pass.
		 */
		bool IsInLayerRange(double wInLambda) const
		{
			size_t layer = WToLayer(wInLambda);
			return layer >= layerRangeStart(_curLayerRangeIndex) && layer < layerRangeStart(_curLayerRangeIndex+1);
		}
		
		
		/**
		 * Determine whether any samples within the specified w-value range
		 * should be gridded in this pass. This can for example be used to
		 * optimize which samples to read from disc; when this function returns
		 * false for a given timestep/baseline, none of the samples are required
		 * in this pass.
		 * This method can only be called after calling @ref StartInversionPass()
		 * or @ref StartPredictionPass(). 
		 * @param wStart W-value of range start in units of number of wavelengths.
		 * @param wEnd W-value of range end in units of number of wavelengths.
		 * @returns true if any of the w-values in the given range are processed in this pass.
		 */
		bool IsInLayerRange(double wStart, double wEnd) const
		{
			size_t
				rangeStart = layerRangeStart(_curLayerRangeIndex),
				rangeEnd = layerRangeStart(_curLayerRangeIndex+1),
				l1 = WToLayer(wStart);
			if(l1 >= rangeStart && l1 < rangeEnd)
				return true;
			size_t l2 = WToLayer(wEnd);
			return ((l2 >= rangeStart && l2 < rangeEnd) // l2 is within the range
			  || (l1 < rangeStart && l2 >= rangeEnd)  // l1 is before, l2 is after range
				|| (l2 < rangeStart && l1 >= rangeEnd) // l2 is before, l1 is after range
			);
		}
		
		/**
		 * Number of passes that are required when not all the w-layers fit in memory at once.
		 * Valid once @ref PrepareWLayers() has been called.
		 * @return Number of passes.
		 */
		size_t NPasses() const
		{
			return _nPasses;
		}
		
#ifndef AVOID_CASACORE
		/**
		 * Grid an array of data values for inversion. The data values should have the same
		 * uvw-values in meters, i.e., they can only differ in frequency. The dataDescId
		 * specifies the frequencies of the array of data. This method requires that the
		 * channel frequencies have been specified beforehand, by calling @ref PrepareBand().
		 * 
		 * This function will internally call @ref AddDataSample() for each visibility, hence
		 * this function is just for convenience, but is not faster than individual calls to
		 * @ref AddDataSample().
		 * 
		 * @param data Array of samples for different channels. The size of this array is given
		 * by the band referred to by dataDescId.
		 * @param dataDescId ID that specifies which band this data is for.
		 * @param uInM U value of UVW coordinate, in meters.
		 * @param vInM V value of UVW coordinate, in meters.
		 * @param wInM W value of UVW coordinate, in meters.
		 */
		void AddData(const std::complex<float>* data, size_t dataDescId, double uInM, double vInM, double wInM);
#endif
		
		/**
		 * Grid a single visibility value for inversion. This method does not require that the
		 * channel frequencies have been specified beforehand.
		 * @param sample Visibility value.
		 * @param uInLambda U value of UVW coordinate, in number of wavelengths.
		 * @param vInLambda V value of UVW coordinate, in number of wavelengths.
		 * @param wInLambda W value of UVW coordinate, in number of wavelengths.
		 */
		void AddDataSample(std::complex<float> sample, double uInLambda, double vInLambda, double wInLambda);
		
		/**
		 * Initialize a new inversion gridding pass. @ref PrepareWLayers() should have been called beforehand.
		 * Each call to @ref StartInversionPass() should be followed by a call to
		 * @ref FinishInversionPass().
		 * @param passIndex The zero-indexed index of this pass, 0 <= @p passIndex < @ref NPasses().
		 */
		void StartInversionPass(size_t passIndex);
		
		/**
		 * Finish an inversion gridding pass. This will perform the Fourier transforms of the currently gridded
		 * w-layers, and add each gridded layer to the final image including w-term corrections.
		 * Therefore, it can take time.
		 * @sa @ref StartInversionPass().
		 */
		void FinishInversionPass();
		
		/**
		 * Finalize inversion once all passes are performed.
		 * @param multiplicationFactor Apply this factor to all pixels. This can be used to normalize
		 * the image for the weighting scheme.
		 * @param correctFFTFactor For normal imaging this should be false. It can be set to true
		 * when no correcting for weighting is performed, but the FFT factor should be taken out.
		 * WSClean uses this for certain values of "-visibility-weighting-mode".
		 */
		void FinalizeImage(double multiplicationFactor, bool correctFFTFactor);
		
		/**
		 * Initialize gridder for prediction and specify image to predict for.
		 * @ref PrepareWLayers() should be called before this method.
		 * Note that this method needs to be called with the same image before
		 * each pass, i.e., before each call to @ref StartPredictionPass().
		 * 
		 * This method is used for prediction of non-complex (IsComplex()==false)
		 * images. Use @ref InitializePrediction(const double*, const double*) for
		 * complex prediction -- see @ref SetIsComplex() for more info.
		 * 
		 * @param image The model image that is to be predicted for. This is an
		 * array of width * height size, index by (x + width*y).
		 */
		void InitializePrediction(const double *image)
		{
			initializePrediction(image, _imageData);
		}
		
		/**
		 * Complex alternative for @ref InitializePrediction(const double*).
		 * 
		 * This method should be used for complex images -- see @ref SetIsComplex() for
		 * more info.
		 * @param real Array of width * height giving the real (model) image values.
		 * @param imaginary Array of width * height giving the imaginary (model) image values.
		 */
		void InitializePrediction(const double *real, const double *imaginary)
		{
			initializePrediction(real, _imageData);
			initializePrediction(imaginary, _imageDataImaginary);
		}

		/**
		 * Start a new prediction pass. One of the @ref InitializePrediction() methods should
		 * be called before each call to this method. This method will perform the fast Fourier
		 * transforms necessary for this pass, so it can take time.
		 * @param passIndex Zero-indexed index of prediction pass, 0 <= @p passIndex < @ref NPasses().
		 */
		void StartPredictionPass(size_t passIndex);
		
#ifndef AVOID_CASACORE
		/**
		 * Predict the values for all channels in a band. This is a convenience function
		 * that will call @ref SampleDataSample() for all data values in the band specified by
		 * @p dataDescId.
		 * 
		 * If a particular value can not be predicted during this pass (because its w-value
		 * belongs to a w-layer of a different pass), its value will be set to NaN.
		 * A call to @ref StartPredictionPass() should have been made before calling this method.
		 * 
		 * @param data Array of samples for different channels. The size of this array is given
		 * by the band referred to by dataDescId.
		 * @param dataDescId ID that specifies which band this data is for.
		 * @param uInM U value of UVW coordinate, in meters.
		 * @param vInM V value of UVW coordinate, in meters.
		 * @param wInM W value of UVW coordinate, in meters.
		 */
		void SampleData(std::complex<float>* data, size_t dataDescId, double uInM, double vInM, double wInM);
#endif
		
		/**
		 * Predict the value of a single visibility with double precision.
		 * If the visibility can not be predicted,
		 * because its w-value belongs to a w-layer that is not processed during this pass,
		 * it will be given a value of NaN.
		 * A call to @ref StartPredictionPass() should have been made before calling this method.
		 * @param value Will be set to the predicted visibility value.
		 * @param uInLambda U value of UVW coordinate, in number of wavelengths.
		 * @param vInLambda V value of UVW coordinate, in number of wavelengths.
		 * @param wInLambda W value of UVW coordinate, in number of wavelengths.
		 */
		void SampleDataSample(std::complex<double>& value, double uInLambda, double vInLambda, double wInLambda);
		
		/**
		 * Predict the value of a single visibility. If the visibility can not be predicted,
		 * because its w-value belongs to a w-layer that is not processed during this pass,
		 * it will be given a value of NaN.
		 * A call to @ref StartPredictionPass() should have been made before calling this method.
		 * @param value Will be set to the predicted visibility value.
		 * @param uInLambda U value of UVW coordinate, in number of wavelengths.
		 * @param vInLambda V value of UVW coordinate, in number of wavelengths.
		 * @param wInLambda W value of UVW coordinate, in number of wavelengths.
		 */
		void SampleDataSample(std::complex<float>& value, double uInLambda, double vInLambda, double wInLambda)
		{
			std::complex<double> doubleValue;
			SampleDataSample(doubleValue, uInLambda, vInLambda, wInLambda);
			value = doubleValue;
		}
		
		/**
		 * Get the image result of inversion. This is an array of size width x height, and can be
		 * indexed with [x + y*width]. It is allowed to change this image, e.g. set the horizon
		 * to zero before saving to fits. This call is only valid once @ref FinalizeImage() has
		 * been called.
		 * 
		 * If a complex image is produced, this image returns the real part. The imaginary part can
		 * be acquired with @ref ImaginaryImage().
		 */
		double *RealImage() { return _imageData[0]; }
		
		/**
		 * Get the imaginary part of a complex image after inversion. Otherwise similar to
		 * @ref RealImage().
		 */
		double *ImaginaryImage() { return _imageDataImaginary[0]; }
		
		/**
		 * Get the number of threads used when performing the FFTs. The w-layers are divided over
		 * the threads such that each thread calculates approximately the same number of w-layers.
		 * Therefore, optimally this value equals the number of (free) CPU cores available.
		 * If the number of w-layers is smaller than the number of threads, some threads will do
		 * nothing. Therefore, when not performing any w-term correction with @ref NWLayers() == 1,
		 * this value has no effect.
		 * 
		 * When @ref NWLayers() == 1 , it is still possible to perform multi-threading by setting
		 * up fftw to do so. The @ref FFTWMultiThreadEnabler class can be used to setup fftw
		 * properly to do so.
		 * 
		 * @returns The number of threads used for the FFTs.
		 */
		size_t NFFTThreads() const { return _nFFTThreads; }
		
		/**
		 * Set the number of threads used to perform the FFTs.
		 * @see @ref NFFTThreads() for an explanation.
		 * @param nfftThreads Number of threads used.
		 */
		void SetNFFTThreads(size_t nfftThreads) { _nFFTThreads = nfftThreads; }
		
		/**
		 * Get the kernel function and its interpolation method used for gridding.
		 * @returns The currently selected gridding mode.
		 */
		enum GridModeEnum GridMode() const { return _gridMode; }
		
		/**
		 * Set the kernel function and its interpolation method used for gridding.
		 * @param mode The new gridding mode.
		 */
		void SetGridMode(enum GridModeEnum mode) {
			if(mode != _gridMode)
			{
				_gridMode = mode;
				if(_gridMode != NearestNeighbourGridding)
					makeKernels();
			}
		}
		
		/**
		 * Whether the image produced by inversion or used by prediction is complex.
		 * In particular, cross-polarized images like XY and YX have complex values,
		 * because the UV grid is no longer Hermitian symmetric.
		 * @returns Whether the gridder is setup for complex images.
		 */
		bool IsComplex() const { return _isComplex; }
		
		/**
		 * Setup the gridder for real or complex images. See @ref IsComplex().
		 * @param isComplex Whether the gridder should produce / predict from
		 * complex images.
		 */
		void SetIsComplex(bool isComplex) { _isComplex = isComplex; }
		
		//void SetImageConjugatePart(bool imageConjugatePart) { _imageConjugatePart = imageConjugatePart; }
		
		/**
		 * Setup the gridder for images with a shifted phase centre. When dl or dm is
		 * non-zero, the centre of the involved images is shifted by the given amount.
		 * This allows inversion/prediction of images that are centred away from the
		 * phase centre, while the projection of the image is kept the same. The benefit
		 * of this can be that the phase centre can be set to zenith, thereby minimizing
		 * w-values, while the actually imaged field is somewhere else, thereby making
		 * it possible to handle non-zenith images with small image sizes and small w-values.
		 * A more extended explanation is given in the WSClean paper by
		 * [Offringa et al. (2014)](http://arxiv.org/abs/1407.1943).
		 * @param dl shift of image centre in 'l' (East) direction
		 * @param dm shift of image centre in 'm' (North) direction
		 */
		void SetDenormalPhaseCentre(double dl, double dm) { _phaseCentreDL = dl; _phaseCentreDM = dm; }
		
		/**
		 * Make an image that contains the effect of the gridding kernel in image space.
		 * The gridder already corrects for the kernel internally, so it is not necessary
		 * to apply this image to the results, but users might be interested in the
		 * response function of the gridding kernel alone.
		 * @param image A width*height array that will be set to the kernel image.
		 */
		void GetGriddingCorrectionImage(double* image) const;
		
		/**
		 * Replace the current result of the imaging with a new image. This is
		 * hardly ever really necessary but can save some memory and/or time when the
		 * WStackingGridder remains in memory with its results, and the results are
		 * passed on but need to be conditionally changed in some way (e.g. as done in
		 * the @ref WSMSGridder class to pass a resized buffer to the WSClean class).
		 * 
		 * After calling this method, @ref RealImage() will return @p newBuffer. The
		 * newBuffer will be cleaned up with the earlier provided ImageBufferAllocator.
		 * The old buffer can no longer be used.
		 * 
		 * This method replaces the real part of the imaging.
		 * @param newBuffer The new buffer that was allocated with the right
		 * ImageBufferAllocator.
		 * @see ReplaceImaginaryImageBuffer().
		 */
		void ReplaceRealImageBuffer(double* newBuffer);
		
		/**
		 * Replace the current imaginary result of the imaging with a new image.
		 * Exactly like @ref ReplaceRealImageBuffer(), except it replaces the imaginary
		 * part.
		 * @param newBuffer The new imaginary buffer part.
		 */
		void ReplaceImaginaryImageBuffer(double* newBuffer);
		
		/**
		 * Retrieve a gridded uv layer. This function can be called after
		 * @ref StartInversionPass() was called, and before @ref FinishInversionPass()
		 * is called.
		 * @param layerIndex Layer index of the grid, with zero being the first
		 * layer of the current pass.
		 * @returns The layer, with the currently gridded samples on it.
		 */
		const std::complex<double>* GetGriddedUVLayer(size_t layerIndex) const
		{
			return _layeredUVData[layerIndex];
		}
		
		/**
		 * Acquire a Kaiser-Bessel kernel. This is mostly a debugging/example function.
		 * @param kernel Array of size @p n
		 * @param n Size of kernel
		 * @param multiplyWithSinc True to multiply the kernel with the sinc function, as
		 * is the case during gridding.
		 */
		void GetKaiserBesselKernel(double* kernel, size_t n, bool multiplyWithSinc);
		
		/**
		 * Get width of image as specified during construction. This is the full width, before
		 * any trimming has been applied.
		 * @returns Width in pixels of image being gridded.
		 */
		size_t Width() const { return _width; }
		
		/**
		 * Get height of image as specified during construction. This is the full height, before
		 * any trimming has been applied.
		 * @returns Height in pixels of image being gridded.
		 */
		size_t Height() const { return _height; }
		
		/**
		 * Get angular width of single pixel in radians. The full image has a width of
		 * @ref Width() * PixelSizeX(), assuming small angle approximation.
		 * @returns Width of a single pixel in radians.
		 */
		double PixelSizeX() const { return _pixelSizeX; }
		
		/**
		 * Get angular height of single pixel in radians. The full image has a height of
		 * @ref Height() * PixelSizeX(), assuming small angle approximation.
		 * @returns Height of a single pixel in radians.
		 */
		double PixelSizeY() const { return _pixelSizeY; }
		
		/**
		 * The @ref ImageBufferAllocator of this gridder.
		 * @returns The image buffer allocator.
		 */
		ImageBufferAllocator* Allocator() const { return _imageBufferAllocator; }
	private:
		size_t layerRangeStart(size_t layerRangeIndex) const
		{
			return (_nWLayers * layerRangeIndex) / _nPasses;
		}
		template<bool IsComplexImpl>
		void projectOnImageAndCorrect(const std::complex<double> *source, double w, size_t threadIndex);
		template<bool IsComplexImpl>
		void copyImageToLayerAndInverseCorrect(std::complex<double> *dest, double w);
		void initializeSqrtLMLookupTable();
		void initializeSqrtLMLookupTableForSampling();
		void initializeLayeredUVData(size_t n);
		void freeLayeredUVData() { initializeLayeredUVData(0); }
		void fftToImageThreadFunction(boost::mutex *mutex, std::stack<size_t> *tasks, size_t threadIndex);
		void fftToUVThreadFunction(boost::mutex *mutex, std::stack<size_t> *tasks);
		void finalizeImage(double multiplicationFactor, std::vector<double*>& dataArray);
		void initializePrediction(const double *image, std::vector<double*>& dataArray);
		
		void makeKernels();
		/**
		 * Make the Kaiser Bessel windowed sinc functions.
		 * Alpha is a parameter of the Kaiser Bessel window function. Values of alpha correspond with the following functions:
		 * -   0 ; the rectangular window
		 * -   5 ; similar to the Hamming window
		 * -   6 ; similar to the Hann window
		 * - 8.6 ; similar to the Blackmann window
		 * -  14 ; similar to prolate spheroidal
		 */
		void makeKaiserBesselKernel(std::vector<double> &kernel, double alpha, size_t overSamplingFactor, bool withSinc=true);
		
		void makeRectangularKernel(std::vector<double> &kernel, size_t overSamplingFactor);
		
		double bessel0(double x, double precision);
		template<bool Inverse>
		void correctImageForKernel(double *image) const;
		
		const size_t _width, _height;
		const double _pixelSizeX, _pixelSizeY;
		size_t _nWLayers, _nPasses, _curLayerRangeIndex;
		double _minW, _maxW, _phaseCentreDL, _phaseCentreDM;
		bool _isComplex, _imageConjugatePart;
#ifndef AVOID_CASACORE
		MultiBandData _bandData;
#endif
		
		enum GridModeEnum _gridMode;
		size_t _overSamplingFactor, _kernelSize;
		std::vector<double> _1dKernel;
		std::vector<std::vector<double>> _griddingKernels;
		
		std::vector<std::complex<double>*> _layeredUVData;
		std::vector<double*> _imageData, _imageDataImaginary;
		std::vector<double> _sqrtLMLookupTable;
		size_t _nFFTThreads;
		ImageBufferAllocator* _imageBufferAllocator;
};

#endif
