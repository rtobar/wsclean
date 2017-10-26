#ifndef MEASUREMENT_SET_GRIDDER_H
#define MEASUREMENT_SET_GRIDDER_H

#include "gridmodeenum.h"

#include "../polarization.h"
#include "../msselection.h"
#include "../weightmode.h"

#include <cmath>
#include <string>
#include <vector>

class MeasurementSetGridder
{
	public:
		/**
		 * This specifies several modi how the visibility weights
		 * (normally stored in the SPECTRUM_WEIGHT column)
		 * are applied to the data.
		 */
		enum VisibilityWeightingMode
		{
			NormalVisibilityWeighting,
			SquaredVisibilityWeighting,
			UnitVisibilityWeighting
		};
		
		MeasurementSetGridder() :
			_imageWidth(0),
			_imageHeight(0),
			_trimWidth(0), _trimHeight(0),
			_nwWidth(0), _nwHeight(0),
			_pixelSizeX((1.0 / 60.0) * M_PI / 180.0),
			_pixelSizeY((1.0 / 60.0) * M_PI / 180.0),
			_wGridSize(0),
			_measurementSets(),
			_dataColumnName("DATA"),
			_doImagePSF(false),
			_doSubtractModel(false),
			_addToModel(false),
			_smallInversion(false),
			_wLimit(0.0),
			_precalculatedWeightInfo(0),
			_polarization(Polarization::StokesI),
			_isComplex(false),
			_weighting(WeightMode::UniformWeighted),
			_verbose(false),
			_antialiasingKernelSize(7),
			_overSamplingFactor(63),
			_normalizeForWeighting(true),
			_visibilityWeightingMode(NormalVisibilityWeighting),
			_gridMode(KaiserBesselKernel),
			_storeImagingWeights(false)
		{
		}
		virtual ~MeasurementSetGridder()
		{
		}
		
		size_t ImageWidth() const { return _imageWidth; }
		size_t ImageHeight() const { return _imageHeight; }
		double PixelSizeX() const { return _pixelSizeX; }
		double PixelSizeY() const { return _pixelSizeY; }
		bool HasWGridSize() const { return _wGridSize != 0; }
		size_t WGridSize() const { return _wGridSize; }
		
		void ClearMeasurementSetList() { _measurementSets.clear(); _selections.clear(); }
		class MSProvider& MeasurementSet(size_t index) const { return *_measurementSets[index]; }
		const MSSelection& Selection(size_t index) const { return _selections[index]; }
		size_t MeasurementSetCount() const { return _measurementSets.size(); }
		void AddMeasurementSet(class MSProvider* msProvider, const MSSelection& selection)
		{
			_measurementSets.push_back(msProvider);
			_selections.push_back(selection);
		}
		
		const std::string &DataColumnName() const { return _dataColumnName; }
		bool DoImagePSF() const { return _doImagePSF; }
		bool DoSubtractModel() const { return _doSubtractModel; }
		bool AddToModel() const { return _addToModel; }
		bool SmallInversion() const { return _smallInversion; }
		PolarizationEnum Polarization() const { return _polarization; }
		WeightMode Weighting() const { return _weighting; }
		class ImageWeights* PrecalculatedWeightInfo() const { return _precalculatedWeightInfo; }
		bool IsComplex() const { return _isComplex; }
		bool Verbose() const { return _verbose; }
		size_t AntialiasingKernelSize() const { return _antialiasingKernelSize; }
		size_t OverSamplingFactor() const { return _overSamplingFactor; }
		bool HasWLimit() const { return _wLimit != 0.0; }
		double WLimit() const { return _wLimit; }
		bool NormalizeForWeighting() const { return _normalizeForWeighting; }
		enum VisibilityWeightingMode VisibilityWeightingMode() const { return _visibilityWeightingMode; }
		bool StoreImagingWeights() const { return _storeImagingWeights; }
		
		void SetImageWidth(size_t imageWidth)
		{
			_imageWidth = imageWidth;
		}
		void SetImageHeight(size_t imageHeight)
		{
			_imageHeight = imageHeight;
		}
		void SetPixelSizeX(double pixelSizeX)
		{
			_pixelSizeX = pixelSizeX;
		}
		void SetPixelSizeY(double pixelSizeY)
		{
			_pixelSizeY = pixelSizeY;
		}
		void SetWGridSize(size_t wGridSize)
		{
			_wGridSize = wGridSize;
		}
		void SetNoWGridSize()
		{
			_wGridSize = 0;
		}
		void SetDataColumnName(const std::string &dataColumnName)
		{
			_dataColumnName = dataColumnName;
		}
		void SetDoImagePSF(bool doImagePSF)
		{
			_doImagePSF = doImagePSF;
		}
		void SetPolarization(PolarizationEnum polarization)
		{
			_polarization = polarization;
		}
		void SetIsComplex(bool isComplex)
		{
			_isComplex = isComplex;
		}
		void SetWeighting(WeightMode weighting)
		{
			_weighting = weighting;
		}
		void SetDoSubtractModel(bool doSubtractModel)
		{
			_doSubtractModel = doSubtractModel;
		}
		void SetAddToModel(bool addToModel)
		{
			_addToModel = addToModel;
		}
		void SetSmallInversion(bool smallInversion)
		{ 
			_smallInversion = smallInversion;
		}
		void SetPrecalculatedWeightInfo(class ImageWeights* precalculatedWeightInfo)
		{ 
			_precalculatedWeightInfo = precalculatedWeightInfo;
		}
		void SetVerbose(bool verbose)
		{
			_verbose = verbose;
		}
		void SetAntialiasingKernelSize(size_t kernelSize)
		{
			_antialiasingKernelSize = kernelSize;
		}
		void SetOverSamplingFactor(size_t factor)
		{
			_overSamplingFactor = factor;
		}
		void SetWLimit(double wLimit)
		{
			_wLimit = wLimit;
		}
		void SetNormalizeForWeighting(bool normalizeForWeighting)
		{
			_normalizeForWeighting = normalizeForWeighting;
		}
		void SetVisibilityWeightingMode(enum VisibilityWeightingMode mode)
		{
			_visibilityWeightingMode = mode;
		}
		void SetStoreImagingWeights(bool storeImagingWeights)
		{
			_storeImagingWeights = storeImagingWeights;
		}
		
		virtual void Invert() = 0;
		
		virtual void Predict(double* image) = 0;
		virtual void Predict(double* real, double* imaginary) = 0;
		
		virtual double *ImageRealResult() = 0;
		virtual double *ImageImaginaryResult() = 0;
		virtual double PhaseCentreRA() const = 0;
		virtual double PhaseCentreDec() const = 0;
		virtual bool HasDenormalPhaseCentre() const { return false; }
		virtual double PhaseCentreDL() const = 0;
		virtual double PhaseCentreDM() const = 0;
		virtual double BeamSize() const = 0;
		virtual double StartTime() const = 0;
		virtual double ImageWeight() const = 0;
		virtual double NormalizationFactor() const = 0;
		
		/**
		 * Deallocate any data that is no longer necessary, but all methods
		 * will still return results from the imaging, with the exception of
		 * ImageReal/ImageResult() and GetGriddingCorrectionImage().
		 */
		virtual void FreeImagingData() {	}
		
		virtual bool HasGriddingCorrectionImage() const = 0;
		virtual void GetGriddingCorrectionImage(double *image) const = 0;
		
		virtual size_t ActualInversionWidth() const { return _imageWidth; }
		virtual size_t ActualInversionHeight() const { return _imageHeight; }
		
		enum GridModeEnum GridMode() const { return _gridMode; }
		void SetGridMode(GridModeEnum gridMode) { _gridMode = gridMode; }
		
		size_t TrimWidth() const { return _trimWidth; }
		size_t TrimHeight() const { return _trimHeight; }
		bool HasTrimSize() const {
			return _trimWidth != 0 || _trimHeight != 0;
		}
		void SetTrimSize(size_t trimWidth, size_t trimHeight) {
			_trimWidth = trimWidth;
			_trimHeight = trimHeight;
		}
		bool HasNWSize() const {
			return _nwWidth!=0 || _nwHeight!=0;
		}
		size_t NWWidth() const { return _nwWidth; }
		size_t NWHeight() const { return _nwHeight; }
		void SetNWSize(size_t nwWidth, size_t nwHeight) {
			_nwWidth = nwWidth;
			_nwHeight = nwHeight;
		}
	private:
		size_t _imageWidth, _imageHeight;
		size_t _trimWidth, _trimHeight;
		size_t _nwWidth, _nwHeight;
		double _pixelSizeX, _pixelSizeY;
		size_t _wGridSize;
		std::vector<MSProvider*> _measurementSets;
		std::string _dataColumnName;
		bool _doImagePSF, _doSubtractModel, _addToModel, _smallInversion;
		double _wLimit;
		class ImageWeights *_precalculatedWeightInfo;
		PolarizationEnum _polarization;
		bool _isComplex;
		WeightMode _weighting;
		bool _verbose;
		std::vector<MSSelection> _selections;
		size_t _antialiasingKernelSize, _overSamplingFactor;
		bool _normalizeForWeighting;
		enum VisibilityWeightingMode _visibilityWeightingMode;
		GridModeEnum _gridMode;
		bool _storeImagingWeights;
};

#endif
