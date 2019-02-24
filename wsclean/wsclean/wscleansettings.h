#ifndef WSCLEAN_SETTINGS_H
#define WSCLEAN_SETTINGS_H

#include "wstackinggridder.h"
#include "measurementsetgridder.h"

#include "../msselection.h"
#include "../system.h"

#include "../deconvolution/deconvolutionalgorithm.h"
#include "../multiscale/multiscaletransforms.h"

enum class DirectFTPrecision { Half, Float, Double, LongDouble };

/**
 * This class describes all settings for a single WSClean run.
 * @sa WSClean
 */
class WSCleanSettings
{
public:
	WSCleanSettings();
	
	void Validate() const;
	
	void Propogate() { RecalculatePaddedDimensions(); }
	
	void RecalculatePaddedDimensions();
	
	std::vector<std::string> filenames;
	enum Mode { ImagingMode, PredictMode, RestoreMode } mode;
	size_t paddedImageWidth, paddedImageHeight;
	size_t trimmedImageWidth, trimmedImageHeight;
	double imagePadding;
	size_t widthForNWCalculation, heightForNWCalculation;
	size_t channelsOut, intervalsOut;
	enum MSSelection::EvenOddSelection evenOddTimesteps;
	bool divideChannelsByGaps;
	double pixelScaleX, pixelScaleY;
	std::string restoreModel, restoreInput, restoreOutput;
	double manualBeamMajorSize, manualBeamMinorSize, manualBeamPA;
	bool fittedBeam, theoreticBeam, circularBeam;
	double beamFittingBoxSize;
	bool continuedRun;
	double memFraction, absMemLimit, minUVWInMeters, maxUVWInMeters, minUVInLambda, maxUVInLambda, wLimit, rankFilterLevel;
	size_t rankFilterSize;
	double gaussianTaperBeamSize, tukeyTaperInLambda, tukeyInnerTaperInLambda, edgeTaperInLambda, edgeTukeyTaperInLambda;
	bool useWeightsAsTaper;
	size_t nWLayers, antialiasingKernelSize, overSamplingFactor, threadCount;
	size_t fieldId;
	size_t startTimestep, endTimestep;
	size_t startChannel, endChannel;
	size_t predictionChannels;
	std::string dataColumnName;
	std::set<PolarizationEnum> polarizations;
	std::set<size_t> spectralWindows;
	WeightMode weightMode;
	std::string prefixName;
	bool joinedPolarizationCleaning, joinedFrequencyCleaning;
	std::set<PolarizationEnum> linkedPolarizations;
	size_t parallelDeconvolutionMaxSize;
	bool smallInversion, makePSF, makePSFOnly, isWeightImageSaved, isUVImageSaved, isDirtySaved, isFirstResidualSaved, isGriddingImageSaved;
	bool writeImagingWeightSpectrumColumn;
	std::string temporaryDirectory;
	bool forceReorder, forceNoReorder, subtractModel, modelUpdateRequired, mfsWeighting;
	size_t fullResOffset, fullResWidth, fullResPad;
	bool applyPrimaryBeam, reusePrimaryBeam, useDifferentialLofarBeam, savePsfPb;
	std::string mwaPath;
	size_t primaryBeamUndersampling;
	bool directFT;
	DirectFTPrecision directFTPrecision;
	bool useIDG;
	std::string atermConfigFilename;
	double atermKernelSize;
	bool gridWithBeam;
	double beamAtermUpdateTime; // in seconds.
	bool saveATerms;
	enum IDGMode { IDG_DEFAULT, IDG_GPU, IDG_CPU, IDG_HYBRID } idgMode;
	enum GridModeEnum gridMode;
	enum MeasurementSetGridder::VisibilityWeightingMode visibilityWeightingMode;
	double baselineDependentAveragingInWavelengths;
	bool simulateNoise;
	double simulatedNoiseStdDev;
	
	/** @{
	 * These settings all relate to the deconvolution.
	 */
	double deconvolutionThreshold, deconvolutionGain, deconvolutionMGain;
	bool autoDeconvolutionThreshold, autoMask;
	double autoDeconvolutionThresholdSigma, autoMaskSigma;
	bool localRMS;
	double localRMSWindow;
	enum LocalRMSMethod { RMSWindow, RMSAndMinimumWindow } localRMSMethod;
	bool saveSourceList;
	size_t deconvolutionIterationCount, majorIterationCount;
	bool allowNegativeComponents, stopOnNegativeComponents;
	bool useMultiscale, useClarkOptimization, squaredJoins;
	//bool forceDynamicJoin;
	bool multiscaleFastSubMinorLoop;
	double multiscaleGain, multiscaleDeconvolutionScaleBias;
	bool multiscaleNormalizeResponse;
	double multiscaleConvolutionPadding;
	ao::uvector<double> multiscaleScaleList;
	MultiScaleTransforms::Shape multiscaleShapeFunction;
	
	double deconvolutionBorderRatio;
	std::string fitsDeconvolutionMask, casaDeconvolutionMask;
	std::string localRMSImage;
	bool useMoreSaneDeconvolution, useIUWTDeconvolution, iuwtSNRTest;
	std::string moreSaneLocation, moreSaneArgs;
	ao::uvector<double> moreSaneSigmaLevels;
	enum SpectralFittingMode spectralFittingMode;
	size_t spectralFittingTerms;
	/**
	 * The number of channels used during deconvolution. This can be used to
	 * image with more channels than deconvolution. Before deconvolution,
	 * channels are averaged, and after deconvolution they are interpolated.
	 * It is 0 when all channels should be used.
	 */
	size_t deconvolutionChannelCount;
	/**
	 * @}
	 */
	
	MSSelection GetMSSelection() const
	{
		MSSelection selection;
		selection.SetInterval(startTimestep, endTimestep);
		selection.SetFieldId(fieldId);
		selection.SetMinUVWInM(minUVWInMeters);
		selection.SetMaxUVWInM(maxUVWInMeters);
		selection.SetEvenOrOddTimesteps(evenOddTimesteps);
		return selection;
	}
	
	bool IsSpectralFittingEnabled() const {
		return spectralFittingMode != NoSpectralFitting;
	}
	
private:
	void checkPolarizations() const;
};

inline WSCleanSettings::WSCleanSettings() :
	filenames(),
	mode(ImagingMode),
	paddedImageWidth(0), paddedImageHeight(0),
	trimmedImageWidth(0), trimmedImageHeight(0),
	imagePadding(1.2),
	widthForNWCalculation(0), heightForNWCalculation(0),
	channelsOut(1), intervalsOut(1),
	evenOddTimesteps(MSSelection::AllTimesteps),
	divideChannelsByGaps(false),
	pixelScaleX(0.0), pixelScaleY(0.0),
	restoreModel(), restoreInput(), restoreOutput(),
	manualBeamMajorSize(0.0), manualBeamMinorSize(0.0),
	manualBeamPA(0.0), fittedBeam(true), theoreticBeam(false), circularBeam(false),
	beamFittingBoxSize(10.0),
	continuedRun(false),
	memFraction(1.0), absMemLimit(0.0),
	minUVWInMeters(0.0), maxUVWInMeters(0.0),
	minUVInLambda(0.0), maxUVInLambda(0.0), wLimit(0.0),
	rankFilterLevel(3.0), rankFilterSize(16),
	gaussianTaperBeamSize(0.0),
	tukeyTaperInLambda(0.0), tukeyInnerTaperInLambda(0.0),
	edgeTaperInLambda(0.0), edgeTukeyTaperInLambda(0.0),
	useWeightsAsTaper(false),
	nWLayers(0), antialiasingKernelSize(7), overSamplingFactor(63),
	threadCount(System::ProcessorCount()),
	fieldId(0),
	startTimestep(0), endTimestep(0),
	startChannel(0), endChannel(0),
	predictionChannels(0),
	dataColumnName(),
	polarizations(),
	weightMode(WeightMode::UniformWeighted),
	prefixName("wsclean"),
	joinedPolarizationCleaning(false), joinedFrequencyCleaning(false),
	linkedPolarizations(),
	parallelDeconvolutionMaxSize(0),
	smallInversion(true), makePSF(false), makePSFOnly(false), isWeightImageSaved(false),
	isUVImageSaved(false), isDirtySaved(true), isFirstResidualSaved(false), isGriddingImageSaved(false),
	writeImagingWeightSpectrumColumn(false),
	temporaryDirectory(),
	forceReorder(false), forceNoReorder(false),
	subtractModel(false),
	modelUpdateRequired(true),
	mfsWeighting(false),
	fullResOffset(0), fullResWidth(0), fullResPad(0),
	applyPrimaryBeam(false), reusePrimaryBeam(false),
	useDifferentialLofarBeam(false),
	savePsfPb(false),
	primaryBeamUndersampling(8),
	directFT(false),
	directFTPrecision(DirectFTPrecision::Double),
	useIDG(false),
	atermConfigFilename(),
	atermKernelSize(5.0),
	gridWithBeam(false),
	beamAtermUpdateTime(300.0),
	saveATerms(false),
	idgMode(IDG_DEFAULT),
	gridMode(KaiserBesselKernel),
	visibilityWeightingMode(MeasurementSetGridder::NormalVisibilityWeighting),
	baselineDependentAveragingInWavelengths(0.0),
	simulateNoise(false),
	simulatedNoiseStdDev(0.0),
// Deconvolution default settings:
	deconvolutionThreshold(0.0),
	deconvolutionGain(0.1),
	deconvolutionMGain(1.0),
	autoDeconvolutionThreshold(false),
	autoMask(false),
	autoDeconvolutionThresholdSigma(0.0),
	autoMaskSigma(0.0),
	localRMS(false),
	localRMSWindow(25.0),
	localRMSMethod(RMSWindow),
	saveSourceList(false),
	deconvolutionIterationCount(0),
	majorIterationCount(20),
	allowNegativeComponents(true), 
	stopOnNegativeComponents(false),
	useMultiscale(false),
	useClarkOptimization(true),
	squaredJoins(false),
	multiscaleFastSubMinorLoop(true),
	multiscaleGain(0.2),
	multiscaleDeconvolutionScaleBias(0.6),
	multiscaleNormalizeResponse(false),
	multiscaleConvolutionPadding(1.1),
	multiscaleScaleList(),
	multiscaleShapeFunction(MultiScaleTransforms::TaperedQuadraticShape),
	deconvolutionBorderRatio(0.0),
	fitsDeconvolutionMask(),
	casaDeconvolutionMask(),
	useMoreSaneDeconvolution(false),
	useIUWTDeconvolution(false),
	iuwtSNRTest(false),
	moreSaneLocation(),
	moreSaneArgs(),
	spectralFittingMode(NoSpectralFitting),
	spectralFittingTerms(0),
	deconvolutionChannelCount(0)
{
	polarizations.insert(Polarization::StokesI);
}


#endif
