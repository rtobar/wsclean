#ifndef WSCLEAN_SETTNGS_H
#define WSCLEAN_SETTNGS_H

#include "wstackinggridder.h"
#include "inversionalgorithm.h"

#include "../msselection.h"
#include "../system.h"

#include "../deconvolution/deconvolutionalgorithm.h"

/**
 * This class describes all settings for a single WSClean run.
 * @sa WSClean
 */
class WSCleanSettings
{
public:
	WSCleanSettings();
	
	void Validate() const;
	
	void Propogate() { setDimensions(); }
	
	std::vector<std::string> filenames;
	enum Mode { ImagingMode, PredictMode, RestoreMode } mode;
	size_t untrimmedImageWidth, untrimmedImageHeight;
	size_t trimmedImageWidth, trimmedImageHeight;
	double imagePadding;
	size_t widthForNWCalculation, heightForNWCalculation;
	size_t channelsOut, intervalsOut;
	double pixelScaleX, pixelScaleY;
	std::string restoreModel, restoreInput, restoreOutput;
	double manualBeamMajorSize, manualBeamMinorSize, manualBeamPA;
	bool fittedBeam, theoreticBeam, circularBeam;
	bool continuedRun;
	double memFraction, absMemLimit, minUVWInMeters, maxUVWInMeters, minUVInLambda, maxUVInLambda, wLimit, rankFilterLevel;
	size_t rankFilterSize;
	double gaussianTaperBeamSize, tukeyTaperInLambda, tukeyInnerTaperInLambda, edgeTaperInLambda, edgeTukeyTaperInLambda;
	size_t nWLayers, antialiasingKernelSize, overSamplingFactor, threadCount;
	size_t fieldId;
	size_t startTimestep, endTimestep;
	size_t startChannel, endChannel;
	bool joinedPolarizationCleaning, joinedFrequencyCleaning;
	size_t predictionChannels;
	std::string dataColumnName;
	std::set<PolarizationEnum> polarizations;
	std::set<size_t> spectralWindows;
	WeightMode weightMode;
	std::string prefixName;
	bool smallInversion, makePSF, makePSFOnly, isWeightImageSaved, isUVImageSaved, isDirtySaved, isGriddingImageSaved, dftPrediction, dftWithBeam;
	std::string temporaryDirectory;
	bool forceReorder, forceNoReorder, subtractModel, modelUpdateRequired, mfsWeighting;
	bool normalizeForWeighting;
	bool applyPrimaryBeam, reusePrimaryBeam, useDifferentialLofarBeam, useIDG;
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
	size_t deconvolutionIterationCount;
	bool allowNegativeComponents, stopOnNegativeComponents;
	bool useMultiscale, useClarkOptimization, squaredJoins, forceDynamicJoin;
	bool multiscaleFastSubMinorLoop;
	double multiscaleGain, multiscaleDeconvolutionScaleBias;
	bool multiscaleNormalizeResponse;
	ao::uvector<double> multiscaleScaleList;
	double deconvolutionBorderRatio;
	std::string fitsDeconvolutionMask, casaDeconvolutionMask;
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
	
	void GetMSSelection(MSSelection& selection)
	{
		selection = MSSelection();
		selection.SetInterval(startTimestep, endTimestep);
		selection.SetFieldId(fieldId);
		selection.SetMinUVWInM(minUVWInMeters);
		selection.SetMaxUVWInM(maxUVWInMeters);
	}
	
	bool IsSpectralFittingEnabled() const {
		return spectralFittingMode != NoSpectralFitting;
	}
	
private:
	void checkPolarizations() const;
	void setDimensions();
};

inline WSCleanSettings::WSCleanSettings() :
	filenames(),
	mode(ImagingMode),
	untrimmedImageWidth(0), untrimmedImageHeight(0),
	trimmedImageWidth(0), trimmedImageHeight(0),
	imagePadding(1.2),
	widthForNWCalculation(0), heightForNWCalculation(0),
	channelsOut(1), intervalsOut(1),
	pixelScaleX(0.0), pixelScaleY(0.0),
	restoreModel(), restoreInput(), restoreOutput(),
	manualBeamMajorSize(0.0), manualBeamMinorSize(0.0),
	manualBeamPA(0.0), fittedBeam(true), theoreticBeam(false), circularBeam(false),
	continuedRun(false),
	memFraction(1.0), absMemLimit(0.0),
	minUVWInMeters(0.0), maxUVWInMeters(0.0),
	minUVInLambda(0.0), maxUVInLambda(0.0), wLimit(0.0),
	rankFilterLevel(3.0), rankFilterSize(16),
	gaussianTaperBeamSize(0.0),
	tukeyTaperInLambda(0.0), tukeyInnerTaperInLambda(0.0),
	edgeTaperInLambda(0.0), edgeTukeyTaperInLambda(0.0),
	nWLayers(0), antialiasingKernelSize(7), overSamplingFactor(63),
	threadCount(System::ProcessorCount()),
	fieldId(0),
	startTimestep(0), endTimestep(0),
	startChannel(0), endChannel(0),
	joinedPolarizationCleaning(false), joinedFrequencyCleaning(false),
	predictionChannels(0),
	dataColumnName(),
	polarizations(),
	weightMode(WeightMode::UniformWeighted),
	prefixName("wsclean"),
	smallInversion(true), makePSF(false), makePSFOnly(false), isWeightImageSaved(false),
	isUVImageSaved(false), isDirtySaved(true), isGriddingImageSaved(false),
	dftPrediction(false), dftWithBeam(false),
	temporaryDirectory(),
	forceReorder(false), forceNoReorder(false),
	subtractModel(false),
	modelUpdateRequired(true),
	mfsWeighting(false),
	normalizeForWeighting(true),
	applyPrimaryBeam(false), reusePrimaryBeam(false),
	useDifferentialLofarBeam(false),
	useIDG(false),
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
	deconvolutionIterationCount(0),
	allowNegativeComponents(true), 
	stopOnNegativeComponents(false),
	useMultiscale(false),
	useClarkOptimization(true),
	squaredJoins(false),
	forceDynamicJoin(false),
	multiscaleFastSubMinorLoop(true),
	multiscaleGain(0.2),
	multiscaleDeconvolutionScaleBias(0.6),
	multiscaleNormalizeResponse(false),
	multiscaleScaleList(),
	deconvolutionBorderRatio(0.05),
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
