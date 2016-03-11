#ifndef WSCLEAN_SETTNGS_H
#define WSCLEAN_SETTNGS_H

#include "wstackinggridder.h"
#include "inversionalgorithm.h"
#include "../msselection.h"

#include "../deconvolution/deconvolutionalgorithm.h"

/**
 * This class describes all settings for a single WSClean run.
 * @sa WSClean
 */
struct WSCleanSettings
{
public:
	WSCleanSettings();
	
	void Validate();
	
	std::vector<std::string> filenames;
	size_t untrimmedImageWidth, untrimmedImageHeight;
	size_t trimmedImageWidth, trimmedImageHeight;
	size_t widthForNWCalculation, heightForNWCalculation;
	size_t channelsOut, intervalsOut;
	double pixelScaleX, pixelScaleY;
	double manualBeamMajorSize, manualBeamMinorSize, manualBeamPA;
	bool fittedBeam, theoreticBeam, circularBeam;
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
	WeightMode weightMode;
	std::string prefixName;
	bool smallInversion, makePSF, makePSFOnly, isWeightImageSaved, isUVImageSaved, isGriddingImageSaved, dftPrediction, dftWithBeam;
	std::string temporaryDirectory;
	bool forceReorder, forceNoReorder, subtractModel, modelUpdateRequired, mfsWeighting;
	bool normalizeForWeighting;
	bool applyPrimaryBeam, reusePrimaryBeam;
	enum WStackingGridder::GridModeEnum gridMode;
	enum InversionAlgorithm::VisibilityWeightingMode visibilityWeightingMode;
	
	/** @{
	 * These settings all relate to the deconvolution.
	 */
	double deconvolutionThreshold, deconvolutionGain, deconvolutionMGain;
	size_t deconvolutionIterationCount;
	bool allowNegativeComponents, stopOnNegativeComponents;
	bool useMultiscale, useFastMultiscale;
	double multiscaleDeconvolutionThresholdBias, multiscaleDeconvolutionScaleBias;
	bool multiscaleNormalizeResponse;
	ao::uvector<double> multiscaleScaleList;
	double deconvolutionBorderRatio;
	std::string fitsDeconvolutionMask, casaDeconvolutionMask;
	bool useMoreSaneDeconvolution, useIUWTDeconvolution;
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
	
	bool IsSpectralFittingEnabled() {
		return spectralFittingMode != NoSpectralFitting;
	}
};

inline WSCleanSettings::WSCleanSettings() :
	filenames(),
	untrimmedImageWidth(2048), untrimmedImageHeight(2048),
	trimmedImageWidth(0), trimmedImageHeight(0),
	widthForNWCalculation(0), heightForNWCalculation(0),
	channelsOut(1), intervalsOut(1),
	pixelScaleX(0.01 * M_PI / 180.0), pixelScaleY(0.01 * M_PI / 180.0),
	manualBeamMajorSize(0.0), manualBeamMinorSize(0.0),
	manualBeamPA(0.0), fittedBeam(true), theoreticBeam(false), circularBeam(false),
	memFraction(1.0), absMemLimit(0.0),
	minUVWInMeters(0.0), maxUVWInMeters(0.0),
	minUVInLambda(0.0), maxUVInLambda(0.0), wLimit(0.0),
	rankFilterLevel(0.0), rankFilterSize(16),
	gaussianTaperBeamSize(0.0),
	tukeyTaperInLambda(0.0), tukeyInnerTaperInLambda(0.0),
	edgeTaperInLambda(0.0), edgeTukeyTaperInLambda(0.0),
	nWLayers(0), antialiasingKernelSize(7), overSamplingFactor(63),
	threadCount(sysconf(_SC_NPROCESSORS_ONLN)),
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
	isUVImageSaved(false), isGriddingImageSaved(false),
	dftPrediction(false), dftWithBeam(false),
	temporaryDirectory(),
	forceReorder(false), forceNoReorder(false),
	subtractModel(false),
	modelUpdateRequired(true),
	mfsWeighting(false),
	normalizeForWeighting(true),
	applyPrimaryBeam(false), reusePrimaryBeam(false),
	gridMode(WStackingGridder::KaiserBesselKernel),
	visibilityWeightingMode(InversionAlgorithm::NormalVisibilityWeighting),
// Deconvolution default settings:
	deconvolutionThreshold(0.0),
	deconvolutionGain(0.1),
	deconvolutionMGain(1.0),
	deconvolutionIterationCount(0),
	allowNegativeComponents(true), 
	stopOnNegativeComponents(false),
	useMultiscale(false),
	useFastMultiscale(false),
	multiscaleDeconvolutionThresholdBias(0.7),
	multiscaleDeconvolutionScaleBias(0.6),
	multiscaleNormalizeResponse(false),
	multiscaleScaleList(),
	deconvolutionBorderRatio(0.05),
	fitsDeconvolutionMask(),
	casaDeconvolutionMask(),
	useMoreSaneDeconvolution(false),
	useIUWTDeconvolution(false),
	moreSaneLocation(),
	moreSaneArgs(),
	spectralFittingMode(NoSpectralFitting),
	spectralFittingTerms(0),
	deconvolutionChannelCount(0)
{
	polarizations.insert(Polarization::StokesI);
}


#endif
