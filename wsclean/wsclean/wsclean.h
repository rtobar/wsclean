#ifndef WSCLEAN_H
#define WSCLEAN_H

#include "../msproviders/msprovider.h"
#include "../msproviders/partitionedms.h"

#include "../msselection.h"
#include "../polarizationenum.h"
#include "../weightmode.h"
#include "../stopwatch.h"

#include "../deconvolution/deconvolution.h"

#include "cachedimageset.h"
#include "wstackinggridder.h"
#include "inversionalgorithm.h"
#include "imagebufferallocator.h"
#include "imagingtable.h"

#include <set>

class WSClean
{
public:
	WSClean();
	~WSClean();
	
	Deconvolution& DeconvolutionInfo() { return _deconvolution; }
	
	void SetImageSize(size_t width, size_t height) { _untrimmedWidth = width; _untrimmedHeight = height; }
	void SetTrimmedImageSize(size_t width, size_t height) { _trimWidth = width; _trimHeight = height; }
	void SetPixelScale(double pixelScale) { _pixelScaleX = pixelScale; _pixelScaleY = pixelScale; }
	void SetNWlayers(size_t nWLayers) { _nWLayers = nWLayers; }
	void SetColumnName(const std::string& columnName) { _columnName = columnName; }
	void SetPolarizations(const std::set<PolarizationEnum>& polarizations) { _polarizations = polarizations; }
	void SetMakePSF(bool makePSF) { _makePSF = makePSF; }
	void SetPrefixName(const std::string& prefixName) { _prefixName = prefixName; }
	void SetGridMode(WStackingGridder::GridModeEnum gridMode) { _gridMode = gridMode; }
	//void SetSmallPSF(bool smallPSF) { _smallPSF = smallPSF; }
	void SetSmallInversion(bool smallInversion) { _smallInversion = smallInversion; }
	void SetIntervalSelection(size_t startTimestep, size_t endTimestep) {
		_globalSelection.SetInterval(startTimestep, endTimestep);
	}
	void SetChannelSelection(size_t startChannel, size_t endChannel) {
		_startChannel = startChannel;
		_endChannel = endChannel;
	}
	void SetFieldSelection(size_t fieldId) {
		_globalSelection.SetFieldId(fieldId);
	}
	void SetChannelsOut(size_t channelsOut) { _channelsOut = channelsOut; }
	void SetIntervalCount(size_t intervalCount) { _intervalCount = intervalCount; }
	void SetMFSWeighting(bool mfsWeighting) { _mfsWeighting = mfsWeighting; }
	void SetWeightMode(enum WeightMode::WeightingEnum weighting) {
		_weightMode.SetMode(WeightMode(weighting));
	}
	void SetBriggsWeighting(double robustness) {
		_weightMode.SetMode(WeightMode::Briggs(robustness));
	}
	void SetSuperWeight(double superWeight) { _weightMode.SetSuperWeight(superWeight); }
	void SetBeamSize(double major, double minor, double positionAngle) {
		_manualBeamMajorSize = major;
		_manualBeamMinorSize = minor;
		_manualBeamPA = positionAngle;
	}
	void SetFittedBeam(bool fittedBeam) { _fittedBeam = fittedBeam; }
	void SetTheoreticBeam(bool theoreticBeam) { _theoreticBeam = theoreticBeam; }
	void SetCircularBeam(bool circularBeam) { _circularBeam = circularBeam; }
	void SetAntialiasingKernelSize(size_t kernelSize) { _antialiasingKernelSize = kernelSize; }
	void SetOversamplingFactor(size_t oversampling) { _overSamplingFactor = oversampling; }
	void SetThreadCount(size_t threadCount) { _threadCount = threadCount; }
	void SetTemporaryDirectory(const std::string& tempDir) { _temporaryDirectory = tempDir; }
	void SetForceReorder(bool forceReorder) { _forceReorder = forceReorder; }
	void SetForceNoReorder(bool forceNoReorder) { _forceNoReorder = forceNoReorder; }
	void SetSubtractModel(bool subtractModel) { _subtractModel = subtractModel; }
	void SetModelUpdateRequired(bool modelUpdateRequired) { _modelUpdateRequired = modelUpdateRequired; }
	void SetMemFraction(double memFraction) { _memFraction = memFraction; }
	void SetMemAbsLimit(double absMemLimit) { _absMemLimit = absMemLimit; }
	void SetMinUVWInM(double minUVW) { _globalSelection.SetMinUVWInM(minUVW); }
	void SetMaxUVWInM(double maxUVW) { _globalSelection.SetMaxUVWInM(maxUVW); }
	void SetMinUVInLambda(double lambda) { _minUVInLambda = lambda; }
	void SetMaxUVInLambda(double lambda) { _maxUVInLambda = lambda; }
	void SetWLimit(double wLimit) { _wLimit = wLimit; }
	void SetCommandLine(const std::string& cmdLine) { _commandLine = cmdLine; }
	void SetSaveWeights(bool saveWeights) { _isWeightImageSaved = saveWeights; }
	void SetSaveUV(bool saveUV) { _isUVImageSaved = saveUV; }
	void SetSaveGriddingImage(bool isGriddingImageSaved) { _isGriddingImageSaved = isGriddingImageSaved; }
	void SetDFTPrediction(bool dftPrediction) { _dftPrediction = dftPrediction; }
	void SetDFTWithBeam(bool applyBeam) { _dftWithBeam = applyBeam; }
	void SetRankFilterLevel(double level) { _rankFilterLevel = level; }
	void SetRankFilterSize(size_t nPixels) { _rankFilterSize = nPixels; }
	void SetJoinPolarizations(bool joinPolarizations) { _joinedPolarizationCleaning = joinPolarizations; }
	bool JoinPolarizations() const { return _joinedPolarizationCleaning; }
	
	void SetJoinChannels(bool joinChannels) { _joinedFrequencyCleaning = joinChannels; }
	bool JoinChannels() const { return _joinedFrequencyCleaning; }
	
	void AddInputMS(const std::string& msPath) { _filenames.push_back(msPath); }
	
	void RunClean();
	
	void RunPredict();
	
	void SetNormalizeForWeighting(bool normalizeForWeighting)
	{
		_normalizeForWeighting = normalizeForWeighting;
	}
	void SetVisibilityWeightingMode(enum InversionAlgorithm::VisibilityWeightingMode mode)
	{
		_visibilityWeightingMode = mode;
	}
private:
	void runIndependentGroup(const ImagingTable& groupTable);
	void predictGroup(const ImagingTable& imagingGroup);
	
	void runFirstInversion(const ImagingTableEntry& entry);
	void prepareInversionAlgorithm(PolarizationEnum polarization);
	
	void validateDimensions();
	void checkPolarizations();
	void performReordering(bool isPredictMode);
	
	void initFitsWriter(class FitsWriter& writer);
	void copyWSCleanKeywords(FitsReader& reader, FitsWriter& writer);
	//void copyDoubleKeywordIfExists(FitsReader& reader, FitsWriter& writer, const char* keywordName);
	void setCleanParameters(class FitsWriter& writer);
	void updateCleanParameters(class FitsWriter& writer, size_t minorIterationNr, size_t majorIterationNr);
	void initializeWeightTapers();
	void initializeImageWeights(const ImagingTableEntry& entry);
	void initializeMFSImageWeights();
	MSProvider* initializeMSProvider(const ImagingTableEntry& entry, const MSSelection& selection, size_t filenameIndex, size_t bandIndex);
	void initializeCurMSProviders(const ImagingTableEntry& entry);
	void clearCurMSProviders();
	void storeAndCombineXYandYX(CachedImageSet& dest, PolarizationEnum polarization, size_t joinedChannelIndex, bool isImaginary, const double* image);
	bool selectChannels(MSSelection& selection, size_t msIndex, size_t bandIndex, const ImagingTableEntry& entry);
	MSSelection selectInterval(MSSelection& fullSelection);
	
	void makeImagingTable();
	void makeImagingTableEntry(const std::vector<double>& channels, size_t outChannelIndex, ImagingTableEntry& entry);
	void addPolarizationsToImagingTable(size_t& joinedGroupIndex, size_t& squaredGroupIndex, size_t outChannelIndex, const ImagingTableEntry& templateEntry);
	
	void imagePSF(size_t currentChannelIndex);
	void imageGridding();
	void imageMainFirst(PolarizationEnum polarization, size_t channelIndex);
	void imageMainNonFirst(PolarizationEnum polarization, size_t channelIndex);
	void predict(PolarizationEnum polarization, size_t channelIndex);
	void dftPredict(const ImagingTable& squaredGroup);
	
	void makeMFSImage(const string& suffix, PolarizationEnum pol, bool isImaginary);
	void writeFits(const string& suffix, const double* image, PolarizationEnum pol, size_t channelIndex, bool isImaginary);
	void saveUVImage(const double* image, PolarizationEnum pol, size_t channelIndex, bool isImaginary, const std::string& prefix);
	void writeFirstResidualImages(const ImagingTable& groupTable);
	void writeModelImages(const ImagingTable& groupTable);
	
	std::string fourDigitStr(size_t val) const
	{
		std::ostringstream str;
		if(val < 1000) str << '0';
		if(val < 100) str << '0';
		if(val < 10) str << '0';
		str << val;
		return str.str();
	}
	
	std::string getPSFPrefix(size_t channelIndex) const
	{
		std::ostringstream partPrefixNameStr;
		partPrefixNameStr << _prefixName;
		if(_intervalCount != 1)
			partPrefixNameStr << "-t" << fourDigitStr(_currentIntervalIndex);
		if(_channelsOut != 1)
			partPrefixNameStr << '-' << fourDigitStr(channelIndex);
		return partPrefixNameStr.str();
	}
	
	std::string getPrefix(PolarizationEnum polarization, size_t channelIndex, bool isImaginary) const
	{
		std::ostringstream partPrefixNameStr;
		partPrefixNameStr << _prefixName;
		if(_intervalCount != 1)
			partPrefixNameStr << "-t" << fourDigitStr(_currentIntervalIndex);
		if(_channelsOut != 1)
			partPrefixNameStr << '-' << fourDigitStr(channelIndex);
		if(_polarizations.size() != 1)
		{
			partPrefixNameStr << '-' << Polarization::TypeToShortString(polarization);
			if(isImaginary)
				partPrefixNameStr << 'i';
		}
		return partPrefixNameStr.str();
	}
	
	std::string getMFSPrefix(PolarizationEnum polarization, bool isImaginary) const
	{
		std::ostringstream partPrefixNameStr;
		partPrefixNameStr << _prefixName;
		if(_intervalCount != 1)
			partPrefixNameStr << "-t" << fourDigitStr(_currentIntervalIndex);
		if(_channelsOut != 1)
			partPrefixNameStr << "-MFS";
		if(_polarizations.size() != 1)
		{
			partPrefixNameStr << '-' << Polarization::TypeToShortString(polarization);
			if(isImaginary)
				partPrefixNameStr << 'i';
		}
		return partPrefixNameStr.str();
	}
	
	bool preferReordering() const
	{
		return (
			(_channelsOut != 1) ||
			(_polarizations.size()>=4) ||
			(_deconvolution.MGain() != 1.0) ||
			_forceReorder
		) && !_forceNoReorder;
	}
	
	size_t _untrimmedWidth, _untrimmedHeight;
	size_t _trimWidth, _trimHeight;
	size_t _channelsOut, _intervalCount;
	double _pixelScaleX, _pixelScaleY;
	double _manualBeamMajorSize, _manualBeamMinorSize, _manualBeamPA;
	bool _fittedBeam, _theoreticBeam, _circularBeam;
	double _memFraction, _absMemLimit, _minUVInLambda, _maxUVInLambda, _wLimit, _rankFilterLevel;
	size_t _rankFilterSize;
	size_t _nWLayers, _antialiasingKernelSize, _overSamplingFactor, _threadCount;
	size_t _startChannel, _endChannel;
	bool _joinedPolarizationCleaning, _joinedFrequencyCleaning;
	MSSelection _globalSelection;
	std::string _columnName;
	std::set<PolarizationEnum> _polarizations;
	WeightMode _weightMode;
	std::string _prefixName;
	bool _smallInversion, _makePSF, _isWeightImageSaved, _isUVImageSaved, _isGriddingImageSaved, _dftPrediction, _dftWithBeam;
	std::string _temporaryDirectory;
	bool _forceReorder, _forceNoReorder, _subtractModel, _modelUpdateRequired, _mfsWeighting;
	enum WStackingGridder::GridModeEnum _gridMode;
	std::vector<std::string> _filenames;
	std::string _commandLine;
	std::vector<double> _inputChannelFrequencies;
	
	struct ChannelInfo {
		ChannelInfo() :
			weight(0.0),
			bandStart(0.0), bandEnd(0.0),
			beamMaj(0.0), beamMin(0.0), beamPA(0.0)
		{ }
		double weight;
		double bandStart, bandEnd;
		double beamMaj, beamMin, beamPA;
	};
	std::vector<ChannelInfo> _infoPerChannel;
	
	std::unique_ptr<class WSMSGridder> _inversionAlgorithm;
	std::unique_ptr<class ImageWeightCache> _imageWeightCache;
	ImageBufferAllocator _imageAllocator;
	Stopwatch _inversionWatch, _predictingWatch, _deconvolutionWatch;
	bool _isFirstInversion, _doReorder;
	size_t _currentIntervalIndex, _majorIterationNr;
	CachedImageSet _psfImages, _modelImages, _residualImages;
	std::vector<PartitionedMS::Handle> _partitionedMSHandles;
	FitsWriter _fitsWriter;
	std::vector<MSProvider*> _currentPolMSes;
	std::vector<MultiBandData> _msBands;
	bool _normalizeForWeighting;
	enum InversionAlgorithm::VisibilityWeightingMode _visibilityWeightingMode;
	Deconvolution _deconvolution;
	ImagingTable _imagingTable;
};

#endif
