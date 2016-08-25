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
#include "imagebufferallocator.h"
#include "imagingtable.h"
#include "wscleansettings.h"

#include <set>

class WSClean
{
public:
	WSClean();
	~WSClean();
	
	WSCleanSettings& Settings() { return _settings; }
	const WSCleanSettings& Settings() const { return _settings; }
	
	void SetCommandLine(const std::string& cmdLine) { _commandLine = cmdLine; }
	
	void RunClean();
	
	void RunPredict();
	
private:
	void runIndependentGroup(const ImagingTable& groupTable);
	void saveRestoredImagesForGroup(const ImagingTableEntry& tableEntry);
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
	void initializeImageWeights(const ImagingTableEntry& entry);
	void initializeMFSImageWeights();
	MSProvider* initializeMSProvider(const ImagingTableEntry& entry, const MSSelection& selection, size_t filenameIndex, size_t dataDescId);
	void initializeCurMSProviders(const ImagingTableEntry& entry);
	void initializeMSProvidersForPB(const ImagingTableEntry& entry, class PrimaryBeam& pb);
	void clearCurMSProviders();
	void storeAndCombineXYandYX(CachedImageSet& dest, PolarizationEnum polarization, size_t joinedChannelIndex, bool isImaginary, const double* image);
	bool selectChannels(MSSelection& selection, size_t msIndex, size_t bandIndex, const ImagingTableEntry& entry);
	MSSelection selectInterval(MSSelection& fullSelection);
	
	void makeImagingTable();
	void makeImagingTableEntry(const std::vector<OrderedChannel>& channels, size_t outChannelIndex, ImagingTableEntry& entry);
	void addPolarizationsToImagingTable(size_t& joinedGroupIndex, size_t& squaredGroupIndex, size_t outChannelIndex, const ImagingTableEntry& templateEntry);
	class ImageWeightCache* createWeightCache();
	
	void multiplyImage(double factor, double* image);
	void imagePSF(const ImagingTableEntry& entry);
	void imageGridding();
	void imageMainFirst(PolarizationEnum polarization, size_t channelIndex);
	void imageMainNonFirst(PolarizationEnum polarization, size_t channelIndex);
	void predict(PolarizationEnum polarization, size_t channelIndex);
	void dftPredict(const ImagingTable& squaredGroup);
	
	void makeMFSImage(const string& suffix, PolarizationEnum pol, bool isImaginary, bool isPSF = false);
	void renderMFSImage(PolarizationEnum pol, bool isImaginary);
	void writeFits(const string& suffix, const double* image, PolarizationEnum pol, const ImagingTableEntry& entry, bool isImaginary);
	void saveUVImage(const double* image, PolarizationEnum pol, const ImagingTableEntry& entry, bool isImaginary, const std::string& prefix);
	void writeFirstResidualImages(const ImagingTable& groupTable);
	void writeModelImages(const ImagingTable& groupTable);
	
	void fitBeamSize(double& bMaj, double& bMin, double& bPA, const double* image, double beamEstimate);
	void determineBeamSize(double& bMaj, double& bMin, double& bPA, const double* image, double theoreticBeam);
	
	void makeBeam();
	
	bool preferReordering() const
	{
		return (
			(_settings.channelsOut != 1) ||
			(_settings.polarizations.size()>=4) ||
			(_settings.deconvolutionMGain != 1.0) ||
			(_settings.baselineDependentAveragingInWavelengths != 0.0) ||
			_settings.forceReorder
		) && !_settings.forceNoReorder;
	}
	
	MSSelection _globalSelection;
	std::string _commandLine;
	std::vector<OrderedChannel> _inputChannelFrequencies;
	
	struct ChannelInfo {
		ChannelInfo() :
			weight(0.0),
			beamMaj(0.0), beamMin(0.0), beamPA(0.0),
			psfNormalizationFactor(1.0)
		{ }
		double weight;
		double beamMaj, beamMin, beamPA;
		double theoreticBeamSize, psfNormalizationFactor;
	};
	WSCleanSettings _settings;
	
	std::vector<ChannelInfo> _infoPerChannel;
	ChannelInfo _infoForMFS;
	
	std::unique_ptr<class MSGridderBase> _gridder;
	std::unique_ptr<class ImageWeightCache> _imageWeightCache;
	std::unique_ptr<class PrimaryBeam> _primaryBeam;
	ImageBufferAllocator _imageAllocator;
	Stopwatch _inversionWatch, _predictingWatch, _deconvolutionWatch;
	bool _isFirstInversion, _doReorder;
	size_t _currentIntervalIndex, _majorIterationNr;
	CachedImageSet _psfImages, _modelImages, _residualImages;
	std::vector<PartitionedMS::Handle> _partitionedMSHandles;
	FitsWriter _fitsWriter;
	std::vector<MSProvider*> _currentPolMSes;
	std::vector<MultiBandData> _msBands;
	Deconvolution _deconvolution;
	ImagingTable _imagingTable;
};

#endif
