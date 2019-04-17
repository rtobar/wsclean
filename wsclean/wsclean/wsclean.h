#ifndef WSCLEAN_H
#define WSCLEAN_H

#include "../msproviders/msprovider.h"
#include "../msproviders/partitionedms.h"

#include "../msselection.h"
#include "../polarization.h"
#include "../weightmode.h"
#include "../stopwatch.h"

#include "../deconvolution/deconvolution.h"

#include "cachedimageset.h"
#include "imagebufferallocator.h"
#include "imagingtable.h"
#include "msgridderbase.h"
#include "outputchannelinfo.h"
#include "wscfitswriter.h"
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
	void runIndependentGroup(ImagingTable& groupTable, std::unique_ptr<class PrimaryBeam>& primaryBeam);
	void saveRestoredImagesForGroup(const ImagingTableEntry& tableEntry, std::unique_ptr<class PrimaryBeam>& primaryBeam) const;
	void predictGroup(const ImagingTable& imagingGroup);
	
	void runFirstInversion(ImagingTableEntry& entry, std::unique_ptr<class PrimaryBeam>& primaryBeam);
	
	void performReordering(bool isPredictMode);
	
	ImageWeights& initializeImageWeights(const ImagingTableEntry& entry, std::vector<std::pair<std::unique_ptr<MSProvider>, MSSelection>>& msList);
	void initializeMFSImageWeights();
	std::unique_ptr<MSProvider> initializeMSProvider(const ImagingTableEntry& entry, const MSSelection& selection, size_t filenameIndex, size_t dataDescId);
	void initializeCurMSProviders(const ImagingTableEntry& entry, class GriddingTask& task);
	void initializeMSProvidersForPB(const ImagingTableEntry& entry, std::vector<std::pair<std::unique_ptr<MSProvider>, MSSelection>>& msList, class PrimaryBeam& pb);
	void storeAndCombineXYandYX(CachedImageSet& dest, PolarizationEnum polarization, size_t joinedChannelIndex, bool isImaginary, const double* image);
	bool selectChannels(MSSelection& selection, size_t msIndex, size_t bandIndex, const ImagingTableEntry& entry);
	MSSelection selectInterval(MSSelection& fullSelection, size_t intervalIndex);
	void readEarlierModelImages(const ImagingTableEntry& entry);
	
	void makeImagingTable(size_t outputIntervalIndex);
	void makeImagingTableEntry(const std::vector<ChannelInfo>& channels, size_t outIntervalIndex, size_t outChannelIndex, ImagingTableEntry& entry);
	void addPolarizationsToImagingTable(size_t& joinedGroupIndex, size_t& squaredGroupIndex, size_t outChannelIndex, const ImagingTableEntry& templateEntry);
	std::unique_ptr<class ImageWeightCache> createWeightCache();
	
	void multiplyImage(double factor, double* image) const;
	void multiplyImage(double factor, ImageBufferAllocator::Ptr& image) const { multiplyImage(factor, image.data()); }
	
	void imagePSF(ImagingTableEntry& entry);
	void imagePSFCallback(ImagingTableEntry& entry, struct GriddingResult& result);
	
	void imageGridding();
	void imageMain(ImagingTableEntry& entry, bool isFirst, bool updateBeamInfo, bool isInitialInversion);
	void imageMainCallback(ImagingTableEntry& entry, struct GriddingResult& result, bool updateBeamInfo, bool isInitialInversion);
	
	void predict(const ImagingTableEntry& entry);
	void predictCallback(const ImagingTableEntry& entry, struct GriddingResult& result);
	
	//void makeMFSImage(const string& suffix, size_t intervalIndex, PolarizationEnum pol, bool isImaginary, bool isPSF = false);
	//void renderMFSImage(size_t intervalIndex, PolarizationEnum pol, bool isImaginary, bool isPBCorrected) const;
	void saveUVImage(const double* image, const ImagingTableEntry& entry, bool isImaginary, const std::string& prefix) const;
	void writeFirstResidualImages(const ImagingTable& groupTable) const;
	void writeModelImages(const ImagingTable& groupTable) const;
	
	double minTheoreticalBeamSize(const ImagingTable& table) const;
	
	void makeBeam();
	
	WSCFitsWriter createWSCFitsWriter(const ImagingTableEntry& entry, bool isImaginary, bool isModel) const;
	
	WSCFitsWriter createWSCFitsWriter(const ImagingTableEntry& entry, PolarizationEnum polarization, bool isImaginary, bool isModel) const;
	
	bool preferReordering() const
	{
		return (
			(_settings.channelsOut != 1) ||
			(_settings.polarizations.size()>=4) ||
			(_settings.deconvolutionMGain != 1.0) ||
			(_settings.baselineDependentAveragingInWavelengths != 0.0) ||
			_settings.simulateNoise ||
			_settings.forceReorder
		) && !_settings.forceNoReorder;
	}
	
	// This must be the first field, because other members might take references
	// of this object and use them in their destructors.
	mutable ImageBufferAllocator _imageAllocator;
	
	MSSelection _globalSelection;
	std::string _commandLine;
	
	WSCleanSettings _settings;
	
	std::vector<OutputChannelInfo> _infoPerChannel;
	OutputChannelInfo _infoForMFS;
	std::map<size_t, MSGridderBase::MetaDataCache> _msGridderMetaCache;
	
	std::unique_ptr<class GriddingTaskManager> _griddingTaskManager;
	std::unique_ptr<class ImageWeightCache> _imageWeightCache;
	Stopwatch _inversionWatch, _predictingWatch, _deconvolutionWatch;
	bool _isFirstInversion, _doReorder;
	size_t _majorIterationNr;
	CachedImageSet _psfImages, _modelImages, _residualImages;
	std::vector<PartitionedMS::Handle> _partitionedMSHandles;
	std::vector<MultiBandData> _msBands;
	Deconvolution _deconvolution;
	ImagingTable _imagingTable;
	double _phaseCentreRA, _phaseCentreDec;
};

#endif
