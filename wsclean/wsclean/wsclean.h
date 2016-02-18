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
	class ImageWeightCache* createWeightCache();
	
	void imagePSF(size_t currentChannelIndex);
	void imageGridding();
	void imageMainFirst(PolarizationEnum polarization, size_t channelIndex);
	void imageMainNonFirst(PolarizationEnum polarization, size_t channelIndex);
	void predict(PolarizationEnum polarization, size_t channelIndex);
	void dftPredict(const ImagingTable& squaredGroup);
	
	void makeMFSImage(const string& suffix, PolarizationEnum pol, bool isImaginary, bool isPSF = false);
	void renderMFSImage(PolarizationEnum pol, bool isImaginary);
	void writeFits(const string& suffix, const double* image, PolarizationEnum pol, size_t channelIndex, bool isImaginary);
	void saveUVImage(const double* image, PolarizationEnum pol, size_t channelIndex, bool isImaginary, const std::string& prefix);
	void writeFirstResidualImages(const ImagingTable& groupTable);
	void writeModelImages(const ImagingTable& groupTable);
	
	void fitBeamSize(double& bMaj, double& bMin, double& bPA, const double* image, double beamEstimate);
	void determineBeamSize(double& bMaj, double& bMin, double& bPA, const double* image, double theoreticBeam);
	
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
		partPrefixNameStr << _settings.prefixName;
		if(_settings.intervalsOut != 1)
			partPrefixNameStr << "-t" << fourDigitStr(_currentIntervalIndex);
		if(_settings.channelsOut != 1)
			partPrefixNameStr << '-' << fourDigitStr(channelIndex);
		return partPrefixNameStr.str();
	}
	
	std::string getPrefix(PolarizationEnum polarization, size_t channelIndex, bool isImaginary) const
	{
		std::ostringstream partPrefixNameStr;
		partPrefixNameStr << _settings.prefixName;
		if(_settings.intervalsOut != 1)
			partPrefixNameStr << "-t" << fourDigitStr(_currentIntervalIndex);
		if(_settings.channelsOut != 1)
			partPrefixNameStr << '-' << fourDigitStr(channelIndex);
		if(_settings.polarizations.size() != 1)
		{
			partPrefixNameStr << '-' << Polarization::TypeToShortString(polarization);
			if(isImaginary)
				partPrefixNameStr << 'i';
		}
		return partPrefixNameStr.str();
	}
	
	std::string getMFSPrefix(PolarizationEnum polarization, bool isImaginary, bool isPSF) const
	{
		std::ostringstream partPrefixNameStr;
		partPrefixNameStr << _settings.prefixName;
		if(_settings.intervalsOut != 1)
			partPrefixNameStr << "-t" << fourDigitStr(_currentIntervalIndex);
		if(_settings.channelsOut != 1)
			partPrefixNameStr << "-MFS";
		if(_settings.polarizations.size() != 1 && !isPSF)
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
			(_settings.channelsOut != 1) ||
			(_settings.polarizations.size()>=4) ||
			(_settings.deconvolutionMGain != 1.0) ||
			_settings.forceReorder
		) && !_settings.forceNoReorder;
	}
	
	MSSelection _globalSelection;
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
		double theoreticBeamSize;
	};
	WSCleanSettings _settings;
	
	std::vector<ChannelInfo> _infoPerChannel;
	ChannelInfo _infoForMFS;
	
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
	Deconvolution _deconvolution;
	ImagingTable _imagingTable;
};

#endif
