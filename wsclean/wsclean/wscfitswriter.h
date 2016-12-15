#ifndef WSC_FITS_WRITER_H
#define WSC_FITS_WRITER_H

#include "imagefilename.h"
#include "imagingtable.h"
#include "outputchannelinfo.h"
#include "wscleansettings.h"

#include "../fitswriter.h"
#include "../polarization.h"

#include <string>

class WSCFitsWriter
{
public:
	WSCFitsWriter(const ImagingTableEntry& entry, bool isImaginary, const WSCleanSettings& settings, const class Deconvolution& deconvolution, size_t majorIterationNr, const class MSGridderBase& gridder, const std::string& commandLine, const OutputChannelInfo& channelInfo);
	
	explicit WSCFitsWriter(FitsReader& templateReader);
	
	FitsWriter& Writer() { return _writer; }
	
	void WriteImage(const std::string& suffix, const double* image);
	
	void WriteUV(const std::string& suffix, const double* image);
	
	void WritePSF(const std::string& fullname, const double* image);
	
	void SetPolarization(PolarizationEnum polarization)
	{
		_writer.SetPolarization(polarization);
	}
	
	/**
	 * Restore an elliptical beam using a FFT deconvolution directly from images.
	 */
	static void Restore(const class WSCleanSettings& settings);
		
private:
	
	void SetSettingsKeywords(const WSCleanSettings& settings, const std::string& commandLine);
	
	void setGridderKeywords(const WSCleanSettings& settings, const class MSGridderBase& gridder);
	
	void setDeconvolutionKeywords(const WSCleanSettings& settings);
	
	void setDeconvolutionResultKeywords(size_t minorIterationNr, size_t majorIterationNr);
	
	void setChannelKeywords(const ImagingTableEntry& entry, const OutputChannelInfo& channelInfo);
	
	void copyWSCleanKeywords(class FitsReader& reader);
	
private:
	FitsWriter _writer;
	std::string _filenamePrefix;
};

#endif

