#include "wscfitsfile.h"

#include "imagefilename.h"
#include "msgridderbase.h"

#include "../fitsreader.h"
#include "../modelrenderer.h"

#include "../deconvolution/deconvolution.h"

WSCFitsWriter::WSCFitsWriter(const ImagingTableEntry& entry, bool isImaginary, const WSCleanSettings& settings, const class Deconvolution& deconvolution, size_t majorIterationNr, const MSGridderBase& gridder, const std::string& commandLine, const OutputChannelInfo& channelInfo)
{
	_filenamePrefix = ImageFilename::GetPrefix(settings, entry.polarization, entry.outputChannelIndex, entry.outputIntervalIndex, isImaginary);
	setGridderKeywords(settings, gridder);
	SetSettingsKeywords(settings, commandLine);
	setChannelKeywords(entry, channelInfo);
	setDeconvolutionKeywords(settings);
	if(deconvolution.IsInitialized())
		setDeconvolutionResultKeywords(deconvolution.GetAlgorithm().IterationNumber(), majorIterationNr);
}

WSCFitsWriter::WSCFitsWriter(FitsReader& templateReader) : _writer(templateReader)
{
	copyWSCleanKeywords(templateReader);
}

void WSCFitsWriter::SetSettingsKeywords(const WSCleanSettings& settings, const std::string& commandLine)
{
	_writer.SetOrigin("WSClean", "W-stacking imager written by Andre Offringa");
	_writer.AddHistory(commandLine);
	if(settings.endChannel != 0)
	{
		_writer.SetExtraKeyword("WSCCHANS", settings.startChannel);
		_writer.SetExtraKeyword("WSCCHANE", settings.endChannel);
	}
	if(settings.endTimestep != 0)
	{
		_writer.SetExtraKeyword("WSCTIMES", settings.startTimestep);
		_writer.SetExtraKeyword("WSCTIMEE", settings.endTimestep);
	}
	_writer.SetExtraKeyword("WSCFIELD", settings.fieldId);
}

void WSCFitsWriter::setGridderKeywords(const WSCleanSettings& settings, const MSGridderBase& gridder)
{
	double
		ra = gridder.PhaseCentreRA(),
		dec = gridder.PhaseCentreDec(),
		pixelScaleX = gridder.PixelSizeX(),
		pixelScaleY = gridder.PixelSizeY(),
		dateObs = gridder.StartTime();
		
	_writer.SetImageDimensions(settings.trimmedImageWidth, settings.trimmedImageHeight, ra, dec, pixelScaleX, pixelScaleY);
	_writer.SetDate(dateObs);
	if(gridder.HasDenormalPhaseCentre())
		_writer.SetPhaseCentreShift(gridder.PhaseCentreDL(), gridder.PhaseCentreDM());
	
	/* This represents the weight of the image when averaging */
	_writer.SetExtraKeyword("WSCIMGWG", gridder.ImageWeight());
	/* This is the normalization factor that was applied. The factor is useful
	 * to undo the normalization for e.g. conversion to Kelvins. */ 
	_writer.SetExtraKeyword("WSCDATAC", gridder.DataColumnName());
	_writer.SetExtraKeyword("WSCWEIGH", gridder.Weighting().ToString());
	_writer.SetExtraKeyword("WSCGKRNL", gridder.AntialiasingKernelSize());
}

void WSCFitsWriter::setDeconvolutionKeywords(const WSCleanSettings& settings)
{
	_writer.SetExtraKeyword("WSCNITER", settings.deconvolutionIterationCount);
	_writer.SetExtraKeyword("WSCTHRES", settings.deconvolutionThreshold);
	_writer.SetExtraKeyword("WSCGAIN", settings.deconvolutionGain);
	_writer.SetExtraKeyword("WSCMGAIN", settings.deconvolutionMGain);
	_writer.SetExtraKeyword("WSCNEGCM", settings.allowNegativeComponents);
	_writer.SetExtraKeyword("WSCNEGST", settings.stopOnNegativeComponents);
}

void WSCFitsWriter::setDeconvolutionResultKeywords(size_t minorIterationNr, size_t majorIterationNr)
{
	_writer.SetExtraKeyword("WSCMINOR", minorIterationNr);
	_writer.SetExtraKeyword("WSCMAJOR", majorIterationNr);
}

void WSCFitsWriter::setChannelKeywords(const ImagingTableEntry& entry, const OutputChannelInfo& channelInfo)
{
	const double
		bandStart = entry.bandStartFrequency,
		bandEnd = entry.bandEndFrequency,
		centreFrequency = 0.5*(bandStart+bandEnd),
		bandwidth = bandEnd-bandStart;
	_writer.SetFrequency(centreFrequency, bandwidth);
	_writer.SetExtraKeyword("WSCIMGWG", channelInfo.weight);
	_writer.SetExtraKeyword("WSCNORMF", channelInfo.normalizationFactor);
	_writer.SetExtraKeyword("WSCNWLAY", channelInfo.wGridSize);
	_writer.SetExtraKeyword("WSCNVIS", channelInfo.visibilityCount);
	_writer.SetExtraKeyword("WSCENVIS", channelInfo.effectiveVisibilityCount);
	_writer.SetBeamInfo(
		channelInfo.beamMaj,
		channelInfo.beamMin,
		channelInfo.beamPA);
	_writer.SetPolarization(entry.polarization);
}

void WSCFitsWriter::copyWSCleanKeywords(FitsReader& reader)
{
	const size_t
		N_STRKEYWORDS=2, N_DBLKEYWORDS=18;
	const char* strKeywords[N_STRKEYWORDS] =
		{ "WSCDATAC", "WSCWEIGH" };
	const char* dblKeywords[N_DBLKEYWORDS] =
		{ "WSCIMGWG", "WSCNWLAY", "WSCGKRNL", "WSCCHANS", "WSCCHANE", "WSCTIMES", "WSCTIMEE", "WSCFIELD",
			"WSCNITER", "WSCTHRES", "WSCGAIN", "WSCMGAIN", "WSCNEGCM", "WSCNEGST",
			"WSCMINOR", "WSCMAJOR", "WSCNVIS", "WSCENVIS"
		};
	for(size_t i=0; i!=N_STRKEYWORDS; ++i)
		_writer.CopyStringKeywordIfExists(reader, strKeywords[i]);
	for(size_t i=0; i!=N_DBLKEYWORDS; ++i)
		_writer.CopyDoubleKeywordIfExists(reader, dblKeywords[i]);
}

void WSCFitsWriter::WriteImage(const std::string& suffix, const double* image)
{
	std::string name = _filenamePrefix + '-' + suffix;
	_writer.Write(name, image);
}

void WSCFitsWriter::WriteUV(const std::string& suffix, const double* image)
{
	std::string name = _filenamePrefix + '-' + suffix;
	FitsWriter::Unit unit = _writer.GetUnit();
	_writer.SetIsUV(true);
	_writer.SetUnit(FitsWriter::Jansky);
	_writer.Write(name, image);
	_writer.SetIsUV(false);
	_writer.SetUnit(unit);
}

void WSCFitsWriter::WritePSF(const std::string& fullname, const double* image)
{
	_writer.Write(fullname, image);
}

void WSCFitsWriter::Restore(const WSCleanSettings& settings)
{
	FitsReader imgReader(settings.restoreInput), modReader(settings.restoreModel);
	if(imgReader.ImageWidth() != modReader.ImageWidth() ||
		imgReader.ImageHeight() != modReader.ImageHeight())
		throw std::runtime_error("Image and model images have different dimensions!");
	ao::uvector<double>
		image(imgReader.ImageWidth() * imgReader.ImageHeight()),
		model(modReader.ImageWidth() * modReader.ImageHeight());
	imgReader.Read(image.data());
	modReader.Read(model.data());
	
	double beamMaj, beamMin, beamPA;
	if(settings.manualBeamMajorSize != 0.0)
	{
		beamMaj = settings.manualBeamMajorSize;
		beamMin = settings.manualBeamMinorSize;
		beamPA = settings.manualBeamPA;
	}
	else {
		beamMaj = imgReader.BeamMajorAxisRad();
		beamMin = imgReader.BeamMinorAxisRad();
		beamPA = imgReader.BeamPositionAngle();
	}
	
	ModelRenderer::Restore(image.data(), model.data(),
					imgReader.ImageWidth(), imgReader.ImageHeight(),
					beamMaj, beamMin, beamPA,
					imgReader.PixelSizeX(),
					imgReader.PixelSizeY());
	
	FitsWriter writer(WSCFitsWriter(imgReader).Writer());
	writer.SetBeamInfo(beamMaj, beamMin, beamPA);
	writer.Write(settings.restoreOutput, image.data());
}
