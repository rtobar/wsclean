#include "angle.h"
#include "wsclean/wsclean.h"
#include "wscversion.h"

#include <boost/algorithm/string.hpp>

#include <iostream>

#ifdef WSCLEAN_NO_MAIN
int wsclean_main(int argc, char *argv[])
#else
int main(int argc, char *argv[])
#endif
{
	std::cout << "\n"
		"WSClean version " WSCLEAN_VERSION_STR " (" WSCLEAN_VERSION_DATE ")\n"
		"This software package is released under the GPL version 3.\n"
	  "Author: André Offringa (offringa@gmail.com).\n\n";
#ifndef NDEBUG
	std::cout << "\n"
		"WARNING: Symbol NDEBUG was not defined; this WSClean version was\n"
		"compiled as a DEBUG version. This can seriously affect performance!\n\n";
#endif
	
	if(argc < 2)
	{
		std::cout << "Syntax: wsclean [options] <input-ms> [<2nd-ms> [..]]\n"
			"Will create cleaned images of the input ms(es).\n"
			"If multiple mses are specified, they need to be phase-rotated to the same point on the sky.\n\n"
			"Options can be:\n\n"
			"  ** GENERAL OPTIONS **\n"
			"-version\n"
			"   Print WSClean's version and exit.\n"
			"-j <threads>\n"
			"   Specify number of computing threads to use, i.e., number of cpu cores that will be used.\n"
			"   Default: use all cpu cores.\n"
			"-mem <percentage>\n"
			"   Limit memory usage to the given fraction of the total system memory. This is an approximate value.\n"
			"   Default: 100.\n"
			"-absmem <memory limit>\n"
			"   Like -mem, but this specifies a fixed amount of memory in gigabytes.\n"
			"-reorder\n"
			"-no-reorder\n"
			"   Force or disable reordering of Measurement Set. This can be faster when the measurement set needs to\n"
			"   be iterated several times, such as with many major iterations or in channel imaging mode.\n"
			"   Default: only reorder when in channel imaging mode.\n"
			"-tempdir <directory>\n"
			"   Set the temporary directory used when reordering files. Default: same directory as input measurement set.\n"
			"-saveweights\n"
			"   Save the gridded weights in the a fits file named <image-prefix>-weights.fits.\n"
			"-saveuv\n"
			"   Save the gridded uv plane, i.e., the FFT of the residual image. The UV plane is complex, hence\n"
			"   two images will be output: <prefix>-uv-real.fits and <prefix>-uv-imag.fits.\n"
			"-update-model-required (default), and\n"
			"-no-update-model-required\n"
			"   These two options specify wether the model data column is required to\n"
			"   contain valid model data after imaging. It can save time to not update\n"
			"   the model data column.\n"
			"\n"
			"  ** INVERSION OPTIONS **\n"
			"-name <image-prefix>\n"
			"   Use image-prefix as prefix for output files. Default is 'wsclean'.\n"
			"-size <width> <height>\n"
			"   Default: 2048 x 2048\n"
			"-trim <width> <height>\n"
			"   After inversion, trim the image to the given size.\n"
			"-scale <pixel-scale>\n"
			"   Scale of a pixel. Default unit is degrees, but can be specificied, e.g. -scale 20asec. Default: 0.01deg.\n"
			"-nwlayers <nwlayers>\n"
			"   Number of w-layers to use. Default: minimum suggested #w-layers for first MS.\n"
			"-channelsout <count>\n"
			"   Splits the bandwidth and makes count nr. of images. Default: 1.\n"
			"-predict <image-prefix>\n"
			"   Only perform a single prediction for an existing image. Doesn't do any imaging or cleaning.\n"
			"-subtract-model\n"
			"   Subtract the model from the data column in the first iteration. This can be used to reimage\n"
			"   an already cleaned image, e.g. at a different resolution.\n"
			"-nosmallinversion and -smallinversion\n"
			"   Perform inversion at the Nyquist resolution and upscale the image to the requested image size afterwards.\n"
			"   This speeds up inversion considerably, but makes aliasing slightly worse. This effect is\n"
			"   in most cases <1%. Default: on.\n"
			"-weight <weightmode>\n"
			"   Weightmode can be: natural, mwa, uniform, briggs. Default: uniform. When using Briggs' weighting,\n"
			"   add the robustness parameter, like: \"-weight briggs 0.5\".\n"
			"-superweight <factor>\n"
			"   Increase the weight gridding box size, similar to Casa's superuniform weighting scheme. Default: 1.0\n"
			"   The factor can be rational and can be less than one for subpixel weighting.\n"
			"-mfsweighting\n"
			"   In spectral mode, calculate the weights as if the image was made using MFS. This makes sure that the sum of\n"
			"   channel images equals the MFS weights. Otherwise, the channel image will become a bit more naturally weighted.\n"
			"   This is only relevant for weighting modes that require gridding (i.e., Uniform, Briggs').\n"
			"   Default: off, unless -joinchannels is specified.\n"
			"-nomfsweighting\n"
			"   Opposite of -mfsweighting; can be used to turn off MFS weighting in -joinchannels mode.\n"
			"-weighting-rank-filter <level>\n"
			"   Filter the weights and set high weights to the local mean. The level parameter specifies\n"
			"   the filter level; any value larger than level*localmean will be set to level*localmean.\n"
			"-weighting-rank-filter-size <size>\n"
			"   Set size of weighting rank filter. Default: 16.\n"
			"-gridmode <nn or kb>\n"
			"   Kernel and mode used for gridding: kb = Kaiser-Bessel (default with 7 pixels), nn = nearest\n"
			"   neighbour (no kernel). Default: kb.\n"
			"-gkernelsize <size>\n"
			"   Gridding antialiasing kernel size. Default: 7.\n"
			"-oversampling <factor>\n"
			"   Oversampling factor used during gridding. Default: 63.\n"
			"-makepsf\n"
			"   Always make the psf, even when no cleaning is performed.\n"
			"-savegridding\n"
			"   Save the gridding correction image. This shows the effect of the antialiasing filter. Default: not saved.\n"
			"-dft-prediction\n"
			"   Predict via a direct Fourier transform. This is slow, but can account for direction-dependent effects. This has\n"
			"   only effect when -mgain is set or -predict is given.\n"
			"-dft-with-beam\n"
			"   Apply the beam during DFT. Currently only works for LOFAR.\n"
			"-visibility-weighting-mode [normal/squared/unit]\n"
			"   Specify visibility weighting modi. Affects how the weights (normally) stored in\n"
			"   WEIGHT_SPECTRUM column are applied. Useful for estimating e.g. EoR power spectra errors.\n"
			"   Normally one would use this in combination with -no-normalize-for-weighting.\n"
			"-no-normalize-for-weighting\n"
			"   Disable the normalization for the weights, which makes the PSF's peak one. See\n"
			"   -visibility-weighting-mode. Only useful with natural weighting.\n"
			"\n"
			"  ** DATA SELECTION OPTIONS **\n"
			"-pol <list>\n"
			"   Default: \'I\'. Possible values: XX, XY, YX, YY, I, Q, U, V, RR, RL, LR or LL (case insensitive).\n"
			"   Multiple values can be separated with commas, e.g.: 'xx,xy,yx,yy'. Two or four polarizations can be\n"
			"   joinedly cleaned (see '-joinpolarizations'), but this is not the default. I, Q, U and V\n"
			"   polarizations will be directly calculated from the visibilities, which is not appropriate for\n"
			"   telescopes with non-orthogonal feeds, such as MWA and LOFAR. The 'xy' polarization will output both\n"
			"   a real and an imaginary image, which allows calculating true Stokes polarizations for those\n"
			"   telescopes.\n"
			"-interval <start-index> <end-index>\n"
			"   Only image the given time interval. Indices specify the timesteps, end index is exclusive.\n"
			"   Default: image all time steps.\n"
			"-intervalsout <count>\n"
			"   Number of intervals to image inside the selected global interval. Default: 1\n"
			"-channelrange <start-channel> <end-channel>\n"
			"   Only image the given channel range. Indices specify channel indices, end index is exclusive.\n"
			"   Default: image all channels.\n"
			"-field <fieldid>\n"
			"   Image the given field id. Default: first field (id 0).\n"
			"-datacolumn <columnname>\n"
			"   Default: CORRECTED_DATA if it exists, otherwise DATA will be used.\n"
			"-maxuvw-m <meters>\n"
			"-minuvw-m <meters>\n"
			"   Set the min/max baseline distance in meters.\n"
			"-maxuv-l <lambda>\n"
			"-minuv-l <lambda>\n"
			"   Set the min/max uv distance in lambda.\n"
			"-maxw <percentage>\n"
			"   Do not grid visibilities with a w-value higher than the given percentage of the max w, to save speed.\n"
			"   Default: grid everything\n"
			"\n"
			"  ** DECONVOLUTION OPTIONS **\n"
			"-niter <niter>\n"
			"   Maximum number of clean iterations to perform. Default: 0\n"
			"-threshold <threshold>\n"
			"   Stopping clean thresholding in Jy. Default: 0.0\n"
			"-gain <gain>\n"
			"   Cleaning gain: Ratio of peak that will be subtracted in each iteration. Default: 0.1\n"
			"-mgain <gain>\n"
			"   Cleaning gain for major iterations: Ratio of peak that will be subtracted in each major\n"
			"   iteration. To use major iterations, 0.85 is a good value. Default: 1.0\n"
			"-joinpolarizations\n"
			"   Perform cleaning by searching for peaks in the sum of squares of the polarizations (either\n"
			"   I^2+Q^2+U^2+V^2, XX^2+real(XY)^2+imag(XY)^2+YY^2, or XX^2+YY^2), but subtract components from\n"
			"   individual channels. Only possible when imaging all Stokes or all linear parameters. Default: off.\n"
			"-joinchannels\n"
			"   Perform cleaning by searching for peaks in the MFS image, but subtract components from individual channels.\n"
			"   This will turn on mfsweighting by default. Default: off.\n"
			"-multiscale\n"
			"   Clean on different scales. This is a new experimental algorithm. Default: off.\n"
			"   This parameter invokes the v1.9 multiscale algorithm, which is slower but more accurate\n"
			"   compared to the older algorithm, and therefore the recommended one to use.\n"
			"   The older algorithm is now invoked with -fast-multiscale.\n"
			"-fast-multiscale\n"
			"   Clean on different scales. This is a new fast experimental algorithm. Default: off.\n"
			"   This method used to be invoked with -multiscale before v1.9, but the newer multiscale\n"
			"   algorithm is somewhat more accurate and therefore recommended.\n"
			"-multiscale-threshold-bias\n"
			"   Parameter to lower the threshold for larger scales. The used threshold for a scale\n"
			"   is threshold(scale)=pointsource_threshold x tbias^scale. A lower bias will clean\n"
			"   larger scales deeper. Default: 0.7\n"
			"-multiscale-scale-bias\n"
			"   Parameter to prevent cleaning small scales in the large-scale iterations. A higher\n"
			"   bias will give more focus to larger scales. Default: 0.6\n"
			"-iuwt\n"
			"   Use the IUWT deconvolution algorithm.\n"
			"-moresane-ext <location>\n"
			"   Use the MoreSane deconvolution algorithm, installed at the specified location.\n"
			"-moresane-arg <arguments>\n"
			"   Pass the specified arguments to moresane. Note that multiple parameters have to be\n"
			"   enclosed in quotes.\n"
			"-moresane-sl <sl1,sl2,...>\n"
			"   MoreSane --sigmalevel setting for each major loop iteration. Useful to start at high\n"
			"   levels and go down with subsequent loops, e.g. 20,10,5\n"
			"-cleanborder <percentage>\n"
			"   Set the border size in which no cleaning is performed, in percentage of the width/height of the image.\n"
			"   With an image size of 1000 and clean border of 1%, each border is 10 pixels. Default: 5 (%).\n"
 			"-fitsmask <mask>\n"
			"   Use the specified fits-file as mask during cleaning.\n"
			"-casamask <mask>\n"
			"   Use the specified CASA mask as mask during cleaning.\n"
			"-smallpsf\n"
			"   Resize the psf to speed up minor clean iterations. Not the default.\n"
			"-nonegative\n"
			"   Do not allow negative components during cleaning. Not the default.\n"
			"-negative\n"
			"   Default on: opposite of -nonegative.\n"
			"-stopnegative\n"
			"   Stop on negative components. Not the default.\n"
			"-fit-spectral-pol <nterms>\n"
			"   Fit a polynomial over frequency to each clean component. This has only effect\n"
			"   when the channels are joined with -joinchannels.\n"
			"-fit-spectral-log-pol <nterms>\n"
			"   Like fit-spectral-pol, but fits a logarithmic polynomial over frequency instead.\n"
			"-deconvolution-channels <nchannels>\n"
			"   Decrease the number of channels as specified by -channelsout to the given number for\n"
			"   deconvolution. Only possible in combination with one of the -fit-spectral options.\n"
			"\n"
			"  ** RESTORATION OPTIONS **\n"
			"-beamsize <arcsec>\n"
			"   Set the FWHM beam size in arcsec for restoring the clean components. Default: longest projected\n"
			"   baseline defines restoring beam.\n"
			"-beamshape <maj in arcsec> <min in arcsec> <position angle in deg>\n"
			"   Set the FWHM beam shape for restoring the clean components. Defaults units for maj and min are arcsec, and\n"
			"   degrees for PA. Can be overriden, e.g. '-beamshape 1amin 1amin 3deg'.\n"
			"-fitbeam\n"
			"   Determine beam shape by fitting the PSF (default if PSF is made).\n"
			"-nofitbeam\n"
			"   Do not determine beam shape from the PSF.\n"
			"-theoreticbeam\n"
			"   Write the beam in output fits files as calculated from the longest projected baseline.\n"
			"   This method results in slightly less accurate integrated fluxes, but in simple imaging provide\n"
			"   a beam size even without making the PSF. Default: off.\n"
			"-circularbeam\n"
			"   Force the beam to be circular: bmin will be set to bmaj.\n"
			"-ellipticalbeam\n"
			"   Allow the beam to be elliptical. Default.\n";
		return -1;
	}
	
	WSClean wsclean;
	int argi = 1;
	bool mfsWeighting = false, noMFSWeighting = false, predictionMode = false;
	while(argi < argc && argv[argi][0] == '-')
	{
		const std::string param = argv[argi][1]=='-' ? (&argv[argi][2]) : (&argv[argi][1]);
		if(param == "version")
		{
			// version already printed: just exit
			return 0;
		}
		else if(param == "tempdir")
		{
			++argi;
			wsclean.SetTemporaryDirectory(argv[argi]);
		}
		else if(param == "saveweights")
		{
			wsclean.SetSaveWeights(true);
		}
		else if(param == "saveuv")
		{
			wsclean.SetSaveUV(true);
		}
		else if(param == "predict")
		{
			predictionMode = true;
		}
		else if(param == "subtract-model")
		{
			wsclean.SetSubtractModel(true);
		}
		else if(param == "size")
		{
			size_t
				width = atoi(argv[argi+1]),
				height = atoi(argv[argi+2]);
			if(width != height)
				throw std::runtime_error("width != height : Can't handle non-square images yet");
			wsclean.SetImageSize(width, height);
			argi += 2;
		}
		else if(param == "trim")
		{
			size_t
				width = atoi(argv[argi+1]),
				height = atoi(argv[argi+2]);
			wsclean.SetTrimmedImageSize(width, height);
			argi += 2;
		}
		else if(param == "scale")
		{
			++argi;
			wsclean.SetPixelScale(Angle::Parse(argv[argi], "scale parameter", Angle::Degrees));
		}
		else if(param == "nwlayers")
		{
			++argi;
			wsclean.SetNWlayers(atoi(argv[argi]));
		}
		else if(param == "gain")
		{
			++argi;
			wsclean.DeconvolutionInfo().SetGain(atof(argv[argi]));
		}
		else if(param == "mgain")
		{
			++argi;
			wsclean.DeconvolutionInfo().SetMGain(atof(argv[argi]));
		}
		else if(param == "niter")
		{
			++argi;
			wsclean.DeconvolutionInfo().SetNIter(atoi(argv[argi]));
		}
		else if(param == "threshold")
		{
			++argi;
			wsclean.DeconvolutionInfo().SetThreshold(atof(argv[argi]));
		}
		else if(param == "datacolumn")
		{
			++argi;
			wsclean.SetColumnName(argv[argi]);
		}
		else if(param == "pol")
		{
			++argi;
			wsclean.SetPolarizations(Polarization::ParseList(argv[argi]));
		}
		else if(param == "negative")
		{
			wsclean.DeconvolutionInfo().SetAllowNegativeComponents(true);
		}
		else if(param == "nonegative")
		{
			wsclean.DeconvolutionInfo().SetAllowNegativeComponents(false);
		}
		else if(param == "stopnegative")
		{
			wsclean.DeconvolutionInfo().SetStopOnNegativeComponents(true);
		}
		else if(param == "iuwt")
		{
			wsclean.DeconvolutionInfo().SetUseIUWT(true);
			// Currently (WSClean 1.9, 2015-08-19) IUWT deconvolution
			// seems not to work when allowing negative components. The algorithm
			// becomes unstable. Hence, turn negative components off.
			wsclean.DeconvolutionInfo().SetAllowNegativeComponents(false);
		}
		else if(param == "moresane-ext")
		{
			++argi;
			wsclean.DeconvolutionInfo().SetUseMoreSane(true);
			wsclean.DeconvolutionInfo().SetMoreSaneLocation(argv[argi]);
		}
		else if(param == "moresane-arg")
		{
			++argi;
			wsclean.DeconvolutionInfo().SetMoreSaneArgs(argv[argi]);
		}
		else if(param == "moresane-sl")
		{
			++argi;
			std::vector<std::string> slevels;
			boost::split(slevels, argv[argi], boost::is_any_of(","));
			wsclean.DeconvolutionInfo().SetMoreSaneSigmaLevels(slevels);
		}
		else if(param == "makepsf")
		{
			wsclean.SetMakePSF(true);
		}
		else if(param == "savegridding")
		{
			wsclean.SetSaveGriddingImage(true);
		}
		else if(param == "dft-prediction")
		{
			wsclean.SetDFTPrediction(true);
		}
		else if(param == "dft-with-beam")
		{
			wsclean.SetDFTWithBeam(true);
		}
		else if(param == "cleanareas")
		{
			throw std::runtime_error("Clean areas is not supported ATM");
			//++argi;
			//wsclean.SetCleanAreasFilename(argv[argi]);
		}
		else if(param == "name")
		{
			++argi;
			wsclean.SetPrefixName(argv[argi]);
			wsclean.DeconvolutionInfo().SetPrefixName(argv[argi]);
		}
		else if(param == "gridmode")
		{
			++argi;
			std::string gridModeStr = argv[argi];
			boost::to_lower(gridModeStr);
			if(gridModeStr == "kb" || gridModeStr == "kaiserbessel" || gridModeStr == "kaiser-bessel")
				wsclean.SetGridMode(WStackingGridder::KaiserBessel);
			else if(gridModeStr == "nn" || gridModeStr == "nearestneighbour")
				wsclean.SetGridMode(WStackingGridder::NearestNeighbour);
			else
				throw std::runtime_error("Invalid gridding mode: should be either kb (Kaiser-Bessel) or nn (NearestNeighbour)");
		}
		else if(param == "smallinversion")
		{
			wsclean.SetSmallInversion(true);
		}
		else if(param == "nosmallinversion")
		{
			wsclean.SetSmallInversion(false);
		}
		else if(param == "interval")
		{
			wsclean.SetIntervalSelection(atoi(argv[argi+1]), atoi(argv[argi+2]));
			argi += 2;
		}
		else if(param == "intervalsout")
		{
			++argi;
			wsclean.SetIntervalCount(atoi(argv[argi]));
		}
		else if(param == "channelrange")
		{
			wsclean.SetChannelSelection(atoi(argv[argi+1]), atoi(argv[argi+2]));
			argi += 2;
		}
		else if(param == "channelsout")
		{
			++argi;
			wsclean.SetChannelsOut(atoi(argv[argi]));
		}
		else if(param == "joinpolarizations")
		{
			wsclean.SetJoinPolarizations(true);
		}
		else if(param == "joinchannels")
		{
			wsclean.SetJoinChannels(true);
		}
		else if(param == "mfsweighting")
		{
			mfsWeighting = true;
		}
		else if(param == "multiscale")
		{
			wsclean.DeconvolutionInfo().SetMultiscale(true);
		}
		else if(param == "fast-multiscale")
		{
			wsclean.DeconvolutionInfo().SetFastMultiscale(true);
		}
		else if(param == "multiscale-threshold-bias")
		{
			++argi;
			wsclean.DeconvolutionInfo().SetMultiscaleThresholdBias(atof(argv[argi]));
		}
		else if(param == "multiscale-scale-bias")
		{
			++argi;
			wsclean.DeconvolutionInfo().SetMultiscaleScaleBias(atof(argv[argi]));
		}
		else if(param == "weighting-rank-filter")
		{
			++argi;
			wsclean.SetRankFilterLevel(atof(argv[argi]));
		}
		else if(param == "weighting-rank-filter-size")
		{
			++argi;
			wsclean.SetRankFilterSize(atoi(argv[argi]));
		}
		else if(param == "cleanborder")
		{
			++argi;
			wsclean.DeconvolutionInfo().SetCleanBorderRatio(atof(argv[argi])*0.01);
		}
		else if(param == "fitsmask")
		{
			++argi;
			wsclean.DeconvolutionInfo().SetFitsMask(argv[argi]);
		}
		else if(param == "casamask")
		{
			++argi;
			wsclean.DeconvolutionInfo().SetCASAMask(argv[argi]);
		}
		else if(param == "nomfsweighting")
		{
			noMFSWeighting = true;
		}
		else if(param == "joinchannels")
		{
			wsclean.SetJoinChannels(true);
		}
		else if(param == "fit-spectral-pol")
		{
			++argi;
			wsclean.DeconvolutionInfo().SetFitSpectralPol(atoi(argv[argi]));
		}
		else if(param == "fit-spectral-log-pol")
		{
			++argi;
			wsclean.DeconvolutionInfo().SetFitSpectralLogPol(atoi(argv[argi]));
		}
		else if(param == "deconvolution-channels")
		{
			++argi;
			wsclean.DeconvolutionInfo().SetDeconvolutionChannels(atoi(argv[argi]));
		}
		else if(param == "field")
		{
			++argi;
			wsclean.SetFieldSelection(atoi(argv[argi]));
		}
		else if(param == "weight")
		{
			++argi;
			std::string weightArg = argv[argi];
			if(weightArg == "natural")
				wsclean.SetWeightMode(WeightMode::NaturalWeighted);
			else if(weightArg == "mwa")
				wsclean.SetWeightMode(WeightMode::DistanceWeighted);
			else if(weightArg == "uniform")
				wsclean.SetWeightMode(WeightMode::UniformWeighted);
			else if(weightArg == "briggs")
			{
				++argi;
				wsclean.SetBriggsWeighting(atof(argv[argi]));
			}
			else throw std::runtime_error("Unknown weighting mode specified");
		}
		else if(param == "superweight")
		{
			++argi;
			wsclean.SetSuperWeight(atof(argv[argi]));
		}
		else if(param == "beamsize")
		{
			++argi;
			double beam = Angle::Parse(argv[argi], "beam size", Angle::Arcseconds);
			wsclean.SetBeamSize(beam, beam, 0.0);
		}
		else if(param == "beamshape")
		{
			double beamMaj = Angle::Parse(argv[argi+1], "beam shape, major axis", Angle::Arcseconds);
			double beamMin = Angle::Parse(argv[argi+2], "beam shape, minor axis", Angle::Arcseconds);
			double beamPA = Angle::Parse(argv[argi+3], "beam shape, position angle", Angle::Degrees);
			argi+=3;
			wsclean.SetBeamSize(beamMaj, beamMin, beamPA);
		}
		else if(param == "fitbeam")
		{
			wsclean.SetFittedBeam(true);
		}
		else if(param == "nofitbeam")
		{
			wsclean.SetFittedBeam(false);
		}
		else if(param == "theoreticbeam")
		{
			wsclean.SetTheoreticBeam(true);
			wsclean.SetFittedBeam(false);
		}
		else if(param == "circularbeam")
		{
			wsclean.SetCircularBeam(true);
		}
		else if(param == "ellipticalbeam")
		{
			wsclean.SetCircularBeam(false);
		}
		else if(param == "gkernelsize")
		{
			++argi;
			wsclean.SetAntialiasingKernelSize(atoi(argv[argi]));
		}
		else if(param == "oversampling")
		{
			++argi;
			wsclean.SetOversamplingFactor(atoi(argv[argi]));
		}
		else if(param == "reorder")
		{
			wsclean.SetForceReorder(true);
			wsclean.SetForceNoReorder(false);
		}
		else if(param == "no-reorder")
		{
			wsclean.SetForceNoReorder(true);
			wsclean.SetForceReorder(false);
		}
		else if(param == "update-model-required")
		{
			wsclean.SetModelUpdateRequired(true);
		}
		else if(param == "no-update-model-required")
		{
			wsclean.SetModelUpdateRequired(false);
		}
		else if(param == "j")
		{
			++argi;
			wsclean.SetThreadCount(atoi(argv[argi]));
		}
		else if(param == "mem")
		{
			++argi;
			wsclean.SetMemFraction(atof(argv[argi]) / 100.0);
		}
		else if(param == "absmem")
		{
			++argi;
			wsclean.SetMemAbsLimit(atof(argv[argi]));
		}
		else if(param == "maxuvw-m")
		{
			++argi;
			wsclean.SetMaxUVWInM(atof(argv[argi]));
		}
		else if(param == "minuvw-m")
		{
			++argi;
			wsclean.SetMinUVWInM(atof(argv[argi]));
		}
		else if(param == "maxuv-l")
		{
			++argi;
			wsclean.SetMaxUVInLambda(atof(argv[argi]));
		}
		else if(param == "minuv-l")
		{
			++argi;
			wsclean.SetMinUVInLambda(atof(argv[argi]));
		}
		else if(param == "maxw")
		{
			// This was to test the optimization suggested in Tasse et al., 2013, Appendix C.
			++argi;
			wsclean.SetWLimit(atof(argv[argi]));
		}
		else if(param == "no-normalize-for-weighting")
		{
			wsclean.SetNormalizeForWeighting(false);
		}
		else if(param == "visibility-weighting-mode")
		{
			++argi;
			std::string modeStr = argv[argi];
			boost::to_lower(modeStr);
			if(modeStr == "normal")
				wsclean.SetVisibilityWeightingMode(InversionAlgorithm::NormalVisibilityWeighting);
			else if(modeStr == "squared")
				wsclean.SetVisibilityWeightingMode(InversionAlgorithm::SquaredVisibilityWeighting);
			else if(modeStr == "unit")
				wsclean.SetVisibilityWeightingMode(InversionAlgorithm::UnitVisibilityWeighting);
			else
				throw std::runtime_error("Unknown weighting mode: " + modeStr);
		}
		else {
			throw std::runtime_error("Unknown parameter: " + param);
		}
		
		++argi;
	}
	
	if(argi == argc)
		throw std::runtime_error("No input measurement sets given.");
	
	wsclean.SetMFSWeighting((wsclean.JoinChannels() && !noMFSWeighting) || mfsWeighting);
	
	for(int i=argi; i != argc; ++i)
		wsclean.AddInputMS(argv[i]);
	
	std::ostringstream commandLineStr;
	commandLineStr << "wsclean";
	for(int i=1; i!=argc; ++i)
		commandLineStr << ' ' << argv[i];
	wsclean.SetCommandLine(commandLineStr.str());
	
	if(predictionMode)
		wsclean.RunPredict();
	else
		wsclean.RunClean();
	return 0;
}
