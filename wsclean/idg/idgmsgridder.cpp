#include "idgmsgridder.h"

#include <idg-cpu.h>

IdgMsGridder::IdgMsGridder()
{
	/**
	 * Comments on idg:
	 * - gridsize; I was not sure this was width*height, or just one dimension. It seemes to be the latter.
	 * Also, dimensions should be specified separately.
	 * - imagesize is a confusing term; is it the real angular size of the image, or is it pixelsize * width().
	 * These become rather different for larger image sizes. Again, x and y pixelsizes/imagesizes should be
	 * independently be specifyable.
	 * - grid_size vs imagesize -> be consistent.
	 * - Proxy should have virtual destructor ~Proxy().
	 * - What's a job size?
	 * - Nr channels can't be a constant, since you can have the situation of different spectral windows, and or
	 * a selection of bands with different (selected) channels, multiple measurement sets with different nr of
	 * (selected) channels, or a "channelsout 2" setting with e.g. 3 measurement sets, which would require gridding
	 * all channels of the first measurement set, and then half the channels of the 2nd measurement set. How do
	 * I do this with this interface?
	 * - Nr stations -- similarly, how do I grid two measurement sets with different amount of stations?
	 * - Nr time / timeslots -- idem
	 * - What's the difference between a time and a timeslot?
	 */
	
	_gridSize = std::max(ImageWidth(), ImageHeight());
	_imageSize = std::max(ImageWidth() * PixelSizeX(), ImageHeight() * PixelSizeY());
	
	idg::Parameters parameters;
	parameters.set_grid_size(_gridSize);
	parameters.set_imagesize(_imageSize);
	//parameters.set_nr_channels(); // impossible to set!
	//parameters.set_nr_time();
	//parameters.get_nr_timeslots();
	//parameters.set_nr_stations();
	//_proxy = new idg::proxy::cpu::Reference(parameters);
	//_proxy->grid_visibilities();
}

IdgMsGridder::~IdgMsGridder()
{
	delete _proxy;
}

void IdgMsGridder::Invert()
{
}

void IdgMsGridder::Predict(double* real)
{
}

void IdgMsGridder::Predict(double* real, double* imaginary)
{
}

double* IdgMsGridder::ImageRealResult() const
{
}

double* IdgMsGridder::ImageImaginaryResult() const
{
}

double IdgMsGridder::ImageWeight() const
{
}

void IdgMsGridder::GetGriddingCorrectionImage(double* image) const
{
}

bool IdgMsGridder::HasGriddingCorrectionImage() const
{
	return false; // For now (TODO)
}

