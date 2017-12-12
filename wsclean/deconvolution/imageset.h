#ifndef IMAGE_SET_H
#define IMAGE_SET_H

#include "../uvector.h"
#include "../wsclean/imagingtable.h"
#include "../wsclean/imagebufferallocator.h"

#include <vector>
#include <map>

class ImageSet
{
public:
	ImageSet(const ImagingTable* table, ImageBufferAllocator& allocator, size_t requestedChannelsInDeconvolution, bool squareJoinedChannels, const std::set<PolarizationEnum>& linkedPolarizations) :
		_images(),
		_imageSize(0),
		_channelsInDeconvolution((requestedChannelsInDeconvolution==0) ? table->SquaredGroupCount() : requestedChannelsInDeconvolution),
		_squareJoinedChannels(squareJoinedChannels),
		_imagingTable(*table),
		_imageIndexToPSFIndex(),
		_linkedPolarizations(linkedPolarizations),
		_allocator(allocator)
	{
		size_t nPol = table->GetSquaredGroup(0).EntryCount();
		size_t nImages = nPol * _channelsInDeconvolution;
		_images.assign(nImages, static_cast<double*>(0));
		_imageIndexToPSFIndex.resize(nImages);
		
		initializePolFactor();
		initializeIndices();
	}
	
	ImageSet(const ImagingTable* table, ImageBufferAllocator& allocator, size_t requestedChannelsInDeconvolution, bool squareJoinedChannels, const std::set<PolarizationEnum>& linkedPolarizations, size_t width, size_t height) :
		_images(),
		_imageSize(width*height),
		_channelsInDeconvolution((requestedChannelsInDeconvolution==0) ? table->SquaredGroupCount() : requestedChannelsInDeconvolution),
		_squareJoinedChannels(squareJoinedChannels),
		_imagingTable(*table),
		_imageIndexToPSFIndex(),
		_linkedPolarizations(linkedPolarizations),
		_allocator(allocator)
	{
		size_t nPol = table->GetSquaredGroup(0).EntryCount();
		size_t nImages = nPol * _channelsInDeconvolution;
		_images.assign(nImages, static_cast<double*>(0));
		_imageIndexToPSFIndex.resize(nImages);
		
		initializePolFactor();
		initializeIndices();
		allocateImages();
	}
	
	~ImageSet()
	{
		free();
	}
	
	void AllocateImages()
	{
		free();
		allocateImages();
	}
	
	void AllocateImages(size_t width, size_t height)
	{
		free();
		_imageSize = width*height;
		allocateImages();
	}
	
	double* Release(size_t imageIndex)
	{
		double* image = _images[imageIndex];
		_images[imageIndex] = 0;
		return image;
	}
	
	void Claim(size_t imageIndex, double* data)
	{
		_allocator.Free(_images[imageIndex]);
		_images[imageIndex] = data;
	}
	
	bool IsAllocated() const
	{
		return _imageSize!=0;
	}
	
	ImageBufferAllocator& Allocator() const
	{
		return _allocator;
	}
	
	void LoadAndAverage(class CachedImageSet& imageSet);
	
	void LoadAndAveragePSFs(class CachedImageSet& psfSet, std::vector<ao::uvector<double>>& psfImages, PolarizationEnum psfPolarization);
	
	void InterpolateAndStore(class CachedImageSet& imageSet, const class SpectralFitter& fitter);
	
	void AssignAndStore(class CachedImageSet& imageSet);
	
	/**
	 * This function will calculate the integration over all images, squaring
	 * images that are in the same squared-image group. For example, with
	 * a squared group of [I, Q, ..] and another group [I2, Q2, ...], this
	 * will calculate:
	 * 
	 * sqrt(I^2 + Q^2 + ..) + sqrt(I2^2 + Q2^2 ..) + ..
	 * ----------------------------------------------
	 *           1          +           1          + ..
	 * 
	 * If the 'squared groups' are of size 1, the average of the groups will be
	 * returned (i.e., without square-rooting the square).
	 * 
	 * If the squared joining option is set in the provided wsclean settings, the
	 * behaviour of this method changes. In that case, it will return the square
	 * root of the average squared value:
	 * 
	 *       I^2 + Q^2 + ..  +  I2^2 + Q2^2 ..  + ..
	 * sqrt( --------------------------------------- )
	 *            1          +        1         + ..
	 * 
	 * These formulae are such that the values will have normal flux values.
	 * @param dest Pre-allocated output array that will be filled with the
	 * integrated image.
	 * @param scratch Pre-allocated scratch space, same size as image.
	 */
	void GetSquareIntegrated(double* dest, double* scratch) const
	{
		if(_squareJoinedChannels)
			getSquareIntegratedWithSquaredChannels(dest);
		else
			getSquareIntegratedWithNormalChannels(dest, scratch);
	}
	
	/**
	 * This function will calculate the 'linear' integration over all images.
	 * This will return the average of all images. Normally, @ref GetSquareIntegrated
	 * should be used for peak finding, but in case negative values should remain
	 * negative, such as with multiscale (otherwise a sidelobe will be fitted with
	 * large scales), this function can be used.
	 * @param dest Pre-allocated output array that will be filled with the average
	 * values.
	 */
	void GetLinearIntegrated(double* dest) const
	{
		if(_squareJoinedChannels)
			getSquareIntegratedWithSquaredChannels(dest);
		else
			getLinearIntegratedWithNormalChannels(dest);
	}

	void GetIntegratedPSF(double* dest, const ao::uvector<const double*>& psfs)
	{
		memcpy(dest, psfs[0], sizeof(double) * _imageSize);
		for(size_t i = 1; i!=PSFCount(); ++i)
		{
			add(dest, psfs[i]);
		}
		multiply(dest, 1.0/double(PSFCount()));
	}
	
	size_t PSFCount() const { return _channelsInDeconvolution; }
	
	size_t ChannelsInDeconvolution() const { return _channelsInDeconvolution; }
	
	ImageSet& operator=(double val)
	{
		for(size_t i=0; i!=size(); ++i)
			assign(_images[i], val);
		return *this;
	}
	
	double* operator[](size_t index)
	{
		return _images[index];
	}
	
	const double* operator[](size_t index) const
	{
		return _images[index];
	}
	
	size_t size() const { return _images.size(); }
	
	size_t PSFIndex(size_t imageIndex) const { return _imageIndexToPSFIndex[imageIndex]; }
	
	const ImagingTable& Table() const { return _imagingTable; }
	
	ImageSet* CreateTrimmed(size_t x1, size_t y1, size_t x2, size_t y2, size_t oldWidth) const
	{
		std::unique_ptr<ImageSet> p(new ImageSet(&_imagingTable, _allocator, _channelsInDeconvolution, _squareJoinedChannels, _linkedPolarizations, x2-x1, y2-y1));
		for(size_t i=0; i!=_images.size(); ++i)
		{
			copySmallerPart(_images[i], p->_images[i], x1, y1, x2, y2, oldWidth);
		}
		return p.release();
	}
	
	ImageSet& operator*=(double factor)
	{
		for(size_t i=0; i!=size(); ++i)
			multiply(_images[i], factor);
		return *this;
	}
	
	ImageSet& operator+=(const ImageSet& other)
	{
		for(size_t i=0; i!=size(); ++i)
			add(_images[i], other._images[i]);
		return *this;
	}
	
	void FactorAdd(ImageSet& rhs, double factor)
	{
		for(size_t i=0; i!=size(); ++i)
			addFactor(_images[i], rhs._images[i], factor);
	}
	
	void Set(size_t index, const double* rhs)
	{
		assign(_images[index], rhs);
	}
	
	bool SquareJoinedChannels() const
	{
		return _squareJoinedChannels; 
	}
	
	const std::set<PolarizationEnum>& LinkedPolarizations() const
	{
		return _linkedPolarizations;
	}
private:
	ImageSet(const ImageSet&) = delete;
	ImageSet& operator=(const ImageSet&) = delete;
		
	void allocateImages()
	{
		for(ao::uvector<double*>::iterator img=_images.begin();
				img!=_images.end(); ++img)
		{
			*img = _allocator.Allocate(_imageSize);
		}
	}
	
	void free()
	{
		for(ao::uvector<double*>::iterator img=_images.begin();
				img!=_images.end(); ++img)
			_allocator.Free(*img);
	}
		
	void assign(double* lhs, const double* rhs) const
	{
		memcpy(lhs, rhs, sizeof(double) * _imageSize);
	}
	
	void assign(double* lhs, const ImageBufferAllocator::Ptr& rhs) const
	{
		memcpy(lhs, rhs.data(), sizeof(double) * _imageSize);
	}
	
	void assignMultiply(double* lhs, const double* rhs, double factor) const
	{
		for(size_t i=0; i!=_imageSize; ++i)
			lhs[i] = rhs[i] * factor;
	}
	
	void assign(double* image, double value) const
	{
		for(size_t i=0; i!=_imageSize; ++i)
			image[i] = value;
	}
	
	void add(double* lhs, const double* rhs) const
	{
		for(size_t i=0; i!=_imageSize; ++i)
			lhs[i] += rhs[i];
	}
	
	void square(double* image) const
	{
		for(size_t i=0; i!=_imageSize; ++i)
			image[i] *= image[i];
	}
	
	void squareRoot(double* image) const
	{
		for(size_t i=0; i!=_imageSize; ++i)
			image[i] = sqrt(image[i]);
	}
	
	void squareRootMultiply(double* image, double factor) const
	{
		for(size_t i=0; i!=_imageSize; ++i)
			image[i] = sqrt(image[i]) * factor;
	}
	
	void addSquared(double* lhs, const double* rhs) const
	{
		for(size_t i=0; i!=_imageSize; ++i)
			lhs[i] += rhs[i]*rhs[i];
	}
	
	void addFactor(double* lhs, const double* rhs, double factor) const
	{
		for(size_t i=0; i!=_imageSize; ++i)
			lhs[i] += rhs[i] * factor;
	}
	
	void multiply(double* image, double fact) const
	{
		if(fact != 1.0)
		{
			for(size_t i=0; i!=_imageSize; ++i)
				image[i] *= fact;
		}
	}
	
	void initializeIndices()
	{
		for(size_t i=0; i!=_imagingTable.EntryCount(); ++i)
		{
			_tableIndexToImageIndex.insert(
				std::make_pair(_imagingTable[i].index, i));
		}
		for(size_t sqIndex = 0; sqIndex!=_channelsInDeconvolution; ++sqIndex)
		{
			ImagingTable subTable = _imagingTable.GetSquaredGroup(sqIndex);
			for(size_t eIndex = 0; eIndex!=subTable.EntryCount(); ++eIndex)
			{
				const ImagingTableEntry& entry = subTable[eIndex];
				size_t imageIndex = _tableIndexToImageIndex.find(entry.index)->second;
				_imageIndexToPSFIndex[imageIndex] = sqIndex;
			}
		}
	}
	
	void initializePolFactor()
	{
		ImagingTable firstChannelGroup = _imagingTable.GetSquaredGroup(0);
		std::set<PolarizationEnum> pols;
		for(size_t i=0; i!=firstChannelGroup.EntryCount(); ++i)
		{
			if(_linkedPolarizations.empty()
				|| _linkedPolarizations.count(firstChannelGroup[i].polarization) !=0 )
			{
				pols.insert(firstChannelGroup[i].polarization);
			}
		}
		bool isDual = pols.size()==2 && Polarization::HasDualPolarization(pols);
		bool isFull = pols.size()==4 && (
			Polarization::HasFullLinearPolarization(pols) ||
			Polarization::HasFullCircularPolarization(pols));
		if(isDual || isFull)
			_polarizationNormalizationFactor = 0.5;
		else
			_polarizationNormalizationFactor = 1.0;
	}
	
	void copySmallerPart(const double* input, double* output, size_t x1, size_t y1, size_t x2, size_t y2, size_t oldWidth) const
	{
		size_t newWidth = x2 - x1;
		for(size_t y=y1; y!=y2; ++y)
		{
			const double* oldPtr = &input[y*oldWidth];
			double* newPtr = &output[(y-y1)*newWidth];
			for(size_t x=x1; x!=x2; ++x)
			{
				newPtr[x - x1] = oldPtr[x];
			}
		}
	}
	
	void directStore(class CachedImageSet& imageSet);
	
	void getSquareIntegratedWithNormalChannels(double* dest, double* scratch) const;
	
	void getSquareIntegratedWithSquaredChannels(double* dest) const;
	
	void getLinearIntegratedWithNormalChannels(double* dest) const;
	
	ao::uvector<double*> _images;
	size_t _imageSize, _channelsInDeconvolution;
	bool _squareJoinedChannels;
	const ImagingTable& _imagingTable;
	std::map<size_t, size_t> _tableIndexToImageIndex;
	ao::uvector<size_t> _imageIndexToPSFIndex;
	double _polarizationNormalizationFactor;
	std::set<PolarizationEnum> _linkedPolarizations;
	ImageBufferAllocator& _allocator;
};

#endif
