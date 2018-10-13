#ifndef PARALLEL_DECONVOLUTION_H
#define PARALLEL_DECONVOLUTION_H

#include "../fftwmanager.h"
#include "../uvector.h"

#include <memory>
#include <mutex>
#include <vector>

class ParallelDeconvolution
{
public:
	ParallelDeconvolution(const class WSCleanSettings& settings) :
		_horImages(0),
		_verImages(0),
		_settings(settings),
		_allocator(nullptr),
		_mask(nullptr),
		_trackPerScaleMasks(false),
		_usePerScaleMasks(false)
	{ }
	
	class DeconvolutionAlgorithm& FirstAlgorithm()
	{
		return *_algorithms.front();
	}
	const class DeconvolutionAlgorithm& FirstAlgorithm() const
	{
		return *_algorithms.front();
	}
	
	void SetAllocator(class ImageBufferAllocator* allocator)
	{
		_allocator = allocator;
	}
	
	void SetAlgorithm(std::unique_ptr<class DeconvolutionAlgorithm> algorithm);
	
	void SetRMSFactorImage(class Image&& image);
	
	void SetThreshold(double threshold);
	
	bool IsInitialized() const
	{
		return !_algorithms.empty();
	}
	
	void SetAutoMaskMode(bool trackPerScaleMasks, bool usePerScaleMasks);
	
	void SetCleanMask(const bool* mask);
	
	void ExecuteMajorIteration(class ImageSet& dataImage, class ImageSet& modelImage, const ao::uvector<const double*>& psfImages, bool& reachedMajorThreshold);
	
	void FreeDeconvolutionAlgorithms()
	{
		_algorithms.clear(); 
		_mask = nullptr;
	}
	
	void SaveSourceList(class CachedImageSet& modelImages, const class ImagingTable& table, long double phaseCentreRA, long double phaseCentreDec);
	
	void SavePBSourceList(class CachedImageSet& modelImages, const class ImagingTable& table, long double phaseCentreRA, long double phaseCentreDec) const;
	
	class FFTWManager& GetFFTWManager() { return _fftwManager; }
	
private:
	struct SubImage {
		size_t index, x, y, width, height;
		ao::uvector<bool> mask;
		double peak;
		bool reachedMajorThreshold;
	};
		
	void runSubImage(SubImage& subImg, ImageSet& dataImage, class ImageSet& modelImage, const ao::uvector<const double*>& psfImages, double majorIterThreshold, bool findPeakOnly, std::mutex* mutex);
	
	void correctChannelForPB(class ComponentList& list, const class ImagingTableEntry& entry) const;
	
	void loadAveragePrimaryBeam(class PrimaryBeamImageSet& beamImages, size_t imageIndex, const class ImagingTable& table) const;
	
	FFTWManager _fftwManager;
	std::vector<std::unique_ptr<class DeconvolutionAlgorithm>> _algorithms;
	size_t _horImages, _verImages;
	const WSCleanSettings& _settings;
	ImageBufferAllocator* _allocator;
	const bool* _mask;
	bool _trackPerScaleMasks, _usePerScaleMasks;
	std::vector<ao::uvector<bool>> _scaleMasks;
};

#endif

