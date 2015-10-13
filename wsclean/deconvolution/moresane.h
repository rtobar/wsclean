#ifndef MORESANE_H
#define MORESANE_H

#include <string>

#include "deconvolutionalgorithm.h"
#include "imageset.h"

class MoreSane : public UntypedDeconvolutionAlgorithm
{
	public:
		MoreSane(const std::string& moreSaneLocation, const std::string& moresaneArguments, 
		         const std::vector<std::string> &moresaneSigmaLevels, const std::string &prefixName,
						 ImageBufferAllocator& allocator ) 
                : _moresaneLocation(moreSaneLocation), _moresaneArguments(moresaneArguments), 
                _moresaneSigmaLevels(moresaneSigmaLevels), _prefixName(prefixName),
                _allocator(&allocator)
		{ }
		
		virtual void ExecuteMajorIteration(DynamicSet& dataImage, DynamicSet& modelImage, const ao::uvector<const double*>& psfImages, size_t width, size_t height, bool& reachedMajorThreshold);
		
		void ExecuteMajorIteration(double* dataImage, double* modelImage, const double* psfImage, size_t width, size_t height, bool& reachedMajorThreshold);
	private:
		const std::string _moresaneLocation, _moresaneArguments;

		const std::vector<std::string> _moresaneSigmaLevels;
		const std::string _prefixName;
		
		ImageBufferAllocator* _allocator;
};

#endif
