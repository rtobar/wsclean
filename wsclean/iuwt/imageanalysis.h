#ifndef IUWT_IMAGE_ANALYSIS_H
#define IUWT_IMAGE_ANALYSIS_H

#include "iuwtdecomposition.h"

class ImageAnalysis
{
public:
	struct Component
	{
		Component(size_t _x, size_t _y, int _scale) : x(_x), y(_y), scale(_scale) { }
		
		std::string ToString() const {
			std::ostringstream str;
			str << x << ',' << y << ", scale " << scale;
			return str.str();
		}
		
		size_t x, y;
		int scale;
	};

	struct Component2D
	{
		Component2D(size_t _x, size_t _y) : x(_x), y(_y) { }
		
		std::string ToString() const {
			std::ostringstream str;
			str << x << ',' << y;
			return str.str();
		}
		
		size_t x, y;
	};

	static void SelectStructures(const IUWTDecomposition& iuwt, IUWTMask& mask, const ao::uvector<double>& thresholds, size_t minScale, size_t endScale, double cleanBorder, const bool* priorMask, size_t& areaSize);
	
	static bool IsHighestOnScale0(const IUWTDecomposition& iuwt, IUWTMask& markedMask, size_t& x, size_t& y, size_t endScale, double& highestScale0);
	
	static void Floodfill(const IUWTDecomposition& iuwt, IUWTMask& mask, const ao::uvector<double>& thresholds, size_t minScale, size_t endScale, const Component& component, double cleanBorder, size_t& areaSize);
	
	static void MaskedFloodfill(const IUWTDecomposition& iuwt, IUWTMask& mask, const ao::uvector<double>& thresholds, size_t minScale, size_t endScale, const Component& component, double cleanBorder, const bool* priorMask, size_t& areaSize);
	
	static void FloodFill2D(const double* image, bool* mask, double threshold, const Component2D& component, size_t width, size_t height, size_t& areaSize);
	
	/**
	 * Exactly like above, but now collecting the components in the
	 * area vector, instead of returning just the area size.
	 */
	static void FloodFill2D(const double* image, bool* mask, double threshold, const Component2D& component, size_t width, size_t height, std::vector<Component2D>& area);
	
private:
	static bool exceedsThreshold(double val, double threshold)
	{
		if(threshold >= 0.0)
			return val > threshold;
		else
			return val < threshold || val > -threshold;
	}
	static bool exceedsThresholdAbs(double val, double threshold)
	{
		return std::fabs(val) > threshold;
	}
};

#endif
