#include "imageanalysis.h"

#include <stack>

bool ImageAnalysis::IsHighestOnScale0(const IUWTDecomposition& iuwt, IUWTMask& markedMask, size_t& x, size_t& y, size_t endScale, double& highestScale0)
{
	const size_t width = iuwt.Width(), height = iuwt.Height();
	Component component(x, y, 0);
	std::stack<Component> todo;
	todo.push(component);
	markedMask[0][x + y*width] = false;
	highestScale0 = iuwt[0][x + y*width];
	double highestOtherScales = 0.0;
	while(!todo.empty())
	{
		Component c = todo.top();
		todo.pop();
		size_t index = c.x + c.y*width;
		if(c.scale == 0)
		{
			if(exceedsThreshold(iuwt[0][index], highestScale0))
			{
				highestScale0 = iuwt[0][index];
				x = c.x;
				y = c.y;
			}
		}
		else {
			if(exceedsThreshold(iuwt[c.scale][index], highestOtherScales))
				highestOtherScales = iuwt[c.scale][index];
		}
		if(c.x > 0)
		{
			if(markedMask[c.scale][index-1])
			{
				markedMask[c.scale][index-1] = false;
				todo.push(Component(c.x-1, c.y, c.scale));
			}
		}
		if(c.x < width-1)
		{
			if(markedMask[c.scale][index+1])
			{
				markedMask[c.scale][index+1] = false;
				todo.push(Component(c.x+1, c.y, c.scale));
			}
		}
		if(c.y > 0)
		{
			if(markedMask[c.scale][index-width])
			{
				markedMask[c.scale][index-width] = false;
				todo.push(Component(c.x, c.y-1, c.scale));
			}
		}
		if(c.y < height-1)
		{
			if(markedMask[c.scale][index+width])
			{
				markedMask[c.scale][index+width] = false;
				todo.push(Component(c.x, c.y+1, c.scale));
			}
		}
		if(c.scale > int(0))
		{
			if(markedMask[c.scale-1][index])
			{
				markedMask[c.scale-1][index] = false;
				todo.push(Component(c.x, c.y, c.scale-1));
			}
		}
		if(c.scale < int(endScale)-1)
		{
			if(markedMask[c.scale+1][index])
			{
				markedMask[c.scale+1][index] = false;
				todo.push(Component(c.x, c.y, c.scale+1));
			}
		}
	}
	return std::fabs(highestScale0) > std::fabs(highestOtherScales);
}

void ImageAnalysis::Floodfill(const IUWTDecomposition& iuwt, IUWTMask& mask, const ao::uvector<double>& thresholds, size_t minScale, size_t endScale, const Component& component, double cleanBorder, size_t& areaSize)
{
	const size_t width = iuwt.Width(), height = iuwt.Height();
	size_t xBorder = cleanBorder*width;
	size_t yBorder = cleanBorder*height;
	size_t minX = xBorder, maxX = width - xBorder;
	size_t minY = yBorder, maxY = height - yBorder;
	
	areaSize = 0;
	endScale = std::min<size_t>(endScale, iuwt.NScales());
	std::stack<Component> todo;
	todo.push(component);
	mask[component.scale][component.x + component.y*width] = true;
	while(!todo.empty())
	{
		Component c = todo.top();
		++areaSize;
		todo.pop();
		size_t index = c.x + c.y*width;
		if(c.x > minX)
		{
			if(exceedsThreshold(iuwt[c.scale][index-1], thresholds[c.scale]) && !mask[c.scale][index-1])
			{
				mask[c.scale][index-1] = true;
				todo.push(Component(c.x-1, c.y, c.scale));
			}
		}
		if(c.x < maxX-1)
		{
			if(exceedsThreshold(iuwt[c.scale][index+1], thresholds[c.scale]) && !mask[c.scale][index+1])
			{
				mask[c.scale][index+1] = true;
				todo.push(Component(c.x+1, c.y, c.scale));
			}
		}
		if(c.y > minY)
		{
			if(exceedsThreshold(iuwt[c.scale][index-width], thresholds[c.scale]) && !mask[c.scale][index-width])
			{
				mask[c.scale][index-width] = true;
				todo.push(Component(c.x, c.y-1, c.scale));
			}
		}
		if(c.y < maxY-1)
		{
			if(exceedsThreshold(iuwt[c.scale][index+width], thresholds[c.scale]) && !mask[c.scale][index+width])
			{
				mask[c.scale][index+width] = true;
				todo.push(Component(c.x, c.y+1, c.scale));
			}
		}
		if(c.scale > int(minScale))
		{
			if(exceedsThreshold(iuwt[c.scale-1][index], thresholds[c.scale-1]) && !mask[c.scale-1][index])
			{
				mask[c.scale-1][index] = true;
				todo.push(Component(c.x, c.y, c.scale-1));
			}
		}
		if(c.scale < int(endScale)-1)
		{
			if(exceedsThreshold(iuwt[c.scale+1][index], thresholds[c.scale+1]) && !mask[c.scale+1][index])
			{
				mask[c.scale+1][index] = true;
				todo.push(Component(c.x, c.y, c.scale+1));
			}
		}
	}
}

void ImageAnalysis::MaskedFloodfill(const IUWTDecomposition& iuwt, IUWTMask& mask, const ao::uvector<double>& thresholds, size_t minScale, size_t endScale, const Component& component, double cleanBorder, const bool* priorMask, size_t& areaSize)
{
	const size_t width = iuwt.Width(), height = iuwt.Height();
	size_t xBorder = cleanBorder*width;
	size_t yBorder = cleanBorder*height;
	size_t minX = xBorder, maxX = width - xBorder;
	size_t minY = yBorder, maxY = height - yBorder;
	
	areaSize = 0;
	endScale = std::min<size_t>(endScale, iuwt.NScales());
	std::stack<Component> todo;
	todo.push(component);
	mask[component.scale][component.x + component.y*width] = true;
	while(!todo.empty())
	{
		Component c = todo.top();
		++areaSize;
		todo.pop();
		size_t index = c.x + c.y*width;
		if(c.x > minX)
		{
			if(exceedsThreshold(iuwt[c.scale][index-1], thresholds[c.scale]) && !mask[c.scale][index-1] && priorMask[index-1])
			{
				mask[c.scale][index-1] = true;
				todo.push(Component(c.x-1, c.y, c.scale));
			}
		}
		if(c.x < maxX-1)
		{
			if(exceedsThreshold(iuwt[c.scale][index+1], thresholds[c.scale]) && !mask[c.scale][index+1] && priorMask[index+1])
			{
				mask[c.scale][index+1] = true;
				todo.push(Component(c.x+1, c.y, c.scale));
			}
		}
		if(c.y > minY)
		{
			if(exceedsThreshold(iuwt[c.scale][index-width], thresholds[c.scale]) && !mask[c.scale][index-width] && priorMask[index-width])
			{
				mask[c.scale][index-width] = true;
				todo.push(Component(c.x, c.y-1, c.scale));
			}
		}
		if(c.y < maxY-1)
		{
			if(exceedsThreshold(iuwt[c.scale][index+width], thresholds[c.scale]) && !mask[c.scale][index+width] && priorMask[index+width])
			{
				mask[c.scale][index+width] = true;
				todo.push(Component(c.x, c.y+1, c.scale));
			}
		}
		if(c.scale > int(minScale))
		{
			if(exceedsThreshold(iuwt[c.scale-1][index], thresholds[c.scale-1]) && !mask[c.scale-1][index] && priorMask[index])
			{
				mask[c.scale-1][index] = true;
				todo.push(Component(c.x, c.y, c.scale-1));
			}
		}
		if(c.scale < int(endScale)-1)
		{
			if(exceedsThreshold(iuwt[c.scale+1][index], thresholds[c.scale+1]) && !mask[c.scale+1][index] && priorMask[index])
			{
				mask[c.scale+1][index] = true;
				todo.push(Component(c.x, c.y, c.scale+1));
			}
		}
	}
}

void ImageAnalysis::SelectStructures(const IUWTDecomposition& iuwt, IUWTMask& mask, const ao::uvector<double>& thresholds, size_t minScale, size_t endScale, double cleanBorder, const bool* priorMask, size_t& areaSize)
{
	const size_t width = iuwt.Width(), height = iuwt.Height();
	const size_t
		xBorder = cleanBorder*width,
		yBorder = cleanBorder*height,
		minX = xBorder, maxX = width - xBorder,
		minY = yBorder, maxY = height - yBorder;
	
	areaSize = 0;
	
	for(size_t scale=minScale; scale!=endScale; ++scale)
	{
		for(size_t y=minY; y!=maxY; ++y)
		{
			for(size_t x=minX; x!=maxX; ++x)
			{
				size_t index = x + y*width;
				bool isInPriorMask = (priorMask==0) || priorMask[index];
				if(exceedsThreshold(iuwt[scale][index], thresholds[scale]) && !mask[scale][index] && isInPriorMask)
				{
					Component component(x, y, scale);
					size_t subAreaSize = 0;
					if(priorMask == 0)
						Floodfill(iuwt, mask, thresholds, minScale, endScale, component, cleanBorder, subAreaSize);
					else
						MaskedFloodfill(iuwt, mask, thresholds, minScale, endScale, component, cleanBorder, priorMask, subAreaSize);
					areaSize += subAreaSize;
				}
			}
		}
	}
}

void ImageAnalysis::FloodFill2D(const double* image, bool* mask, double threshold, const ImageAnalysis::Component2D& component, size_t width, size_t height, size_t& areaSize)
{
	areaSize = 0;
	std::stack<Component2D> todo;
	todo.push(component);
	mask[component.x + component.y*width] = true;
	while(!todo.empty())
	{
		Component2D c = todo.top();
		++areaSize;
		todo.pop();
		size_t index = c.x + c.y*width;
		if(c.x > 0)
		{
			if(exceedsThreshold(image[index-1], threshold) && !mask[index-1])
			{
				mask[index-1] = true;
				todo.push(Component2D(c.x-1, c.y));
			}
		}
		if(c.x < width-1)
		{
			if(exceedsThreshold(image[index+1], threshold) && !mask[index+1])
			{
				mask[index+1] = true;
				todo.push(Component2D(c.x+1, c.y));
			}
		}
		if(c.y > 0)
		{
			if(exceedsThreshold(image[index-width], threshold) && !mask[index-width])
			{
				mask[index-width] = true;
				todo.push(Component2D(c.x, c.y-1));
			}
		}
		if(c.y < height-1)
		{
			if(exceedsThreshold(image[index+width], threshold) && !mask[index+width])
			{
				mask[index+width] = true;
				todo.push(Component2D(c.x, c.y+1));
			}
		}
	}
}

void ImageAnalysis::FloodFill2D(const double* image, bool* mask, double threshold, const ImageAnalysis::Component2D& component, size_t width, size_t height, std::vector<Component2D>& area)
{
	area.clear();
	std::stack<Component2D> todo;
	todo.push(component);
	mask[component.x + component.y*width] = true;
	while(!todo.empty())
	{
		Component2D c = todo.top();
		area.push_back(c);
		todo.pop();
		size_t index = c.x + c.y*width;
		if(c.x > 0)
		{
			if(exceedsThresholdAbs(image[index-1], threshold) && !mask[index-1])
			{
				mask[index-1] = true;
				todo.push(Component2D(c.x-1, c.y));
			}
		}
		if(c.x < width-1)
		{
			if(exceedsThresholdAbs(image[index+1], threshold) && !mask[index+1])
			{
				mask[index+1] = true;
				todo.push(Component2D(c.x+1, c.y));
			}
		}
		if(c.y > 0)
		{
			if(exceedsThresholdAbs(image[index-width], threshold) && !mask[index-width])
			{
				mask[index-width] = true;
				todo.push(Component2D(c.x, c.y-1));
			}
		}
		if(c.y < height-1)
		{
			if(exceedsThresholdAbs(image[index+width], threshold) && !mask[index+width])
			{
				mask[index+width] = true;
				todo.push(Component2D(c.x, c.y+1));
			}
		}
	}
}
