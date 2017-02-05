#ifndef COMPONENT_LIST_H
#define COMPONENT_LIST_H

#include "../image.h"
#include "../uvector.h"

#include "../wsclean/imagebufferallocator.h"

#include <vector>

class ComponentList
{
public:
	ComponentList(size_t width, size_t height, size_t nScales, size_t nFrequencies, ImageBufferAllocator& allocator) :
		_width(width), _height(height),
		_nScales(nScales), _nFrequencies(nFrequencies),
		_componentsAddedSinceLastMerge(0),
		_maxComponentsBeforeMerge(100000),
		_listPerScale(nScales),
		_allocator(allocator)
	{
	}
	
	void Add(size_t x, size_t y, size_t scaleIndex, const double* values)
	{
		_listPerScale[scaleIndex].values.push_back(values, values + _nFrequencies);
		_listPerScale[scaleIndex].positions.emplace_back(x, y);
		++_componentsAddedSinceLastMerge;
		if(_componentsAddedSinceLastMerge >= _maxComponentsBeforeMerge)
			MergeDuplicates();
	}
	
	void Write(const class MultiScaleAlgorithm& multiscale, const class WSCleanSettings& settings, long double pixelScaleX, long double pixelScaleY, long double phaseCentreRA, long double phaseCentreDec);
	
	void MergeDuplicates()
	{
		for(size_t scaleIndex=0; scaleIndex!=_nScales; ++scaleIndex)
		{
			mergeDuplicates(scaleIndex);
		}
		_componentsAddedSinceLastMerge = 0;
	}
	
	size_t ComponentCount(size_t scaleIndex) const
	{ return _listPerScale[scaleIndex].positions.size(); }
	
	void GetComponent(size_t scaleIndex, size_t index, size_t& x, size_t& y, double* values) const
	{
		x = _listPerScale[scaleIndex].positions[index].x;
		y = _listPerScale[scaleIndex].positions[index].y;
		for(size_t f=0; f!=_nFrequencies; ++f)
			values[f] = _listPerScale[scaleIndex].values[index * _nFrequencies + f];
	}
private:
	struct Position {
		Position(size_t _x, size_t _y) :
			x(_x), y(_y)
		{ }
		size_t x, y;
	};
	struct ScaleList
	{
	/**
	 * This list contains nFrequencies values for each
	 * component, such that _positions[i] corresponds with the values
	 * starting at _values[i * _nFrequencies].
	 */
		ao::uvector<double> values;
		ao::uvector<Position> positions;
	};
	
	void mergeDuplicates(size_t scaleIndex)
	{
		ScaleList& list = _listPerScale[scaleIndex];
		ao::uvector<double> newValues;
		ao::uvector<Position> newPositions;
		
		std::vector<Image> images(_nFrequencies);
		for(Image& image : images)
			image = Image(_width, _height, 0.0, _allocator);
		size_t valueIndex = 0;
		for(size_t index=0; index!=list.positions.size(); ++index)
		{
			size_t position = list.positions[index].x + list.positions[index].y * _width;
			for(size_t frequency = 0; frequency != _nFrequencies; ++frequency)
			{
				images[frequency][position] += list.values[valueIndex];
				valueIndex++;
			}
		}
		
		list.values.clear();
		list.positions.clear();
		
		for(size_t imageIndex = 0; imageIndex != images.size(); ++imageIndex)
		{
			Image& image = images[imageIndex];
			size_t posIndex = 0;
			for(size_t y=0; y!=_height; ++y)
			{
				for(size_t x=0; x!=_width; ++x)
				{
					if(image[posIndex] != 0.0)
					{
						for(size_t i=0; i!=images.size(); ++i)
						{
							newValues.push_back(images[i][posIndex]);
							images[i][posIndex] = 0.0;
						}
						newPositions.emplace_back(x, y);
					}
					++posIndex;
				}
			}
		}
		std::swap(_listPerScale[scaleIndex].values, newValues);
		std::swap(_listPerScale[scaleIndex].positions, newPositions);
	}
	const size_t _width, _height;
	size_t _nScales, _nFrequencies;
	size_t _componentsAddedSinceLastMerge;
	size_t _maxComponentsBeforeMerge;
	std::vector<ScaleList> _listPerScale;
	ImageBufferAllocator& _allocator;
};

#endif
