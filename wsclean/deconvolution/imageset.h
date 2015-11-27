#ifndef CLEANABLE_IMAGE_SET_H
#define CLEANABLE_IMAGE_SET_H

#include "../wsclean/cachedimageset.h"
#include "../wsclean/imagebufferallocator.h"

#include "dynamicset.h"

#include <vector>

namespace deconvolution {
		
	class SingleImageSet {
	public:
		struct Value {
			double value;
			Value() { }
			Value(double _value) : value(_value) { }
			double GetValue(size_t i) { 
				return value;
			}
			static Value Zero() { return Value(0.0); }
		};
		
		SingleImageSet(size_t size, SingleImageSet& prototype) :
			image(prototype._allocator->Allocate(size)),
			_allocator(prototype._allocator)
		{
		}
		
		SingleImageSet(size_t size, ImageBufferAllocator& allocator) :
			image(allocator.Allocate(size)),
			_allocator(&allocator)
		{
		}
		
		SingleImageSet(double* allocatedData, ImageBufferAllocator& allocator) :
			image(allocatedData),
			_allocator(&allocator)
		{
		}
		
		~SingleImageSet()
		{
			_allocator->Free(image);
		}
		
		void Load(CachedImageSet& set, PolarizationEnum pol, size_t freqIndex)
		{
			set.Load(image, pol, freqIndex, false);
		}
		
		void Store(CachedImageSet& set, PolarizationEnum pol, size_t freqIndex) const
		{
			set.Store(image, pol, freqIndex, false);
		}
		
		void Transfer(DynamicSet& destination)
		{
			destination.Claim(0, image);
			image = 0;
		}
		
		Value Get(size_t pixelIndex) const
		{
			return Value(image[pixelIndex]);
		}
		
		double JoinedValue(size_t pixelIndex) const
		{
			return image[pixelIndex];
		}
		
		double JoinedValueNormalized(size_t pixelIndex) const
		{
			return image[pixelIndex];
		}
		
		double AbsJoinedValue(size_t pixelIndex) const
		{
			return fabs(image[pixelIndex]);
		}
		
		bool IsComponentNegative(size_t pixelIndex) const
		{
			return image[pixelIndex]<0.0;
		}
		
		void AddComponent(size_t pixelIndex, const double* values)
		{
			image[pixelIndex] += values[0];
		}
		
		size_t ImageCount() const { return 1; }
		
		static size_t StaticImageCount() { return 1; }
		
		double* GetImage(size_t imageIndex)
		{
			return image;
		}
		static size_t PSFIndex(size_t imageIndex)
		{
			return 0;
		}
		ImageBufferAllocator* Allocator() { return _allocator; }
		
		double* Data() { return image; }
	private:
		double *image;
		
		ImageBufferAllocator* _allocator;
	};
	
	template<size_t PolCount>
	class PolarizedImageSet {
	public:
		struct Value {
			double data[PolCount];
			double GetValue(size_t i) { 
				return data[i];
			}
			static Value Zero() {
				Value zero;
				for(size_t i=0; i!=PolCount; ++i)
					zero.data[i] = 0.0;
				return zero;
			}
		};
		
		PolarizedImageSet(size_t size, PolarizedImageSet<PolCount>& prototype) :
			_allocator(prototype._allocator)
		{
			for(size_t i=0; i!=PolCount; ++i)
				images[i] = _allocator->Allocate(size);
		}
		
		PolarizedImageSet(size_t size, ImageBufferAllocator& allocator) :
			_allocator(&allocator)
		{
			for(size_t i=0; i!=PolCount; ++i)
				images[i] = _allocator->Allocate(size);
		}
		
		PolarizedImageSet(DynamicSet& source, ImageBufferAllocator& allocator) :
			_allocator(&allocator)
		{
			for(size_t i=0; i!=PolCount; ++i)
				images[i] = source.Release(i);
		}
		
		PolarizedImageSet(DynamicSet& source, size_t offset, ImageBufferAllocator& allocator) :
			_allocator(&allocator)
		{
			for(size_t i=0; i!=PolCount; ++i)
				images[i] = source.Release(i+offset);
		}
		
		~PolarizedImageSet()
		{
			for(size_t i=0; i!=PolCount; ++i)
				_allocator->Free(images[i]);
		}
		
		void Transfer(DynamicSet& destination, size_t offset)
		{
			for(size_t i=0; i!=PolCount; ++i)
			{
				destination.Claim(offset+i, images[i]);
				images[i] = 0;
			}
		}
		
		void Load(CachedImageSet& set, PolarizationEnum polarization, size_t freqIndex)
		{
			set.Load(images[0], polarization, freqIndex, false);
		}
		
		void Load(CachedImageSet& set, const std::set<PolarizationEnum>& polarizations, size_t freqIndex)
		{
			std::set<PolarizationEnum>::const_iterator p=polarizations.begin();
			for(size_t i = 0; i!=PolCount; ++i, ++p)
			{
				if(*p == Polarization::YX)
					set.Load(images[i], PolarizationEnum::XY, freqIndex, true);
				else
					set.Load(images[i], *p, freqIndex, false);
			}
		}
		
		void Store(CachedImageSet& set, PolarizationEnum polarization, size_t freqIndex)
		{
			set.Store(images[0], polarization, freqIndex, false);
		}
		
		void Store(CachedImageSet& set, const std::set<PolarizationEnum>& polarizations, size_t freqIndex) const
		{
			std::set<PolarizationEnum>::const_iterator p=polarizations.begin();
			for(size_t i = 0; i!=PolCount; ++i, ++p)
			{
				if(*p == Polarization::YX)
					set.Store(images[i], PolarizationEnum::XY, freqIndex, true);
				else
					set.Store(images[i], *p, freqIndex, false);
			}
		}
		
		Value Get(size_t index) const
		{
			Value v;
			for(size_t i=0; i!=PolCount; ++i)
				v.data[i] = images[i][index];
			return v;
		}
		
		double JoinedValue(size_t index) const
		{
			return SquaredSum(index);
		}
		
		double JoinedValueNormalized(size_t index) const
		{
			return sqrt(SquaredSum(index)*0.5);
		}
		
		double AbsJoinedValue(size_t index) const
		{
			return SquaredSum(index);
		}
		
		double SquaredSum(size_t index) const
		{
			if(PolCount == 4)
			{
				return
					images[0][index]*images[0][index] +
					images[1][index]*images[1][index] + images[2][index]*images[2][index] +
					images[3][index]*images[3][index];
			}
			else {
				double sum = 0.0;
				for(size_t i=0; i!=PolCount; ++i)
					sum += images[i][index]*images[i][index];
				return sum;
			}
		}
		
		bool IsComponentNegative(size_t index) const
		{
			if(PolCount == 4)
				return images[0][index]<0.0 || images[3][index]<0.0;
			else if(PolCount == 2)
				return images[0][index]<0.0 || images[1][index]<0.0;
			else if(PolCount == 1)
				return images[0][index]<0.0;
			else
				return false;
		}
		
		void AddComponent(size_t index, const double* values)
		{
			for(size_t i=0; i!=PolCount; ++i)
				images[i][index] += values[i];
		}
		
		size_t ImageCount() const { return PolCount; }
		
		static size_t StaticImageCount() { return PolCount; }
		
		double* GetImage(size_t imageIndex)
		{
			return images[imageIndex];
		}
		static size_t PSFIndex(size_t imageIndex)
		{
			return 0;
		}
		ImageBufferAllocator* Allocator() { return _allocator; }
	private:
		double *images[PolCount];
		
		ImageBufferAllocator* _allocator;
	};
	
	template<typename SingleImageSetType>
	class MultiImageSet {
	public:
		struct Value {
			std::vector<typename SingleImageSetType::Value> values;
			double GetValue(size_t i)
			{
				return values[i/SingleImageSetType::StaticImageCount()].GetValue(i%SingleImageSetType::StaticImageCount());
			}
			static Value Zero() { return Value(); }
		};
		
		MultiImageSet(size_t imageSize, MultiImageSet& prototype)
		{
			for(size_t i=0; i!=prototype._sets.size(); ++i)
			{
				_sets.push_back(new SingleImageSetType(imageSize, *prototype.Allocator()));
			}
		}
		
		MultiImageSet(size_t imageSize, size_t count, ImageBufferAllocator& allocator)
		{
			for(size_t i=0; i!=count; ++i)
			{
				_sets.push_back(new SingleImageSetType(imageSize, allocator));
			}
		}
		
		MultiImageSet(DynamicSet& source, size_t count, ImageBufferAllocator& allocator)
		{
			size_t imgIndex = 0;
			for(size_t i=0; i!=count; ++i)
			{
				SingleImageSetType* img = new SingleImageSetType(source, imgIndex, allocator);
				size_t imagesInGroup = img->ImageCount();
				imgIndex += imagesInGroup;
				_sets.push_back(img);
			}
		}
		
		~MultiImageSet()
		{
			for(typename std::vector<SingleImageSetType*>::const_iterator i=_sets.begin(); i!=_sets.end(); ++i)
			{
				delete *i;
			}
		}
		
		void Transfer(DynamicSet& destination)
		{
			size_t imgIndex = 0;
			for(size_t i=0; i!=_sets.size(); ++i)
			{
				size_t imagesInGroup = _sets[i]->ImageCount();
				_sets[i]->Transfer(destination, imgIndex);
				imgIndex += imagesInGroup;
			}
		}
		
		void Load(CachedImageSet& set, PolarizationEnum polarization, size_t freqIndex)
		{
			_sets[freqIndex]->Load(set, polarization, freqIndex);
		}
		
		void Load(CachedImageSet& set, const std::set<PolarizationEnum>& polarizations, size_t i)
		{
			_sets[i]->Load(set, polarizations, i);
		}
		
		void Store(CachedImageSet& set, PolarizationEnum polarization, size_t freqIndex)
		{
			_sets[freqIndex]->Store(set, polarization, freqIndex);
		}
		
		void Store(CachedImageSet& set, const std::set<PolarizationEnum>& polarizations, size_t i) const
		{
			_sets[i]->Store(set, polarizations, i);
		}
		
		double JoinedValue(size_t index) const
		{
			double val = 0.0;
			for(typename std::vector<SingleImageSetType*>::const_iterator i=_sets.begin(); i!=_sets.end(); ++i)
			{
				val += (*i)->JoinedValueNormalized(index);
			}
			return val;
		}
		
		double JoinedValueNormalized(size_t index) const
		{
			return JoinedValue(index) / _sets.size();
		}
		
		double AbsJoinedValue(size_t index) const
		{
			return std::fabs(JoinedValue(index));
		}
		
		bool IsComponentNegative(size_t index) const
		{
			for(typename std::vector<SingleImageSetType*>::const_iterator i=_sets.begin(); i!=_sets.end(); ++i)
			{
				if((*i)->IsComponentNegative(index)) return true;
			}
			return false;
		}
		
		void AddComponent(size_t index, const double* values)
		{
			size_t vi = 0;
			for(size_t i=0; i!=_sets.size(); ++i)
			{
				_sets[i]->AddComponent(index, &values[vi]);
				vi += _sets[i]->ImageCount();
			}
		}
		
		Value Get(const size_t index)
		{
			Value v;
			v.values.resize(_sets.size());
			for(size_t i=0; i!=_sets.size(); ++i)
				v.values[i] = _sets[i]->Get(index);
			return v;
		}
		
		size_t ImageCount() const { return SingleImageSetType::StaticImageCount() * _sets.size(); }
		
		double* GetImage(size_t imageIndex)
		{
			return _sets[imageIndex/SingleImageSetType::StaticImageCount()]->GetImage(imageIndex%SingleImageSetType::StaticImageCount());
		}
		
		double* GetImage(size_t polIndex, size_t freqIndex)
		{
			return _sets[freqIndex]->GetImage(polIndex);
		}
		
		static size_t PSFIndex(size_t imageIndex)
		{
			return imageIndex/SingleImageSetType::StaticImageCount();
		}
		ImageBufferAllocator* Allocator()
		{ 
			return _sets.front()->Allocator();
		}
	private:
		std::vector<SingleImageSetType*> _sets;
	};	
}

#endif
