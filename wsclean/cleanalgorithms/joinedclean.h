#ifndef JOINED_CLEAN_H
#define JOINED_CLEAN_H

#include "cleanalgorithm.h"
#include "simpleclean.h"

#include "../imagebufferallocator.h"
#include "../cachedimageset.h"

namespace ao {
	template<typename T> class lane;
}

namespace joined_pol_clean {
		
	class SingleImageSet {
	public:
		struct Value {
			double xx, xyr, xyi, yy;
			double GetValue(size_t i) { 
				switch(i) {
					default:
					case 0: return xx;
					case 1: return xyr;
					case 2: return xyi;
					case 3: return yy;
				}
			}
			static Value Zero() {
				Value zero;
				zero.xx = 0.0; zero.xyr = 0.0;
				zero.xyi = 0.0; zero.yy = 0.0;
				return zero;
			}
		};
		
		double *xx, *xyr, *xyi, *yy;
		
		SingleImageSet(size_t size, ImageBufferAllocator<double>& allocator) :
			xx(allocator.Allocate(size)),
			xyr(allocator.Allocate(size)),
			xyi(allocator.Allocate(size)),
			yy(allocator.Allocate(size)),
			_allocator(&allocator)
		{
		}
		
		~SingleImageSet()
		{
			_allocator->Free(xx);
			_allocator->Free(xyr);
			_allocator->Free(xyi);
			_allocator->Free(yy);
		}
		
		void LoadLinear(CachedImageSet& set, size_t freqIndex)
		{
			set.Load(xx, PolarizationEnum::XX, freqIndex, false);
			set.Load(xyr, PolarizationEnum::XY, freqIndex, false);
			set.Load(xyi, PolarizationEnum::XY, freqIndex, true);
			set.Load(yy, PolarizationEnum::YY, freqIndex, false);
		}
		
		void LoadStokes(CachedImageSet& set, size_t freqIndex)
		{
			set.Load(xx, PolarizationEnum::StokesI, freqIndex, false);
			set.Load(xyr, PolarizationEnum::StokesQ, freqIndex, false);
			set.Load(xyi, PolarizationEnum::StokesU, freqIndex, false);
			set.Load(yy, PolarizationEnum::StokesV, freqIndex, false);
		}
		
		void StoreLinear(CachedImageSet& set, size_t freqIndex) const
		{
			set.Store(xx, PolarizationEnum::XX, freqIndex, false);
			set.Store(xyr, PolarizationEnum::XY, freqIndex, false);
			set.Store(xyi, PolarizationEnum::XY, freqIndex, true);
			set.Store(yy, PolarizationEnum::YY, freqIndex, false);
		}
		
		void StoreStokes(CachedImageSet& set, size_t freqIndex) const
		{
			set.Store(xx, PolarizationEnum::StokesI, freqIndex, false);
			set.Store(xyr, PolarizationEnum::StokesQ, freqIndex, false);
			set.Store(xyi, PolarizationEnum::StokesU, freqIndex, false);
			set.Store(yy, PolarizationEnum::StokesV, freqIndex, false);
		}
		
		Value Get(size_t index) const
		{
			Value v;
			v.xx = xx[index];
			v.xyr = xyr[index];
			v.xyi = xyi[index];
			v.yy = yy[index];
			return v;
		}
		
		double JoinedValue(size_t index) const
		{
			return SquaredSum(index);
		}
		
		double JoinedValueNormalized(size_t index) const
		{
			return sqrt(SquaredSum(index));
		}
		
		double SquaredSum(size_t index) const
		{
			return
				xx[index]*xx[index] +
				xyr[index]*xyr[index] + xyi[index]*xyi[index] +
				yy[index]*yy[index];
		}
		
		bool IsComponentNegative(size_t index) const
		{
			return xx[index]<0.0 || yy[index]<0.0;
		}
		
		void AddComponent(const SingleImageSet& source, size_t index, double factor)
		{
			xx[index] += source.xx[index] * factor;
			xyr[index] += source.xyr[index] * factor;
			xyi[index] += source.xyi[index] * factor;
			yy[index] += source.yy[index] * factor;
		}
		
		size_t ImageCount() const { return 4; }
		
		static size_t StaticImageCount() { return 4; }
		
		double* GetImage(size_t imageIndex)
		{
			double* vals[4] = { xx, xyr, xyi, yy };
			return vals[imageIndex];
		}
		static size_t PSFIndex(size_t imageIndex)
		{
			return 0;
		}
	private:
		ImageBufferAllocator<double> *_allocator;
	};
	
	class MultiImageSet {
	public:
		struct Value {
			std::vector<SingleImageSet::Value> values;
			double GetValue(size_t i)
			{
				return values[i/4].GetValue(i%4);
			}
			static Value Zero() { return Value(); }
		};
		MultiImageSet(size_t imageSize, size_t count, ImageBufferAllocator<double>& allocator)
		{
			for(size_t i=0; i!=count; ++i)
			{
				_sets.push_back(new SingleImageSet(imageSize, allocator));
			}
		}
		
		~MultiImageSet()
		{
			for(std::vector<SingleImageSet*>::const_iterator i=_sets.begin(); i!=_sets.end(); ++i)
			{
				delete *i;
			}
		}
		
		void LoadLinear(CachedImageSet& set, size_t i)
		{
			_sets[i]->LoadLinear(set, i);
		}
		
		void StoreLinear(CachedImageSet& set, size_t i) const
		{
			_sets[i]->StoreLinear(set, i);
		}
		
		void LoadStokes(CachedImageSet& set, size_t i)
		{
			_sets[i]->LoadStokes(set, i);
		}
		
		void StoreStokes(CachedImageSet& set, size_t i) const
		{
			_sets[i]->StoreStokes(set, i);
		}
		
		double JoinedValue(size_t index) const
		{
			double val = 0.0;
			for(std::vector<SingleImageSet*>::const_iterator i=_sets.begin(); i!=_sets.end(); ++i)
			{
				val += (*i)->JoinedValueNormalized(index);
			}
			return val;
		}
		
		double JoinedValueNormalized(size_t index) const
		{
			return JoinedValue(index) / _sets.size();
		}
		
		bool IsComponentNegative(size_t index) const
		{
			for(std::vector<SingleImageSet*>::const_iterator i=_sets.begin(); i!=_sets.end(); ++i)
			{
				if((*i)->IsComponentNegative(index)) return true;
			}
			return false;
		}
		
		void AddComponent(const MultiImageSet& source, size_t index, double factor)
		{
			for(size_t i=0; i!=_sets.size(); ++i)
			{
				_sets[i]->AddComponent(*source._sets[i], index, factor);
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
		
		size_t ImageCount() const { return SingleImageSet::StaticImageCount() * _sets.size(); }
		
		double* GetImage(size_t imageIndex)
		{
			return _sets[imageIndex/4]->GetImage(imageIndex%4);
		}
		
		double* GetImage(size_t polIndex, size_t freqIndex)
		{
			return _sets[freqIndex]->GetImage(polIndex);
		}
		
		static size_t PSFIndex(size_t imageIndex)
		{
			return imageIndex/4;
		}
	private:
		std::vector<SingleImageSet*> _sets;
	};	
}

template<typename ImageSetType = joined_pol_clean::SingleImageSet>
class JoinedClean : public CleanAlgorithm
{
public:
	void ExecuteMajorIteration(ImageSetType& dataImage, ImageSetType& modelImage, std::vector<double*> psfImages, size_t width, size_t height, bool& reachedStopGain);
	typedef ImageSetType ImageSet;
private:
	size_t _width, _height;
	
	struct CleanTask
	{
		size_t cleanCompX, cleanCompY;
		typename ImageSetType::Value peak;
	};
	struct CleanResult
	{
		CleanResult() : nextPeakX(0), nextPeakY(0), peakLevelUnnormalized(0.0)
		{ }
		size_t nextPeakX, nextPeakY;
		double peakLevelUnnormalized;
	};
	struct CleanThreadData
	{
		size_t startY, endY;
		ImageSetType* dataImage;
		std::vector<double*> psfImages;
	};

	void findPeak(const ImageSetType& image, size_t& x, size_t& y) const
	{
		findPeak(image, x, y, 0, _height);
	}
	void findPeak(const ImageSetType& image, size_t& x, size_t& y, size_t startY, size_t stopY) const;
	
	std::string peakDescription(const ImageSetType& image, size_t& x, size_t& y);
	void cleanThreadFunc(ao::lane<CleanTask>* taskLane, ao::lane<CleanResult>* resultLane, CleanThreadData cleanData);
	
	void subtractImage(double *image, const double *psf, size_t x, size_t y, double factor, size_t startY, size_t endY) const
	{
		SimpleClean::PartialSubtractImage(image, _width, _height, psf, _width, _height, x, y, factor, startY, endY);
	}
};

#endif