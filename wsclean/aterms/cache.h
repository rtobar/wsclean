#ifndef ATERM_CACHE_H
#define ATERM_CACHE_H

#include <algorithm>
#include <complex>
#include <memory>
#include <vector>

/**
 * A simple a-term cache that stores aterms for all frequencies in a single timestep.
 * Since each timestep will have the same frequency values, the cache maintains buffers
 * for each frequency value, and once all frequencies have been stored at least once,
 * no more allocations are performed.
 * 
 * Used by the FitsATerm class.
 */
class Cache
{
public:
	/**
	 * Construct the cache.
	 */
	Cache(size_t atermSize) : _atermSize(atermSize) { }
	
	Cache(const Cache&) = delete;
	Cache(Cache&&) = delete;
	Cache& operator=(const Cache&) = delete;
	Cache& operator=(Cache&&) = delete;
	
	static const size_t NOT_FOUND;
	
	/**
	 * Clears all stored values, such that e.g. the cache is ready for a new
	 * timestep. After this call, @ref Find() will always return @ref NOT_FOUND .
	 */
	void Reset()
	{
		// The cache is 'emptied', but we don't entirely clear the data vector at
		// the top level, as we can reuse the space from the individual entries
		// and prevent reallocation.
		for(Entry& entry : _entries)
		{
			entry.isValid = false;
		}
	}
	
	/**
	 * Returns the index of the frequency, or @ref NOT_FOUND if not found.
	 */
	size_t Find(double frequency) const
	{
		auto iter = std::lower_bound(_frequencies.begin(), _frequencies.end(), frequency);
		if(iter == _frequencies.end() || *iter != frequency)
			return NOT_FOUND;
		else {
			size_t index = iter - _frequencies.begin();
			if(_entries[index].isValid)
				return index;
			else
				return NOT_FOUND;
		}
	}
	
	/**
	 * Retrieve a buffer given a frequency index (as returned by @ref Find() ).
	 * @param destination Array with @ref ATermSize() elements.
	 */
	void Get(size_t index, std::complex<float>* destination) const
	{
		std::copy_n(_entries[index].ptr.get(), _atermSize, destination);
	}
	
	/**
	 * Store the data for the given frequency in the cache. The data array should
	 * be an array with @ref ATermSize() elements.
	 */
	void Store(double frequency, const std::complex<float>* data)
	{
		auto iter = std::lower_bound(_frequencies.begin(), _frequencies.end(), frequency);
		if(iter!=_frequencies.end() && *iter == frequency)
		{
			size_t index = iter - _frequencies.begin();
			std::copy_n(data, _atermSize, _entries[index].ptr.get());
			_entries[index].isValid = true;
		}
		else {
			size_t index = iter - _frequencies.begin();
			_frequencies.emplace(iter, frequency);
			Entry newEntry;
			newEntry.ptr.reset( new std::complex<float>[_atermSize] );
			newEntry.isValid = true;
			std::copy_n(data, _atermSize, newEntry.ptr.get());
			_entries.emplace(_entries.begin() + index, std::move(newEntry));
		}
	}
	
	/**
	 * Size of one aterm buffer (in number of complex float values).
	 */
	size_t ATermSize() const { return _atermSize; }
	
private:
	// This array is always kept sorted
	std::vector<double> _frequencies;
	size_t _atermSize;
	
	struct Entry {
		std::unique_ptr<std::complex<float>[]> ptr;
		bool isValid;
	};
	
	// _entries[index] corresponds with _frequencies[index]
	std::vector<Entry> _entries;
};

#endif

