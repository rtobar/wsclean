#ifndef WSCLEAN_IMAGING_TABLE_H
#define WSCLEAN_IMAGING_TABLE_H

#include <string>
#include <vector>

#include "../uvector.h"
#include "../polarization.h"

#include "../msproviders/partitionedms.h"

class ImagingTableEntry
{
public:
	struct MSBandInfo
	{
		size_t bandIndex;
		size_t partIndex;
	};
	
	struct MSInfo
	{
		std::vector<MSBandInfo> bands;
	};
	
	ImagingTableEntry();
	
	size_t index;
	
	/**
	 * Note that mses might have overlapping frequencies.
	 */
	double lowestFrequency, highestFrequency;
	double bandStartFrequency, bandEndFrequency;
	size_t inputChannelCount;
	
	PolarizationEnum polarization;
	
	size_t outputChannelIndex;
	
	size_t outputIntervalIndex;
	
	/**
	 * This vector links a filename index to MS data
	 */
	std::vector<MSInfo> msData;
	
	/**
	 * The group of entries with equal squaredDeconvolutionIndex
	 * should be 'joinedly' deconvolved by adding their squared powers
	 * together. Normally, all the polarizations from a single
	 * (output)channel / timestep form such a group.
	 */
	size_t squaredDeconvolutionIndex;
	
	/**
	 * Entries with equal joinedGroupIndex are joinedly deconvolved.
	 * Such a group of entries can be further split up in 'squared'
	 * deconvolution groups.
	 */
	size_t joinedGroupIndex;
	
	/**
	 * A normal inversion results in '1' image. However, an XY
	 * imaging run results in 2 (real and imaginary), while an
	 * YX imaging run results in 0, as it is added to XY.
	 */
	size_t imageCount;
	
	std::string tmpFilePrefix;
	
	std::string ToString();
	
	double CentralFrequency() const
	{
		return 0.5 * (bandStartFrequency + bandEndFrequency);
	}
	
	/**
	 * A number that scales with the estimated inverse-variance of the image. It can be used when
	 * averaging images or fitting functions through the images to get the optimal sensitivity.
	 * It is set after the first inversion.
	 */
	double imageWeight;
};

class ImagingTable
{
public:
	size_t IndependentGroupCount() const
	{
		return _independentGroupLookup.size();
	}
	
	ImagingTable GetIndependentGroup(size_t index) const;
	
	size_t SquaredGroupCount() const
	{
		return _squaredGroupLookup.size();
	}
	
	ImagingTable GetSquaredGroup(size_t index) const;
	
	size_t EntryCount() const
	{
		return _entries.size();
	}
	
	ImagingTableEntry& operator[](size_t index)
	{
		return *_entries[index];
	}
	const ImagingTableEntry& operator[](size_t index) const
	{
		return *_entries[index];
	}
	
	size_t ImageCount() const
	{
		return _imageLookup.size();
	}
	
	void GetImageInfo(size_t imageIndex, bool& isImaginary);
	
	void Clear() { _entries.clear(); }
	
	ImagingTableEntry& AddEntry()
	{
		_entries.emplace_back(new ImagingTableEntry());
		return *_entries.back();
	}
	
	void Update()
	{
		updateIndependentGroupLookup();
		updateSquaredGroupLookup();
		updateImageLookup();
	}
	
	void Print();
	
	ImagingTableEntry& Front() { return *_entries.front(); }
	const ImagingTableEntry& Front() const { return *_entries.front(); }
	
	const ImagingTableEntry* FirstWithHigherFrequency(double frequency) const
	{
		double currentDistance = std::numeric_limits<double>::max();
		ImagingTableEntry* entry = nullptr;
		
		for(auto& e : _entries)
		{
			if(e->CentralFrequency() > frequency &&
				e->CentralFrequency() - frequency < currentDistance)
			{
				currentDistance = e->CentralFrequency() - frequency;
				entry = &*e;
			}
		}
		return entry;
	}
	
	const ImagingTableEntry* FirstWithLowerFrequency(double frequency) const
	{
		double currentDistance = std::numeric_limits<double>::max();
		ImagingTableEntry* entry = nullptr;
		
		for(auto& e : _entries)
		{
			if(e->CentralFrequency() < frequency &&
				frequency - e->CentralFrequency() < currentDistance)
			{
				currentDistance = frequency - e->CentralFrequency();
				entry = &*e;
			}
		}
		return entry;
	}
	
private:
	void printIndependentGroup(bool isFinal);
	void updateIndependentGroupLookup();
	void updateSquaredGroupLookup();
	void updateImageLookup();
	
	typedef std::shared_ptr<ImagingTableEntry> ImagingTableEntryPtr;
	std::vector<ImagingTableEntryPtr> _entries;
	
	std::vector<std::vector<ImagingTableEntryPtr>> _independentGroupLookup;
	std::vector<std::vector<ImagingTableEntryPtr>> _squaredGroupLookup;
	std::vector<std::pair<ImagingTableEntryPtr,bool>> _imageLookup;
};

#endif
