#include "../deconvolution/dynamicset.h"
#include "../deconvolution/spectralfitter.h"

#include "../wsclean/cachedimageset.h"

#include <boost/test/unit_test.hpp>

struct DynSetFixture
{
	DynSetFixture()
	{
		addToDynamicSet(table, 0, 0, 0, 0, Polarization::XX, 100);
		addToDynamicSet(table, 1, 0, 0, 0, Polarization::YY, 100);
		addToDynamicSet(table, 2, 0, 1, 1, Polarization::XX, 200);
		addToDynamicSet(table, 3, 0, 1, 1, Polarization::YY, 200);
		table.Update();
	}
	
	void addToDynamicSet(ImagingTable& table, size_t index, size_t j, size_t c, size_t s, PolarizationEnum p, size_t frequencyMHz)
	{
		ImagingTableEntry& e = table.AddEntry();
		e.index = index;
		e.joinedGroupIndex = j;
		e.outputChannelIndex = c;
		e.squaredDeconvolutionIndex = s;
		e.polarization = p;
		e.lowestFrequency = frequencyMHz;
		e.highestFrequency = frequencyMHz;
		e.bandStartFrequency = frequencyMHz;
		e.bandEndFrequency = frequencyMHz;
		e.imageCount = 1;
	}

	ImagingTable table;
	ImageBufferAllocator allocator;
};

BOOST_FIXTURE_TEST_SUITE(dynamicset, DynSetFixture)

BOOST_AUTO_TEST_CASE( squaredGroupCount )
{
	BOOST_CHECK_EQUAL(table.SquaredGroupCount(), 2);
}

BOOST_AUTO_TEST_CASE( entryCount )
{
	BOOST_CHECK_EQUAL(table.EntryCount(), 4);
}

BOOST_AUTO_TEST_CASE( entriesInGroup )
{
	BOOST_CHECK_EQUAL(table.GetSquaredGroup(0).EntryCount(), 2);
}

BOOST_AUTO_TEST_CASE( psfCount1 )
{
	DynamicSet dset(&table, allocator, 1, false, 2, 2);
	BOOST_CHECK_EQUAL(dset.PSFCount(), 1);
}

BOOST_AUTO_TEST_CASE( psfCount2 )
{
	DynamicSet dset(&table, allocator, 2, false, 2, 2);
	BOOST_CHECK_EQUAL(dset.PSFCount(), 2);
}

struct AdvDynamicSetFixture : public DynSetFixture
{
	FitsWriter writer;
	CachedImageSet cSet;
	ao::uvector<double> image;
	
	AdvDynamicSetFixture() :
		image(4, 0.0)
	{
		writer.SetImageDimensions(2, 2);
		cSet.Initialize(writer, 2, 2, "wsctest", allocator);
		image[0] = 2.0;
		cSet.Store(image.data(), Polarization::XX, 0, false);
		image[0] = -1.0;
		cSet.Store(image.data(), Polarization::YY, 0, false);
		image[0] = 20.0;
		cSet.Store(image.data(), Polarization::XX, 1, false);
		image[0] = -10.0;
		cSet.Store(image.data(), Polarization::YY, 1, false);
	}
};

BOOST_FIXTURE_TEST_CASE( load , AdvDynamicSetFixture)
{
	cSet.Load(image.data(), Polarization::XX, 1, false);
	BOOST_CHECK_EQUAL(image[0], 20.0);
}
	
BOOST_FIXTURE_TEST_CASE( loadAndAverage , AdvDynamicSetFixture)
{
	DynamicSet dset(&table, allocator, 1, false, 2, 2);
	dset.LoadAndAverage(cSet);
	BOOST_CHECK_CLOSE_FRACTION(dset[0][0], 0.5*( 2.0 + 20.0), 1e-8);
	BOOST_CHECK_CLOSE_FRACTION(dset[1][0], 0.5*(-1.0 - 10.0), 1e-8);
}

BOOST_FIXTURE_TEST_CASE( interpolateAndStore , AdvDynamicSetFixture)
{
	DynamicSet dset(&table, allocator, 2, false, 2, 2);
	SpectralFitter fitter(NoSpectralFitting, 2);
	dset.LoadAndAverage(cSet);
	dset.InterpolateAndStore(cSet, fitter);
	BOOST_CHECK_CLOSE_FRACTION(dset[0][0], 2.0, 1e-8);
	BOOST_CHECK_CLOSE_FRACTION(dset[1][0],-1.0, 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
