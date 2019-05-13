#include "../deconvolution/imageset.h"
#include "../deconvolution/spectralfitter.h"

#include "../wsclean/cachedimageset.h"
#include "../image.h"

#include <boost/test/unit_test.hpp>

struct ImageSetFixtureBase
{
	void addToImageSet(ImagingTable& table, size_t index, size_t joinedGroup, size_t outChannel, size_t squaredIndex, PolarizationEnum pol, size_t frequencyMHz, double imageWeight = 1.0)
	{
		ImagingTableEntry& e = table.AddEntry();
		e.index = index;
		e.joinedGroupIndex = joinedGroup;
		e.outputChannelIndex = outChannel;
		e.squaredDeconvolutionIndex = squaredIndex;
		e.polarization = pol;
		e.lowestFrequency = frequencyMHz;
		e.highestFrequency = frequencyMHz;
		e.bandStartFrequency = frequencyMHz;
		e.bandEndFrequency = frequencyMHz;
		e.imageCount = 1;
		e.imageWeight = imageWeight;
	}
	
	void checkLinearValue(size_t index, double value, const ImageSet& dset)
	{
		Image dest(2, 2, 1.0, allocator);
		dset.GetLinearIntegrated(dest.data());
		BOOST_CHECK_CLOSE_FRACTION(dest[index], value, 1e-6);
	}
	
	void checkSquaredValue(size_t index, double value, const ImageSet& dset)
	{
		Image dest(2, 2, 1.0, allocator), scratch(2, 2, allocator);
		dset.GetSquareIntegrated(dest.data(), scratch.data());
		BOOST_CHECK_CLOSE_FRACTION(dest[index], value, 1e-6);
	}
	
	ImagingTable table;
	ImageBufferAllocator allocator;
	WSCleanSettings settings;
};

struct ImageSetFixture : public ImageSetFixtureBase
{
	ImageSetFixture()
	{
		settings.deconvolutionChannelCount = 1;
		settings.squaredJoins = false;
		settings.linkedPolarizations = std::set<PolarizationEnum>();
	
		addToImageSet(table, 0, 0, 0, 0, Polarization::XX, 100);
		addToImageSet(table, 1, 0, 0, 0, Polarization::YY, 100);
		addToImageSet(table, 2, 0, 1, 1, Polarization::XX, 200);
		addToImageSet(table, 3, 0, 1, 1, Polarization::YY, 200);
		table.Update();
	}
	
};

BOOST_FIXTURE_TEST_SUITE(imageset, ImageSetFixture)

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
	WSCleanSettings settings;
	settings.deconvolutionChannelCount = 1;
	settings.squaredJoins = false;
	settings.linkedPolarizations = std::set<PolarizationEnum>();
	ImageSet dset(&table, allocator, settings, 2, 2);
	BOOST_CHECK_EQUAL(dset.PSFCount(), 1);
}

BOOST_AUTO_TEST_CASE( psfCount2 )
{
	WSCleanSettings settings;
	settings.deconvolutionChannelCount = 2;
	settings.squaredJoins = false;
	settings.linkedPolarizations = std::set<PolarizationEnum>();
	ImageSet dset(&table, allocator, settings, 2, 2);
	BOOST_CHECK_EQUAL(dset.PSFCount(), 2);
}

struct AdvImageSetFixture : public ImageSetFixture
{
	FitsWriter writer;
	CachedImageSet cSet;
	ao::uvector<double> image;
	
	AdvImageSetFixture() :
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

BOOST_FIXTURE_TEST_CASE( load , AdvImageSetFixture)
{
	cSet.Load(image.data(), Polarization::XX, 1, false);
	BOOST_CHECK_EQUAL(image[0], 20.0);
}
	
BOOST_FIXTURE_TEST_CASE( loadAndAverage , AdvImageSetFixture)
{
	ImageSet dset(&table, allocator, settings, 2, 2);
	dset.LoadAndAverage(cSet);
	BOOST_CHECK_CLOSE_FRACTION(dset[0][0], 0.5*( 2.0 + 20.0), 1e-8);
	BOOST_CHECK_CLOSE_FRACTION(dset[1][0], 0.5*(-1.0 - 10.0), 1e-8);
}

BOOST_FIXTURE_TEST_CASE( interpolateAndStore , AdvImageSetFixture)
{
	settings.deconvolutionChannelCount = 2;
	ImageSet dset(&table, allocator, settings, 2, 2);
	SpectralFitter fitter(NoSpectralFitting, 2);
	dset.LoadAndAverage(cSet);
	dset.InterpolateAndStore(cSet, fitter);
	BOOST_CHECK_CLOSE_FRACTION(dset[0][0], 2.0, 1e-8);
	BOOST_CHECK_CLOSE_FRACTION(dset[1][0],-1.0, 1e-8);
}

BOOST_FIXTURE_TEST_CASE( xxNormalization , ImageSetFixtureBase )
{
	settings.linkedPolarizations = std::set<PolarizationEnum>{ Polarization::XX };
	addToImageSet(table, 0, 0, 0, 0, Polarization::XX, 100);
	table.Update();
	ImageSet dset(&table, allocator, settings, 2, 2);
	dset[0][1] = 5.0;
	checkLinearValue(1, 5.0, dset);
	checkSquaredValue(1, 5.0, dset);
}

BOOST_FIXTURE_TEST_CASE( iNormalization , ImageSetFixtureBase )
{
	addToImageSet(table, 0, 0, 0, 0, Polarization::StokesI, 100);
	table.Update();
	settings.linkedPolarizations = std::set<PolarizationEnum>{ Polarization::StokesI };
	ImageSet dset(&table, allocator, settings, 2, 2);
	dset[0][2] = 6.0;
	checkLinearValue(2, 6.0, dset);
	checkSquaredValue(2, 6.0, dset);
}

BOOST_FIXTURE_TEST_CASE( i_2channel_Normalization , ImageSetFixtureBase )
{
	addToImageSet(table, 0, 0, 0, 0, Polarization::StokesI, 100);
	addToImageSet(table, 1, 0, 1, 1, Polarization::StokesI, 200);
	table.Update();
	settings.deconvolutionChannelCount = 2;
	settings.linkedPolarizations = std::set<PolarizationEnum>{ Polarization::StokesI };
	std::set<PolarizationEnum> pols{ Polarization::StokesI };
	ImageSet dset(&table, allocator, settings, 2, 2);
	dset[0][0] = 12.0;
	dset[1][0] = 13.0;
	checkLinearValue(0, 12.5, dset);
	checkSquaredValue(0, 12.5, dset);
}

BOOST_FIXTURE_TEST_CASE( xxyyNormalization , ImageSetFixtureBase )
{
	addToImageSet(table, 0, 0, 0, 0, Polarization::XX, 100);
	addToImageSet(table, 1, 0, 0, 0, Polarization::YY, 100);
	table.Update();
	settings.linkedPolarizations = std::set<PolarizationEnum>{ Polarization::XX, Polarization::YY };
	ImageSet dset(&table, allocator, settings, 2, 2);
	dset[0][3] = 7.0;
	dset[1][3] = 8.0;
	checkLinearValue(3, 7.5, dset);
	dset[0][3] = -7.0;
	checkSquaredValue(3, sqrt((7.0*7.0 + 8.0*8.0) * 0.5), dset);
}

BOOST_FIXTURE_TEST_CASE( iqNormalization , ImageSetFixtureBase )
{
	addToImageSet(table, 0, 0, 0, 0, Polarization::StokesI, 100);
	addToImageSet(table, 1, 0, 0, 0, Polarization::StokesQ, 100);
	table.Update();
	settings.linkedPolarizations = std::set<PolarizationEnum>{ Polarization::StokesI, Polarization::StokesQ };
	ImageSet dset(&table, allocator, settings, 2, 2);
	dset[0][0] = 6.0;
	dset[1][0] = -1.0;
	checkLinearValue(0, 5.0, dset);
	checkSquaredValue(0, sqrt(6.0*6.0 + -1.0*-1.0), dset);
}

BOOST_FIXTURE_TEST_CASE( linkedINormalization , ImageSetFixtureBase )
{
	addToImageSet(table, 0, 0, 0, 0, Polarization::StokesI, 100);
	addToImageSet(table, 1, 0, 0, 0, Polarization::StokesQ, 100);
	table.Update();
	settings.linkedPolarizations = std::set<PolarizationEnum>{ Polarization::StokesI };
	ImageSet dset(&table, allocator, settings, 2, 2);
	dset[0][0] = 3.0;
	dset[1][0] = -1.0;
	checkLinearValue(0, 3.0, dset);
	checkSquaredValue(0, 3.0, dset);
}

BOOST_FIXTURE_TEST_CASE( iquvNormalization , ImageSetFixtureBase )
{
	addToImageSet(table, 0, 0, 0, 0, Polarization::StokesI, 100);
	addToImageSet(table, 1, 0, 0, 0, Polarization::StokesQ, 100);
	addToImageSet(table, 2, 0, 0, 0, Polarization::StokesU, 100);
	addToImageSet(table, 3, 0, 0, 0, Polarization::StokesV, 100);
	table.Update();
	settings.linkedPolarizations = std::set<PolarizationEnum>{
		Polarization::StokesI, Polarization::StokesQ,
		Polarization::StokesU, Polarization::StokesV
	};
	ImageSet dset(&table, allocator, settings, 2, 2);
	dset[0][0] = 9.0;
	dset[1][0] = 0.2;
	dset[2][0] = 0.2;
	dset[3][0] = 0.2;
	checkLinearValue(0, 9.6, dset);
	checkSquaredValue(0, sqrt(9.0*9.0 + 3.0*0.2*0.2), dset);
}

BOOST_FIXTURE_TEST_CASE( xx_xy_yx_yyNormalization , ImageSetFixtureBase )
{
	addToImageSet(table, 0, 0, 0, 0, Polarization::XX, 100);
	addToImageSet(table, 1, 0, 0, 0, Polarization::XY, 100);
	addToImageSet(table, 2, 0, 0, 0, Polarization::YX, 100);
	addToImageSet(table, 3, 0, 0, 0, Polarization::YY, 100);
	table.Update();
	settings.linkedPolarizations = std::set<PolarizationEnum>{
		Polarization::XX, Polarization::XY,
		Polarization::YX, Polarization::YY
	};
	ImageSet dset(&table, allocator, settings, 2, 2);
	dset[0][1] = 10.0;
	dset[1][1] = 0.25;
	dset[2][1] = 0.25;
	dset[3][1] = 10.0;
	checkLinearValue(1, 10.25, dset);
	checkSquaredValue(1, sqrt((10.0*10.0*2.0 + 0.25*0.25*2.0)*0.5), dset);
}

BOOST_FIXTURE_TEST_CASE( xx_xy_yx_yy_2channel_Normalization , ImageSetFixtureBase )
{
	addToImageSet(table, 0, 0, 0, 0, Polarization::XX, 100);
	addToImageSet(table, 1, 0, 0, 0, Polarization::XY, 100);
	addToImageSet(table, 2, 0, 0, 0, Polarization::YX, 100);
	addToImageSet(table, 3, 0, 0, 0, Polarization::YY, 100);
	addToImageSet(table, 4, 0, 1, 1, Polarization::XX, 200);
	addToImageSet(table, 5, 0, 1, 1, Polarization::XY, 200);
	addToImageSet(table, 6, 0, 1, 1, Polarization::YX, 200);
	addToImageSet(table, 7, 0, 1, 1, Polarization::YY, 200);
	table.Update();
	settings.deconvolutionChannelCount = 2;
	settings.linkedPolarizations = std::set<PolarizationEnum>{
		Polarization::XX, Polarization::XY,
		Polarization::YX, Polarization::YY
	};
	ImageSet dset(&table, allocator, settings, 2, 2);
	dset[0][2] = 5.0;
	dset[1][2] = 0.1;
	dset[2][2] = 0.2;
	dset[3][2] = 6.0;
	dset[4][2] = 7.0;
	dset[5][2] = 0.3;
	dset[6][2] = 0.4;
	dset[7][2] = 8.0;
	double
		sqVal1 = 0.0,
		sqVal2 = 0.0;
	for(size_t i=0; i!=4; ++i)
	{
		sqVal1 += dset[i][2] * dset[i][2];
		sqVal2 += dset[i+4][2] * dset[i+4][2];
	}
	checkLinearValue(2, 27.0*0.25, dset);
	checkSquaredValue(2, (sqrt(sqVal1*0.5) + sqrt(sqVal2*0.5))*0.5, dset);
}

BOOST_FIXTURE_TEST_CASE( linked_xx_yy_2channel_Normalization , ImageSetFixtureBase )
{
	addToImageSet(table, 0, 0, 0, 0, Polarization::XX, 100);
	addToImageSet(table, 1, 0, 0, 0, Polarization::XY, 100);
	addToImageSet(table, 2, 0, 0, 0, Polarization::YX, 100);
	addToImageSet(table, 3, 0, 0, 0, Polarization::YY, 100);
	addToImageSet(table, 4, 0, 1, 1, Polarization::XX, 200);
	addToImageSet(table, 5, 0, 1, 1, Polarization::XY, 200);
	addToImageSet(table, 6, 0, 1, 1, Polarization::YX, 200);
	addToImageSet(table, 7, 0, 1, 1, Polarization::YY, 200);
	table.Update();
	settings.deconvolutionChannelCount = 2;
	settings.linkedPolarizations = std::set<PolarizationEnum>{
		Polarization::XX, Polarization::YY
	};
	ImageSet dset(&table, allocator, settings, 2, 2);
	dset[0][2] = 7.5;
	dset[1][2] = 0.1;
	dset[2][2] = -0.2;
	dset[3][2] = 6.5;
	dset[4][2] = 8.5;
	dset[5][2] = 0.3;
	dset[6][2] = -0.4;
	dset[7][2] = 9.5;
	double
		sqVal1 = dset[0][2] * dset[0][2] + dset[3][2] * dset[3][2],
		sqVal2 = dset[4][2] * dset[4][2] + dset[7][2] * dset[7][2];
	checkLinearValue(2, 32.0*0.25, dset);
	checkSquaredValue(2, (sqrt(sqVal1*0.5) + sqrt(sqVal2*0.5))*0.5, dset);
}

BOOST_FIXTURE_TEST_CASE( linked_xx_2channel_Normalization , ImageSetFixtureBase )
{
	addToImageSet(table, 0, 0, 0, 0, Polarization::XX, 100);
	addToImageSet(table, 1, 0, 0, 0, Polarization::XY, 100);
	addToImageSet(table, 2, 0, 0, 0, Polarization::YX, 100);
	addToImageSet(table, 3, 0, 0, 0, Polarization::YY, 100);
	addToImageSet(table, 4, 0, 1, 1, Polarization::XX, 200);
	addToImageSet(table, 5, 0, 1, 1, Polarization::XY, 200);
	addToImageSet(table, 6, 0, 1, 1, Polarization::YX, 200);
	addToImageSet(table, 7, 0, 1, 1, Polarization::YY, 200);
	table.Update();
	settings.deconvolutionChannelCount = 2;
	settings.linkedPolarizations = std::set<PolarizationEnum>{ Polarization::XX };
	ImageSet dset(&table, allocator, settings, 2, 2);
	dset[0][2] = 7.5;
	dset[1][2] = 0.1;
	dset[2][2] = -0.2;
	dset[3][2] = 6.5;
	dset[4][2] = 8.5;
	dset[5][2] = 0.3;
	dset[6][2] = -0.4;
	dset[7][2] = 9.5;
	double
		sqVal1 = dset[0][2] * dset[0][2],
		sqVal2 = dset[4][2] * dset[4][2];
	checkLinearValue(2, 32.0*0.25, dset);
	checkSquaredValue(2, (sqrt(sqVal1) + sqrt(sqVal2))*0.5, dset);
}

BOOST_FIXTURE_TEST_CASE( deconvchannels_normalization , ImageSetFixtureBase )
{
	addToImageSet(table, 0, 0, 0, 0, Polarization::StokesI, 100, 1);
	addToImageSet(table, 1, 0, 1, 1, Polarization::StokesI, 200, 1);
	addToImageSet(table, 2, 0, 2, 2, Polarization::StokesI, 300, 2);
	addToImageSet(table, 3, 0, 3, 3, Polarization::StokesI, 400, 2);
	table.Update();
	settings.deconvolutionChannelCount = 2;
	settings.linkedPolarizations = std::set<PolarizationEnum>{ Polarization::StokesI };
	ImageSet dset(&table, allocator, settings, 2, 2);
	dset[0][0] = 10.0;
	dset[1][0] = 13.0;
	checkLinearValue(0, 12.0, dset);
	checkSquaredValue(0, 12.0, dset);
}

BOOST_FIXTURE_TEST_CASE( deconvchannels_zeroweight , ImageSetFixtureBase )
{
	addToImageSet(table, 0, 0, 0, 0, Polarization::StokesI, 100, 1);
	addToImageSet(table, 1, 0, 1, 1, Polarization::StokesI, 200, 0);
	addToImageSet(table, 2, 0, 2, 2, Polarization::StokesI, 300, 2);
	addToImageSet(table, 3, 0, 3, 3, Polarization::StokesI, 400, 2);
	table.Update();
	settings.deconvolutionChannelCount = 2;
	settings.linkedPolarizations = std::set<PolarizationEnum>{ Polarization::StokesI };
	ImageSet dset(&table, allocator, settings, 2, 2);
	dset[0][0] = 10.0;
	dset[1][0] = 5.0;
	checkLinearValue(0, 6.0, dset);
	checkSquaredValue(0, 6.0, dset);
}

BOOST_FIXTURE_TEST_CASE( psfindex , ImageSetFixtureBase )
{
	for(size_t ch=0; ch!=4; ++ch)
	{
		addToImageSet(table, ch*4+0, 0, ch, ch, Polarization::XX, 100, 1);
		addToImageSet(table, ch*4+1, 0, ch, ch, Polarization::XY, 200, 0);
		addToImageSet(table, ch*4+2, 0, ch, ch, Polarization::YX, 300, 2);
		addToImageSet(table, ch*4+3, 0, ch, ch, Polarization::YY, 400, 2);
	}
	table.Update();
	settings.deconvolutionChannelCount = 2;
	settings.linkedPolarizations = std::set<PolarizationEnum>{
		Polarization::XX, Polarization::XY,
		Polarization::YX, Polarization::YY
	};
	ImageSet dset(&table, allocator, settings, 2, 2);
	
	BOOST_CHECK_EQUAL(dset.PSFIndex(0), 0);
	BOOST_CHECK_EQUAL(dset.PSFIndex(1), 0);
	BOOST_CHECK_EQUAL(dset.PSFIndex(2), 0);
	BOOST_CHECK_EQUAL(dset.PSFIndex(3), 0);
	BOOST_CHECK_EQUAL(dset.PSFIndex(4), 1);
	BOOST_CHECK_EQUAL(dset.PSFIndex(5), 1);
	BOOST_CHECK_EQUAL(dset.PSFIndex(6), 1);
	BOOST_CHECK_EQUAL(dset.PSFIndex(7), 1);
}

BOOST_AUTO_TEST_SUITE_END()
