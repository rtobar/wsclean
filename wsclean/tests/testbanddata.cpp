#include <boost/test/unit_test.hpp>

#include "../banddata.h"

BOOST_AUTO_TEST_SUITE(banddata)

BOOST_AUTO_TEST_CASE( construction )
{
	BandData bandData1;
	BOOST_CHECK_EQUAL(bandData1.ChannelCount(), 0);
	
	// Check that these methods don't crash:
	bandData1.CentreFrequency();
	bandData1.CentreWavelength();
	bandData1.FrequencyStep();
	bandData1.HighestFrequency();
	bandData1.LongestWavelength();
	bandData1.LowestFrequency();
	bandData1.SmallestWavelength();
	
	std::vector<ChannelInfo> channels(3, ChannelInfo(150e6, 10e6));
	channels[1] = ChannelInfo(160e6, 10e6);
	channels[2] = ChannelInfo(170e6, 10e6);
	BandData bandData2(channels);
	BOOST_CHECK_EQUAL(bandData2.ChannelCount(), 3);
	BOOST_CHECK_CLOSE_FRACTION(bandData2.ChannelFrequency(0), 150e6, 1e-5);
	BOOST_CHECK_CLOSE_FRACTION(bandData2.ChannelFrequency(1), 160e6, 1e-5);
	BOOST_CHECK_CLOSE_FRACTION(bandData2.ChannelFrequency(2), 170e6, 1e-5);
	BOOST_CHECK_CLOSE_FRACTION(bandData2.FrequencyStep(), 10e6, 1e-5);
	bandData2 = bandData1;
	BOOST_CHECK_EQUAL(bandData2.ChannelCount(), 0);
	bandData2 = BandData(channels);
	BOOST_CHECK_EQUAL(bandData2.ChannelCount(), 3);
	BandData bandData3(std::move(bandData2));
	BOOST_CHECK_EQUAL(bandData3.ChannelCount(), 3);
	BandData bandData4(bandData3);
	BOOST_CHECK_EQUAL(bandData4.ChannelCount(), 3);
	BOOST_CHECK_CLOSE_FRACTION(bandData4.ChannelFrequency(0), 150e6, 1e-5);
	BOOST_CHECK_CLOSE_FRACTION(bandData4.ChannelFrequency(1), 160e6, 1e-5);
	BOOST_CHECK_CLOSE_FRACTION(bandData4.ChannelFrequency(2), 170e6, 1e-5);
}

BOOST_AUTO_TEST_SUITE_END()

