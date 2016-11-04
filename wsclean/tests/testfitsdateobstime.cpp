#include "../fitsreader.h"
#include "../fitswriter.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include <string>

BOOST_AUTO_TEST_SUITE(fits_date_obs_time)

const std::pair<std::string, double> timeValues[] = {
	{"2013-08-22T00:00:00.0", 0.0},
	{"2013-08-22T00:00:00.1", 0.1/60/60/24},
	{"2013-08-22T00:00:00.9", 0.9/60/60/24},
	{"2013-08-22T00:00:59.9", 59.9/60/60/24},
	{"2013-08-22T00:59:59.9", (59.9/60 + 59)/60/24},
	{"2013-08-22T23:59:59.9", ((59.9/60 + 59)/60 + 23)/24},
	{"2013-08-22T23:59:59.9", ((59.99/60 + 59)/60 + 23)/24},
	{"2013-08-22T00:00:00.0", 1.0},
	{"2013-08-22T14:00:16.0", (((16 - 0.00000000000069)/60 + 00)/60 + 14)/24}
};

BOOST_DATA_TEST_CASE( functionMJDToHMS, boost::unit_test::data::xrange(std::end(timeValues)-std::begin(timeValues)) )
{
	const std::string& str = timeValues[sample].first;
	double mjd = timeValues[sample].second;

	int hour, min, sec, dsec;
	FitsWriter::MJDToHMS(mjd, hour, min, sec, dsec);
	char dateStrA[40], dateStrB[40];
  std::sprintf(dateStrA, "%02d:%02d:%02d.%01d", hour, min, sec, dsec);
	BOOST_CHECK_EQUAL(str.substr(11), std::string(dateStrA));
	
	double fMJD = FitsReader::ParseFitsDateToMJD(str.c_str());
	FitsWriter::MJDToHMS(fMJD, hour, min, sec, dsec);
  std::sprintf(dateStrB, "%02d:%02d:%02d.%01d", hour, min, sec, dsec);
	BOOST_CHECK_EQUAL(str.substr(11), std::string(dateStrB));
}

BOOST_AUTO_TEST_SUITE_END()

