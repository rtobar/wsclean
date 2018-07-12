#include <boost/test/unit_test.hpp>

#include "../parsetreader.h"

BOOST_AUTO_TEST_SUITE(parsetreader)

BOOST_AUTO_TEST_CASE( normal_parset )
{
	const char* testparset =
		"# This is a test parset\n"
		"	aterms = [ tec, rm, amplitude ] # of phase, diagonal, fulljones, beam, ...\n"
		"	tec.type = tec\n"
		"	tec.images = [ tecimages.fits ]\n"
		"	rm.type = rm\n"
		"	rm.images = [ rmimages.fits ]\n"
		"	amplitude.type = amplitude\n"
		"	rm.images = [ amplitudesXX.fits amplitudesYY.fits ]\n"
		"# Comments are allowed\n"
		"	beam.differential = true\n";
	std::istringstream stream(testparset);
	ParsetReader reader(stream);
	
	BOOST_CHECK_EQUAL(reader.GetStringList("aterms").size(), 3);
	BOOST_CHECK_EQUAL(reader.GetStringList("aterms")[0], "tec");
	BOOST_CHECK_EQUAL(reader.GetStringList("aterms")[1], "rm");
	BOOST_CHECK_EQUAL(reader.GetStringList("aterms")[2], "amplitude");
	BOOST_CHECK_EQUAL(reader.GetString("tec.type"), "tec");
	BOOST_CHECK_EQUAL(reader.GetStringList("tec.images").size(), 1);
	BOOST_CHECK_EQUAL(reader.GetStringList("tec.images")[0], "tecimages.fits");
	BOOST_CHECK_EQUAL(reader.GetString("rm.type"), "rm");
	BOOST_CHECK_EQUAL(reader.GetString("amplitude.type"), "amplitude");
	BOOST_CHECK_EQUAL(reader.GetString("beam.differential"), "true");
	BOOST_CHECK_EQUAL(reader.GetBool("beam.differential"), true);
	
	// Or values
	BOOST_CHECK_EQUAL(reader.GetBoolOr("beam.differential", false), true);
	BOOST_CHECK_EQUAL(reader.GetBoolOr("notexisting", true), true);
	BOOST_CHECK_EQUAL(reader.GetBoolOr("notexisting", false), false);
	
	BOOST_CHECK_EQUAL(reader.GetStringOr("rm.type", "-"), "rm");
	BOOST_CHECK_EQUAL(reader.GetStringOr("notexisting", "!"), "!");
}

BOOST_AUTO_TEST_CASE( key_not_found)
{
	const char* testparset =
		"# This is a test parset\n"
		"	aterms = [ tec, rm, amplitude ] # of phase, diagonal, fulljones, beam, ...\n"
		"	tec.type = tec\n";
	std::istringstream stream(testparset);
	ParsetReader reader(stream);
	
	BOOST_CHECK_THROW(reader.GetString("notakey"), std::runtime_error);
	BOOST_CHECK_THROW(reader.GetStringList("notakey"), std::runtime_error);
	BOOST_CHECK_THROW(reader.GetBool("notakey"), std::runtime_error);
}

BOOST_AUTO_TEST_CASE( invalid_value_type )
{
	const char* testparset =
		"# This is a test parset\n"
		"	aterms = [ tec, rm, amplitude ] # of phase, diagonal, fulljones, beam, ...\n"
		"	tec.type = tec\n";
	std::istringstream stream(testparset);
	ParsetReader reader(stream);
	
	BOOST_CHECK_THROW(reader.GetString("aterms"), std::runtime_error);
	BOOST_CHECK_THROW(reader.GetStringList("tec.type"), std::runtime_error);
	BOOST_CHECK_THROW(reader.GetBool("tec.type"), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
