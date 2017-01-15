#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "../msproviders/averagingmsrowprovider.h"
#include "../msproviders/directmsrowprovider.h"

#include <boost/filesystem/operations.hpp>

#include "../wsclean/logger.h"

BOOST_AUTO_TEST_SUITE(baseline_dependent_averaging)

BOOST_AUTO_TEST_CASE( noAveraging )
{
	Logger::SetVerbosity(Logger::QuietVerbosity);
	std::string filename("/home/anoko/Data/3C196-TestSet/L258627_SAP000_SB120_uv.MS.1ch12s.dppp");
	if(boost::filesystem::exists(filename))
	{
		size_t nAnt=70;
		MSSelection selection;
		std::map<size_t, size_t> dataDescIds;
		dataDescIds.insert(std::make_pair(0, 0));
		AveragingMSRowProvider avgProvider(1e-8, filename, selection, dataDescIds, "DATA", false);
		DirectMSRowProvider directProvider(filename, selection, dataDescIds, "DATA", false);
		size_t nRow = 0, nFinite = 0;
		casacore::IPosition shape(2, 4, 1);
		MSRowProvider::DataArray dataArrayAvg(shape), dataArrayDirect(shape);
		MSRowProvider::FlagArray flagArrayAvg(shape), flagArrayDirect(shape);
		MSRowProvider::WeightArray weightsAvg(shape), weightsDirect(shape);
		while(!avgProvider.AtEnd() && !directProvider.AtEnd() && nRow<nAnt*(nAnt-1)/2 * 10)
		{
			double uA, vA, wA, uD, vD, wD;
			uint32_t a1A=1337, a2A=1338, a1D=1339, a2D=1340;
			uint32_t descA=1341, descD=1342;
			avgProvider.ReadData(dataArrayAvg, flagArrayAvg, weightsAvg, uA, vA, wA, descA, a1A, a2A);
			directProvider.ReadData(dataArrayDirect, flagArrayDirect, weightsDirect, uD, vD, wD, descD, a1D, a2D);
					
			for(size_t p=0; p!=4; ++p)
			{
				std::complex<float> valA = dataArrayAvg.data()[p];
				std::complex<float> valD = dataArrayDirect.data()[p];
				bool isAFinite = std::isfinite(valA.real()) && std::isfinite(valA.imag());
				bool isDFinite = std::isfinite(valD.real()) && std::isfinite(valD.imag());
				
				BOOST_REQUIRE_EQUAL(isAFinite, isDFinite);
				if(isAFinite) {
					BOOST_REQUIRE_CLOSE_FRACTION(valA, valD, 1e-8);
					++nFinite;
				}
				
				BOOST_REQUIRE_EQUAL(flagArrayAvg.data()[p], flagArrayDirect.data()[p]);
				BOOST_REQUIRE_CLOSE_FRACTION(weightsAvg.data()[p], weightsDirect.data()[p], 1e-8);
			}
			
			BOOST_REQUIRE_CLOSE_FRACTION(uA, uD, 1e-8);
			BOOST_REQUIRE_CLOSE_FRACTION(vA, vD, 1e-8);
			BOOST_REQUIRE_CLOSE_FRACTION(wA, wD, 1e-8);
			
			BOOST_REQUIRE_EQUAL(descA, descD);
			
			BOOST_REQUIRE_EQUAL(a1A, a1D);
			BOOST_REQUIRE_EQUAL(a2A, a2D);
			
			avgProvider.NextRow();
			directProvider.NextRow();
			++nRow;
		}
		
		BOOST_CHECK_GT(nRow, 0);
		BOOST_CHECK_GT(nFinite, 0);
		BOOST_CHECK_EQUAL(avgProvider.AtEnd(), directProvider.AtEnd());
	}
	else {
		BOOST_TEST_WARN(false, "Test set not found -- skipping test averagingRowProvider");
	}
}

/*
BOOST_AUTO_TEST_CASE( extremeAveraging )
{
	std::string filename("/home/anoko/Data/3C196-TestSet/L258627_SAP000_SB120_uv.MS.1ch12s.dppp");
	size_t nAnt=70;
	if(boost::filesystem::exists(filename))
	{
		MSSelection selection;
		std::map<size_t, size_t> dataDescIds;
		dataDescIds.insert(std::make_pair(0, 0));
		AveragingMSRowProvider avgProvider(1e8, filename, selection, dataDescIds, "DATA", false);
		size_t nRow = 0;
		casacore::IPosition shape(2, 4, 1);
		MSRowProvider::DataArray dataArrayAvg(shape), dataArrayDirect(shape);
		MSRowProvider::FlagArray flagArrayAvg(shape), flagArrayDirect(shape);
		MSRowProvider::WeightArray weightsAvg(shape), weightsDirect(shape);
		std::set<std::pair<size_t,size_t>> baselines;
		while(!avgProvider.AtEnd())
		{
			double u, v, w;
			uint32_t a1=1337, a2=1338;
			uint32_t desc=1341;
			avgProvider.ReadData(dataArrayAvg, flagArrayAvg, weightsAvg, u, v, w, desc, a1, a2);
			std::pair<size_t, size_t> b(a1, a2);
			BOOST_CHECK(baselines.find(b) == baselines.end());
			baselines.insert(b);
					
			avgProvider.NextRow();
			++nRow;
		}
		
		BOOST_CHECK_EQUAL(nRow, nAnt*(nAnt-1)/2);
		BOOST_CHECK(avgProvider.AtEnd());
	}
	else {
		BOOST_TEST_WARN(false, "Test set not found -- skipping test averagingRowProvider");
	}
}*/

BOOST_AUTO_TEST_SUITE_END()
