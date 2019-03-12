#ifndef TELESCOPE_H
#define TELESCOPE_H

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/tables/Tables/ScalarColumn.h>

#include <boost/algorithm/string/case_conv.hpp>

class Telescope {
public:
	enum TelescopeType {
		AARTFAAC, LOFAR, MWA
	};
	
	static TelescopeType GetType(casacore::MeasurementSet& ms)
	{
		return GetType(GetTelescopeName(ms));
	}
	
	static TelescopeType GetType(casacore::MeasurementSet&& ms)
	{
		return GetType(GetTelescopeName(ms));
	}
	
	static TelescopeType GetType(const std::string& telescopeName)
	{
		std::string up = boost::to_upper_copy(telescopeName);
		if(up == "LOFAR")
			return LOFAR;
		else if(up == "AARTFAAC")
			return AARTFAAC;
		else if(up == "MWA")
			return MWA;
		else
			throw std::runtime_error("Telescope name not recognized: " + telescopeName);
	}
	
	template<typename MS>
	static std::string GetTelescopeName(MS ms)
	{
		casacore::MSObservation obsTable = ms.observation();
		casacore::ScalarColumn<casacore::String> telescopeNameCol(obsTable, obsTable.columnName(casacore::MSObservationEnums::TELESCOPE_NAME));
		if(obsTable.nrow() != 0)
			return telescopeNameCol(0);
		else
			return std::string();
	}
};

#endif
