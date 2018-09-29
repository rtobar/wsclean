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
		std::string telescopeName = GetTelescopeName(ms);
		boost::to_upper(telescopeName);
		
		if(telescopeName == "LOFAR")
			return LOFAR;
		else if(telescopeName == "AARTFAAC")
			return AARTFAAC;
		else if(telescopeName == "MWA")
			return MWA;
		else
			throw std::runtime_error("Telescope name not recognized: " + telescopeName);
	}
	
	static std::string GetTelescopeName(casacore::MeasurementSet& ms)
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
