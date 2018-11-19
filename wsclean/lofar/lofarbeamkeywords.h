#ifndef LOFAR_BEAM_KEYWORDS
#define LOFAR_BEAM_KEYWORDS

#include "../wsclean/logger.h"
#include "../units/radeccoord.h"

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/TableRecord.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MeasureHolder.h>

class LOFARBeamKeywords
{
public:
	static bool GetPreappliedBeamDirection(casacore::MeasurementSet& ms, const std::string& dataColumnName, bool useDifferentialBeam, casacore::MDirection& preappliedBeamDir)
	{
		casacore::ScalarMeasColumn<casacore::MDirection> referenceDirColumn(ms.field(), casacore::MSField::columnName(casacore::MSFieldEnums::REFERENCE_DIR));
		preappliedBeamDir = referenceDirColumn(0);
		
		// Read beam keywords of input datacolumn
		casacore::ArrayColumn<std::complex<float>> dataCol(ms, dataColumnName);
		bool wasBeamApplied = false;
		if(dataCol.keywordSet().isDefined("LOFAR_APPLIED_BEAM_MODE"))
		{
			std::string mode = dataCol.keywordSet().asString("LOFAR_APPLIED_BEAM_MODE");
			if(mode == "None")
				wasBeamApplied = false;
			else {
				if(mode == "Element" || mode == "ArrayFactor")
					throw std::runtime_error("This observation was corrected for the " + mode + " beam. WSClean can only handle a full pre-applied beam (both arrayfactor + element).");
				else if(mode == "Full")
				{
					wasBeamApplied = true;
					casacore::String error;
					casacore::MeasureHolder mHolder;
					if(!mHolder.fromRecord(error, dataCol.keywordSet().asRecord("LOFAR_APPLIED_BEAM_DIR")))
						throw std::runtime_error("Error while reading LOFAR_APPLIED_BEAM_DIR keyword: " + error);
					preappliedBeamDir = mHolder.asMDirection();
				}
				else
					throw std::runtime_error("Measurement set specifies an unknown beam correction: " + mode);
				
			}
		}
		if(useDifferentialBeam)
		{
			if(wasBeamApplied)
				Logger::Warn <<
					"* This measurement set has keywords specifying the pre-applied beam mode & direction:"
					"  it is therefore unnecessary to specify '-use-differential-lofar-beam'.\n";
			else
				Logger::Warn <<
					"!!!\n"
					"!!! The differential beam is being applied, but the measurement specifies that no\n"
					"!!! beam was pre-applied. Either an old version of DP3 was used, or you are applying\n"
					"!!! the wrong beam!\n"
					"!!!\n";
		}
		if(wasBeamApplied || useDifferentialBeam)
		{
			useDifferentialBeam = true;
			double ra = preappliedBeamDir.getAngle().getValue()[0];
			double dec = preappliedBeamDir.getAngle().getValue()[1];
			Logger::Debug << "Applying differential beam from direction " << RaDecCoord::RaDecToString(ra, dec) << ".\n";
		}
		return useDifferentialBeam;
	}
};

#endif

