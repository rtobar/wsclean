#include "msgridderbase.h"
#include "logger.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/tables/Tables/ArrColDesc.h>

MSGridderBase::MSGridderBase() :
	MeasurementSetGridder(),
	_hasFrequencies(false),
	_freqHigh(0.0), _freqLow(0.0),
	_bandStart(0.0), _bandEnd(0.0),
	_startTime(0.0),
	_phaseCentreRA(0.0), _phaseCentreDec(0.0),
	_phaseCentreDL(0.0), _phaseCentreDM(0.0),
	_denormalPhaseCentre(false)
{
}

void MSGridderBase::initializePhaseCentre(casacore::MeasurementSet& ms, size_t fieldId)
{
	casacore::MSAntenna aTable = ms.antenna();
	size_t antennaCount = aTable.nrow();
	if(antennaCount == 0) throw std::runtime_error("No antennae in set");
	casacore::MPosition::ROScalarColumn antPosColumn(aTable, aTable.columnName(casacore::MSAntennaEnums::POSITION));
	casacore::MPosition ant1Pos = antPosColumn(0);
	
	casacore::MEpoch::ROScalarColumn timeColumn(ms, ms.columnName(casacore::MSMainEnums::TIME));
	casacore::MSField fTable(ms.field());
	casacore::MDirection::ROScalarColumn phaseDirColumn(fTable, fTable.columnName(casacore::MSFieldEnums::PHASE_DIR));
	casacore::MDirection phaseDir = phaseDirColumn(fieldId);
	casacore::MEpoch curtime = timeColumn(0);
	casacore::MeasFrame frame(ant1Pos, curtime);
	casacore::MDirection::Ref j2000Ref(casacore::MDirection::J2000, frame);
	casacore::MDirection j2000 = casacore::MDirection::Convert(phaseDir, j2000Ref)();
	casacore::Vector<casacore::Double> j2000Val = j2000.getValue().get();
	_phaseCentreRA = j2000Val[0];
	_phaseCentreDec = j2000Val[1];
	if(fTable.keywordSet().isDefined("WSCLEAN_DL"))
		_phaseCentreDL = fTable.keywordSet().asDouble(casacore::RecordFieldId("WSCLEAN_DL"));
	else _phaseCentreDL = 0.0;
	if(fTable.keywordSet().isDefined("WSCLEAN_DM"))
		_phaseCentreDM = fTable.keywordSet().asDouble(casacore::RecordFieldId("WSCLEAN_DM"));
	else _phaseCentreDM = 0.0;

	_denormalPhaseCentre = _phaseCentreDL != 0.0 || _phaseCentreDM != 0.0;
	if(_denormalPhaseCentre)
		Logger::Info << "Set has denormal phase centre: dl=" << _phaseCentreDL << ", dm=" << _phaseCentreDM << '\n';
}
