#include "msgridderbase.h"
#include "logger.h"
#include "smallinversionoptimization.h"

#include "../msproviders/msprovider.h"

#include "../imageweights.h"

#include "../units/angle.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/tables/Tables/ArrColDesc.h>

MSGridderBase::MSData::MSData() : msIndex(0), matchingRows(0), totalRowsProcessed(0)
{ }

MSGridderBase::MSData::~MSData()
{ }

MSGridderBase::MSGridderBase() :
	MeasurementSetGridder(),
	_theoreticalBeamSize(0.0),
	_actualInversionWidth(0), _actualInversionHeight(0),
	_actualPixelSizeX(0), _actualPixelSizeY(0),
	_metaDataCache(nullptr),
	_hasFrequencies(false),
	_freqHigh(0.0), _freqLow(0.0),
	_bandStart(0.0), _bandEnd(0.0),
	_startTime(0.0),
	_phaseCentreRA(0.0), _phaseCentreDec(0.0),
	_phaseCentreDL(0.0), _phaseCentreDM(0.0),
	_denormalPhaseCentre(false),
	_griddedVisibilityCount(0),
	_totalWeight(0.0),
	_maxGriddedWeight(0.0),
	_visibilityWeightSum(0.0)
{ }

MSGridderBase::~MSGridderBase()
{ }

void MSGridderBase::GetPhaseCentreInfo(casacore::MeasurementSet& ms, size_t fieldId, double& ra, double& dec, double& dl, double& dm)
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
	ra = j2000Val[0];
	dec = j2000Val[1];
	if(fTable.keywordSet().isDefined("WSCLEAN_DL"))
		dl = fTable.keywordSet().asDouble(casacore::RecordFieldId("WSCLEAN_DL"));
	else dl = 0.0;
	if(fTable.keywordSet().isDefined("WSCLEAN_DM"))
		dm = fTable.keywordSet().asDouble(casacore::RecordFieldId("WSCLEAN_DM"));
	else dm = 0.0;
}

void MSGridderBase::initializePhaseCentre(casacore::MeasurementSet& ms, size_t fieldId)
{
	GetPhaseCentreInfo(ms, fieldId, _phaseCentreRA, _phaseCentreDec, _phaseCentreDL, _phaseCentreDM);
	
	_denormalPhaseCentre = _phaseCentreDL != 0.0 || _phaseCentreDM != 0.0;
	if(_denormalPhaseCentre)
		Logger::Debug << "Set has denormal phase centre: dl=" << _phaseCentreDL << ", dm=" << _phaseCentreDM << '\n';
}

void MSGridderBase::initializeBandData(casacore::MeasurementSet& ms, MSGridderBase::MSData& msData)
{
	msData.bandData = MultiBandData(ms.spectralWindow(), ms.dataDescription());
	if(Selection(msData.msIndex).HasChannelRange())
	{
		msData.startChannel = Selection(msData.msIndex).ChannelRangeStart();
		msData.endChannel = Selection(msData.msIndex).ChannelRangeEnd();
		Logger::Debug << "Selected channels: " << msData.startChannel << '-' << msData.endChannel << '\n';
		const BandData& firstBand = msData.bandData.FirstBand();
		if(msData.startChannel >= firstBand.ChannelCount() || msData.endChannel > firstBand.ChannelCount()
			|| msData.startChannel == msData.endChannel)
		{
			std::ostringstream str;
			str << "An invalid channel range was specified! Measurement set only has " << firstBand.ChannelCount() << " channels, requested imaging range is " << msData.startChannel << " -- " << msData.endChannel << '.';
			throw std::runtime_error(str.str());
		}
	}
	else {
		msData.startChannel = 0;
		msData.endChannel = msData.bandData.FirstBand().ChannelCount();
	}
}

template<size_t NPolInMSProvider>
void MSGridderBase::calculateWLimits(MSGridderBase::MSData& msData)
{
	Logger::Info << "Determining min and max w & theoretical beam size... ";
	Logger::Info.Flush();
	msData.maxW = 0.0;
	msData.maxWWithFlags = 0.0;
	msData.minW = 1e100;
	msData.maxBaselineUVW = 0.0;
	msData.maxBaselineInM = 0.0;
	MultiBandData selectedBand = msData.SelectedBand();
	std::vector<float> weightArray(selectedBand.MaxChannels() * NPolInMSProvider);
	msData.msProvider->Reset();
	while(msData.msProvider->CurrentRowAvailable())
	{
		size_t dataDescId;
		double uInM, vInM, wInM;
		msData.msProvider->ReadMeta(uInM, vInM, wInM, dataDescId);
		const BandData& curBand = selectedBand[dataDescId];
		double wHi = fabs(wInM / curBand.SmallestWavelength());
		double wLo = fabs(wInM / curBand.LongestWavelength());
		double baselineInM = sqrt(uInM*uInM + vInM*vInM + wInM*wInM);
		double halfWidth = 0.5*ImageWidth(), halfHeight = 0.5*ImageHeight();
		if(wHi > msData.maxW || wLo < msData.minW || baselineInM / curBand.SmallestWavelength() > msData.maxBaselineUVW)
		{
			msData.msProvider->ReadWeights(weightArray.data());
			const float* weightPtr = weightArray.data();
			for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
			{
				const double wavelength = curBand.ChannelWavelength(ch);
				double wInL = wInM/wavelength;
				msData.maxWWithFlags = std::max(msData.maxWWithFlags, fabs(wInL));
				if(*weightPtr != 0.0)
				{
					double
						uInL = uInM/wavelength, vInL = vInM/wavelength,
						x = uInL * PixelSizeX() * ImageWidth(),
						y = vInL * PixelSizeY() * ImageHeight(),
						imagingWeight = PrecalculatedWeightInfo()->GetWeight(uInL, vInL);
					if(imagingWeight != 0.0)
					{
						if(floor(x) > -halfWidth  && ceil(x) < halfWidth &&
							floor(y) > -halfHeight && ceil(y) < halfHeight)
						{
							msData.maxW = std::max(msData.maxW, fabs(wInL));
							msData.minW = std::min(msData.minW, fabs(wInL));
							msData.maxBaselineUVW = std::max(msData.maxBaselineUVW, baselineInM / wavelength);
							msData.maxBaselineInM = std::max(msData.maxBaselineInM, baselineInM);
						}
					}
				}
				weightPtr += NPolInMSProvider;
			}
		}
		
		msData.msProvider->NextRow();
	}
	
	if(msData.minW == 1e100)
	{
		msData.minW = 0.0;
		msData.maxWWithFlags = 0.0;
		msData.maxW = 0.0;
	}
	
	Logger::Info << "DONE (w=[" << msData.minW << ":" << msData.maxW << "] lambdas, maxuvw=" << msData.maxBaselineUVW << " lambda)\n";
	if(msData.maxWWithFlags != msData.maxW)
	{
		Logger::Debug << "Discarded data has higher w value of " << msData.maxWWithFlags << " lambda.\n";
	}
}

template void MSGridderBase::calculateWLimits<1>(MSGridderBase::MSData& msData);
template void MSGridderBase::calculateWLimits<4>(MSGridderBase::MSData& msData);

void MSGridderBase::initializeMSDataVector(std::vector<MSGridderBase::MSData>& msDataVector)
{
	if(MeasurementSetCount() == 0)
		throw std::runtime_error("Something is wrong during inversion: no measurement sets given to inversion algorithm");
	msDataVector = std::vector<MSGridderBase::MSData>(MeasurementSetCount());
	
	resetMetaData();
	
	bool hasCache = !_metaDataCache->msDataVector.empty();
	if(!hasCache)
		_metaDataCache->msDataVector.resize(MeasurementSetCount());
	for(size_t i=0; i!=MeasurementSetCount(); ++i)
	{
		msDataVector[i].msIndex = i;
		initializeMeasurementSet(msDataVector[i], _metaDataCache->msDataVector[i], hasCache);
	}
	
	calculateOverallMetaData(msDataVector.data());
}

void MSGridderBase::initializeMetaData(casacore::MeasurementSet& ms, size_t fieldId)
{
	casacore::MSObservation oTable = ms.observation();
	size_t obsCount = oTable.nrow();
	if(obsCount == 0) throw std::runtime_error("No observations in set");
	casacore::ROScalarColumn<casacore::String> telescopeNameColumn(oTable, oTable.columnName(casacore::MSObservation::TELESCOPE_NAME));
	casacore::ROScalarColumn<casacore::String> observerColumn(oTable, oTable.columnName(casacore::MSObservation::OBSERVER));
	_telescopeName = telescopeNameColumn(0);
	_observer = observerColumn(0);
	
	casacore::MSField fTable = ms.field();
	casacore::ROScalarColumn<casacore::String> fieldNameColumn(fTable, fTable.columnName(casacore::MSField::NAME));
	_fieldName = fieldNameColumn(fieldId);
}

void MSGridderBase::initializeMeasurementSet(MSGridderBase::MSData& msData, MetaDataCache::Entry& cacheEntry, bool isCacheInitialized)
{
	MSProvider& msProvider = MeasurementSet(msData.msIndex);
	msData.msProvider = &msProvider;
	SynchronizedMS ms(msProvider.MS());
	if(ms->nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	initializeBandData(*ms, msData);
	
	calculateMSLimits(msData.SelectedBand(), msProvider.StartTime());
	
	initializePhaseCentre(*ms, Selection(msData.msIndex).FieldId());
	
	initializeMetaData(*ms, Selection(msData.msIndex).FieldId());
	
	if(isCacheInitialized)
	{
		msData.maxW = cacheEntry.maxW;
		msData.maxWWithFlags = cacheEntry.maxWWithFlags;
		msData.minW = cacheEntry.minW;
		msData.maxBaselineUVW = cacheEntry.maxBaselineUVW;
		msData.maxBaselineInM = cacheEntry.maxBaselineInM;
	}
	else {
		if (msProvider.Polarization() == Polarization::Instrumental)
			calculateWLimits<4>(msData);
		else
			calculateWLimits<1>(msData);
		cacheEntry.maxW = msData.maxW;
		cacheEntry.maxWWithFlags = msData.maxWWithFlags;
		cacheEntry.minW = msData.minW;
		cacheEntry.maxBaselineUVW = msData.maxBaselineUVW;
		cacheEntry.maxBaselineInM = msData.maxBaselineInM;
	}
}

void MSGridderBase::calculateOverallMetaData(const MSData* msDataVector)
{
	_maxW = 0.0;
	_minW = std::numeric_limits<double>::max();
	double maxBaseline = 0.0;
	
	for(size_t i=0; i!=MeasurementSetCount(); ++i)
	{
		const MSData& msData = msDataVector[i];
		
		maxBaseline = std::max(maxBaseline, msData.maxBaselineUVW);
		_maxW = std::max(_maxW, msData.maxW);
		_minW = std::min(_minW, msData.minW);
	}
	if(_minW > _maxW)
	{
		_minW = _maxW;
		Logger::Error << "*** Error! ***\n"
			"*** Calculating maximum and minimum w values failed! Make sure the data selection and scale settings are correct!\n"
			"***\n";
	}
	
	_theoreticalBeamSize = 1.0 / maxBaseline;
	Logger::Info << "Theoretic beam = " << Angle::ToNiceString(_theoreticalBeamSize) << "\n";
	if(HasWLimit()) {
		_maxW *= (1.0 - WLimit());
		if(_maxW < _minW) _maxW = _minW;
	}

	if(!HasTrimSize())
		SetTrimSize(ImageWidth(), ImageHeight());
	
	_actualInversionWidth = ImageWidth();
	_actualInversionHeight = ImageHeight();
	_actualPixelSizeX = PixelSizeX();
	_actualPixelSizeY = PixelSizeY();
	
	if(SmallInversion())
	{
		size_t optWidth, optHeight, minWidth, minHeight;
		SmallInversionOptimization::DetermineOptimalSize(_actualInversionWidth, _actualPixelSizeX, _theoreticalBeamSize, minWidth, optWidth);
		SmallInversionOptimization::DetermineOptimalSize(_actualInversionHeight, _actualPixelSizeY, _theoreticalBeamSize, minHeight, optHeight);
		if(optWidth < _actualInversionWidth || optHeight < _actualInversionHeight)
		{
			size_t newWidth = std::max(std::min(optWidth, _actualInversionWidth), size_t(32));
			size_t newHeight = std::max(std::min(optHeight, _actualInversionHeight), size_t(32));
			Logger::Info << "Minimal inversion size: " << minWidth << " x " << minHeight << ", using optimal: " << newWidth << " x " << newHeight << "\n";
			_actualPixelSizeX = (double(_actualInversionWidth) * _actualPixelSizeX) / double(newWidth);
			_actualPixelSizeY = (double(_actualInversionHeight) * _actualPixelSizeY) / double(newHeight);
			_actualInversionWidth = newWidth;
			_actualInversionHeight = newHeight;
		}
		else {
			Logger::Info << "Small inversion enabled, but inversion resolution already smaller than beam size: not using optimization.\n";
		}
	}
	
	if(Verbose() || !HasWGridSize())
	{
		size_t suggestedGridSize = getSuggestedWGridSize();
		if(!HasWGridSize())
			SetActualWGridSize(suggestedGridSize);
		else
			SetActualWGridSize(WGridSize());
	}
	else 
		SetActualWGridSize(WGridSize());
}

template<size_t PolarizationCount>
void MSGridderBase::readAndWeightVisibilities(MSProvider& msProvider, InversionRow& rowData, const BandData& curBand, float* weightBuffer, std::complex<float>* modelBuffer, const bool* isSelected)
{
	if(DoImagePSF())
	{
		for(size_t chp=0; chp!=curBand.ChannelCount() * PolarizationCount; ++chp)
			rowData.data[chp] = 1.0;
		if(HasDenormalPhaseCentre())
		{
			double lmsqrt = sqrt(1.0-PhaseCentreDL()*PhaseCentreDL()- PhaseCentreDM()*PhaseCentreDM());
			double shiftFactor = 2.0*M_PI* (rowData.uvw[2] * (lmsqrt-1.0));
			rotateVisibilities<PolarizationCount>(curBand, shiftFactor, rowData.data);
		}
	}
	else {
		msProvider.ReadData(rowData.data);
	}
	rowData.rowId = msProvider.RowId();
	
	if(DoSubtractModel())
	{
		msProvider.ReadModel(modelBuffer);
		std::complex<float>* modelIter = modelBuffer;
		for(std::complex<float>* iter = rowData.data; iter!=rowData.data+(curBand.ChannelCount()*PolarizationCount); ++iter)
		{
			*iter -= *modelIter;
			modelIter++;
		}
	}
	
	msProvider.ReadWeights(weightBuffer);
	
	// Any visibilities that are not gridded in this pass
	// should not contribute to the weight sum, so set these
	// to have zero weight.
	for(size_t ch=0; ch!=curBand.ChannelCount()*PolarizationCount; ++ch)
	{
		if(!isSelected[ch])
			weightBuffer[ch] = 0.0;
	}
	
	switch(VisibilityWeightingMode())
	{
	case NormalVisibilityWeighting:
		// The weight buffer already contains the visibility weights: do nothing
		break;
	case SquaredVisibilityWeighting:
		// Square the visibility weights
		for(size_t chp=0; chp!=curBand.ChannelCount() * PolarizationCount; ++chp)
			weightBuffer[chp] *= weightBuffer[chp];
		break;
	case UnitVisibilityWeighting:
		// Set the visibility weights to one
		for(size_t chp=0; chp!=curBand.ChannelCount() * PolarizationCount; ++chp)
		{
			if(weightBuffer[chp] != 0.0)
				weightBuffer[chp] = 1.0f;
		}
		break;
	}
	
	// Calculate imaging weights
	std::complex<float>* dataIter = rowData.data;
	float* weightIter = weightBuffer;
	_scratchWeights.resize(curBand.ChannelCount());
	for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
	{
		double
			u = rowData.uvw[0] / curBand.ChannelWavelength(ch),
			v = rowData.uvw[1] / curBand.ChannelWavelength(ch),
			imageWeight = PrecalculatedWeightInfo()->GetWeight(u, v);
		_scratchWeights[ch] = imageWeight;
		
		for(size_t p=0; p!=PolarizationCount; ++p)
		{
			double cumWeight = *weightIter * imageWeight;
			if(p == 0 && cumWeight != 0.0) {
				// Visibility weight sum is the sum of weights excluding imaging weights
				_visibilityWeightSum += *weightIter;
				_maxGriddedWeight = std::max(cumWeight, _maxGriddedWeight);
				++_griddedVisibilityCount;
				// Total weight includes imaging weights
				_totalWeight += cumWeight;
			}
			*weightIter = cumWeight;
			*dataIter *= *weightIter;
			++dataIter;
			++weightIter;
		}
	}
	if(StoreImagingWeights())
		msProvider.WriteImagingWeights(rowData.rowId, _scratchWeights.data());
}

template void MSGridderBase::readAndWeightVisibilities<1>(MSProvider& msProvider, InversionRow& newItem, const BandData& curBand, float* weightBuffer, std::complex<float>* modelBuffer, const bool* isSelected);

template void MSGridderBase::readAndWeightVisibilities<4>(MSProvider& msProvider, InversionRow& newItem, const BandData& curBand, float* weightBuffer, std::complex<float>* modelBuffer, const bool* isSelected);

template<size_t PolarizationCount>
void MSGridderBase::rotateVisibilities(const BandData& bandData, double shiftFactor, std::complex<float>* dataIter)
{
	for(unsigned ch=0; ch!=bandData.ChannelCount(); ++ch)
	{
		const double wShiftRad = shiftFactor / bandData.ChannelWavelength(ch);
		float
			rotSin = sin(wShiftRad),
			rotCos = cos(wShiftRad);
		std::complex<float> v = *dataIter;
		*dataIter = std::complex<float>(
			v.real() * rotCos  -  v.imag() * rotSin,
			v.real() * rotSin  +  v.imag() * rotCos);
		++dataIter;
	}
}

template void MSGridderBase::rotateVisibilities<1>(const BandData& bandData, double shiftFactor, std::complex<float>* dataIter);
template void MSGridderBase::rotateVisibilities<4>(const BandData& bandData, double shiftFactor, std::complex<float>* dataIter);
