#include "msprovider.h"

#include "../wsclean/logger.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrColDesc.h>

#include "../msselection.h"

void MSProvider::copyWeightedData(std::complex<float>* dest, size_t startChannel, size_t endChannel, const std::vector<PolarizationEnum>& polsIn, const casacore::Array<std::complex<float>>& data, const casacore::Array<float>& weights, const casacore::Array<bool>& flags, PolarizationEnum polOut)
{
	const size_t polCount = polsIn.size();
	casacore::Array<std::complex<float> >::const_contiter inPtr = data.cbegin() + startChannel * polCount;
	casacore::Array<float>::const_contiter weightPtr = weights.cbegin() + startChannel * polCount;
	casacore::Array<bool>::const_contiter flagPtr = flags.cbegin() + startChannel * polCount;
	const size_t selectedChannelCount = endChannel - startChannel;
		
	size_t polIndex;
	if(polOut == Polarization::Instrumental)
	{
		for(size_t ch=0; ch!=selectedChannelCount*polsIn.size(); ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] = *inPtr * (*weightPtr);
			}
			else {
				dest[ch] = 0;
			}
			weightPtr++;
			inPtr++;
			flagPtr++;
		}
	}
	else if(Polarization::TypeToIndex(polOut, polsIn, polIndex)) {
		inPtr += polIndex;
		weightPtr += polIndex;
		flagPtr += polIndex;
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] = *inPtr * (*weightPtr);
			}
			else {
				dest[ch] = 0;
			}
			weightPtr += polCount;
			inPtr += polCount;
			flagPtr += polCount;
		}
	}
	else {
		// Copy the right visibilities with conversion if necessary.
		// Note that many conversions require dividing by two, e.g.
		// I = (XX + YY)/2. This division is done with weighting, so
		// weighted I = (w1 XX + w2 YY),
		// which is given a weight (w1 + w2),
		// and hence unweighted I = (w1 XX + w2 YY) / (w1 + w2), which if
		// w1 = w2 results in I = (XX + YY) / 2.
		switch(polOut)
		{
		case Polarization::StokesI: {
			size_t polIndexA=0, polIndexB=0;
			bool hasXX = Polarization::TypeToIndex(Polarization::XX, polsIn, polIndexA);
			bool hasYY = Polarization::TypeToIndex(Polarization::YY, polsIn, polIndexB);
			if(!hasXX || !hasYY)
			{
				bool hasRR = Polarization::TypeToIndex(Polarization::RR, polsIn, polIndexA);
				bool hasLL = Polarization::TypeToIndex(Polarization::LL, polsIn, polIndexB);
				if(!hasRR || !hasLL)
					throw std::runtime_error("Can not form requested polarization (Stokes I) from available polarizations");
			}
			
			for(size_t ch=0; ch!=selectedChannelCount; ++ch)
			{
				weightPtr += polIndexA;
				inPtr += polIndexA;
				flagPtr += polIndexA;
				
				bool flagA = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
				casacore::Complex valA = *inPtr * (*weightPtr);
				
				weightPtr += polIndexB - polIndexA;
				inPtr += polIndexB - polIndexA;
				flagPtr += polIndexB - polIndexA;
				
				bool flagB = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
				if(flagA || flagB)
					dest[ch] = 0.0;
				else {
					// I = XX + YY
					dest[ch] = (*inPtr * (*weightPtr)) + valA;
				}
				
				weightPtr += polCount - polIndexB;
				inPtr += polCount - polIndexB;
				flagPtr += polCount - polIndexB;
			}
		} break;
		case Polarization::StokesQ: {
			size_t polIndexA=0, polIndexB=0;
			bool hasXX = Polarization::TypeToIndex(Polarization::XX, polsIn, polIndexA);
			bool hasYY = Polarization::TypeToIndex(Polarization::YY, polsIn, polIndexB);
			if(hasXX && hasYY)
			{
				// Convert to StokesQ from XX and YY
				for(size_t ch=0; ch!=selectedChannelCount; ++ch)
				{
					weightPtr += polIndexA;
					inPtr += polIndexA;
					flagPtr += polIndexA;
					
					bool flagA = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					casacore::Complex valA = *inPtr * (*weightPtr);
					
					weightPtr += polIndexB - polIndexA;
					inPtr += polIndexB - polIndexA;
					flagPtr += polIndexB - polIndexA;
					
					bool flagB = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					if(flagA || flagB)
						dest[ch] = 0.0;
					else {
						// Q = (XX - YY)/2
						dest[ch] = valA - (*inPtr * (*weightPtr));
					}
					
					weightPtr += polCount - polIndexB;
					inPtr += polCount - polIndexB;
					flagPtr += polCount - polIndexB;
				}
			}
			else {
				// Convert to StokesQ from RR and LL
				bool hasRL = Polarization::TypeToIndex(Polarization::RL, polsIn, polIndexA);
				bool hasLR = Polarization::TypeToIndex(Polarization::LR, polsIn, polIndexB);
				if(!hasRL || !hasLR)
					throw std::runtime_error("Can not form requested polarization (Stokes Q) from available polarizations");
				for(size_t ch=0; ch!=selectedChannelCount; ++ch)
				{
					weightPtr += polIndexA;
					inPtr += polIndexA;
					flagPtr += polIndexA;
					
					bool flagA = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					casacore::Complex valA = *inPtr * (*weightPtr);
					
					weightPtr += polIndexB - polIndexA;
					inPtr += polIndexB - polIndexA;
					flagPtr += polIndexB - polIndexA;
					
					bool flagB = *flagPtr || !std::isfinite(inPtr->real()) || !std::isfinite(inPtr->imag());
					if(flagA || flagB)
						dest[ch] = 0.0;
					else {
						// Q = (RL + LR)/2
						dest[ch] = (*inPtr * (*weightPtr)) + valA;
					}
					
					weightPtr += polCount - polIndexB;
					inPtr += polCount - polIndexB;
					flagPtr += polCount - polIndexB;
				}
			}
		} break;
		case Polarization::StokesU: {
			size_t polIndexA=0, polIndexB=0;
			bool hasXY = Polarization::TypeToIndex(Polarization::XY, polsIn, polIndexA);
			bool hasYX = Polarization::TypeToIndex(Polarization::YX, polsIn, polIndexB);
			if(hasXY && hasYX)
			{
				// Convert to StokesU from XY and YX
				for(size_t ch=0; ch!=selectedChannelCount; ++ch)
				{
					weightPtr += polIndexA;
					inPtr += polIndexA;
					flagPtr += polIndexA;
					
					bool flagA = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					casacore::Complex valA = *inPtr * (*weightPtr);
					
					weightPtr += polIndexB - polIndexA;
					inPtr += polIndexB - polIndexA;
					flagPtr += polIndexB - polIndexA;
					
					bool flagB = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					if(flagA || flagB)
						dest[ch] = 0.0;
					else
						dest[ch] = valA + (*inPtr * (*weightPtr)); // U = (XY + YX)/2
					
					weightPtr += polCount - polIndexB;
					inPtr += polCount - polIndexB;
					flagPtr += polCount - polIndexB;
				}
			}
			else {
				// Convert to StokesU from RR and LL
				bool hasRL = Polarization::TypeToIndex(Polarization::RL, polsIn, polIndexA);
				bool hasLR = Polarization::TypeToIndex(Polarization::LR, polsIn, polIndexB);
				if(!hasRL || !hasLR)
					throw std::runtime_error("Can not form requested polarization (Stokes U) from available polarizations");
				for(size_t ch=0; ch!=selectedChannelCount; ++ch)
				{
					weightPtr += polIndexA;
					inPtr += polIndexA;
					flagPtr += polIndexA;
					
					bool flagA = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					casacore::Complex valA = *inPtr * (*weightPtr);
					
					weightPtr += polIndexB - polIndexA;
					inPtr += polIndexB - polIndexA;
					flagPtr += polIndexB - polIndexA;
					
					bool flagB = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					if(flagA || flagB)
						dest[ch] = 0.0;
					else {
						casacore::Complex diff = (valA - *inPtr * (*weightPtr));
						// U = -i (RL - LR)/2
						dest[ch] = casacore::Complex(diff.imag(), -diff.real());
					}
					
					weightPtr += polCount - polIndexB;
					inPtr += polCount - polIndexB;
					flagPtr += polCount - polIndexB;
				}
			}
		} break;
		case Polarization::StokesV: {
			size_t polIndexA=0, polIndexB=0;
			bool hasXY = Polarization::TypeToIndex(Polarization::XY, polsIn, polIndexA);
			bool hasYX = Polarization::TypeToIndex(Polarization::YX, polsIn, polIndexB);
			if(hasXY && hasYX)
			{
				// Convert to StokesV from XX and YY
				for(size_t ch=0; ch!=selectedChannelCount; ++ch)
				{
					weightPtr += polIndexA;
					inPtr += polIndexA;
					flagPtr += polIndexA;
					
					bool flagA = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					casacore::Complex valA = *inPtr * (*weightPtr);
					
					weightPtr += polIndexB - polIndexA;
					inPtr += polIndexB - polIndexA;
					flagPtr += polIndexB - polIndexA;
					
					bool flagB = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					if(flagA || flagB)
						dest[ch] = 0.0;
					else {
						casacore::Complex diff = valA - (*inPtr * (*weightPtr));
						// V = -i(XY - YX)/2
						dest[ch] = casacore::Complex(diff.imag(), -diff.real());
					}
					
					weightPtr += polCount - polIndexB;
					inPtr += polCount - polIndexB;
					flagPtr += polCount - polIndexB;
				}
			}
			else {
				// Convert to StokesV from RR and LL
				bool hasRL = Polarization::TypeToIndex(Polarization::RR, polsIn, polIndexA);
				bool hasLR = Polarization::TypeToIndex(Polarization::LL, polsIn, polIndexB);
				if(!hasRL || !hasLR)
					throw std::runtime_error("Can not form requested polarization (Stokes V) from available polarizations");
				for(size_t ch=0; ch!=selectedChannelCount; ++ch)
				{
					weightPtr += polIndexA;
					inPtr += polIndexA;
					flagPtr += polIndexA;
					
					bool flagA = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					casacore::Complex valA = *inPtr * (*weightPtr);
					
					weightPtr += polIndexB - polIndexA;
					inPtr += polIndexB - polIndexA;
					flagPtr += polIndexB - polIndexA;
					
					bool flagB = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					if(flagA || flagB)
						dest[ch] = 0.0;
					else {
						// V = (RR - LL)/2
						dest[ch] = valA - *inPtr * (*weightPtr);
					}
					
					weightPtr += polCount - polIndexB;
					inPtr += polCount - polIndexB;
					flagPtr += polCount - polIndexB;
				}
			}
		} break;
		default:
			throw std::runtime_error("Could not convert ms polarizations to requested polarization");
		}
	}
}

template<typename NumType>
void MSProvider::copyWeights(NumType* dest, size_t startChannel, size_t endChannel, const std::vector<PolarizationEnum>& polsIn, const casacore::Array<std::complex<float>>& data, const casacore::Array<float>& weights, const casacore::Array<bool>& flags, PolarizationEnum polOut)
{
	const size_t polCount = polsIn.size();
	casacore::Array<std::complex<float> >::const_contiter inPtr = data.cbegin() + startChannel * polCount;
	casacore::Array<float>::const_contiter weightPtr = weights.cbegin() + startChannel * polCount;
	casacore::Array<bool>::const_contiter flagPtr = flags.cbegin() + startChannel * polCount;
	const size_t selectedChannelCount = endChannel - startChannel;
		
	size_t polIndex;
	if(polOut == Polarization::Instrumental)
	{
		for(size_t ch=0; ch!=selectedChannelCount * polsIn.size(); ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
				dest[ch] = *weightPtr;
			else
				dest[ch] = 0;
			inPtr++;
			weightPtr++;
			flagPtr++;
		}
	}
	else if(Polarization::TypeToIndex(polOut, polsIn, polIndex)) {
		inPtr += polIndex;
		weightPtr += polIndex;
		flagPtr += polIndex;
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
				dest[ch] = *weightPtr;
			else
				dest[ch] = 0;
			inPtr += polCount;
			weightPtr += polCount;
			flagPtr += polCount;
		}
	}
	else {
		size_t polIndexA=0, polIndexB=0;
		switch(polOut) {
			case Polarization::StokesI: {
				bool hasXY = Polarization::TypeToIndex(Polarization::XX, polsIn, polIndexA);
				bool hasYX = Polarization::TypeToIndex(Polarization::YY, polsIn, polIndexB);
				if(!hasXY || !hasYX) {
					Polarization::TypeToIndex(Polarization::RR, polsIn, polIndexA);
					Polarization::TypeToIndex(Polarization::LL, polsIn, polIndexB);
				}
			}
			break;
			case Polarization::StokesQ: {
				bool hasXX = Polarization::TypeToIndex(Polarization::XX, polsIn, polIndexA);
				bool hasYY = Polarization::TypeToIndex(Polarization::YY, polsIn, polIndexB);
				if(!hasXX || !hasYY) {
					Polarization::TypeToIndex(Polarization::RL, polsIn, polIndexA);
					Polarization::TypeToIndex(Polarization::LR, polsIn, polIndexB);
				}
			}
			break;
			case Polarization::StokesU: {
				bool hasXY = Polarization::TypeToIndex(Polarization::XY, polsIn, polIndexA);
				bool hasYX = Polarization::TypeToIndex(Polarization::YX, polsIn, polIndexB);
				if(!hasXY || !hasYX) {
					Polarization::TypeToIndex(Polarization::RL, polsIn, polIndexA);
					Polarization::TypeToIndex(Polarization::LR, polsIn, polIndexB);
				}
			}
			break;
			case Polarization::StokesV: {
				bool hasXY = Polarization::TypeToIndex(Polarization::XY, polsIn, polIndexA);
				bool hasYX = Polarization::TypeToIndex(Polarization::YX, polsIn, polIndexB);
				if(!hasXY || !hasYX) {
					Polarization::TypeToIndex(Polarization::RR, polsIn, polIndexA);
					Polarization::TypeToIndex(Polarization::LL, polsIn, polIndexB);
				}
			}
			break;
			default:
				throw std::runtime_error("Could not convert ms polarizations to requested polarization");
			break;
		}
		
		weightPtr += polIndexA;
		inPtr += polIndexA;
		flagPtr += polIndexA;
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
				dest[ch] = *weightPtr;
			else
				dest[ch] = 0;
			inPtr += polIndexB-polIndexA;
			weightPtr += polIndexB-polIndexA;
			flagPtr += polIndexB-polIndexA;
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
				dest[ch] += *weightPtr;
			else
				dest[ch] = 0.0;
			weightPtr += polCount - polIndexB + polIndexA;
			inPtr += polCount - polIndexB + polIndexA;
			flagPtr += polCount - polIndexB + polIndexA;
		}
	}
}

template
void MSProvider::copyWeights<float>(float* dest, size_t startChannel, size_t endChannel, const std::vector<PolarizationEnum>& polsIn, const casacore::Array<std::complex<float>>& data, const casacore::Array<float>& weights, const casacore::Array<bool>& flags, PolarizationEnum polOut);

template
void MSProvider::copyWeights<std::complex<float>>(std::complex<float>* dest, size_t startChannel, size_t endChannel, const std::vector<PolarizationEnum>& polsIn, const casacore::Array<std::complex<float>>& data, const casacore::Array<float>& weights, const casacore::Array<bool>& flags, PolarizationEnum polOut);

void MSProvider::reverseCopyData(casacore::Array<std::complex<float>>& dest, size_t startChannel, size_t endChannel, const std::vector<PolarizationEnum> &polsDest, const std::complex<float>* source, PolarizationEnum polSource)
{
	size_t polCount = polsDest.size();
	const size_t selectedChannelCount = endChannel - startChannel;
	casacore::Array<std::complex<float>>::contiter dataIter = dest.cbegin() + startChannel * polCount;
	
	size_t polIndex;
	if(polSource == Polarization::Instrumental)
	{
		for(size_t chp=0; chp!=selectedChannelCount * polsDest.size(); ++chp)
		{
			if(std::isfinite(source[chp].real()))
			{
				*dataIter = source[chp];
			}
			dataIter++;
		}
	}
	else if(Polarization::TypeToIndex(polSource, polsDest, polIndex)) {
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(std::isfinite(source[ch].real()))
			{
				*(dataIter+polIndex) = source[ch];
			}
			dataIter += polCount;
		}
	}
	else {
		switch(polSource) {
			case Polarization::StokesI: {
				size_t polIndexA=0, polIndexB=0;
				bool hasXX = Polarization::TypeToIndex(Polarization::XX, polsDest, polIndexA);
				bool hasYY = Polarization::TypeToIndex(Polarization::YY, polsDest, polIndexB);
				if(!hasXX || !hasYY) {
					Polarization::TypeToIndex(Polarization::RR, polsDest, polIndexA);
					Polarization::TypeToIndex(Polarization::LL, polsDest, polIndexB);
				}
				for(size_t ch=0; ch!=selectedChannelCount; ++ch)
				{
					if(std::isfinite(source[ch].real()))
					{
						*(dataIter + polIndexA) = source[ch]; // XX = I (or rr = I)
						*(dataIter + polIndexB) = source[ch]; // YY = I (or ll = I)
					}
					dataIter += polCount;
				}
			}
			break;
			case Polarization::StokesQ: {
				size_t polIndexA=0, polIndexB=0;
				bool hasXX = Polarization::TypeToIndex(Polarization::XX, polsDest, polIndexA);
				bool hasYY = Polarization::TypeToIndex(Polarization::YY, polsDest, polIndexB);
				if(hasXX && hasYY) {
					// StokesQ to linear
					for(size_t ch=0; ch!=selectedChannelCount; ++ch)
					{
						if(std::isfinite(source[ch].real()))
						{
							casacore::Complex stokesI = casacore::Complex::value_type(0.5) * (*(dataIter + polIndexB) + *(dataIter + polIndexA));
							*(dataIter + polIndexA) = stokesI + source[ch]; // XX = I + Q
							*(dataIter + polIndexB) = stokesI - source[ch]; // YY = I - Q
						}
						dataIter += polCount;
					}
				}
				else {
					// StokesQ to circular
					Polarization::TypeToIndex(Polarization::RL, polsDest, polIndexA);
					Polarization::TypeToIndex(Polarization::LR, polsDest, polIndexB);
					for(size_t ch=0; ch!=selectedChannelCount; ++ch)
					{
						if(std::isfinite(source[ch].real()))
						{
							*(dataIter + polIndexA) = source[ch]; // rl = Q + iU (with U still zero)
							*(dataIter + polIndexB) = source[ch]; // lr = Q - iU (with U still zero)
						}
						dataIter += polCount;
					}
				}
			}
			break;
			case Polarization::StokesU: {
				size_t polIndexA=0, polIndexB=0;
				bool hasXY = Polarization::TypeToIndex(Polarization::XY, polsDest, polIndexA);
				bool hasYX = Polarization::TypeToIndex(Polarization::YX, polsDest, polIndexB);
				if(hasXY && hasYX) {
					// StokesU to linear
					for(size_t ch=0; ch!=selectedChannelCount; ++ch)
					{
						if(std::isfinite(source[ch].real()))
						{
							*(dataIter + polIndexA) = source[ch]; // XY = (U + iV), V still zero
							*(dataIter + polIndexB) = source[ch]; // YX = (U - iV), V still zero
						}
						dataIter += polCount;
					}
				}
				else {
					// StokesU to circular
					Polarization::TypeToIndex(Polarization::RL, polsDest, polIndexA);
					Polarization::TypeToIndex(Polarization::LR, polsDest, polIndexB);
					for(size_t ch=0; ch!=selectedChannelCount; ++ch)
					{
						if(std::isfinite(source[ch].real()))
						{
							// Q = (RL + LR) / 2
							casacore::Complex stokesQ = casacore::Complex::value_type(0.5) * (*(dataIter + polIndexA) + *(dataIter + polIndexB));
							casacore::Complex iTimesStokesU = casacore::Complex(-source[ch].imag(), source[ch].real());
							*(dataIter + polIndexA) = stokesQ + iTimesStokesU; // rl = Q + iU
							*(dataIter + polIndexB) = stokesQ - iTimesStokesU; // lr = Q - iU
						}
						dataIter += polCount;
					}
				}
			}
			break;
			case Polarization::StokesV: {
				size_t polIndexA=0, polIndexB=0;
				bool hasXY = Polarization::TypeToIndex(Polarization::XY, polsDest, polIndexA);
				bool hasYX = Polarization::TypeToIndex(Polarization::YX, polsDest, polIndexB);
				if(hasXY && hasYX) {
					// StokesV to linear
					for(size_t ch=0; ch!=selectedChannelCount; ++ch)
					{
						if(std::isfinite(source[ch].real()))
						{
							// U = (YX + XY)/2
							casacore::Complex stokesU = casacore::Complex::value_type(0.5) * (*(dataIter + polIndexB) + *(dataIter + polIndexA));
							casacore::Complex iTimesStokesV = casacore::Complex(-source[ch].imag(), source[ch].real());
							*(dataIter + polIndexA) = stokesU + iTimesStokesV; // XY = (U + iV)
							*(dataIter + polIndexB) = stokesU - iTimesStokesV; // YX = (U - iV)
						}
						dataIter += polCount;
					}
				}
				else {
					// StokesV to circular
					Polarization::TypeToIndex(Polarization::RR, polsDest, polIndexA);
					Polarization::TypeToIndex(Polarization::LL, polsDest, polIndexB);
					for(size_t ch=0; ch!=selectedChannelCount; ++ch)
					{
						if(std::isfinite(source[ch].real()))
						{
							// I = (RR + LL)/2
							casacore::Complex stokesI = casacore::Complex::value_type(0.5) * (*(dataIter + polIndexA) + *(dataIter + polIndexB));
							*(dataIter + polIndexA) = stokesI + source[ch]; // RR = I + V
							*(dataIter + polIndexB) = stokesI - source[ch]; // LL = I - V
						}
						dataIter += polCount;
					}
				}
			}
			break;
			default:
				throw std::runtime_error("Can't store polarization in set (not implemented or conversion not possible)");
		}
	}
}

void MSProvider::getRowRange(casacore::MeasurementSet& ms, const MSSelection& selection, size_t& startRow, size_t& endRow)
{
	startRow = 0;
	endRow = ms.nrow();
	if(selection.HasInterval())
	{
		Logger::Info << "Determining first and last row index... ";
		Logger::Info.Flush();
		casacore::ROScalarColumn<double> timeColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::TIME));
		double time = timeColumn(0);
		size_t timestepIndex = 0;
		for(size_t row = 0; row!=ms.nrow(); ++row)
		{
			if(time != timeColumn(row))
			{
				++timestepIndex;
				if(timestepIndex == selection.IntervalStart())
					startRow = row;
				if(timestepIndex == selection.IntervalEnd())
				{
					endRow = row;
					break;
				}
				time = timeColumn(row);
			}
		}
		Logger::Info << "DONE (" << startRow << '-' << endRow << ")\n";
	}
}

void MSProvider::getRowRangeAndIDMap(casacore::MeasurementSet& ms, const MSSelection& selection, size_t& startRow, size_t& endRow, const std::set<size_t>& dataDescIds, vector<size_t>& idToMSRow)
{
	startRow = 0;
	endRow = ms.nrow();
	
	Logger::Info << "Mapping measurement set rows... ";
	Logger::Info.Flush();
	casacore::ROArrayColumn<double> uvwColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::UVW));
	casacore::ROScalarColumn<int> antenna1Column(ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA1));
	casacore::ROScalarColumn<int> antenna2Column(ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA2));
	casacore::ROScalarColumn<int> fieldIdColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::FIELD_ID));
	casacore::ROScalarColumn<double> timeColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::TIME));
	casacore::ROScalarColumn<int> dataDescIdColumn(ms, ms.columnName(casacore::MSMainEnums::DATA_DESC_ID));
	double time = timeColumn(0);
	size_t timestepIndex = 0;
	bool timeStepSelected = !selection.HasInterval() || timestepIndex == selection.IntervalStart();
	for(size_t row = 0; row!=ms.nrow(); ++row)
	{
		if(time != timeColumn(row))
		{
			++timestepIndex;
			if(selection.HasInterval() && timestepIndex == selection.IntervalStart())
			{
				startRow = row;
				timeStepSelected = true;
			}
			if(timestepIndex == selection.IntervalEnd())
			{
				if(selection.HasInterval())
					endRow = row;
				break;
			}
			time = timeColumn(row);
		}
		if(timeStepSelected)
		{
			const int
				a1 = antenna1Column(row), a2 = antenna2Column(row),
				fieldId = fieldIdColumn(row), dataDescId = dataDescIdColumn(row);
			casacore::Vector<double> uvw = uvwColumn(row);
			std::set<size_t>::const_iterator dataDescIdIter = dataDescIds.find(dataDescId);
			if(selection.IsSelected(fieldId, timestepIndex, a1, a2, uvw) && dataDescIdIter != dataDescIds.end())
				idToMSRow.push_back(row);
		}
	}
	Logger::Info << "DONE (" << startRow << '-' << endRow << "; " << idToMSRow.size() << " rows)\n";
}

void MSProvider::initializeModelColumn(casacore::MeasurementSet& ms)
{
	casacore::ROArrayColumn<casacore::Complex> dataColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::DATA));
	if(ms.isColumn(casacore::MSMainEnums::MODEL_DATA))
	{
		casacore::ArrayColumn<casacore::Complex> modelColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::MODEL_DATA));
		casacore::IPosition dataShape = dataColumn.shape(0);
		bool isDefined = modelColumn.isDefined(0);
		bool isSameShape = false;
		if(isDefined)
		{
			casacore::IPosition modelShape = modelColumn.shape(0);
			isSameShape = modelShape == dataShape;
		}
		if(!isDefined || !isSameShape)
		{
			Logger::Warn << "WARNING: Your model column does not have the same shape as your data column: resetting MODEL column.\n";
			casacore::Array<casacore::Complex> zeroArray(dataShape);
			for(casacore::Array<casacore::Complex>::contiter i=zeroArray.cbegin(); i!=zeroArray.cend(); ++i)
				*i = std::complex<float>(0.0, 0.0);
			for(size_t row=0; row!=ms.nrow(); ++row)
				modelColumn.put(row, zeroArray);
		}
	}
	else { //if(!_ms.isColumn(casacore::MSMainEnums::MODEL_DATA))
		Logger::Info << "Adding model data column... ";
		Logger::Info.Flush();
		casacore::IPosition shape = dataColumn.shape(0);
		casacore::ArrayColumnDesc<casacore::Complex> modelColumnDesc(ms.columnName(casacore::MSMainEnums::MODEL_DATA), shape);
		try {
			ms.addColumn(modelColumnDesc, "StandardStMan", true, true);
		} catch(std::exception& e)
		{
			ms.addColumn(modelColumnDesc, "StandardStMan", false, true);
		}
		
		casacore::Array<casacore::Complex> zeroArray(shape);
		for(casacore::Array<casacore::Complex>::contiter i=zeroArray.cbegin(); i!=zeroArray.cend(); ++i)
			*i = std::complex<float>(0.0, 0.0);
		
		casacore::ArrayColumn<casacore::Complex> modelColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::MODEL_DATA));
		
		for(size_t row=0; row!=ms.nrow(); ++row)
			modelColumn.put(row, zeroArray);
		
		Logger::Info << "DONE\n";
	}
}

vector<PolarizationEnum> MSProvider::GetMSPolarizations(casacore::MeasurementSet& ms)
{
	std::vector<PolarizationEnum> pols;
	casacore::MSPolarization polTable(ms.polarization());
	casacore::ROArrayColumn<int> corrTypeColumn(polTable, casacore::MSPolarization::columnName(casacore::MSPolarizationEnums::CORR_TYPE));
	casacore::Array<int> corrTypeVec(corrTypeColumn(0));
	for(casacore::Array<int>::const_contiter p=corrTypeVec.cbegin(); p!=corrTypeVec.cend(); ++p)
		pols.push_back(Polarization::AipsIndexToEnum(*p));
	return pols;
}

bool MSProvider::openWeightSpectrumColumn(casacore::MeasurementSet& ms, std::unique_ptr<casacore::ROArrayColumn<float>>& weightColumn, const casacore::IPosition& dataColumnShape)
{
	bool isWeightDefined;
	if(ms.isColumn(casacore::MSMainEnums::WEIGHT_SPECTRUM))
	{
		weightColumn.reset(new casacore::ROArrayColumn<float>(ms, casacore::MS::columnName(casacore::MSMainEnums::WEIGHT_SPECTRUM)));
		isWeightDefined = weightColumn->isDefined(0);
	} else {
		isWeightDefined = false;
	}
	casacore::Array<float> weightArray(dataColumnShape);
	if(isWeightDefined)
	{
		casacore::IPosition weightShape = weightColumn->shape(0);
		isWeightDefined = (weightShape == dataColumnShape);
	}
	if(!isWeightDefined)
	{
		Logger::Warn << "WARNING: This measurement set has no or an invalid WEIGHT_SPECTRUM column; will use less informative WEIGHT column.\n";
		weightColumn.reset();
	}
	return isWeightDefined;
}
