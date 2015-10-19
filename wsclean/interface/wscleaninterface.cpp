#define WSCLEAN_NO_MAIN
#include "../wscleanmain.cpp"
#undef WSCLEAN_NO_MAIN

#include "wscleaninterface.h"

#include "../banddata.h"

#include <string>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

struct WSCleanUserData
{
	std::string msPath;
	unsigned int width;
	unsigned int height;
	double pixelScaleX;
	double pixelScaleY;
	std::string extraParameters;
	
	std::string dataColumn;
	bool hasAImage, hasAtImage;
	size_t nACalls, nAtCalls;
};

template<typename T>
std::string str(T i)
{
	std::ostringstream s;
	s << i;
	return s.str();
}

void wsclean_main(const std::vector<std::string>& parms)
{
	char** argv = new char*[parms.size()];
	for(size_t i=0; i!=parms.size(); ++i)
	{
		std::cout << parms[i] << ' ';
		argv[i] = new char[parms[i].size()+1];
		memcpy(argv[i], parms[i].c_str(), parms[i].size()+1);
	}

	wsclean_main(parms.size(), argv);
	
	for(size_t i=0; i!=parms.size(); ++i)
		delete[] argv[i];
	delete[] argv;
}

void wsclean_initialize(
	void** userData,
	const imaging_parameters* domain_info,
	imaging_data* data_info
)
{
	WSCleanUserData* wscUserData = new WSCleanUserData();
	wscUserData->msPath = domain_info->msPath;
	wscUserData->width = domain_info->imageWidth;
	wscUserData->height = domain_info->imageHeight;
	wscUserData->pixelScaleX = domain_info->pixelScaleX;
	wscUserData->pixelScaleY = domain_info->pixelScaleY;
	wscUserData->extraParameters = domain_info->extraParameters;
	wscUserData->hasAImage = false;
	wscUserData->hasAtImage = false;
	wscUserData->nACalls = 0;
	wscUserData->nAtCalls = 0;
	(*userData) = static_cast<void*>(wscUserData);
	
	// Number of vis is nchannels x selected nrows; calculate both.
	// (Assuming Stokes I polarization for now)
	casacore::MeasurementSet ms(wscUserData->msPath);
	casacore::ROScalarColumn<int> a1Col(ms, casacore::MeasurementSet::columnName(casacore::MSMainEnums::ANTENNA1));
	casacore::ROScalarColumn<int> a2Col(ms, casacore::MeasurementSet::columnName(casacore::MSMainEnums::ANTENNA2));
	BandData bandData(ms.spectralWindow());
	size_t nChannel = bandData.ChannelCount();
	size_t selectedRows = 0;
	for(size_t row=0; row!=ms.nrow(); ++row)
	{
		if(a1Col(row) != a2Col(row))
			++selectedRows;
	}
	
	data_info->dataSize = selectedRows * nChannel;
	data_info->lhs_data_type = imaging_data::DATA_TYPE_COMPLEX_DOUBLE;
	data_info->rhs_data_type = imaging_data::DATA_TYPE_DOUBLE;
	//data_info->deinitialize_function = wsclean_deinitialize;
	//data_info->read_function = wsclean_read;
	//data_info->write_function = wsclean_write;
	//data_info->operator_A_function = wsclean_operator_A;
	//data_info->operator_At_function = wsclean_operator_At;
	
	bool hasCorrected = ms.tableDesc().isColumn("CORRECTED_DATA");
	if(hasCorrected) {
		std::cout << "First measurement set has corrected data: tasks will be applied on the corrected data column.\n";
		wscUserData->dataColumn = "CORRECTED_DATA";
	} else {
		std::cout << "No corrected data in first measurement set: tasks will be applied on the data column.\n";
		wscUserData->dataColumn = "DATA";
	}
}

void wsclean_deinitialize(void* userData)
{
	WSCleanUserData* wscUserData = static_cast<WSCleanUserData*>(userData);
	delete wscUserData;
}

void wsclean_read(void* userData, DCOMPLEX* data, double* weights)
{
	WSCleanUserData* wscUserData = static_cast<WSCleanUserData*>(userData);
	casacore::MeasurementSet ms(wscUserData->msPath);
	BandData bandData(ms.spectralWindow());
	size_t nChannels = bandData.ChannelCount();
	
	casacore::ROScalarColumn<int> a1Col(ms, casacore::MeasurementSet::columnName(casacore::MSMainEnums::ANTENNA1));
	casacore::ROScalarColumn<int> a2Col(ms, casacore::MeasurementSet::columnName(casacore::MSMainEnums::ANTENNA2));
	
	casacore::ROArrayColumn<casacore::Complex> dataCol(ms, wscUserData->dataColumn);
	casacore::ROArrayColumn<float> weightCol(ms, casacore::MeasurementSet::columnName(casacore::MSMainEnums::WEIGHT_SPECTRUM));
	casacore::ROArrayColumn<bool> flagCol(ms, casacore::MeasurementSet::columnName(casacore::MSMainEnums::FLAG));
	
	DCOMPLEX* dataPtr = data;
	double* weightPtr = weights;
	casacore::IPosition shape = dataCol.shape(0);
	size_t polarizationCount = shape[0];
	casacore::Array<casacore::Complex> dataArr(shape);
	casacore::Array<bool> flagArr(shape);
	casacore::Array<float> weightArr(shape);
	for(size_t row=0; row!=ms.nrow(); ++row)
	{
		if(a1Col(row) != a2Col(row))
		{
			dataCol.get(row, dataArr);
			flagCol.get(row, flagArr);
			weightCol.get(row, weightArr);
			
			casacore::Array<casacore::Complex>::const_contiter di = dataArr.cbegin();
			casacore::Array<bool>::const_contiter fi = flagArr.cbegin();
			casacore::Array<float>::const_contiter wi = weightArr.cbegin();
			
			for(size_t ch=0; ch!=nChannels; ++ch)
			{
				// TODO this only works for XX/YY and LL/RR pol, but not if
				// MS contains IQUV
				std::complex<double> val = 0.5*(std::complex<double>(*di) + std::complex<double>(*(di+polarizationCount-1)));
				double weight = 0.5*(double(*wi) + (*wi+polarizationCount-1));
				bool flag = *fi || *(fi+polarizationCount-1);
				if(!std::isfinite(val.real()) || !std::isfinite(val.imag()))
				{
					val = 0.0;
					weight = 0.0;
				}
				
				dataPtr[ch] = val;
				weightPtr[ch] = flag ? 0.0 : weight;
				
				di += polarizationCount;
				fi += polarizationCount;
				wi += polarizationCount;
			}
			dataPtr += nChannels;
			weightPtr += nChannels;
		}
	}
}

void wsclean_write(void* userData, const double* image)
{
	WSCleanUserData* wscUserData = static_cast<WSCleanUserData*>(userData);
	FitsWriter writer;
	writer.SetImageDimensions(wscUserData->width, wscUserData->height, wscUserData->pixelScaleX, wscUserData->pixelScaleY);
	if(wscUserData->hasAtImage)
	{
		FitsReader reader("tmp-operator-At-0-image.fits");
		writer = FitsWriter(reader);
	}
	writer.Write("wsclean-model.fits", image);
}

void getCommandLine(std::vector<std::string>& commandline, const WSCleanUserData& userData)
{
	commandline.push_back("wsclean");
	commandline.push_back("-size");
	commandline.push_back(str(userData.width));
	commandline.push_back(str(userData.height));
	commandline.push_back("-scale");
	commandline.push_back(Angle::ToNiceString(userData.pixelScaleX));
	if(!userData.extraParameters.empty())
	{
		size_t pos = 0;
		size_t nextPos = userData.extraParameters.find(' ', 0);
		while(nextPos!=std::string::npos)
		{
			commandline.push_back(userData.extraParameters.substr(pos, nextPos-pos));
			pos = nextPos+1;
			nextPos = userData.extraParameters.find(' ', pos);
		}
		commandline.push_back(userData.extraParameters.substr(pos));
	}
	if(userData.pixelScaleX != userData.pixelScaleY)
		throw std::runtime_error("pixelscaleX should be equal to pixelscaleY for WSClean");
}

// Go from image to visibilities
// dataIn :  double[] of size width*height
// dataOut : complex double[] of size nvis: nchannels x nbaselines x ntimesteps
void wsclean_operator_A(void* userData, DCOMPLEX* dataOut, const double* dataIn)
{
	WSCleanUserData* wscUserData = static_cast<WSCleanUserData*>(userData);
	std::cout << "------ wsclean_operator_A(), image: " << wscUserData->width << " x " << wscUserData->height << ", pixelscale=" << Angle::ToNiceString(wscUserData->pixelScaleX) << "," << Angle::ToNiceString(wscUserData->pixelScaleY) << '\n';
	
	// Remove non-finite values
	// TODO skipped now -- can't change input array without making a copy
	size_t nonFiniteValues = 0;
	double imageSum = 0.0;
	for(size_t i=0; i!=wscUserData->width * wscUserData->height; ++i)
	{
		if(!std::isfinite(dataIn[i]))
		{
			//dataIn[i] = 0.0;
			++nonFiniteValues;
		}
		else {
			imageSum += dataIn[i];
		}
	}
	if(nonFiniteValues != 0)
		std::cout << "Warning: input image contains " << nonFiniteValues << " non-finite values!\n";
	std::cout << "Mean value in image: " << imageSum/(wscUserData->width*wscUserData->height-nonFiniteValues) << '\n';

	std::ostringstream filenameStr;
	filenameStr << "tmp-operator-A-" << wscUserData->nACalls;
	
	// Write dataIn to a fits file
	FitsWriter writer;
	writer.SetImageDimensions(wscUserData->width, wscUserData->height, wscUserData->pixelScaleX, wscUserData->pixelScaleY);
	writer.Write(filenameStr.str() + "-model.fits", dataIn);
	wscUserData->hasAImage = true;
	
	// Run WSClean -predict (creates/fills new column MODEL_DATA)
	std::vector<std::string> commandline;
	getCommandLine(commandline, *wscUserData);
	commandline.push_back("-name");
	commandline.push_back(filenameStr.str());
	commandline.push_back("-predict");
	commandline.push_back(wscUserData->msPath);
	wsclean_main(commandline);
	
	// Read MODEL_DATA into dataOut
	casacore::MeasurementSet ms(wscUserData->msPath);
	BandData bandData(ms.spectralWindow());
	size_t nChannels = bandData.ChannelCount();
	
	casacore::ROScalarColumn<int> a1Col(ms, casacore::MeasurementSet::columnName(casacore::MSMainEnums::ANTENNA1));
	casacore::ROScalarColumn<int> a2Col(ms, casacore::MeasurementSet::columnName(casacore::MSMainEnums::ANTENNA2));
	
	casacore::ROArrayColumn<casacore::Complex> dataCol(ms,  casacore::MeasurementSet::columnName(casacore::MSMainEnums::MODEL_DATA));
	
	DCOMPLEX* dataPtr = (DCOMPLEX*) dataOut;
	casacore::IPosition shape = dataCol.shape(0);
	size_t polarizationCount = shape[0];
	casacore::Array<casacore::Complex> dataArr(shape);
	for(size_t row=0; row!=ms.nrow(); ++row)
	{
		if(a1Col(row) != a2Col(row))
		{
			dataCol.get(row, dataArr);
			casacore::Array<casacore::Complex>::contiter di = dataArr.cbegin();
			
			for(size_t ch=0; ch!=nChannels; ++ch)
			{
				dataPtr[ch] = 0.5 * (std::complex<double>(*di) + std::complex<double>(*(di+polarizationCount-1)));
				// This *might* change the weighting; but if a value is not finite, it should already
				// have received zero weight during the initial read out.
				if(!std::isfinite(dataPtr[ch].real()) || !std::isfinite(dataPtr[ch].imag()))
					dataPtr[ch] = 0.0;
				di += polarizationCount;
			}
			
			dataPtr += nChannels;
		}
	}
	++(wscUserData->nACalls);
	std::cout << "------ end of wsclean_operator_A()\n";
}

// Go from visibilities to image
void wsclean_operator_At(void* userData, double* dataOut, const DCOMPLEX* dataIn)
{
	// Write dataIn to the MODEL_DATA column
	WSCleanUserData* wscUserData = static_cast<WSCleanUserData*>(userData);
	std::cout << "------ wsclean_operator_At(), image: " << wscUserData->width << " x " << wscUserData->height << ", pixelscale=" << Angle::ToNiceString(wscUserData->pixelScaleX) << "," << Angle::ToNiceString(wscUserData->pixelScaleY) << '\n';
	casacore::MeasurementSet ms(wscUserData->msPath, casacore::Table::Update);
	BandData bandData(ms.spectralWindow());
	size_t nChannels = bandData.ChannelCount();
	
	casacore::ROScalarColumn<int> a1Col(ms, casacore::MeasurementSet::columnName(casacore::MSMainEnums::ANTENNA1));
	casacore::ROScalarColumn<int> a2Col(ms, casacore::MeasurementSet::columnName(casacore::MSMainEnums::ANTENNA2));
	
	casacore::ArrayColumn<casacore::Complex> dataCol(ms,  casacore::MeasurementSet::columnName(casacore::MSMainEnums::MODEL_DATA));
	
	DCOMPLEX* dataPtr = (DCOMPLEX*) dataIn;
	casacore::IPosition shape = dataCol.shape(0);
	size_t polarizationCount = shape[0];
	casacore::Array<casacore::Complex> dataArr(shape);
	for(size_t row=0; row!=ms.nrow(); ++row)
	{
		if(a1Col(row) != a2Col(row))
		{
			dataCol.get(row, dataArr);
			casacore::Array<casacore::Complex>::contiter di = dataArr.cbegin();
			
			for(size_t ch=0; ch!=nChannels; ++ch)
			{
				*di = std::complex<float>(dataPtr[ch]);
				*(di+polarizationCount-1) = std::complex<float>(dataPtr[ch]);
				di += polarizationCount;
			}
			
			dataCol.put(row, dataArr);
			dataPtr += nChannels;
		}
	}
	
	std::ostringstream prefixName;
	prefixName << "tmp-operator-At-" << wscUserData->nAtCalls;
	
	// Run WSClean to create dirty image
	std::vector<std::string> commandline;
	getCommandLine(commandline, *wscUserData);
	commandline.push_back("-name");
	commandline.push_back(prefixName.str());
	commandline.push_back("-datacolumn");
	commandline.push_back("MODEL_DATA");
	commandline.push_back(wscUserData->msPath);
	wsclean_main(commandline);
	wscUserData->hasAtImage = true;
	
	// Read dirty image and store in dataOut
	FitsReader reader(prefixName.str() + "-image.fits");
	reader.Read(dataOut);
	++(wscUserData->nAtCalls);
	std::cout << "------ end of wsclean_operator_At()\n";
}

double wsclean_parse_angle(const char* angle)
{
	return Angle::Parse(angle, "angle", Angle::Degrees);
}
