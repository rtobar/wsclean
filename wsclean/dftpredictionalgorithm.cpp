#include "dftpredictionalgorithm.h"
#include "units/imagecoordinates.h"
#include "matrix2x2.h"
#include "model/model.h"
#include "progressbar.h"

#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>

DFTPredictionImage::DFTPredictionImage(size_t width, size_t height, ImageBufferAllocator& allocator) :
	_width(width),
	_height(height),
	_allocator(&allocator)
{
	const size_t s = width*height;
	for(size_t p=0; p!=4; ++p)
	{
		_allocator->Allocate(s, _images[p]);
		for(size_t i=0; i!=s; ++i)
			_images[p][i] = 0.0;
	}
}

void DFTPredictionImage::Add(PolarizationEnum polarization, const double* image)
{
	size_t size = _width * _height;
	for(size_t i=0; i!=size; ++i)
	{
		std::complex<double> temp[4];
		double stokes[4];
		switch(polarization)
		{
		case Polarization::XX:
			temp[0] = image[i];
			Polarization::LinearToStokes(temp, stokes);
			_images[0][i] += stokes[0];
			_images[1][i] += stokes[1];
			break;
		case Polarization::XY:
		case Polarization::YX:
		case Polarization::Instrumental:
			throw std::runtime_error("Invalid call to Add()");
		case Polarization::YY:
			temp[3] = image[i];
			Polarization::LinearToStokes(temp, stokes);
			_images[0][i] += stokes[0];
			_images[1][i] += stokes[1];
			break;
		case Polarization::StokesI:
			_images[0][i] += image[i];
			break;
		case Polarization::StokesQ:
			_images[1][i] += image[i];
			break;
		case Polarization::StokesU:
			_images[2][i] += image[i];
			break;
		case Polarization::StokesV:
			_images[3][i] += image[i];
			break;
		case Polarization::RR:
			temp[0] = image[i];
			Polarization::CircularToStokes(temp, stokes);
			_images[0][i] += stokes[0];
			_images[3][i] += stokes[3];
			break;
		case Polarization::RL:
		case Polarization::LR:
			throw std::runtime_error("Invalid call to Add()");
		case Polarization::LL:
			temp[0] = image[i];
			Polarization::CircularToStokes(temp, stokes);
			_images[0][i] += stokes[0];
			_images[3][i] += stokes[3];
			break;
		}
	}
	_pols.push_back(polarization);
}

void DFTPredictionImage::Add(PolarizationEnum polarization, const double* real, const double* imaginary)
{
	size_t size = _width * _height;
	for(size_t i=0; i!=size; ++i)
	{
		std::complex<double> temp[4];
		double stokes[4];
		switch(polarization)
		{
		case Polarization::XX:
		case Polarization::YY:
		case Polarization::StokesI:
		case Polarization::StokesQ:
		case Polarization::StokesU:
		case Polarization::StokesV:
		case Polarization::RR:
		case Polarization::LL:
		case Polarization::Instrumental:
			throw std::runtime_error("Invalid call to Add()");
		case Polarization::XY:
			temp[1] = std::complex<double>(real[i], imaginary[i]);
			Polarization::LinearToStokes(temp, stokes);
			_images[2][i] += stokes[2];
			_images[3][i] += stokes[3];
			break;
		case Polarization::YX:
			temp[2] = std::complex<double>(real[i], imaginary[i]);
			Polarization::LinearToStokes(temp, stokes);
			_images[2][i] += stokes[2];
			_images[3][i] += stokes[3];
			break;
		case Polarization::RL:
			temp[1] = std::complex<double>(real[i], imaginary[i]);
			Polarization::CircularToStokes(temp, stokes);
			_images[1][i] += stokes[1];
			_images[2][i] += stokes[2];
			break;
		case Polarization::LR:
			temp[2] = std::complex<double>(real[i], imaginary[i]);
			Polarization::CircularToStokes(temp, stokes);
			_images[1][i] += stokes[1];
			_images[2][i] += stokes[2];
			break;
		}
	}
	_pols.push_back(polarization);
}

void DFTPredictionImage::FindComponents(DFTPredictionInput& destination, double phaseCentreRA, double phaseCentreDec, double pixelSizeX, double pixelSizeY, double dl, double dm, size_t channelCount)
{
	size_t index = 0;
	for(size_t y=0; y!=_height; ++y)
	{
		for(size_t x=0; x!=_width; ++x)
		{
			if(_images[0][index] != 0.0 || _images[1][index] != 0.0 ||
				_images[2][index] != 0.0 || _images[3][index] != 0.0)
			{
				double l, m;
				ImageCoordinates::XYToLM(x, y, pixelSizeX, pixelSizeY, _width, _height, l, m);
				l += dl; m += dm;
				double ra, dec;
				ImageCoordinates::LMToRaDec(l, m, phaseCentreRA, phaseCentreDec, ra, dec);
				double stokes[4] = { _images[0][index], _images[1][index],
					_images[2][index], _images[3][index] };
				std::complex<double> linear[4];
				Polarization::StokesToLinear(stokes, linear);
				destination.AddComponent(DFTPredictionComponent(ra, dec, l, m, linear, channelCount));
			}
			++index;
		}
	}
}

struct ComponentInfo
{
	std::vector<size_t> count;
	std::vector<MC2x2> beamValues;
};
	
void DFTPredictionInput::ConvertApparentToAbsolute(casacore::MeasurementSet& ms)
{
	std::vector<ComponentInfo> compInfos(_components.size());
	
	const BandData band(ms.spectralWindow());
	LBeamEvaluator evaluator(ms);
	casacore::MEpoch::ROScalarColumn timeColumn(ms, ms.columnName(casacore::MSMainEnums::TIME));
	size_t nrow = ms.nrow();

	for(std::vector<ComponentInfo>::iterator cInfo=compInfos.begin(); cInfo!=compInfos.end(); ++cInfo)
	{
		cInfo->beamValues.assign(band.ChannelCount(), MC2x2::Zero());
		cInfo->count.assign(band.ChannelCount(), 0);
	}
	
	ProgressBar progress("Evaluating beam");
	for(size_t row=0; row!=nrow; ++row)
	{
		casacore::MEpoch time = timeColumn(row);
		if(time.getValue().get() != evaluator.Time().getValue().get())
		{
			evaluator.SetTime(time);
			LBeamEvaluator::PrecalcPosInfo posInfo;
			for(size_t i=0; i!=_components.size(); ++i)
			{
				const DFTPredictionComponent& c = _components[i];
				ComponentInfo& cInfo = compInfos[i];
				evaluator.PrecalculatePositionInfo(posInfo, c.RA(), c.Dec());
				for(size_t ch=0; ch!=band.ChannelCount(); ++ch)
				{
					MC2x2 timeStepValue;
					evaluator.EvaluateFullArray(posInfo, band.ChannelFrequency(ch), timeStepValue);
					cInfo.beamValues[ch] += timeStepValue;
					++cInfo.count[ch];
				}
			}
		}
		progress.SetProgress(row+1,nrow);
	}
	
	for(size_t i=0; i!=_components.size(); ++i)
	{
		DFTPredictionComponent& c = _components[i];
		ComponentInfo& cInfo = compInfos[i];
		for(size_t ch=0; ch!=band.ChannelCount(); ++ch)
		{
			cInfo.beamValues[ch] /= double(cInfo.count[ch]);
			cInfo.beamValues[ch].Invert();
			if(ch==band.ChannelCount()/2)
			{
				std::cout << RaDecCoord::RAToString(c.RA()) << " " << RaDecCoord::DecToString(c.Dec()) << " :";
				for(size_t p=0; p!=4; ++p)
					std::cout << " " << cInfo.beamValues[ch][p];
				std::cout << " -> ";
			}
			MC2x2 temp;
			MC2x2::ATimesB(temp, cInfo.beamValues[ch], c.LinearFlux(ch));
			MC2x2::ATimesHermB(c.LinearFlux(ch), temp, cInfo.beamValues[ch]);
			if(ch==band.ChannelCount()/2)
				std::cout << c.LinearFlux(ch).ToString() << " (" << c.L() << "," << c.M() << ")\n";
		}
	}
	
	
}

void DFTPredictionInput::InitializeFromModel(const Model& model, long double phaseCentreRA, long double phaseCentreDec, const BandData& band)
{
	for(Model::const_iterator s=model.begin(); s!=model.end(); ++s)
	{
		for(ModelSource::const_iterator c=s->begin(); c!=s->end(); ++c)
		{
			const ModelComponent& comp = *c;
			long double l, m;
			DFTPredictionComponent& component = AddComponent();
			ImageCoordinates::RaDecToLM(comp.PosRA(), comp.PosDec(), phaseCentreRA, phaseCentreDec, l, m);
			component.SetPosition(comp.PosRA(), comp.PosDec(), l, m);
			if(comp.Type() == ModelComponent::GaussianSource)
			{
				component.SetGaussianInfo(comp.PositionAngle(), comp.MajorAxis(), comp.MinorAxis());
			}
			component.SetChannelCount(band.ChannelCount());
			for(size_t ch=0; ch!=band.ChannelCount(); ++ch)
			{
				MC2x2& flux = component.LinearFlux(ch);
				double stokes[4];
				for(size_t p=0; p!=4; ++p)
					stokes[p] = comp.SED().FluxAtFrequency(band.ChannelFrequency(ch), Polarization::IndexToStokes(p));
				Polarization::StokesToLinear(stokes, flux.Data());
			}
		}
	}
}

void DFTPredictionAlgorithm::Predict(MC2x2& dest, double u, double v, double w, size_t channelIndex, size_t a1, size_t a2)
{
	dest = MC2x2::Zero();
	MC2x2 single;
	for(DFTPredictionInput::const_iterator c=_input.begin(); c!=_input.end(); ++c)
	{
		predict(single, u, v, w, channelIndex, a1, a2, *c);
		dest += single;
	}
}

void DFTPredictionAlgorithm::predict(MC2x2& dest, double u, double v, double w, size_t channelIndex, size_t a1, size_t a2, const DFTPredictionComponent& component)
{
	double l = component.L(), m = component.M(), lmsqrt = component.LMSqrt();
	double angle = 2.0*M_PI*(u*l + v*m + w*(lmsqrt-1.0));
	double
		sinangleOverLMS = sin(angle),
		cosangleOverLMS = cos(angle);
	MC2x2 appFlux;
	if(_hasBeam)
	{
		MC2x2 temp;
		MC2x2::ATimesB(temp, component.AntennaInfo(a1).BeamValue(channelIndex), component.LinearFlux(channelIndex));
		MC2x2::ATimesHermB(appFlux, temp, component.AntennaInfo(a2).BeamValue(channelIndex));
	}
	else {
		appFlux = component.LinearFlux(channelIndex);
	}
	if(component.IsGaussian())
	{
		const double* gausTrans = component.GausTransformationMatrix();
		double uTemp = u*gausTrans[0] + v*gausTrans[1];
		v = u*gausTrans[2] + v*gausTrans[3];
		u = uTemp;
		double gaus = exp(-u*u - v*v);
		for(size_t p=0; p!=4; ++p)
		{
			std::complex<double> val = appFlux[p] * gaus;
			dest[p] = std::complex<double>(
				val.real() * cosangleOverLMS - val.imag() * sinangleOverLMS,
				val.real() * sinangleOverLMS + val.imag() * cosangleOverLMS);
		}
	}
	else {
		for(size_t p=0; p!=4; ++p)
		{
			std::complex<double> val = appFlux[p];
			dest[p] = std::complex<double>(
				val.real() * cosangleOverLMS - val.imag() * sinangleOverLMS,
				val.real() * sinangleOverLMS + val.imag() * cosangleOverLMS);
		}
	}
}

void DFTPredictionAlgorithm::UpdateBeam(LBeamEvaluator& beamEvaluator)
{
	if(!_hasBeam)
	{
		_hasBeam = true;
	}
	for(DFTPredictionInput::iterator component=_input.begin(); component!=_input.end(); ++component)
	{
		LBeamEvaluator::PrecalcPosInfo posInfo;
		beamEvaluator.PrecalculatePositionInfo(posInfo, component->RA(), component->Dec());
		
		for(size_t antenna=0; antenna!=component->AntennaCount(); ++antenna)
		{
			DFTAntennaInfo& antennaInfo = component->AntennaInfo(antenna);
			
			for(size_t channel=0; channel!=antennaInfo.ChannelCount(); ++channel)
			{
				double freq = _band.ChannelFrequency(channel);
				beamEvaluator.Evaluate(posInfo, freq, antenna, antennaInfo.BeamValue(channel));
				/*if(channel == antennaInfo.ChannelCount()/2 && antenna==0)
				{
					std::cout << antennaInfo.BeamValue(channel).ToString() << ' ';
				}*/
			}
		}
	}
}
