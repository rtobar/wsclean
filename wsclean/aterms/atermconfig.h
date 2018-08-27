#ifndef ATERM_CONFIG_H
#define ATERM_CONFIG_H

#include <string>

#include "atermbase.h"
#include "fitsaterm.h"
#include "lofarbeamterm.h"

#include "../matrix2x2.h"
#include "../parsetreader.h"

#include "../wsclean/wscleansettings.h"

class ATermConfig : public ATermBase
{
public:
	ATermConfig(casacore::MeasurementSet& ms, size_t nAntenna, size_t width, size_t height, double dl, double dm, double phaseCentreDL, double phaseCentreDM,
		const WSCleanSettings& settings
	) :
		_ms(ms),
		_nAntenna(nAntenna),
		_width(width),
		_height(height),
		_dl(dl), _dm(dm),
		_phaseCentreDL(phaseCentreDL),
		_phaseCentreDM(phaseCentreDM),
		_settings(settings)
	{ }
	
	void Read(const std::string& parset)
	{
		ParsetReader reader(parset);
		std::vector<std::string> aterms = reader.GetStringList("aterms");
		if(aterms.empty())
			throw std::runtime_error("No a-term correction given in parset (aterms key is an empty list)");
		
		for(const std::string atermName : aterms)
		{
			std::string atermType = reader.GetStringOr(atermName + ".type", atermName);
			if(atermType == "tec")
			{
				std::vector<std::string> tecFiles = reader.GetStringList(atermName + ".images");
				if(tecFiles.size() != 1)
					throw std::runtime_error("A TEC aterm should consist of only one image");
				std::unique_ptr<FitsATerm> f(new FitsATerm(_nAntenna, _width, _height, _dl, _dm, _phaseCentreDL, _phaseCentreDM));
				f->OpenTECFile(tecFiles.front());
				_aterms.emplace_back(std::move(f));
			}
			else if(atermType == "beam")
			{
				// TODO this option is currently specific to lofar, but should
				// determine telescope from measurement set, etc.
				bool differential = reader.GetBoolOr("beam.differential", false);
				double updateInterval = reader.GetDoubleOr("beam.update_interval", 120.0);
				std::unique_ptr<LofarBeamTerm> beam(new LofarBeamTerm(_ms, _width, _height, _dl, _dm, _phaseCentreDL, _phaseCentreDM, updateInterval, differential));
				beam->SetSaveATerms(_settings.saveATerms);
				_aterms.emplace_back(std::move(beam));
			}
		}
		Logger::Debug << "Constructed an a-term configuration with " << _aterms.size() << " terms.\n";
		if(_aterms.empty())
		{
			throw std::runtime_error("The specified a-term configuration does not define any terms to apply");
		}
		if(_aterms.size() > 1)
		{
			_previousAterms.resize(_aterms.size());
			for(ao::uvector<std::complex<float>>& buf : _previousAterms)
				buf.resize(_width * _height * _nAntenna * 4);
		}
	}
	
	virtual bool Calculate(std::complex<float>* buffer, double time, double frequency) final override
	{
		if(_aterms.size() == 1)
			return _aterms.front()->Calculate(buffer, time, frequency);
		else {
			bool isUpdated = false;
			for(size_t i=0; i!=_aterms.size(); ++i)
			{
				bool atermUpdated = _aterms[i]->Calculate(_previousAterms[i].data(), time, frequency);
				isUpdated = isUpdated || atermUpdated;
			}
			
			if(isUpdated)
			{
				std::copy(_previousAterms[0].begin(), _previousAterms[0].end(), buffer);
				for(size_t i=1; i!=_aterms.size(); ++i)
				{
					for(size_t j=0; j!=_width*_height*_nAntenna*4; j+=4)
					{
						std::complex<float> scratch[4];
						Matrix2x2::ATimesB(scratch, &_previousAterms[i][j], &buffer[j]);
						Matrix2x2::Assign(&buffer[j], scratch);
					}
				}
			}
			
			return isUpdated;
		}
	}
private:
	casacore::MeasurementSet& _ms;
	size_t _nAntenna, _width, _height;
	double _dl, _dm, _phaseCentreDL, _phaseCentreDM;
	std::vector<std::unique_ptr<ATermBase>> _aterms;
	std::vector<ao::uvector<std::complex<float>>> _previousAterms;
	const WSCleanSettings& _settings;
};

#endif
