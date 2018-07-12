#ifndef ATERM_CONFIG_H
#define ATERM_CONFIG_H

#include <string>

#include "atermbase.h"
#include "fitsaterm.h"
#include "lofarbeamterm.h"

#include "../parsetreader.h"

class ATermConfig : public ATermBase
{
public:
	ATermConfig(casacore::MeasurementSet& ms, size_t nAntenna, size_t width, size_t height, double dl, double dm, double phaseCentreDL, double phaseCentreDM) :
		_ms(ms),
		_nAntenna(nAntenna),
		_width(width),
		_height(height),
		_dl(dl), _dm(dm),
		_phaseCentreDL(phaseCentreDL),
		_phaseCentreDM(phaseCentreDM)
	{ }
	
	void Read(const std::string& parset)
	{
		ParsetReader reader(parset);
		std::vector<std::string> aterms = reader.GetStringList("aterms");
		if(aterms.size() != 1)
			throw std::runtime_error("At the moment, WSClean can only handle one a-term");
		
		std::string atermName = aterms[0];
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
			// TODO this whole option is currently specific to lofar...
			bool differential = reader.GetBoolOr("beam.differential", false);
			std::unique_ptr<LofarBeamTerm> beam(new LofarBeamTerm(_ms, _width, _height, _dl, _dm, _phaseCentreDL, _phaseCentreDM, differential));
			_aterms.emplace_back(std::move(beam));
		}
		Logger::Debug << "Constructed an a-term configuration with " << _aterms.size() << " terms.\n";
		if(_aterms.empty())
		{
			throw std::runtime_error("The specified a-term configuration does not define any terms to apply");
		}
	}
	
	virtual void Calculate(std::complex<float>* buffer, double time, double frequency) final override
	{
		// TODO iterate and multiple
		_aterms.front()->Calculate(buffer, time, frequency);
	}
private:
	casacore::MeasurementSet& _ms;
	size_t _nAntenna, _width, _height;
	double _dl, _dm, _phaseCentreDL, _phaseCentreDM;
	std::vector<std::unique_ptr<ATermBase>> _aterms;
};

#endif
