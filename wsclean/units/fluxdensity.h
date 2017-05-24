#ifndef FLUX_DENSITY_H
#define FLUX_DENSITY_H

#include <string>
#include <sstream>
#include <stdexcept>
#include <cmath>

#include <boost/algorithm/string/case_conv.hpp>

class FluxDensity
{
public:
	enum Unit { Jansky, MilliJansky, MicroJansky, NanoJansky, KiloJansky };
	/**
	 * Parse the string to a flux density, possibly with unit specification, and return in Jansky.
	 * @return The flux density
	 */
	static double Parse(const std::string& s, const std::string& valueDescription, Unit defaultUnit = Jansky);
	
	static std::string ToNiceString(double angleRad);
	
private:
	static size_t findNumberEnd(const std::string& s);
	static bool isDigit(const char c) { return c>='0' && c<='9'; }
	static bool isWhitespace(const char c) { return c==' ' || c=='\t'; }
};

inline std::string FluxDensity::ToNiceString(double valueJansky)
{
	std::ostringstream str;
	if(valueJansky == 0.0)
		return "0 Jy";
	else {
		if(valueJansky < 0.0)
		{
			str << "-";
			valueJansky = -valueJansky;
		}
		if(valueJansky >= 1000.0)
			str << round(valueJansky*0.1)/100.0 << " KJy";
		else if(valueJansky >= 1.0)
			str << round(valueJansky*1e2)/100.0 << " Jy";
		else if(valueJansky >= 1e-3)
			str << round(valueJansky*1e5)/100.0 << " mJy";
		else if(valueJansky >= 1e-6)
			str << round(valueJansky*1e8)/100.0 << " µJy";
		else if(valueJansky >= 1e-9)
			str << round (valueJansky*1e11)/100.0 << " nJy";
		return str.str();
	}
}

inline double FluxDensity::Parse(const std::string& s, const std::string& valueDescription, Unit defaultUnit)
{
  size_t end = findNumberEnd(s);
	if(end == 0)
		throw std::runtime_error("Error parsing " + valueDescription);
	std::string number = s.substr(0, end);
	double val = atof(number.c_str());
	// Skip whitespace after number
	const char *c = s.c_str();
	while(isWhitespace(c[end]))
		++end;
	std::string unitStr = std::string(&c[end]);
	std::string lUnitStr = boost::to_lower_copy(unitStr);

	// Unit string empty? Than use default unit.
	if(unitStr.empty())
	{
		switch(defaultUnit)
		{
			case KiloJansky:
				return val * 1000.0;
			case Jansky:
				return val;
			case MilliJansky:
				return val / 1e3;
			case MicroJansky:
				return val / 1e6;
			case NanoJansky:
				return val / 1e9;
		}
	}
	else if(lUnitStr=="jy" || lUnitStr=="jansky")
		return val;
	else if(unitStr=="mjy" || unitStr=="mJy" || lUnitStr == "millijansky")
		return val*1e-3;
	else if(unitStr=="KJy" || unitStr=="kjy" || lUnitStr=="kilojansky")
		return val*1e3;
	else if(unitStr=="µJy" || unitStr=="µjy" || lUnitStr=="microjansky")
		return val*1e-6;
	else if(unitStr=="njy" || unitStr=="nJy" || lUnitStr=="nanojansky")
		return val*1e-9;
	
	throw std::runtime_error("Invalid unit specification in flux density given for " + valueDescription);
}

inline size_t FluxDensity::findNumberEnd(const std::string& s)
{
	const char* c = s.c_str();
	size_t pos = 0;
	while(isWhitespace(c[pos]))
		++pos;
	while(isDigit(c[pos]))
		++pos;
	if(c[pos]=='.')
		++pos;
	while(isDigit(c[pos]))
		++pos;
	if(c[pos]=='e' || c[pos]=='E')
	{
		++pos;
		if(c[pos]=='-' || c[pos]=='+')
			++pos;
		while(isDigit(c[pos]))
			++pos;
	}
	return pos;
}

#endif

