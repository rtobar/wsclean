#include "parsetreader.h"

#include <fstream>

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <boost/tokenizer.hpp>

ParsetReader::ParsetReader(const std::string& filename)
{
	std::ifstream stream(filename);
	if(!stream)
		throw std::runtime_error("Could not open parset file '" + filename + "'");
	read(stream);
}

ParsetReader::ParsetReader(std::istream& stream)
{
	read(stream);
}

void ParsetReader::read(std::istream& stream)
{
	while(stream)
	{
		std::string line;
		std::getline(stream, line);
		if(stream)
		{
			size_t hash = line.find('#');
			if(hash != line.npos)
				line = line.substr(0, hash);
			if(!line.empty())
			{
				ParsetEntry entry(line);
				_entries.emplace(entry.Key(), std::move(entry));
			}
		}
	}
}

ParsetReader::ParsetEntry::ParsetEntry(const std::string& line)
{
	size_t eqsym = line.find('=');
	if(eqsym == line.npos)
		throw std::runtime_error("Parset error: expecting equals sign ('=') in line:\n" + line);
	_key = line.substr(0, eqsym);
	boost::algorithm::trim(_key);
	if(_key.empty())
		throw std::runtime_error("Parset error: no key given before equals sign ('=')");
	std::string valStr = line.substr(eqsym+1);
	boost::algorithm::trim(valStr);
	if(valStr.empty() || valStr[0] != '[')
	{
		std::unique_ptr<StringValue> value(new StringValue());
		value->_value = valStr;
		_value = std::move(value);
	}
	else {
		// It's a list
		boost::char_separator<char> listSep("[,] ");
		boost::tokenizer<boost::char_separator<char>> listTokenizer(valStr, listSep);
		std::unique_ptr<StringListValue> value(new StringListValue());
		for(auto item : listTokenizer)
			value->_value.push_back(item);
		_value = std::move(value);
	}
}

const std::string& ParsetReader::ParsetEntry::GetStringValue() const {
	StringValue* value = dynamic_cast<StringValue*>(_value.get());
	if(value == nullptr)
		throw std::runtime_error("Value is not a string");
	else
		return value->_value;
}

const std::vector<std::string>& ParsetReader::ParsetEntry::GetStringListValue() const {
	StringListValue* value = dynamic_cast<StringListValue*>(_value.get());
	if(value == nullptr)
		throw std::runtime_error("Value is not a string list");
	else
		return value->_value;
}

const std::string& ParsetReader::GetString(const std::string& key) const
{
	auto iter = _entries.find(key);
	if(iter == _entries.end())
		throw std::runtime_error("Key not found");
	return iter->second.GetStringValue();
}

const std::string& ParsetReader::GetStringOr(const std::string& key, const std::string& orValue) const
{
	auto iter = _entries.find(key);
	if(iter == _entries.end())
		return orValue;
	else
		return iter->second.GetStringValue();
}

const std::vector<std::string>& ParsetReader::GetStringList(const std::string& key) const
{
	auto iter = _entries.find(key);
	if(iter == _entries.end())
		throw std::runtime_error("Key not found");
	return iter->second.GetStringListValue();
}

bool ParsetReader::GetBool(const std::string& key) const
{
	std::string v(GetString(key));
	boost::to_lower(v);
	if(v == "false")
		return false;
	else if(v == "true")
		return true;
	else
		throw std::runtime_error("Bad boolean value for key " + key + ": '" + v + "'");
}

bool ParsetReader::GetBoolOr(const std::string& key, bool orValue) const
{
	auto iter = _entries.find(key);
	if(iter == _entries.end())
		return orValue;
	else {
		std::string v(iter->second.GetStringValue());
		boost::to_lower(v);
		if(v == "false")
			return false;
		else if(v == "true")
			return true;
		else
			throw std::runtime_error("Bad boolean value for key " + key + ": '" + v + "'");
	}
}

double ParsetReader::GetDouble(const std::string& key) const
{
	return atof(GetString(key).c_str());
}

double ParsetReader::GetDoubleOr(const std::string& key, double orValue) const
{
	auto iter = _entries.find(key);
	if(iter == _entries.end())
		return orValue;
	else
		return atof(iter->second.GetStringValue().c_str());
}
