#ifndef PARSET_READER_H
#define PARSET_READER_H

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

class ParsetReader
{
public:
	ParsetReader(const std::string& filename);
	ParsetReader(std::istream& stream);
	
	class ParsetEntry
	{
	public:
		ParsetEntry(ParsetEntry&& source) = default;
		enum Type { String, StringList };
		ParsetEntry(const std::string& line);
		bool operator<(const ParsetEntry& rhs) const
		{
			return _key < rhs._key;
		}
		const std::string& Key() const { return _key; }
		const std::string& GetStringValue() const;
		const std::vector<std::string>& GetStringListValue() const;
	private:
		struct Value { virtual ~Value() { } };
		struct StringValue : public Value { std::string _value; } ;
		struct StringListValue : public Value { std::vector<std::string> _value; } ;
		
		std::string _key;
		std::unique_ptr<Value> _value;
	};
	
	const std::string& GetString(const std::string& key) const;
	const std::string& GetStringOr(const std::string& key, const std::string& orValue) const;
	const std::vector<std::string>& GetStringList(const std::string& key) const;
	bool GetBool(const std::string& key) const;
	bool GetBoolOr(const std::string& key, bool orValue) const;
	
private:
	void read(std::istream& stream);
	
	std::map<std::string, ParsetEntry> _entries;
};

#endif

