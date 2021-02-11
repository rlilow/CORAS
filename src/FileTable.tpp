#include <sstream>

////////////////////////////////////////////////////////////////////////////////////////////////////
// public

template <typename T>
std::vector<T> FileTable::column(const std::size_t columnIndex) const
{
	std::vector<T> parsedColum;

	for (const std::string &entry : TableValues[columnIndex])
	{
		std::stringstream parser(entry);
		T parsedEntry;

		parser >> parsedEntry;

		parsedColum.push_back(parsedEntry);
	}

	return parsedColum;
}