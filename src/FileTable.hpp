#ifndef CORAS_FILE_TABLE_H
#define CORAS_FILE_TABLE_H

#include <cstddef>
#include <string>
#include <vector>

/**
 * \addtogroup AUXILIARY
 * @{
 */

/**
 * \brief Class implementing a routine to read tabulated values from a file.
 */
class FileTable
{
public:
	/**
	 * Constructor instantiating a table that reads the values in the file specified by the path \a fileName, containing
	 * \a columnNumber columns. Columns are separated by \a delimiterCharacter (default: white space), lines starting
	 * with \a commentCharacter (default: #) are ignored, and the first \a skippedColumnNumber (default: 0) lines are
	 * skipped.
	 */
	FileTable(const std::string &fileName, std::size_t columnNumber, char delimiterCharacter = ' ', char commentCharacter = '#', std::size_t skippedColumnNumber = 0);

	/**
	 * Return a reference to the vector of strings contained in the column with index \a columnIndex (starting at index
	 * 0).
	 */
	std::vector<std::string> &column(std::size_t columnIndex);

	/**
	 * Return a reference to the \c const vector of strings contained in the column with index \a columnIndex (starting
	 * at index 0).
	 */
	const std::vector<std::string> &column(std::size_t columnIndex) const;

	/**
	 * Return the vector of values of type \a T contained in the column with index \a columnIndex (starting at index 0).
	 */
	template <typename T>
	std::vector<T> column(std::size_t columnIndex) const;

	/**
	 * Return the number of columns.
	 */
	std::size_t number_of_columns() const;

	/**
	 * Return the number of rows.
	 */
	std::size_t number_of_rows() const;

private:
	/**
	 * Read the tabulated values in the file specified by the path \a fileName, and write those to the array
	 * FileTable::TableValues. Columns are separated by \a delimiterCharacter, lines starting with \a commentCharacter
	 * are ignored, columns are separated by \a delimiterCharacter, and the first \a skippedColumnNumber lines are
	 * skipped.
	 */
	void read_values_from_file(const std::string &fileName, char delimiterCharacter, char commentCharacter, std::size_t skippedColumnNumber);

	/**
	 * Number of columns.
	 */
	std::size_t ColumnNumber;

	/**
	 * Array containing the values read from file.
	 */
	std::vector<std::vector<std::string>> TableValues;
};

/** @} */

#include "FileTable.tpp"

#endif