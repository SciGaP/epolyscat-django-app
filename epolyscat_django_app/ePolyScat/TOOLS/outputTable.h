#ifndef OUTPUTTABLE_H
#define OUTPUTTABLE_H

#include <string>
#include <vector>

/// table-formated data for human-readable print, adjust column width to fit all
/// (used in PrintOutput)
class OutputTable
{
    // _row is certainly redundand
    size_t _row,_col;
    std::vector<std::vector<std::string>> _table;
    std::vector<size_t> _columnWidths;
public:
    OutputTable(); ///< emtpy output table
    OutputTable(const std::vector<std::string> &Title);
    /// (re-)calculate column widths and return rows starting FromLine
    std::vector<std::string> format(size_t FromLine=0);

    size_t rows(){return _table.size()?_table[0].size():0;} ///< number of rows in table
    void newRow();                   ///<open a new row in table
    void addRow(const std::vector<std::string> &Row);
    void rowItem(double Value, size_t Width=4); ///< append double item to current row
    void rowItem(   int Value);  ///< append int item to current row
    void rowItem(size_t Value);///< append unsigned int item to current row
    void rowItem(std::string Value); ///< append string item to current row
};

#endif // OUTPUTTABLE_H
