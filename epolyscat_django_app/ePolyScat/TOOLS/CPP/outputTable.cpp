#include "../outputTable.h"

#include "abort.h"
#include "tools.h"
#include "str.h"

#include <sstream>
#include "printOutput.h"

OutputTable::OutputTable():_row(0),_col(0){}

OutputTable::OutputTable(const std::vector<std::string> &Title){
    for(auto t: Title){
        _table.push_back({});
        _table.back().push_back(t);
    }
    _row=1;
}
//void newRow();                   ///<open a new row in table

//void rowItem(double Value, size_t Width=4); ///< append double item to current row
//void rowItem(       int Value);  ///< append int item to current row
//void rowItem(size_t Value);///< append size_t item to current row
//void rowItem(std::string Value); ///< append string item to current row
//void flush(); ///< write current output, without terminating
//void updateTable(); ///< append new outputs to table
//void end(); ///< flush and explicitly terminate any set of outputs (line,table...)
void OutputTable::newRow(){_row++;_col=1;}
void OutputTable::addRow(const std::vector<std::string> &Row){
    newRow();
    for(auto i: Row)rowItem(i);
}
void OutputTable::rowItem(double Value,size_t Width){rowItem(tools::str(Value,Width,DBL_MAX/2));}
void OutputTable::rowItem(   int Value){rowItem(tools::str(Value));}
void OutputTable::rowItem(size_t Value){rowItem(tools::str(Value));}
void OutputTable::rowItem(std::string Value){
    if(_row==0 or _col==0)ABORT("open newRow befor inserting rowItem");
    for (size_t j=_table.size();j<_col;j++)_table.push_back(std::vector<std::string>(0));
    for (size_t i=_table[_col-1].size();i<_row;i++)_table[_col-1].push_back("");
    _table[_col-1][_row-1]=Value;
    _col++;
}
static std::string PrintOutput_columnSeparator=" ";
static std::string flushRight(size_t Width, std::string Item){
    if(Width<Item.length())return Item.substr(Item.length()-Width);
    return std::string(Width-Item.length(),' ')+Item;
}

std::vector<std::string> OutputTable::format(size_t FromLine){
    std::vector<std::string> rows;
    // values in table
    if(_columnWidths.size()>0 or _table.size()>0){
        bool newTable=_columnWidths.size()==0;
        size_t length=0;
        for (size_t n=0;n<_table.size();n++){
            if(newTable){
                _columnWidths.push_back(_table[n][0].length());
                for (size_t k=1;k<_table[n].size();k++)
                    _columnWidths.back()=std::max(_columnWidths.back(),(size_t)(_table[n][k].length()));
            }
            length=std::max(length,(size_t)(_table[n].size()));
        }


        for(size_t l=FromLine;l<length;l++){
            rows.push_back("");
            for(size_t n=0;n<_table.size();n++){
                if(n==0)rows.back()+=PrintOutput::indent();
                if(_table[n].size()>l)rows.back()+=flushRight(_columnWidths[n],_table[n][l]);
                else rows.back()+=std::string(_columnWidths[n],' ');
                if(n<_table.size()-1)rows.back()+=PrintOutput_columnSeparator;
            }

            if(newTable and l==0){
                // first line in new table is header
                rows.push_back(PrintOutput::indent());
                for(size_t n=0;n<_table.size();n++){
                    if(n>0)rows.back()+=" ";
                    rows.back()+="  ";
                    rows.back()+=std::string(_columnWidths[n]>3?_columnWidths[n]-2:2,'-');
                }
            }
        }
        _columnWidths.clear();
        _col=0;
    }

    return rows;
}
