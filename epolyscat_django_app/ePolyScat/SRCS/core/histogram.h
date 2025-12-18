// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BINCOEFFS_H
#define BINCOEFFS_H

#include "coordinateTrans.h"
#include "tree.h"
#include <memory>

class Coefficients;
class Index;
class AsciiFile;
class CoorSystem;

/// \ingroup Plot

/// \brief histogram of modulus-squared of grid values wrt. to one coordinate for each coordinate sub-system
class Histogram: public Tree<Histogram>{
    struct Info{
        std::shared_ptr<CoorSystem> _from,_to;
        int _binI,_size; // which coordinate to use for binning, number of 1d bins
        double _low,_width; // beginning of axis, bin width

        Info(std:: string FromCoor, std::string ToCoor, const std::string BinCoor, double Low, double Up, int NBins);
        Info(std::string From_To_Which_Low_Up_Size);
        int size() const {return _size;}
        std::string str() const;

    };

public:
    /// given (i0,i1,...,iN), map for _subLevel 1,3,6 is tree (i1,i3,i6)
    /// leaf's contain Data, i.e. bin-number, weight, and original coordinates
    struct Map:public Tree<Map>{
        // NOTE: this Map class has a more general appeal to it...
        std::vector<unsigned int> _subLevel;
        struct Data{
            int _bin;        /// add into _bin of Histogram
            double _weight;  /// multiply value by _weight (Jacobian determinant for transformation)
            std ::vector<double> _fromCoor; /// original coordinates (for cross-checking)
            Data(int Bin,double Weight,std::vector<double> Coor)
                :_bin(Bin),_weight(Weight),_fromCoor(Coor){}
            std::string str() const;
        };
        std::shared_ptr<Data> _dat;

        Map(const Index* FromIndex, Info Hist, std::vector<double> Coor=std::vector<double>());
        const std::shared_ptr<Data> data(const std::vector<unsigned int> Idx){return nodeAt(Idx)->_dat;}
        std::string strNode(int Level) const;
    };
private:
    std::shared_ptr<Map> _map;
    double _value;
    double _center;

    int bin(std::vector<unsigned int> Sub){return _map->nodeAt(Sub)->_dat->_bin;}
    double weight(std::vector<unsigned int> Sub){return _map->nodeAt(Sub)->_dat->_weight;}
    void fill(const Coefficients *GridVals, std::vector<unsigned int> Idx);
    void add(std::vector<unsigned int> Idx, double Val);
    void write(AsciiFile & File, std::vector<double> &Row);
    bool empty(){return _map.use_count()==0;}

public:
    Histogram(const Index* FromIndex, std::vector<std::string> InfStr=std::vector<std::string>(),int Depth=0);
    void fill(Coefficients *GridVals);
    void write(std::string FileName);
    std::string strNode(int Precision=0) const;
};

#endif // BINCOEFFS_H
