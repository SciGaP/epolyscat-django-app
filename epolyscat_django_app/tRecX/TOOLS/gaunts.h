// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef GAUNTS_H
#define GAUNTS_H

#include <stdlib.h>
#include <map>
#include <vector>

///@brief efficiently compute/retrieve and compactly store Gaunt factors as needed
class Gaunts
{
    std::map<int,std::map<int,std::vector<std::vector<std::vector<double>>> > > _store;
public:
    Gaunts(){}
    struct M{
        const std::vector<std::vector<std::vector<double>>> * ll;
        bool swapped;
        const std::vector<double> * retrieve(int L1, int L2){
            if(swapped)return &(*ll)[L2][L1];
            else       return &(*ll)[L1][L2];
        }
        M(int M1,int M2){
            swapped=std::abs(M1)<std::abs(M2);
            if(swapped)std::swap(M1,M2);
            if(M2<0){
                M1=-M1;
                M2=-M2;
            }
        }
    };

    ///@brief all non-zero <Y[l1,m1]Y[l2,m2]|Y[l3,m1+m2]>
    ///
    /// vals(...)[ll]=<Y[l1,m1]Y[l2,m2]|Y[Gaunts::l3Min(...)+2*ll,m1+m2]>
    ///
    /// only (l1+l2+l3)%2=0 are non-zero (parity conservation)
    /// <br>retrieves from table, extends table as needed
    /// <br>SIGNS: for Y[lm] the tRecX internal convention is used, i.e. Y[l,m]^* = Y[l,-m]
    /// <br>to connect to the convention used for GSL m-dependent sign factors must be included
    std::vector<double> & vals(int M1, int M2, int L1, int L2, int &L3Min /** lowest L3 with non-zero Gaunt */);

    ///@brief lowest non-zero L3 for given M1,M2,L1,L2
    static int l3Min(int M1, int M2, int L1, int L2){
        int lmin=std::abs(L1-L2);
        while(lmin<std::abs(M1+M2))lmin+=2;
        return lmin;
    }
};

#endif // GAUNTS_H
