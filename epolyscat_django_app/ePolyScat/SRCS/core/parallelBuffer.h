// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef PARALLELBUFFER_H
#define PARALLELBUFFER_H

#include "mpiWrapper.h"

class Coefficients;
class ParallelBuffer:public MPIwrapper::Buffer
{
    void add(const std::complex<double>* Data,int Size){
        val.insert(val.end(),Data,Data+Size);
    }
    int extract(int Pos,std::complex<double>* Data, int Size){
        for(int k=0;k<Size;k++)Data[k]=val[Pos+k];
        return Pos+Size;
    }

public:
    ParallelBuffer(){}
    void resize(int Size){val.resize(Size);}

    void add(const std::vector<std::complex<double> > & Data){add(Data.data(),Data.size());}
    int extract(int Pos, std::vector<std::complex<double> >&Data){return extract(Pos,Data.data(),Data.size());}

    void add(const Coefficients* C);
    int extract(int Pos, Coefficients* C);

    void bcast(int From){MPIwrapper::Bcast(val.data(),val.size(),From);}
    void allGather();
};

#endif // PARALLELBUFFER_H
