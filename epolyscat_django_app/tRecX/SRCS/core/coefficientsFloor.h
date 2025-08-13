// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef __COEFFICIENTS_FLOOR_
#define __COEFFICIENTS_FLOOR_

#include <complex>
#include <vector>


//forward declarations
class Coefficients;
class Index;
class UseMatrix;

class CoefficientsFloor {
    friend class Coefficients;
    friend class OperatorSingle;
    friend class SpectralDiscretization; //must be able to access neighbours to produce continuous eigenvectors without weight in complex scaled area
    friend class MatrixData;
public:

    CoefficientsFloor(){}

    /// floor view (does not own data)
    CoefficientsFloor(const Index *Idx, std::complex<double>* CData){}

    /// floor view (does not own data)
    CoefficientsFloor(const CoefficientsFloor *Floor, std::complex<double>* CData){}

    // ====== functions =====================================================
//    std::complex<double> innerProduct(const Index * Idx, CoefficientsFloor& ket, bool pseudoScalar=false); ///< sum_i conjg(c_i) * ket.c_i
//    void setToZero(const Index *Idx, std::complex<double> *Data); ///< set all coefficients =0
//    bool isZero(const Index *Idx, std::complex<double> *Data, double eps=0.); ///< true if coefficients are = 0

//    void print(const Index *Idx, std::complex<double> *Data, std::ofstream & stream) const;///< ascii write
//    void write(const Index *Idx, std::complex<double> *Data, std::ofstream & stream);///< binary write
//    void  read(const Index *Idx, std::complex<double> *Data, std::ifstream & stream);///< binary read

};

#endif
