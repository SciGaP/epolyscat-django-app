// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISMATMATRIX_H
#define BASISMATMATRIX_H

#include <memory>
#include <map>
#include <vector>

#include "qtEigenDense.h"
#include "eigenTools.h"

class BasisAbstract;
class ReadInput;
class UseMatrix;

///@brief matrix for pair of bases that need not be integrable
class BasisMatMatrix
{
    static std::map<std::string,std::shared_ptr<BasisMatMatrix> >_list;
protected:
    Eigen::MatrixXcd _mat;
    static std::complex<double> preFactor(std::string Op);
    BasisMatMatrix(){}
    BasisMatMatrix(std::string Op, const BasisAbstract* IBas, const BasisAbstract * JBas);
public:
    virtual ~BasisMatMatrix(){}
    static void read(ReadInput & Inp);
    static const BasisMatMatrix* factory(std::string Op, const BasisAbstract* IBas, const BasisAbstract * JBas);
    static void add(std::string Def,const Eigen::MatrixXcd & Mat);
    std::string str(int Digits=2) const {return EigenTools::str(_mat,Digits);}

    bool isEmpty() const {return _mat.size()==0;}
    const UseMatrix useMat() const;
    const Eigen::MatrixXcd & mat() const{return _mat;}
    const std::vector<const Eigen::MatrixXcd*> mats() const{return std::vector<const Eigen::MatrixXcd*>(1,&_mat);}
};

#endif // BASISMATMATRIX_H
