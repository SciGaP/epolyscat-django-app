// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORMAP_H
#define OPERATORMAP_H


#define EIGEN_MATRIX_PLUGIN "EigenAddonMatrix.h"
#include "qtEigenDense.h"
#include "operatorAbstract.h"
#include "tree.h"
#include <memory>

class Index;
class BasisAbstract;
class CoefficientsViewDeep;

/** \ingroup Structures */
/// map between hierarchies of equal structure
class OperatorMap : public OperatorAbstract, public Tree<OperatorMap>
{
    // for creating contigous and deep views, as needed
    CoefficientsViewDeep *_deepX,*_deepY;
    Coefficients *_highX,*_highY,*_viewX,*_viewY;
    Index *_iX,*_iY;

    Coefficients * _tmp;
    const Eigen::MatrixXcd * _tensor;
    int _iVec,_jVec;

    class Storage{
        std::vector<std::unique_ptr<Eigen::MatrixXcd>> _matrices;
    public:
        static Storage main;
        const Eigen::MatrixXcd* get(Eigen::MatrixXcd&& from);
    };


    void axpy3(std::complex<double> A, const Coefficients *X, Coefficients *Y) const;
    OperatorMap(const Index* IIndex, const Index* JIndex, std::string Derivative, const std::complex<double> Multiplier, bool EntryLevel);
    bool isIdentity(double Eps=0.,bool Stochastic=true) const override;
public:
    using Tree::str; // there is also OperatorAbstract::str
    ~OperatorMap();
    OperatorMap(const Index* IIndex, const Index* JIndex, std::string Derivative="")
        :OperatorMap(IIndex,JIndex,Derivative,1.,true){}
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const override;
    void update(double Time, const Coefficients* C=nullptr) override {if(Time==0 and C)return;}//shut up warning about unused
    static Eigen::MatrixXcd basisMap(const BasisAbstract *IBas, const BasisAbstract *JBas, int JDerivative=0);

    std::string strNode(int Level) const override;
};

#endif // OPERATORMAP_H
