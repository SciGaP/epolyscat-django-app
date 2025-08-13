// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISDVR_H
#define BASISDVR_H

#include <vector>
#include <complex>
#include <memory>

#include "basisIntegrable.h"

class BasisSetDef;
class UseMatrix;
class OrthogonalPolynomial;
/** \ingroup Basissets */

///@brief general DVR basis (see tsurff.pdf for details of the definition)
class BasisDVR: public BasisIntegrable
{
protected:
    // should go to BasisIntegrable
    static double snapBoundary(double X, double Low, double Up){
        double x(X);
        double eps=std::min(std::abs(Up*1.e-14-Low*1.e-14),1.e-10);
        if(std::abs(x-Low)<eps)x=Low;
        if(std::abs(x- Up)<eps)x=Up;
        return x;
    }

    // the next few lines are not DVR specific and could go into a class "BasisOrthoPol
    std::shared_ptr<OrthogonalPolynomial> _opol;
    double _scale,_shiftQ,_shiftX;
    inline double qFromX(const double X) const {return snapBoundary((X-_shiftX)/_scale+_shiftQ,_opol->lowerBoundary(),_opol->upperBoundary());}
    inline double xFromQ(const double Q) const {return snapBoundary((Q-_shiftQ)*_scale+_shiftX,lowBound(),upBound());}

    std::vector<double> _dvrX,_dvrW;
    int _nBeg; ///< first Lagrange polylonmial to use
    int _size; ///< size of basis (after possible removal of boundary functions)

    bool _isDVR; ///< treat as DVR in standard matrix calculations
public:
    static bool femDVR;
    void setDVR(bool Dvr){_isDVR=Dvr;} ///< set the DVR integration status
    unsigned int size() const {return _size;}
    bool isDVR() const;
    unsigned int order() const {return _dvrX.size();}
    void quadRule(int N, std::vector<double> & QuadX, std::vector<double> & QuadW) const;

    //    static const BasisAbstract* factory(const BasisSetDef & Def);
    BasisDVR():BasisIntegrable(0.,0.),_isDVR(true){}
    BasisDVR(const BasisSetDef & Def);
    BasisDVR(const std::string & Def);

    void dvrRule(std::vector<double> & QuadX, std::vector<double> & QuadW) const {QuadX=_dvrX;QuadW=_dvrW;}
    void valDer(const std::vector<std::complex<double> >&X,std::vector<std::complex<double> >&Val,std::vector<std::complex<double> >&Der, bool ZeroOutside) const;
    virtual std::vector<double> valNodes() const; ///< values of basis functions at DVR nodes
    const std::vector<double> nodes() const{ return std::vector<double>(_dvrX.data()+_nBeg,_dvrX.data()+_nBeg+size());}   ///< nodes where basis function = 1
    const std::vector<double> weights() const{ return std::vector<double>(_dvrW.data()+_nBeg,_dvrW.data()+_nBeg+size());} ///< weights where basis function = 1
    void dvrRule(UseMatrix & QuadX, UseMatrix & QuadW) const; ///< old style interface (obsolescent)

    std::string str(int Level) const;
    std::string strDefinition() const;
    bool operator==(const BasisAbstract& Other) const;
    int nBeg() const {return _nBeg;}
};

#endif // BASISDVR_H
