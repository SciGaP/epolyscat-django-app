// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "basisBesselCoulomb.h"

#ifdef _USE_GSL_
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_sf_bessel.h>
#endif

#include <cmath>
#include <map>
#include "readInput.h"
#include "printOutput.h"
#include "qtAlglib.h"
#include "qtAlglib.h"
#include "useMatrix.h"
#include "basisSetDef.h"
#include "quadratureRule.h"

static std::map<bool,std::map<bool,std::map<int,std::map<double,std::vector<Eigen::MatrixXcd> > > > > valDerMap;

BasisBesselCoulomb::BasisBesselCoulomb(double Rc, double Rx, unsigned int Langular, std::vector<double> KGrid)
    : BasisIntegrable(0.,Rx),besselRadius(Rc),lAngular(Langular),kGrid(KGrid),_untransformed(0)
{
    if(besselRadius>_upBound) DEVABORT("Cannot calculate Bessel-Coulomb function values for Rc < Rx, is: "+tools::str(besselRadius)+", "+tools::str(_upBound));

    _invNorms.assign(KGrid.size(),1.);

    // compute norms
    std::vector<std::complex<double> > ovr=QuadratureRule::integralsBasisBasis(this);
    for(int k=0;k<size();k++)_invNorms[k]=1./sqrt(ovr[k+k*size()]).real();
}

BasisBesselCoulomb::BasisBesselCoulomb(const BasisSetDef &Def)
    : BasisIntegrable(0.,Def.upBound()),lAngular(Def.par[0]),_untransformed(0)
{
    // construct kGrid from Def
    std::string par=tools::stringInBetween(Def.funcs,"[","]");
    std::vector<std::string> elem;
    if(par!=Def.funcs)elem=tools::splitString(par,',');
    if(elem.size()>2)ABORT("specify at most 2 parameters, found: "+Def.funcs);
    std::vector<double> radii;
    if(ReadInput::main.file()=="FLAG_ONLY.inp"){
        radii.push_back(tools::string_to_double(elem[0]));
    }
    else{
        ReadInput::main.read("Surface","points",radii,"","save values and derivatives at surface points (blank-separated list)");
        if(elem.size()>0 and radii[0]!=tools::string_to_double(elem[0]))
            PrintOutput::warning("Surface: points and Rc in besselCoulomb[Rc,...] must match, is: "+tools::str(radii[0])+", "+Def.funcs);
    }

    // lowest k-grid point defaults to 0;
    if(elem.size()<2)
        kGrid.assign(1,0.);
    else
        kGrid.assign(1,tools::string_to_double(elem[1]));

    while(kGrid.size()<Def.order)kGrid.push_back(kGrid.back()+math::pi/(0.5*_upBound));

    _invNorms.assign(size(),1.);

    // use Bessel functions up to old surface radius
    besselRadius=radii[0];
    if(besselRadius>_upBound)DEVABORT("Cannot calculate Bessel-Coulomb function values for Rc < Rx, is: "+tools::str(besselRadius)+", "+tools::str(_upBound));
    if(Def.shift!=0.) DEVABORT("Not yet implemented for shift!=0., is: "+tools::str(Def.shift));
    if(Def.coor.name()!="Rn")ABORT("BesselCoulomb can only be used on Coordinate Rn, is: "+Def.coor.name());

    // compute norms
    std::vector<std::complex<double> > ovr=QuadratureRule::integralsBasisBasis(this);
    for(int k=0;k<size();k++)_invNorms[k]=1./sqrt(ovr[k+k*size()]).real();
}

void BasisBesselCoulomb::initializeTransformation(){
    // evaluate Bessel-Coulomb functions at Rx and move largest value to the end
    bVec.resize(kGrid.size());
    std::vector<std::complex<double> > rXVec,bValVec,bDerVec;
    rXVec.push_back(_upBound);
    _mIdx = 0;
    double bMax = -DBL_MAX;
    for(unsigned int k=0;k<kGrid.size();k++){
        valDer(rXVec,kGrid[k],bValVec,bDerVec);
        bVec[k] = bValVec[0];
        if(bVec[k].real()>bMax){
            _mIdx = k;
            bMax = bVec[k].real();
        }
    }
    // permute the largest entry of bVec to the bottom
    for(unsigned int k=_mIdx;k<kGrid.size()-1;k++){
        std::swap(bVec[k],bVec[k+1]);
    }

    // calculate final b-vector
    for(unsigned int k=0;k<kGrid.size()-1;k++){
        bVec[k] /= (-bMax);
    }
    bVec[kGrid.size()-1] = 1./bMax;
}

// "vector" version of valDer
void BasisBesselCoulomb::valDer(const std::vector<std::complex<double> > &X, std::vector<std::complex<double> > &Val, std::vector<std::complex<double> > &Der, bool ZeroOutside) const{
    Val.clear();
    Der.clear();
    for(unsigned int j=0;j<kGrid.size();j++){
        std::vector<std::complex<double> > val,der;
        valDer(X,kGrid[j],val,der);
        for(unsigned int i=0;i<X.size();i++){
            Val.push_back(val[i]*_invNorms[j]);
            Der.push_back(der[i]*_invNorms[j]);
        }
    }
}

// the actual computations are done here
void BasisBesselCoulomb::valDer(const std::vector<std::complex<double> > &X, const double &K, std::vector<std::complex<double> > &Val, std::vector<std::complex<double> > &Der) const
{
    Val=X;
    Der.assign(X.size(),1.);
    if(K==0.)return; // include k=0

#ifndef _USE_GSL_
    ABORT("for using BasisBesselCoulomb, need to compile with option -D_USE_GSL");
#else
    double kBessel;
    if(besselRadius==_upBound) kBessel = K;
    else kBessel = sqrt(K*K+2./besselRadius);
    double eta = -1./K; // Z/k, Sommerfeld parameter
    double k = 0; // shift, not needed
    std::complex<double> norm = 1.;
    std::complex<double> a_L, b_L; // continuity coefficients
    if(not besselRadius==upBound()){
        // calculate the coefficients a_l(k) and b_l(k) to ensure continuity of value and derivative:
        gsl_sf_result f2,f2p,g2,g2p;
        double exp_f2,exp_g2;
        double kRc = K*besselRadius;
        int status0 = gsl_sf_coulomb_wave_FG_e(eta,kRc,lAngular,k,&f2,&f2p,&g2,&g2p,&exp_f2,&exp_g2);
        if(status0==17){
            PrintOutput::warning("GSL_ELOSS - loss of accuracy in gsl_sf_coulomb_wave_FG_e",5);
        }
        double f2_L = f2.val*exp(exp_f2);
        double g2_L = g2.val*exp(exp_g2);
        double f2p_L = f2p.val*exp(exp_f2);
        double g2p_L = g2p.val*exp(exp_g2);

        std::complex<double> besselVal,besselDer;
        gsl_sf_result besselL,besselLp1;
        int status1 = gsl_sf_bessel_jl_e(lAngular,kBessel*besselRadius,&besselL);
        int status2 = gsl_sf_bessel_jl_e(lAngular+1,kBessel*besselRadius,&besselLp1);
        besselVal = besselL.val;
        besselDer = lAngular*besselL.val/besselRadius-kBessel*besselLp1.val;

        b_L = (besselDer*kRc*kRc*f2_L-besselVal*K*kRc*(f2p_L*kRc-f2_L))/(K*(g2p_L*kRc-g2_L)*f2_L-g2_L*(f2p_L*kRc-f2_L)*K);
        a_L = (besselVal*kRc-b_L*g2_L)/f2_L;
        norm = sqrt(a_L*a_L + b_L*b_L);
    }

    for(int i=X.size()-1;i>=0;i--){
        if(X[i].real()<=besselRadius){
            gsl_sf_result besselL,besselLp1;
            int status1 = gsl_sf_bessel_jl_e(lAngular,kBessel*X[i].real(),&besselL);
            int status2 = gsl_sf_bessel_jl_e(lAngular+1,kBessel*X[i].real(),&besselLp1);
            Val[i] = X[i].real()*besselL.val / norm;
            Der[i] = lAngular*besselL.val-X[i].real()*kBessel*besselLp1.val / norm;
        }
        else{
            double exp_f,exp_g;
            gsl_sf_result f,fp,g,gp; // factor *r already included in Coulomb waves
            int status3 = gsl_sf_coulomb_wave_FG_e(eta,K*X[i].real(),lAngular,k,&f,&fp,&g,&gp,&exp_f,&exp_g);

            if(status3==17){
                PrintOutput::warning("GSL_ELOSS - loss of accuracy in gsl_sf_coulomb_wave_FG_e",5);
            }
            // set values below first NaN to zero
            if(std::isnan(f.val) or std::isnan(g.val) or std::isnan(fp.val) or std::isnan(gp.val)){
                std::cout << "X[i] = " << X[i] << std::endl;
                PrintOutput::warning("NaN values in gsl_sf_coulomb_wave_FG_e",5);
                break;
            }

            double f_L = f.val*exp(exp_f);
            double g_L = g.val*exp(exp_g);
            double fp_L = fp.val*exp(exp_f);
            double gp_L = gp.val*exp(exp_g);

            // final solutions: (absorb r^2)
            Val[i] = (a_L * f_L + b_L * g_L) / (K * norm);
            Der[i] = (a_L * (fp_L - f_L/(K*X[i].real()) ) + b_L * (gp_L - g_L/(K*X[i].real())) ) / norm;
        }
    }
#endif
}

bool BasisBesselCoulomb::operator==(const BasisAbstract &Other) const{
    if(dynamic_cast<const BasisBesselCoulomb*>(&Other)!=0){
        // compare kGrid and lAngular (probably more checks should be added here)
        if(this->lAngular!=dynamic_cast<const BasisBesselCoulomb*>(&Other)->lAngular) return false;
        else if(this->kGrid!=dynamic_cast<const BasisBesselCoulomb*>(&Other)->kGrid) return false;
        else return true;
    }
    else return false;
}

std::string BasisBesselCoulomb::str(int Level) const {
    std::string s=name()+" ["+tools::str(lowBound(),3)+","+tools::str(upBound(),3)+"] "+tools::str(order())+"["+tools::str(size())+"]";
    return s;
}

void BasisBesselCoulomb::quadRule(int Npoints, std::vector<double> &X, std::vector<double> &W) const{

    //HACK
    if(Npoints==size()){
        if(Npoints!=size())PrintOutput::DEVwarning(Str("points>basis size")+Npoints+size());
        QuadratureRule::pointsAndWeights(this,X,W);
        return;
    }

    DEVABORT("reimplment if needed");

//    // piecewise Gauss-Legendre quadrature
//    // divide quadrature points onto the two regions [0,Rc] and [Rc,Rx]
//    int nPointsRc,nPointsRx;
//    if(besselRadius==_upBound) nPointsRc = Npoints;
//    else{
//        nPointsRc = (int) besselRadius/upBound()*Npoints;
//        nPointsRx = Npoints - nPointsRc;
//        // make sure there is at least one quadrature point per interval
//        if(nPointsRc==0){nPointsRc++;nPointsRx--;}
//        else if(nPointsRx==0){nPointsRc--;nPointsRx++;}
//    }
//    int maxOrder = 30; // maximum number of quadrature points per interval
//    std::vector<int> quadRc,quadRx;

//    // distribute the quadrature points over the corresponding intervals
//    distributeQuadPoints(nPointsRc,maxOrder,quadRc);
//    if(besselRadius!=_upBound) distributeQuadPoints(nPointsRx,maxOrder,quadRx);

//    // calculate the quadrature rules and write them into the matrices
//    X = UseMatrix(Npoints,1);
//    W = UseMatrix(Npoints,1);
//    int pos = 0;
//    calculateQuadRules(quadRc,pos,X,W,0.,besselRadius);
//    if(besselRadius!=_upBound) calculateQuadRules(quadRx,pos,X,W,besselRadius,_upBound);

//    // equidistant quadrature points
//    //    X=UseMatrix(Npoints,1);
//    //    for(unsigned int k=0;k<Npoints;k++)X(k)=(k+1)/double(Npoints)*upBound();
//    //    W=UseMatrix::Constant(Npoints,1,upBound()/double(Npoints));
}

void BasisBesselCoulomb::distributeQuadPoints(unsigned int NPoints, int MaxOrder, std::vector<int> &QuadPoints) const {
    QuadPoints.resize(0);
    if(NPoints%MaxOrder==0){
        QuadPoints.assign(NPoints/MaxOrder,MaxOrder);
    }
    else{
        int ints, order; // number of quadrature intervals (ints) and points per interval (order)
        for(unsigned int i=1;;i++){
            if(NPoints/i>(MaxOrder-1)) continue;
            ints = i;
            order = NPoints/i;
            break;
        }
        QuadPoints.assign(ints,order);
        for(unsigned int i=0;i<NPoints-ints*order;i++){
            QuadPoints[i]++;
        }
    }
}

void BasisBesselCoulomb::calculateQuadRules(const std::vector<int> &QuadPoints, int &Pos, UseMatrix &X, UseMatrix &W, const double Rmin, const double Rmax) const {
    // calculate the quadrature rules for the standard interval [-1,+1] and transform them
    alglib::ae_int_t n=QuadPoints[0];
    alglib::ae_int_t info;
    alglib::real_1d_array xq;
    alglib::real_1d_array wq;
    alglib::gqgenerategausslegendre(n,info,xq,wq);
    int ints = QuadPoints.size(); // number of quadrature intervals
    for(unsigned int i=0; i<ints; i++){
        if(QuadPoints[i]!=n){
            n=QuadPoints[i];
            alglib::gqgenerategausslegendre(n,info,xq,wq);
        }
        for(unsigned int j=0; j<QuadPoints[i]; j++){
            X(Pos) = Rmin+i*(Rmax-Rmin)/ints+(Rmax-Rmin)/(2.*ints)*(1-xq[QuadPoints[i]-j-1]);
            W(Pos) = (Rmax-Rmin)/(2.*ints)*wq[QuadPoints[i]-j-1];
            Pos++;
        }
    }
}



