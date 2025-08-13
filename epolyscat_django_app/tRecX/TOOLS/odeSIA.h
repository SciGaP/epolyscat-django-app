// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef ODESIA_H
#define ODESIA_H

#include "odeStep.h"

#include <complex>
#include <vector>
#include <complex>
#include "arnoldi.h"
#include "abort.h"

/// Short Iterative Arnoldi
template<class Der,class V>
class OdeSIA: public OdeStep<Der,V>
{
    /// propagate first Krylov vector in present Krylov subspace by Tstep
    void krylovExtrapolate(double Tstep,std::vector<std::complex<double> > &Coefs){
        Coefs.assign(eVal.size(),0.);
        for(unsigned int k=0,ik=0,k0=0;k<eVal.size();k++,k0++)
            for(unsigned int i=0;i<eVal.size();i++,ik++)
                Coefs[i]+=rVec[ik]*std::exp(Tstep*eVal[k])*lVec[k0];
    }

    unsigned int nStep; // increase Krylov by nStep, base error estimate on this, nStep>=2 recommended
    Arnoldi<Der,V> * arn;// Arnoldi scheme
    double epsSqu;      // L2-norm accuracy estimate
    unsigned int maxK;  // maximal Krylov dimension
    std::vector<std::complex<double> > eVal,lVec,rVec,cExpIt; // keep for external use of krylovExtrapolate

public:
    OdeSIA(Der *Op, unsigned int MaxKrylov, double EpsSquared)
        :OdeStep<Der,V>("SIA",Op),arn(new Arnoldi<Der,V>(*Op,Op->lhsVector())),
          maxK(std::min(MaxKrylov,(unsigned int)Op->lhsVector().size())),
          epsSqu(EpsSquared),nStep(2){}

    V & step(V &Vec, double Tstart, double Tstep)
    {
        OdeStep<Der,V>::derOde->update(Tstart+Tstep/2.);

        std::complex<double> nrm=sqrt(Vec.scalarProduct(Vec));
        arn->reset(Vec);
        arn->extend(nStep);
        eVal.clear();
        lVec.clear();
        rVec.clear();
        cExpIt.clear();

        bool invariantSubspace=false;
        unsigned int dim=arn->krylovDim();
        double errSqu=1.;
        std::vector<std::complex<double> > cOld;
        while(not invariantSubspace and dim<maxK and errSqu>epsSqu){
            dim=std::min(dim+nStep,maxK);

            // extend the Krylov space by Arnoldi iteration
            invariantSubspace=arn->extend(dim);

            // solve eigenproblem for Arnoldi matrix
            arn->eigen(eVal,lVec,rVec);

            // get the new time-propagated coefficients
            cOld=cExpIt;
            krylovExtrapolate(Tstep,cExpIt);

            // check convergence of L2-norm
            if(cOld.size()>2)errSqu=0.;
            for(unsigned int k=0;k<cOld.size();k++)errSqu+=std::norm(cExpIt[k]-cOld[k]);
        }
        OdeStep<Der,V>::nCallsStep=arn->krylovDim();
        OdeStep<Der,V>::nCalls+=OdeStep<Der,V>::nApplyStep();


        if(invariantSubspace)std::cout<<"!!! WARNING: Arnoldi invariant subspace, dimension="<<dim<<std::endl;
        if(dim==maxK)std::cout<<"!!! WARNING maximal of N="<<arn->krylovDim()<<" Arnoldi iteration dimensions reached"<<std::endl;

        Vec.axpy(cExpIt[0]*nrm,arn->krylovVectors()[0],0.);
        for(unsigned int k=1;k<cExpIt.size();k++)Vec.axpy(cExpIt[k]*nrm,arn->krylovVectors()[k],1.);
        return Vec;
    }

    V & stepErrorSia(V &Vec, double Tstart, double Tstep, V &Err){
        ABORT("disabled");
        double fCenter=1./7.; // this ratio determines factor in front (and whether it is realistic)
        double tLeft =Tstart+(1.-fCenter)/2.*Tstep;
        double tRight=tLeft+fCenter*Tstep;

        // propagate begin to tLeft
        step(Vec,Tstart,tLeft-Tstart);

        // extrapolate left Krylov to tCenter using left interval operator
        std::vector<std::complex<double> > bErr;
        krylovExtrapolate(Tstep/2.,bErr);
        Err.axpy(bErr[0],arn->krylovVectors()[0],0.);
        for(unsigned int k=1;k<bErr.size();k++)Err.axpy(bErr[k],arn->krylovVectors()[k],1.);

        // propagate tLeft to tRight
        step(Vec,tLeft,tRight-tLeft);

        // propagate tRight to end
        step(Vec,tRight,Tstart+Tstep-tRight);

        // back-extrapolate right Krylov to tCenter using right interval operator
        krylovExtrapolate(-fCenter*Tstep/2.,bErr);
        for(unsigned int k=0;k<bErr.size();k++)Err.axpy(-bErr[k],arn->krylovVectors()[k],1.);

        // (later, the step size control will be administred here directly)
    }

    std::string name() const{return "SIA";}
    unsigned int consistencyOrder() const {return 2;}
    unsigned int currentKrylov(){return arn->krylovDim();}
    unsigned int maxKrylov(){return arn->krylovVectors().size()-1;}
    static void test();
protected:
    double safetyFactor() const {return 0.8;}
};

#endif // ODESIA_H
