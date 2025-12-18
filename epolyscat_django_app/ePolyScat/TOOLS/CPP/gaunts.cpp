// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "gaunts.h"
#include "orthopol.h"

#include "qtEigenDense.h"

#include "timer.h"
#include "printOutput.h"

#ifdef _USE_GSL_
#include <gsl/gsl_sf_coupling.h>
#endif

TIMERSAMPLE(gaunts,)
TIMERSAMPLE(gaunts1,)
TIMER(gauntsCalc,)
std::vector<double> &Gaunts::vals(int M1, int M2, int L1, int L2, int & L3Min){
    // sort by increasing magnitude
    if(abs(M2)<abs(M1)){
        std::swap(M1,M2);
        std::swap(L1,L2);
    }
    // flip signs of M to have negative first
    if(M2<0){
        M1=-M1;
        M2=-M2;
    }

    size_t lm1=L1-abs(M1);
    size_t lm2=L2-abs(M2);
    L3Min=l3Min(M1,M2,L1,L2);

    int lMax;
    std::vector<std::vector<std::vector<double>>> & storeMM=_store[M1][M2];
    if(storeMM.size()<=lm1
            or storeMM[lm1].size()<=lm2
            or storeMM[lm1][lm2].size()<size_t((L1+L2-L3Min)/2)){
        START(gauntsCalc);
        if(lm1<0 or lm2<0)ABORT(Sstr+"illegal angular momenta: l,m"+M1+L1+"or"+M2+L2);

        // recalculate table for M1,M2, and max(L1,L2)
        lMax=std::max(L1,L2);
        _store[M1][M2].clear();

        std::vector<double> x,w;
        OrthogonalNassocLegendre assoc0(0);
        assoc0.quadratureGauss(2*lMax+3,x,w); // quadrature

        for(size_t k=0;k<x.size();k++){

            std::vector<double> v1(OrthogonalNassocLegendre(abs(M1)   ).val(lMax-abs(M1)+1     ,x[k]));
            std::vector<double> v2(OrthogonalNassocLegendre(abs(M2)   ).val(lMax-abs(M2)+1     ,x[k]));
            std::vector<double> v3(OrthogonalNassocLegendre(abs(M1+M2)).val(2*lMax-abs(M1+M2)+1,x[k]));

            if(_store[M1][M2].size()==0)_store[M1][M2].resize(v1.size());
            for(size_t ll1=0;ll1<v1.size();ll1++){
                double w1=w[k]*v1[ll1]*sqrt(2./M_PI);
                if(_store[M1][M2][ll1].size()==0)_store[M1][M2][ll1].resize(v2.size());

                for(size_t ll2=0;ll2<v2.size();ll2++){
                    int l1=ll1+abs(M1);
                    int l2=ll2+abs(M2);

                    double w12=w1*v2[ll2];
                    int llmin=l3Min(M1,M2,l1,l2);
                    if(_store[M1][M2][ll1][ll2].size()==0)_store[M1][M2][ll1][ll2].assign((l1+l2-llmin)/2+1,0.);
                    for(int ll=0;llmin+2*ll<=l1+l2;ll++)
                        _store[M1][M2][ll1][ll2][ll]+=w12*v3[llmin+2*ll-abs(M1+M2)];

                }
            }
        }
#ifdef _USE_GSL_
        if(lMax<40){
            // WARNING: tests fail at large angular mommenta - very likely GSL fails
            // see Johannson, SIAM J. SCI. COMPUT. Vol. 38, No. 1, pp. A376–A384

            // connect to the standard sign convention

            int sig=1;
            if(M1<0 and M1%2)sig=-sig;
            if(-M1-M2<0 and (M1+M2)%2)sig=-sig;

            int nonZero=0,totalSize=0,nSucc=0;
            for(int lm1=0;lm1<=lMax-abs(M1);lm1++){
                for(int lm2=0;lm2<=lMax-abs(M2);lm2++){
                    int l1=lm1+abs(M1);
                    int l2=lm2+abs(M2);
                    int llmin=l3Min(M1,M2,l1,l2);
                    for(int l3=llmin,ll=0;ll<_store[M1][M2][lm1][lm2].size();ll++,l3+=2){
                        totalSize++;
                        if(std::abs(_store[M1][M2][lm1][lm2][ll])>1.e-12){
                            nonZero++;
                            //check
                            double rat=sqrt(double((2*l1+1)*(2*l2+1)*(2*l3+1))/(4*M_PI))
                                    *gsl_sf_coupling_3j(2*l1,2*l3,2*l2,2*M1,-2*(M1+M2),2*M2)
                                    *gsl_sf_coupling_3j(2*l1,2*l3,2*l2,0,0,0)
                                    /(sig*storeMM[lm1][lm2][ll]);
                            if(std::abs(rat-1.)*std::abs(_store[M1][M2][lm1][lm2][ll])>1.e-9){
                                Sstr+"failure"+(rat)+"[success"+nSucc+"] m1,m2,l1,l2,l3"+M1+M2+l1+l1+l3+_store[M1][M2][lm1][lm2][ll]+Sendl;
                                Sstr+"parts"
                                        +sqrt(double((2*l1+1)*(2*l2+1)*(2*l3+1))/(4*M_PI))
                                        +gsl_sf_coupling_3j(2*l1,2*l3,2*l2,0,0,0)
                                        +gsl_sf_coupling_3j(2*l1,2*l3,2*l2,2*M1,-2*(M1+M2),2*M2)
                                        +storeMM[lm1][lm2][ll]
                                        +Sendl;
                            }
                            else nSucc++;
                        }
                        else if(abs(gsl_sf_coupling_3j(2*l1,2*l3,2*l2,2*M1,-2*(M1+M2),2*M2))>1.e-10){
                            Sstr+"zero failure"+"[success"+nSucc+"] m1,m2,l1,l2,l3"+M1+M2+l1+l2+l3+_store[M1][M2][lm1][lm2][ll]+Sendl;
                            Sstr+"parts"
                                    +sqrt(double((2*l1+1)*(2*l2+1)*(2*l3+1))/(4*M_PI))
                                    +gsl_sf_coupling_3j(2*l1,2*l3,2*l2,0,0,0)
                                    +gsl_sf_coupling_3j(2*l1,2*l3,2*l2,2*M1,-2*(M1+M2),2*M2)
                                    +storeMM[lm1][lm2][ll]
                                    +Sendl;
                            COUNTDOWN("zeros",20);
                        }
                    }
                }
            }
            if(abs(totalSize-nonZero)>200)PrintOutput::DEVwarning(Sstr+"large number of zeros in Gaunt: total vs nonZero"+totalSize+nonZero);
        }
#endif
        STOP(gauntsCalc);
    }
    return storeMM[lm1][lm2];
}
