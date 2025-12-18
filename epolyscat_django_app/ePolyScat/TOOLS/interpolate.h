// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <complex>
#include <vector>
#include <complex>
#include "abort.h"
#include "tools.h"

/// coefficients for Newton interpolation by diveded differences
template<class Bas,class Val>
void dividedDifference(const std::vector<Bas> & B /** base points */ ,
                       const std::vector<Val> & V /** values */,
                       std::vector<Val> & Coeff){
    // Newton interpolation coefficients by divided differences
    // this algorithm is data-local to the level of Val and minimizes divisions
    // it may under- or overflow at high orders
    // algorithm works on Val from any linear space
    // for longer vectors Val it is not optimal in terms of data locality
    Coeff=V;
    std::vector<Bas> fac(B.size(),1.); // temporary factor on divided differences
    for(unsigned int n=1;n<V.size();n++){
        for(unsigned int i=B.size()-1;i>=n;i--){
            // standard algorithm
            Coeff[i]-=Coeff[i-1];
            Coeff[i]*=(1./(B[i]-B[i-n]));
            // smart(?) algorithm
            //            Coeff[i].axpy(-fac[i],Coeff[i-1]);
            //            fac[i]*=B[i]-B[i-n];
        }
        //        Coeff[n]*=1./fac[n]; // recursion for n terminated, remove factor
    }
}

/// Newton interpolation template class
/// Bas .....element of algebra
/// Val .....linear space over Bas
template<class Bas, class Val>
class Interpolate
{
protected:
    typename std::vector<Bas>::iterator iBas;
    std::vector<Bas> bas;
    std::vector<std::vector<Val> > coeff;
public:
    virtual Val& val(const Bas &X, Val & Res)=0;

    unsigned int points(){return coeff[0].size();}

    /// evaluate a polynomial p(X) using Horner's scheme
    ///
    /// not optimally data-local for larger result vectors
    Val horner(const Bas & X,
               const std::vector<Bas> & B /** base points */ ,
               const std::vector<Val> & C /** coefficients */)
    {
        // evaluate by Horner scheme
        // for longer vectors res it is not optimal in terms of data locality
        Val res(C.back());
        for(unsigned int k=C.size();k>0;k--){
            res*=(X-B[k-1]);
            res+=C[k-1];
        }
        return res;
    }

    /// check for mononitonic support points, optionally remove duplicates
    /// returns higher index of non-increasing pair
    /// return value = 0 indicates success
    static size_t notIncreasing(std::vector<Bas>& B, std::vector<Val> V, bool RemoveDuplicate=false){
        double eps=(B.back()-B.front())/B.size()*1.e-7;
        for(unsigned int k=B.size()-1;k>0;k--){
            if(B.size()!=V.size())DEVABORT("bas and values do not have equal size");
            if(B[k]-B[k-1]<eps){
                if(std::abs(B[k]-B[k-1])<=eps
                        and (V[k]-V[k-1]).maxAbsVal()<1.e-7*(V[k]).maxAbsVal()
                        and RemoveDuplicate)
                {
                    std::cout<<"removing "<<k<<std::endl;
                    V.erase(V.begin()+k);
                    B.erase(B.begin()+k);
                }
                else {
                    std::cout<<"diffs "<<B[k]-B[k-1]<<" "<<(V[k]-V[k-1]).maxAbsVal()<<" "<<1.e-7*(V[k]).maxAbsVal()<<std::endl;
                    return k;
                }
            }
        }
        return 0;
    }
};

template<class Bas, class Val>
class InterpolateNewton: public Interpolate<Bas,Val>{
    using Interpolate<Bas,Val>::bas;    // C++? surprising it does not know where it is derived from
    using Interpolate<Bas,Val>::coeff;  // C++? surprising it does not know where it is derived from
    using Interpolate<Bas,Val>::iBas;   // C++? surprising it does not know where it is derived from
    using Interpolate<Bas,Val>::horner; // C++? surprising it does not know where it is derived from
    using Interpolate<Bas,Val>::points; // C++? surprising it does not know where it is derived from
public:
    InterpolateNewton(const std::vector<Bas> &BasX, const std::vector<Val> &V, unsigned int Points=4)
    {
        bas=BasX;
        // checks
        if(BasX.size()<Points)ABORT("too few bas points");
        // base points must be strictly in increasing order and cannot bee too narrowly spaced
        double eps=(bas.back()-bas.front())/bas.size()*1.e-7;
        for(unsigned int k=1;k<bas.size();k++)
            if(bas[k-1]>=bas[k]-eps){
                std::cout<<"k, x[k-1], x[k] "<<k<<" "<<bas[k-1]<<" "<<bas[k]<<std::endl;
                ABORT("non-increasing or near-conincident points");
            }

        // local coefficients by divided differences
        coeff.clear();
        for(iBas=bas.begin();iBas!=bas.end()-Points;iBas++){
            coeff.push_back(std::vector<Val>());
            dividedDifference(std::vector<Bas>(iBas,iBas+Points),
                              std::vector<Val>(V.begin()+(iBas-bas.begin()),V.begin()+(iBas-bas.begin())+Points),
                              coeff.back());

            //            // check
//            for(typename std::vector<Bas>::iterator it=iBas;it!=iBas+Points;it++){
//                if((V[it-bas.begin()]-horner(*it,std::vector<Bas>(iBas,iBas+Points),coeff.back())).maxAbsVal()>1.e-10){
//                    std::cout<<"base: "<<tools::str(V[it-bas.begin()])<<std::endl;
//                    std::cout<<"      "<<tools::str(horner(*it,std::vector<Bas>(iBas,iBas+Points),coeff.back()))<<std::endl;
//                    ABORT("failed");
//                }
//            }
        }
        iBas=bas.begin();
    }

    Val & val(const Bas &X, Val &Res){

        // get first iBas above X (avoid margins)
        if(X>*iBas)iBas=std::upper_bound(iBas+1,bas.end(),X);
        else       iBas=std::lower_bound(bas.begin(),iBas,X);

        // check for margins
        if(iBas-bas.begin()<points()/2){
            if(X<bas.front())ABORT("X below range, cannot extrapolate");
            iBas=bas.begin();
        } else if(bas.end()-iBas<points()+1){
            if(X>bas.back())ABORT("X above range, cannot extrapolate");
            iBas=bas.end()-points()-1;
        } else {
            iBas-=points()/2;
        }

        // evaluate by Horner scheme
        return Res=horner(X, std::vector<Bas>(iBas,iBas+coeff[iBas-bas.begin()].size()),coeff[iBas-bas.begin()]);
    }

};


#endif // INTERPOLATE_H
