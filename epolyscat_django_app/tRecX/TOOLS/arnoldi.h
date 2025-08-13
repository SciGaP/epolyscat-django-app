// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef ARNOLDI_H
#define ARNOLDI_H

#include <complex>
#include <vector>
#include "abort.h"
#ifdef _HAS_LAPACKE_
#include "lapacke.h"
#endif

/// Arnoldi reduction to upper Hessenberg form
///
/// M... linear map, as specified in LinSpaceMap <br>
/// V... Hilbert space vector, as specified in LinSpaceHilbert
template<class M, class V>
class Arnoldi
{
    const M * map;
    std::vector<std::complex<double> > mat; // upper Hessenberg matrix in packed storage
    std::vector<V> kry;                     // storage for Krylov vectors, size() = maximal dim+1 during object lifetime
public:
    Arnoldi(const M& Map,const V & Vec):map(&Map){kry.assign(1,Vec);reset(Vec);}

    void reset(const V & Vec){
        mat.assign(1,std::sqrt(Vec.scalarProduct(Vec)));
        kry[0]=Vec;
        if(mat[0]!=0.)kry[0]*=1./mat[0];
    }

    /// extend Krylov space and Arnoldi matrix (true if invariant subspace)
    bool extend(unsigned int & DimExt){

        if(DimExt==krylovDim())return false; // no extension needed
        if(DimExt<krylovDim())ABORT("cannot shrink Arnoldi matrix, use reset(..) to for new starting vector");
        if(DimExt>kry[0].size())ABORT("cannot extend beyond linear space dimension");

        // algorithm (see Wikipedia):
        // kry[krylovDim()]...present residual vector (=rk)
        // map ...............upper Hessenberg matrix in packed storage, last element norm of residual

        unsigned int k=krylovDim();
        unsigned int jk=mat.size()-1; // last element contains old residual norm
        mat.resize((DimExt*(DimExt+1))/2+DimExt+1,0.);
        if(kry.size()<DimExt+1)kry.resize(DimExt+1,kry[0]);
        for(;k<DimExt;k++,jk++){

            // new vector
            map->apply(1.,kry[k],0.,kry[k+1]);

            // orthogonalize
            for (unsigned int j=0;j<k+1;j++,jk++){
                mat[jk]=kry[j].scalarProduct(kry[k+1]);
                kry[k+1].axpy(-mat[jk],kry[j],1.);
            }
            // residual norm
            mat[jk]=std::sqrt(kry[k+1].scalarProduct(kry[k+1]));

            if(k+1<kry[0].size() and std::abs(mat[jk])<1.e-14){
                // invariant subspace
                DimExt=k+1;
                mat.resize((DimExt*(DimExt+1))/2+DimExt+1,0.);
                return true;
            }
            // normalize
            kry[k+1]*=1./mat[jk];
        }
        return false;
    }

    /// return Arnoldi matrix in full storage (column major)
    std::vector<std::complex<double> > matrix() const {
        unsigned int dim=krylovDim();
        std::vector<std::complex<double> > m(dim*dim,0.);
        unsigned int ij=0;
        for(unsigned int j=0;j<dim;j++)
            for(unsigned int i=0;i<std::min(j+2,dim);i++,ij++)
                m[i+dim*j]=mat[ij];
        return m;
    }

    /// all eigenvalues and -vectors of Arnoldi matrix
    void eigen(std::vector<std::complex<double> > & Val,
               std::vector<std::complex<double> > & LVec,
               std::vector<std::complex<double> > & RVec)
    {
        std::vector<std::complex<double> > h(matrix());
        Val.resize(krylovDim());
        RVec.resize(krylovDim()*krylovDim());
        LVec.resize(krylovDim()*krylovDim());

        // special case: krylovDim=1
        if(krylovDim()<2){
            Val[0]=h[0];
            RVec[0]=1.;
            LVec[0]=1.;
            return;
        }

        DEVABORT("eliminating lapacke from code - disabled");
//        lapack_int one=1,dim=krylovDim();
//        double abnrm;
//        std::vector<double> scale(krylovDim()),rconde(krylovDim()),rcondv(krylovDim());

//        int info=LAPACKE_zgeevx(LAPACK_COL_MAJOR,'B','V','V','N',krylovDim(),h.data(),krylovDim(),
//                                Val.data(),LVec.data(),krylovDim(),RVec.data(),krylovDim(),
//                                &one,&dim,scale.data(),&abnrm,rconde.data(),rcondv.data());
//        if(info!=0)ABORT("eigenproblem for Arnoldi matrix failed (zgeevx)");

//        // LU factorize right eigenvectors
//        std::vector<lapack_int> ipiv(dim);
//        LVec=RVec;
//        info=LAPACKE_zgetrf(LAPACK_COL_MAJOR,dim,dim,LVec.data(),dim,ipiv.data());
//        if(info!=0)ABORT("eigenvectors not linearly independent, LU failed (zgetrf)");

//        // get inverse
//        info=LAPACKE_zgetri(LAPACK_COL_MAJOR,dim,LVec.data(),dim,ipiv.data());
//        if(info!=0)ABORT("inverse failed (zgetri)");

    }

    /// return current Krylov vector storage (CAUTION: krylovVector.size() != krylovDim())
    const std::vector<V> & krylovVectors() const {return kry;}

    /// actual dimension of Krylov space (CAUTION: krylovVector.size() != krylovDim())
    unsigned int krylovDim() const {return (unsigned int)(std::sqrt(double(2*mat.size()+0.25000001))-1.49999);}

};

#endif // ARNOLDI_H
