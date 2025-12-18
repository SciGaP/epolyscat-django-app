// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
//#include "../multiPole.h"
//#include "useMatrix.h"

//#include "basisAbstract.h"

//using namespace std;

//vector<const BasisSet*> MultiPole::Ints::basis;
//vector<complex<double> > MultiPole::Ints::funcVal;
//int MultiPole::Ints::test=0;

//MultiPole::MultiPole(const BasisSet &Ibas,const BasisSet &Jbas,  unsigned int Lmax, double Epsilon)
//{

//    // return if already in table
//    vector<vector<UseMatrix> > * tab=&svdVl[hash(Ibas,Jbas)];
//    if(tab->size()>Lmax)return;

//    // set up the integrator
//    vector<vector<unsigned int> > nquad(1,vector<unsigned int>(1,Ibas.size()+2));
//    nquad.push_back(vector<unsigned int>(1,Jbas.size()+2));
//    Ints integ(Lmax+1,Ibas,Jbas,Epsilon,Epsilon);

//    // recursively compute the integrals
//    vector<vector<double> > vol(2,vector<double>(2));
//    for(unsigned int i=0;i<2;i++){
//        vol[i][0]=Ints::basis[i]->lowBound();
//        vol[i][1]=Ints::basis[i]->upBound();
//    }
//    vector<complex<double> > mpInts=integ.recursive(vol,Ints::multiPot);
//    UseMatrix::UseMap(mpInts.data(),mpInts.size()/(Lmax+1),Lmax+1).show("matrices");

//    // perform basic test for non-overlapping regions
//    vector<const BasisSet*> bas(Ints::basis);
//    if(not basisSet()->hasOverlap(*bas[1])){
//        vector<vector<complex<double> > > mpFactor;
//        for(unsigned int i=0;i<2;i++){
//            if(vol[i][0]<vol[(i+1)%2][0])Ints::test=+1; // current is at smaller radius
//            else                         Ints::test=-1;
//            vector<vector<double> > vol1(1,vol[i]);
//            Ints intTest(Lmax+1,*bas[i],*bas[i],Epsilon,Epsilon);
//            mpFactor.push_back(intTest.recursive(vol1,Ints::multiPot));
//        }
//        unsigned int ij=0,isiz=basisSet()->size(),jsiz=bas[1]->size();
//        vector<complex<double> > prod;
//        for(unsigned int l=0;l<Lmax+1;l++)
//            for(unsigned int j=0;j<mpFactor[1].size()/(Lmax+1);j++)
//                for(unsigned int i=0;i<mpFactor[0].size()/(Lmax+1);i++){
//                    prod.push_back(mpFactor[0][i+l*isiz]*mpFactor[1][j+l*jsiz]);
//                }
//        double eps=max(Epsilon,1.e-5);
//        UseMatrix diff=UseMatrix::UseMap(prod.data(),prod.size()/(Lmax+1),Lmax+1).relativeDifference(
//                    UseMatrix::UseMap(mpInts.data(),mpInts.size()/(Lmax+1),Lmax+1),eps);
//        if(diff.maxAbsVal()>Epsilon)diff.show("relative differences exceed "+tools::str(eps));
//        else PrintOutput::message("OK: product form matches 2d integration");

//        Ints::test=0;
//    }

//    // singular value decompose and put into table
//    unsigned int lMin=tab->size();
//    tab->resize(Lmax+1,vector<UseMatrix>(2));
//    for(unsigned int l=lMin;l<=Lmax;l++){
//        UseMatrix::UseMap(mpInts.data()+l*(mpInts.size()/(Lmax+1)),Ibas.size(),Jbas.size()).svd(tab->at(l)[0],tab->at(l)[1]);
//    }
//}

//MultiPole::Ints::Ints(const unsigned int Nintegrands, const BasisSet &IBas, const BasisSet &JBas,
//                      double AccRel, double AccAbs, std::vector<std::vector<unsigned int> > NQuad,const string Kind, string KindInf)
//    :Tools(AccRel,AccAbs,NQuad,Kind,KindInf, 20)
//{

//    currentDim=0;
//    funcVal.resize(Nintegrands);

//    // currentBasis sets for integration
//    basis.clear();
//    basis.push_back(&IBas);
//    if(test==0)basis.push_back(&JBas); // test!=0 use only first set of functions

//    // default quadrature points
//    if(nQuad.size()==0){
//        for(unsigned int k=0;k<basis.size();k++)nQuad.push_back(vector<unsigned int>(1,8));
//    }
//}

//string MultiPole::hash(const BasisSet &Ibas, const BasisSet &Jbas) const{
//    return Ibas.str()+"|"+Jbas.str();
//}

//vector<complex<double> > MultiPole::Ints::multiPot(const vector<double>&Q){
//    std::vector<std::complex<double> >zQ;
//    for(unsigned int k=0;k<basis.size();k++)
//        zQ.push_back(basis[k]->comsca().xScaled(Q[k]));

//    complex<double> r;
//    switch (test){
//    case 0:
//        if(Q[0]<Q[1]){funcVal[0]=zQ[0]*zQ[0]*zQ[1];r=zQ[0]/zQ[1];}
//        else         {funcVal[0]=zQ[0]*zQ[1]*zQ[1];r=zQ[1]/zQ[0];}
//        break;
//    case 1:
//        // calculate multiple moment
//        if(basis.size()!=1)ABORT("invalid call of multipole test");
//        funcVal[0]=zQ[0]*zQ[0];r=zQ[0];
//        break;
//    case -1:
//        // integrate over multipole potential
//        if(basis.size()!=1)ABORT("invalid call of multipole test");
//        funcVal[0]=zQ[0];r=1./zQ[0];
//        break;
//    default: ABORT("illegal test flag="+tools::str(test));
//    }
//    for(unsigned int i=1;i<funcVal.size();i++)funcVal[i]=funcVal[i-1]*r;

//    return funcVal;
//}

//vector<complex<double> > MultiPole::Ints::nDim(const vector<vector<double> > Vol,  tr1::function<vector<complex<double> >(std::vector<double>)> Func, const vector<double> Params)
//{
//    // n-dimensional integration, recursive procedure starting from the last variable
//    // adapted for multiple integrals for the form I[k0,...,kN;n] = int dx0... int dxN G[n](x0,...,xN) Prod[i] f[ki](xi)
//    // result is stored in a vector, where k0 is the fasted running index, n is slowest

//    if(currentDim==0){
//        xvals.resize(Vol.size());
//    }
//    if(currentDim==Vol.size()){
//        funcCount++;
//        return multiPot(xvals);
//    }

//    // return matrix corresponding to lower level
//    vector<double> qPoin(nQuad[currentDim][0]);
//    vector<double> qWeig(qPoin.size());
//    quadRule(Vol[currentDim],qPoin,qWeig);

//    // next lower result storage
//    unsigned int ldBlock=1;
//    for(unsigned int k=currentDim+1;k<basis.size();k++)ldBlock*=basis[k]->size();
//    vector<complex<double > > block(ldBlock*funcVal.size());

//    // result storage
//    unsigned int ldResult=ldBlock*basis[currentDim]->size();
//    vector<complex<double> > result(ldResult*funcVal.size(),0.);

//    // loop through quadrature points
//    for (unsigned int q=0;q<qPoin.size();q++){
//        xvals[currentDim]=qPoin[q];
//        UseMatrix val=basis[currentDim]->val(UseMatrix::Constant(1,1,qPoin[q]));
//        currentDim++;
//        block=nDim(Vol,Func,Params);
//        currentDim--;
//        for(unsigned int i=0;i<block.size();i++)block[i]*=(qWeig[q]*basis[currentDim]->jacobian(qPoin[q])*basis[currentDim]->eta());

//        // threshold for skipping function values
//        double epsVal=1.e-10*val.maxAbsVal();

//        // loop through integrands
//        for (unsigned int n=0;n<funcVal.size();n++){

//            // loop through factor functions
//            for (unsigned int k=0;k<val.cols();k++){
//                if(abs(val(k).complex())>epsVal){
//                    // current factor fuction value (UseMatrix element access is slow)
//                    complex<double> valk=val(k).complex();

//                    // add into section of column
//                    // "column-wise" sorting: lowest function index runs fastest
//                    unsigned int cBlock =n*ldBlock;         // start index in block
//                    unsigned int cResult=n*ldResult+k;      // start index in result
//                    for(unsigned int i=0;i<ldBlock;i++)result[i*ldBlock+cResult]+=block[i+cBlock]*valk;
//                }

//            }
//        }
//    }
//    return result;
//}
