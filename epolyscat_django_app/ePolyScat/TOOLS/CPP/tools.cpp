// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "tools.h"
#include "time.h"

#include <thread>
#include <chrono>

#include "useMatrix.h"
#include "qtEigenDense.h"
#include "qtEigenDense.h"
#ifdef _HAS_LAPACKE_
#include "lapacke.h"
#endif
#include "qtAlglib.h"    // get the quadrature rules

using namespace std;

namespace tools {

char zero(const complex<double> & val) {
    if (abs(val)<1.e-13) {
        return '.';
    } else if (abs(val)<1.e-7) {
        return 'o';
    } else if (abs(val)<1.e-1) {
        return 'x';
    } else {
        return 'X';
    }
}

vector<int> fractionSmallPrimes(double Value, double Eps, double MaxDen)
{
    if(Value-int(Value)==0)return {int(Value),1};

    // the following can be done better....
    vector<int> prim={2,3,5,7};
    long long denom=1;
    for (int p: prim){
        int fac=1;
        while (fac*p<MaxDen)fac*=p;
        denom*=fac;
    }

    double dValNum=Value*denom;
    long long num=(long long)(dValNum*(1.+1e-12));

    if(abs(dValNum-num)>Eps*abs(dValNum))return {0,0};

    for(int p: prim)
        while(num%p==0 and denom%p==0){
            num/=p;
            denom/=p;
        }

    return {int(num),int(denom)};
}


bool isEquidistant(const std::vector<double> & Vec)
{
    for(size_t k=2;k<Vec.size();k++){
        if(abs((Vec[k]-Vec[k-1])-(Vec[1]-Vec[0]))>abs(Vec.back()-Vec[0])*1.e-12)return false;
    }
    return true;
}

void sleepSeconds(double Secs){
    std::this_thread::sleep_for(std::chrono::microseconds(int(Secs*1.e6)));
}

string newFile(const string File){
    ifstream file(File.c_str(),ifstream::in);
    if(not file.is_open())return File;

    int i=-1;
    do {
        i++;
        file.close();
        file.open((File+tools::str(i)).c_str(),ifstream::in);
    } while (file.is_open());
    return File+tools::str(i);
}

/// Cholesky decomposition for a (complex) symmetric tri-diagonal matrix
void pseudoCholeskyTri(vector<complex<double> > & Diag, vector<complex<double> > & Sub){

    if(Diag.size()!=Sub.size()+1)ABORT("diagonal must be exactly 1 more than sub-diagaonal");

    vector<complex<double> >d(Diag),s(Sub); // save for checking

    Diag[0]=sqrt(Diag[0]);
    for(unsigned int i=1;i<Diag.size();i++){
        Sub[i-1]/=Diag[i-1];
        Diag[i]=sqrt(Diag[i]-Sub[i-1]*Sub[i-1]);
    }

    // check decomposition
    for(unsigned int i=1;i<Diag.size();i++){
        if(abs(d[i]-Diag[i]*Diag[i]-Sub[i-1]*Sub[i-1])>1.e-12*abs(d[i]))
            ABORT("pseudo-Cholesky failed, diag "+str(i));
        if(abs(s[i-1]-Diag[i-1]*Sub[i-1])>1.e-12*abs(d[i]))
            ABORT("pseudo-Cholesky failed, sub "+str(i));
    }
}

/// Gram-Schmidt transformation Trans^T Metric Trans = 1
/// algorithm: see tsurff.pdf
void gramSchmidtTrans(const vector<vector<double> > &Metric, vector<vector<double> > &Trans){

    // sanity checks for Metric
    if(Trans.size()==0){
        for(unsigned int i=0;i<Metric.size();i++){
            if(Metric[i].size()!=Metric.size())ABORT("Metric matrix is not square");
            for(unsigned int j=0;j<i;j++)
                if(pow(Metric[i][j]-Metric[j][i],2)>1.e-12*Metric[i][i]*Metric[j][j]){
                    for(unsigned int i=0;i<Metric.size();i++)cout<<str(Metric[i],5,",")<<endl;
                    ABORT("Metric is not symmetric");
                }
            if(Metric[i][i]<=0.)ABORT("Metric diagonal is non-positive");
        }
    }

    // add row to transformation matrix
    Trans.push_back(vector<double>(Trans.size()+1,0.));
    for(unsigned int n=0;n<Trans.back().size()-1;n++){
        double tmn=0;
        for(unsigned int j=0;j<Trans[n].size();j++)
            tmn+=Trans[n][j]*Metric[Trans.size()-1][j];
        for(unsigned int k=0;k<Trans[n].size();k++)
            Trans.back()[k]-=Trans[n][k]*tmn;
    }
    Trans.back().back()=1.;

    // compute norm of new row
    double norm=0;
    for(unsigned int k=0;k<Trans.size();k++){
        double msk=0;
        for(unsigned int l=0;l<Trans.back().size();l++)
            msk+=Metric[k][l]*Trans.back()[l];
        norm+=Trans.back()[k]*msk;
    }

    // normalize new row of Trans
    if(abs(norm)<1.e-28)ABORT("Metric ill-conditioned");
    norm=1/sqrt(norm);
    for(unsigned int k=0;k<Trans.back().size();k++)Trans.back()[k]*=norm;

    // descend while needed
    if(Trans.size()<Metric.size()){
        gramSchmidtTrans(Metric,Trans);
    }

    // consistency check: Trans^T Metric Trans = 1
    else {
        vector<double> maxT;
        for(unsigned int n=0;n<Trans.size();n++){
            maxT.push_back(0.);
            for(unsigned int k=0;k<Metric.size();k++)maxT.back()=max(maxT.back(),abs(Metric[n][k]));
            maxT.back()=sqrt(maxT.back());
            vector<double> mtn(Trans.size(),0);
            for(unsigned int i=0;i<Trans.size();i++)
                for(unsigned int j=0;j<Trans[n].size();j++)
                    mtn[i]+=Metric[i][j]*Trans[n][j];

            for(unsigned int m=0;m<=n;m++){
                double dmn=0;
                for(unsigned int i=0;i<Trans[m].size();i++)
                    dmn+=Trans[m][i]*mtn[i];
                if(m==n)dmn-=1.;
                if(abs(dmn)>1.e-12*maxT[m]*maxT[n])ABORT("transformation failed");
            }
        }
    }
}


vector<double> insertInterval(double A,double B, vector<double> Split){
    if(A>B)ABORT("need interval A<B");
    vector<double> ints;
    ints.push_back(A);
    for(unsigned int k=0;k<Split.size();k++){
        if(k>0 and Split[k-1]>Split[k])ABORT("split points must be in ascending order, is: "+tools::str(Split,8,", "));
        if(Split[k]>=B)break;
        if(Split[k]>A)ints.push_back(Split[k]);

    }
    ints.push_back(B);
    return ints;
}

// various comparisons
bool lessReal(const complex<double> &A, const complex<double> &B){return A.real()<B.real();}
bool lessImag(const complex<double> &A, const complex<double> &B){return A.imag()<B.imag();}
bool lessAbs( const complex<double> &A, const complex<double> &B){return abs(A)  <abs(B);}

/// @brief elementary interface to lapack_zggev general complex eigensolver (NOTE: lapack destroys  both matrices)
//void lapack_zggev(UseMatrix &A, UseMatrix &B, vector<complex<double> > & Eval, UseMatrix & vr, bool Vectors){
//    UseMatrix alpha(A.rows(),1);
//    UseMatrix beta(A.rows(),1);
//    UseMatrix vl;
//    vr=UseMatrix(A.rows(),A.cols());
//    UseMatrix Bsave(B); // copy for normalization, as zggev destroys the matrices
//    ///////////////////////////////

//    //check whether A and B are hermitian
//    if(A.isHermitian(1.e-10) and B.isHermitian(1.e-10)) {
//        Eigen::MatrixXcd evecs;
//        Eigen::MatrixXcd evals;
//        if(A.isReal() and B.isReal()){
//            Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(Eigen::Map<Eigen::MatrixXcd>(A.data(),A.rows(),A.cols()).real(),
//                                                                         Eigen::Map<Eigen::MatrixXcd>(B.data(),B.rows(),B.cols()).real());
//            if(Vectors)evecs = es.eigenvectors().cast<complex<double> >();
//            evals = es.eigenvalues().cast<complex<double> >();
//        } else {
//            Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> es(Eigen::Map<Eigen::MatrixXcd>(A.data(),A.rows(),A.cols()),
//                                                                          Eigen::Map<Eigen::MatrixXcd>(B.data(),B.rows(),B.cols()));
//            if(Vectors)evecs = es.eigenvectors().cast<complex<double> >();
//            evals = es.eigenvalues().cast<complex<double> >();
//        }
//        cout<<"evals: "<<evals.size()<<endl;
//        if(Vectors)vr = MapMatrix(evecs.data(),evecs.rows(),evecs.cols());
//        alpha = MapMatrix(evals.data(),evals.rows(),evals.cols());
//        beta = UseMatrix::Constant(alpha.rows(),1,1.0);
//    } else {
//        char vectors='N';
//        if(Vectors)vectors='V';
//        int info=LAPACKE_zggev(LAPACK_COL_MAJOR,'N',vectors,A.rows(),A.data(),A.rows(),B.data(),B.rows(),
//                               alpha.data(),beta.data(),vl.data(),A.rows(),vr.data(),A.rows());
//        if (info!=0){
//            cerr << endl;
//            cerr << "--------------------------------------------------------------------------" << endl;
//            cerr << " WARNING: Stop will be called in eigenvalues_zggev: lapacke error code " << info << endl;
//            cerr << endl;
//            cerr << " Note, for your problem: N = " << A.cols() << endl;
//            cerr << endl;
//            cerr << " ERROR CODE EXPLANATION:" << endl;
//            cerr << "\t = 0:  successful exit " << endl <<
//                    "\t < 0:  if INFO = -i, the i-th argument had an illegal value." << endl <<
//                    "\t =1,...,N:" << endl <<
//                    "\t\t The QZ iteration failed.  No eigenvectors have been " << endl <<
//                    "\t\t calculated, but ALPHA(j) and BETA(j) should be" << endl <<
//                    "\t\t correct for j=INFO+1,...,N." << endl <<
//                    "\t > N:\t =N+1: other then QZ iteration failed in DHGEQZ," << endl <<
//                    "\t\t =N+2: error return from DTGEVC." << endl;
//            cerr << "--------------------------------------------------------------------------" << endl;
//            //            ABORT("Stop!");
//        }
//    }
//    ///////////////////////////////////
//    // (pseudo-)normalize to c^T B c=1
//    for(unsigned int j=0;j<vr.cols();j++){
//        UseMatrix Bvrj=Bsave*vr.col(j);
//        complex<double> norm=0;
//        for (unsigned int i=0;i<vr.rows();i++)norm+=vr(i,j)*Bvrj(i);
//        vr.col(j)/=sqrt(norm);
//    }

//    Eval.clear();
//    vector<vector<complex<double> > > Evec;
//    for(unsigned int i=0;i<vr.cols();i++){
//        if(abs(beta(i))<1e-12 and abs(alpha(i))>1e-12){
//            cerr <<"special case encountered: create fake eigenvalue on imaginary axis" << endl;
//            beta(i)=1.;
//            alpha(i)=complex<double>(0.,-10.);
//        }
//        Eval.push_back(alpha(i)/beta(i));
//        Evec.push_back(vector<complex<double> >(0));
//        if(Vectors)for (unsigned int j=0;j<vr.rows();j++)Evec.back().push_back(vr(j,i));
//}

//UseMatrix::UseMap(Eval.data(),1,std::min(A.rows(),(unsigned int) 20)).print("Eval");


//// write back into UseMatrix
//if(Vectors)
//for (unsigned int i=0;i<Evec.size();i++)
//for (unsigned int j=0;j<Evec[i].size();j++)
//vr(j,i)=Evec[i][j];
//}

void get_quadrature_rule(int order, quadrature_kind kind, alglib::real_1d_array &points, alglib::real_1d_array &weights, const std::vector<double> &endpoints)
{
    static alglib::ae_int_t info; // not used, yet mandatory for alglib quadrature generator

    // memoization for performance
    static int previous_order;
    static quadrature_kind previous_kind;
    static alglib::real_1d_array previous_points, previous_weights;

    if (order==previous_order and kind==previous_kind) { // memoization
        points =previous_points ;
        weights=previous_weights;
    }
    else { // get new quadrature
        switch (kind)
        {
        case gq_legendre:
        {
            alglib::gqgenerategausslegendre(order, info, points, weights);
            break;
        }
        case gq_laguerre:
        {
            alglib::gqgenerategausslaguerre(order, 0., info, points, weights);
            break;
        }
        case gq_hermite:
        {
            alglib::gqgenerategausshermite(order, info, points, weights);
            break;
        }
        case q_equidistant: //standard interval = [-1, 1]
        {
            points.setlength(order);
            weights.setlength(order);
            for (int i=0; i!= order; ++i) {
                points[i] = -1.+2.*(0.5+i)/order;
                weights[i]= 2./order;
            }
            break;
        }
        default:
            ABORT("Unknown quadrature desired.");
        }
        // reset for memoization
        previous_points  = points ;
        previous_weights = weights;
        previous_order   = order;
        previous_kind    = kind;
    }

    if (endpoints.empty()) {return;} // standard interval (depends on quadrature rule)

    switch (kind)
    {
    case gq_legendre:
    case q_equidistant:
    {
        if (endpoints.size()!=2) {
            std::cout << "Don't know how to do quadrature " << kind << " with #" << endpoints.size() << "endpoints!\nUse standard quadrature interval\n";
        } else {
            for (int i=0; i!= points.length(); ++i) {
                points[i]  =(endpoints[0]*(1.-points[i])+endpoints[1]*(1.+points[i]))/2.;
                weights[i]*=(endpoints[1]-endpoints[0])/2.;
            }
            //            std::for_each( points.begin(),  points.end(), _1=(endpoints[0]*(1.-_1)+endpoints[1]*(1.+_1))/2.);
            //            std::for_each(weights.begin(), weights.end(), _1=_1*(endpoints[1]-endpoints[0])/2.);
        }
        break;
    }
    case gq_laguerre:
    {
        if (endpoints.size()!=1) {
            cout << "Don't know how to do quadrature " << kind << " with #" << endpoints.size() << "endpoints!\nUse standard quadrature interval\n";
        } else {
            for (int i=0; i!=points.length(); ++i) {
                points[i]+=endpoints[0];
            }
        }
        break;
    }
        // Gauss-Hermite quadrature has intervals (-\infty, \infty), so there is no adaption
    default:
    {
        cout << "Something unexpected happened in quadrature " << kind << endl;
        cout << "Did you try to adjust a standard interval (-infty, infty)? Will nevertheless continue with standard interval";
    }
    }
}

} // end tools namespace




