// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "basisFunction.h"
#include "orthopol.h"
#include "tools.h"
#include "constants.h"
#include "printOutput.h"
#include "readInput.h"

#include "qtAlglib.h"
#ifdef _USE_HACC_
#include "basisfunctionciion.h"
#include "basisfunctioncineutral.h"
#endif

using namespace std;
using namespace physics;
using namespace math;

map<string,const BasisFunction*> BasisFunction::tableNew;

bool BasisFunction::asympZero(const string Name){
    if(Name=="grid")return false;
    if(Name.find("assocLegendre")==0)return false;
    if(Name.find("custom")==0)return false;
    if(Name.find("CI")!=string::npos)return false;
//    if(Name.find("sqrt*")==0)return get(Name.substr(5),3,vector<double>())->asympZero();
    if(Name.find("sqrt")==0)return get(Name.substr(Name.find("*")+1),3,vector<double>())->asympZero();
    if(Name.find("Rpow")==0)return get(Name.substr(Name.find("*")+1),3,vector<double>())->asympZero();
    if(Name.find("besselCoulomb")==0)return false;
    if(Name.find("Eigenbasis[")==0)return false;
    if(Name.find("Orbital[")==0)return false;
    if(Name.find("Channel")==0)return false;
    if(Name.find("PolarOff")==0)return false;
    const BasisFunction* b=get(Name,3,vector<double>());
    return b->asympZero();
}

BasisFunction::BasisFunction(std::string NameStr,unsigned int Order):order(Order),leftInfty(false),rightInfty(false),nameStr(NameStr){}

const BasisFunction * BasisFunction::get(string NameStr, unsigned int Order, std::vector<double> Par){
//    if(not ReadInput::main.flag("DEBUGallowGet","allow use of BasisFunction get"))DEVABORT("for using, specify -DEBUGallowGet");
    string hash=NameStr;
    hash+=tools::str(Order)+"|"+tools::str(Par);
    if(tableNew.count(hash)==1)return tableNew[hash];

    // small factory...
    //OBSOLESCENT only needed to get some paramters...
    const BasisFunction* cur(0);
    if(NameStr.find("Rl*")==0)
        cur=(new BasisFunctionRL(get(NameStr.substr(3),Order,Par),Par));
    if(NameStr.find("sqrt*")==0)
    {DEVABORT("sqrt*-type functions out of service: "+NameStr);}
    else if(NameStr=="useIndex")
        cur=(new BasisFunctionIndex(Order));
#ifdef _USE_HACC_
    else if(NameStr.find("CIion")==0)
        cur=(new BasisFunctionCIion(Order));
    else if(NameStr.find("CIneut")==0)
        cur=(new BasisFunctionCINeutral(Order));
#endif
    else if(NameStr.find("polExp[")==0)
        cur=(new BasisFunctionDvr(new OrthogonalLaguerre,Order));
    else if(NameStr=="polynomial")
        cur=(new BasisFunctionDvr(new OrthogonalLegendre,Order));
    else if(NameStr=="cosSin")
        cur=(new BasisFunctionCosSin(Order));
    else if(NameStr.find("expIm")==0)
        cur=(new BasisFunctionPlaneWave(Order,Par));
    else if(NameStr.find("assocLegendre")==0){
        cur=(new BasisFunctionOpolWeig(new OrthogonalNassocLegendre((int)Par[0]),Order));
        const_cast<BasisFunction*>(cur)->par=vector<double>(1,Par[0]);
    }
    else
        ABORT("unknown function name: "+NameStr+", known: Rl*, sqrt*, polExp, polynomial, cosSin, expIm, assocLegendre, useIndex");


    if(cur)tableNew[hash]=cur;
    return cur;
}

UseMatrix BasisFunction::val(const UseMatrix &X) const {
    UseMatrix V,D;
    valDer(X,V,D);
    return V;
}
UseMatrix BasisFunction::der(const UseMatrix &X) const {
    UseMatrix V,D;
    valDer(X,V,D);
    return D;
}

void BasisFunction::valDer(const UseMatrix &X, UseMatrix &Val, UseMatrix &Der, bool ZeroOutside) const{
    vector<complex<double> > val,der;
    Val=UseMatrix(X.size(),order);
    Der=UseMatrix(X.size(),order);
    for(unsigned int i=0;i<X.size();i++){
        valDer(X(i).complex(),val,der,ZeroOutside);
        for(unsigned int j=0;j<val.size();j++){
            Val(i,j)=val[j];
            Der(i,j)=der[j];
        }
    }
}

void BasisFunctionCosSin::valDer(const std::complex<double> &X, std::vector<std::complex<double> > &Val, std::vector<std::complex<double> > &Der, bool ZeroOutside) const
{
    Val.assign(order,0.);
    Der.assign(order,0.);
    checkStandardInterval(std::real(X),ZeroOutside);

    for(unsigned int j=0;j<order;j++){
        if(j==0){
            Val[j]=1.;
            Der[j]=0.;
        }
        else if(j%2==0){
            Val[j]=             cos(2*pi*par[j]*std::real(X))*sqrt(2.);
            Der[j]=-2*pi*par[j]*sin(2*pi*par[j]*std::real(X))*sqrt(2.);
        }
        else {
            Val[j]=             sin(2*pi*par[j]*std::real(X))*sqrt(2.);
            Der[j]= 2*pi*par[j]*cos(2*pi*par[j]*std::real(X))*sqrt(2.);
        }
    }
}

void BasisFunctionPlaneWave::valDer(const std::complex<double> &X, std::vector<std::complex<double> > &Val, std::vector<std::complex<double> > &Der, bool ZeroOutside) const
{
    Val.assign(order,0.);
    Der.assign(order,0.);
    checkStandardInterval(std::real(X),ZeroOutside);
    complex<double> expIx=complex<double>(cos(deltaK*real(X)),sin(deltaK*real(X)));
    for(unsigned int j=0;j<order;j++){
        Val[j]=std::pow(expIx,(int)par[j]);
        Der[j]=complex<double>(0,deltaK*par[j])*Val[j];
    }
}


int mAssocLeg=0;
double wAssocLeg(double X){return 1.;}
double dAssocLeg(double X){if(mAssocLeg==0) return 0.; return mAssocLeg*X*pow(1-X*X,mAssocLeg-1);}

void BasisFunction::Test(bool Print){
    BasisFunctionOpolWeig(new OrthogonalLaguerre(),6).test(Print);
    BasisFunctionDvr(new OrthogonalLegendre(),6).test(Print);
    BasisFunctionDvr(new OrthogonalLaguerre,10).test(Print);
    BasisFunctionCosSin(9).test(Print);
    BasisFunctionPlaneWave(9).test(Print);
    BasisFunctionOpolWeig(new OrthogonalLaguerre,10).test(Print);
    BasisFunctionOpolWeig(new OrthogonalLegendre(),5).test(Print);
    BasisFunctionOpolWeig(new OrthogonalNassocLegendre(0),5).test(Print);
}

void BasisFunction::test(bool Print) const{
    UseMatrix x,w;
    unsigned int nMin=order;
    // for periodic functions, need one more point
    if((name()=="CosSin" or name()=="PlaneWave") and nMin%2==0)nMin++;
    for (unsigned int o=nMin;o<nMin+7;o+=3){

        if(name().find("DVR:")==0){
            if(o!=order)continue;
            dvrRule(x,w);
        } else {
            quadRule(o,x,w);
        }

        UseMatrix vals=val(x);
        UseMatrix weig=UseMatrix(w.rows(),w.rows());
        for (unsigned int i=0;i<w.rows();i++) weig(i,i)=w(i);
        UseMatrix ovr;
        ovr=vals.adjoint()*weig*vals;

        if(not ovr.isIdentity(1.e-12)){
            ovr.print("BAD BasisFunction "+name()+" test failed, see plot BAD:....");
            plot("BAD:");

        } else {
            cout<<"OK "+name()+" order="<<order<<" quadPts="<<o;
            cout<<", normsq: ";
            for(unsigned int k=0;k<ovr.rows();k++)cout<<" "<<ovr(k,k).real();
            cout<<endl;
        }

    }
}

void BasisFunction::plot(string Dir) const{
    UseMatrix q(101,1),val,der;
    double scal=1.;
    if(rightInfty)scal=2.;
    for(unsigned int k=0;k<q.rows();k++)q(k)=scal*k/double(q.rows()-1);
    valDer(q,val,der);
    ofstream out;
    string file=Dir+"Basis:"+name();
    out.open(file.c_str());
    if (not out.is_open())ABORT("could not open file "+file);
    for(unsigned int k=0;k<val.rows();k++){
        out<<q(k).real();
        for (unsigned int l=0;l<val.cols();l++)
            out<<", "<<val(k,l).real()<<", "<<der(k,l).real();
        out<<endl;
    }
}

BasisFunctionRL::BasisFunctionRL(const BasisFunction * BasFun,std::vector<double>Par):BasisFunction("rL*"+BasFun->name(),BasFun->order),
    basFun(BasFun)
{
    if(nameStr.find("Rl*")!=0)ABORT("illegal name "+nameStr+" must start with \"Rl*\"");
    par.push_back(Par.back());

}

void BasisFunctionRL::valDer(const UseMatrix &X, UseMatrix &Val, UseMatrix &Der, bool ZeroOutside) const{
    basFun->valDer(X,Val,Der,ZeroOutside);
    for(unsigned int k=0;k<X.size();k++){
        std::complex<double> rl=1.,rm=0.;
        for(unsigned int l=0;l<(unsigned int)(par[0]);l++){
            rm=rl;
            rl*=X(k).complex();
        }
        Der.row(k)=Val.row(k)*rm*par[0]+Der.row(k)*rl;
        Val.row(k)*=rl;
    }
}

BasisFunctionSqrt::BasisFunctionSqrt(const BasisFunction * BasFun):BasisFunction("sqrt*"+BasFun->name(),BasFun->order),
    basFun(BasFun)
{
    if(nameStr.find("sqrt*")!=0)ABORT("illegal name "+nameStr+" must start with \"sqrt*\"");
}

void BasisFunctionSqrt::valDer(const UseMatrix &X, UseMatrix &Val, UseMatrix &Der, bool ZeroOutside) const{
    basFun->valDer(X,Val,Der,ZeroOutside);
    for(unsigned int k=0;k<X.size();k++){
        if(abs(X(k).real())<1.e-12){
            Der.row(k)=0.;
            Val.row(k)=0.;
        } else {
            double sqrtX=sqrt(X(k).real());
            Der.row(k)=Val.row(k)*0.5/sqrtX+Der.row(k)*sqrtX;
            Val.row(k)*=sqrtX;
        }
        X.print("X");
        Val.row(k).print("valrhow");
    }
}
void BasisFunctionSqrt::valDer(const std::complex<double> & X, std::vector<std::complex<double> >  & Val, std::vector<std::complex<double> >  & Der, bool ZeroOutside) const{
    UseMatrix val,der;
    valDer(UseMatrix::UseMap(const_cast<complex<double>*>(&X),1,1),val,der,ZeroOutside);
    Val=vector<complex<double> >(val.data(),val.data()+val.size());
    Der=vector<complex<double> >(der.data(),der.data()+der.size());
}

void BasisFunctionOpolWeig::quadRule(const unsigned int &N, UseMatrix &QuadX, UseMatrix &QuadW) const{

    vector<double>x,w;
    oPol->quadrature(N,x,w);

    // multiply by weight function and scale into standard interval
    double scale=1.;
    if(oPol->lowerBoundary()<-DBL_MAX)
        ABORT("not implemented for "+oPol->name())
    else if (oPol->upperBoundary()<DBL_MAX/2)
        scale=1./(oPol->upperBoundary()-oPol->lowerBoundary());

    QuadX=UseMatrix(x.size(),1);
    QuadW=UseMatrix(w.size(),1);
    for( int i = 0; i<x.size(); i++){
        QuadW(i)=scale* w[i]/oPol->weight(x[i]);
        QuadX(i)=scale*(x[i]-oPol->lowerBoundary());
    }
}

void BasisFunctionOpolWeig::valDer(const std::complex<double> &X, std::vector<std::complex<double> > &Val, std::vector<std::complex<double> > &Der, bool ZeroOutside) const{

    double x=X.real();
    Val.assign(order,0);
    Der.assign(order,0);
    if(not checkStandardInterval(x,ZeroOutside))return;

    // shift and scale to oPol interval
    double shift=0,scale=1.;
    if(oPol->lowerBoundary()>-DBL_MAX/2)shift=oPol->lowerBoundary();
    if(oPol->upperBoundary()<DBL_MAX/2)scale=oPol->upperBoundary()-oPol->lowerBoundary();
    x=shift+x*scale;

    vector<double> v,d;
    oPol->valDer(order,x,v,d);
    double sqrtW=sqrt(oPol->weight(x));
    for(unsigned int j=0;j<v.size();j++){
        Val[j] =v[j]*sqrtW;
        Der[j] =d[j]*sqrtW+v[j]*(0.5*oPol->derWeight(x)/sqrtW);
        Val[j]*=sqrt(scale/oPol->normsq(j));
        Der[j]*=sqrt(scale/oPol->normsq(j))*scale;
    }
}

void BasisFunctionOpolWeig::dvrRule(std::vector<double> &Node, std::vector<double> &Weig) const{
    oPol->quadratureWithEnds(order,Node,Weig);
    for(unsigned int k=0;k<Weig.size();k++)Weig[k]/=oPol->weight(Node[k]);

    // scale into standard interval
    if(oPol->lowerBoundary()<-DBL_MAX/2){
        ABORT("not implemented for lowerBoundary =-Infty");
    } else if(oPol->upperBoundary()>DBL_MAX/2){
        if(oPol->lowerBoundary()!=0.)ABORT("only implmented for [0,infty)");
    } else {
        for(unsigned int k=0;k<Node.size();k++){
            Node[k]=(Node[k]-oPol->lowerBoundary())/(oPol->upperBoundary()-oPol->lowerBoundary());
            Weig[k]= Weig[k]                       /(oPol->upperBoundary()-oPol->lowerBoundary());
        }
    }
}
void BasisFunctionOpolWeig::dvrRule(UseMatrix &QuadX, UseMatrix &QuadW) const{
    vector<double> nodes,weights;
    dvrRule(nodes,weights);
    QuadX=UseMatrix(nodes.size(),1);
    QuadW=UseMatrix(nodes.size(),1);
    for(unsigned int k=0;k<nodes.size();k++){
        QuadX(k)=nodes[k];
        QuadW(k)=weights[k];
    }
}

BasisFunctionDvr::BasisFunctionDvr(const OrthogonalPolynomial *OPol, unsigned int Order)
    :BasisFunctionOpolWeig(OPol,Order)
{
    BasisFunctionOpolWeig::dvrRule(dvrNodes,qNorms);
    for(unsigned int k=0;k<dvrNodes.size();k++)qNorms[k]=1./sqrt(qNorms[k]*OPol->weight(dvrNodes[k]));
}

void BasisFunctionDvr::valDer(const std::complex<double> &X, std::vector<std::complex<double> > &Val, std::vector<std::complex<double> > &Der, bool ZeroOutside) const
{
    // lagrange polynomials at dvrNodes
    Val.assign(order,1);
    Der.assign(order,0);
    for(unsigned int j=0;j<dvrNodes.size();j++){
        for(unsigned int k=0;k<dvrNodes.size();k++){
            if(k!=j) {
                double r=real(X)-dvrNodes[k];
                double q=1./(dvrNodes[j]-dvrNodes[k]);
                Der[j]=(Der[j]*r+Val[j])*q;
                Val[j]=(Val[j]*r       )*q;
            }
        }
    }

    vector<complex<double> > sqrtWeig,derWeig;
    BasisFunctionOpolWeig::valDer(X,sqrtWeig,derWeig,ZeroOutside);

    // normalize functions such that norms approximated by lobatto/radau quadratures are = 1 (i.e. overlap = identity)
    for(unsigned int j=0;j<Val.size();j++){
        Der[j]=qNorms[j]*(Der[j]*sqrtWeig[0]+Val[j]*derWeig[0]);
        Val[j]=qNorms[j]* Val[j]*sqrtWeig[0];
    }

    if(Val.size()>2){
        std::swap(Val[1],Val.back());
        std::swap(Der[1],Der.back());
    }
}

void BasisFunctionGrid::valDer(const std::complex<double> &X, std::vector<std::complex<double> > &Val, std::vector<std::complex<double> > &Der, bool ZeroOutside) const{
    if(imag(X)!=0)ABORT("grid must be strictly real. X = "+tools::str(X));
    double eps=(points.back()-points[0])*1.e-14;
    Val.assign(points.size(),0);
    Der.assign(points.size(),0);
    for(unsigned int k=0;k<points.size();k++){
        if(abs(real(X)-points[k])<eps){
            Val[k]=1.;
            return;
        }
    }
    ABORT("point not in grid: X="+tools::str(real(X))+", grid="+tools::str(points,8,","));
}

BasisFunctionIndex::BasisFunctionIndex(unsigned int Order)
    :BasisFunction("useIndex",Order)
{
    for(unsigned int k=0;k<Order;k++){
        vals.push_back(vector<complex<double> >(Order,0.));
        vals[k][k]=1.;
    }
}

BasisFunctionIndex::BasisFunctionIndex(const vector<int>& IVals)
    :BasisFunction("useIndex",IVals.size())
{
    for(unsigned int k=0;k<IVals.size();k++){
        vals.push_back(vector<complex<double> >(IVals.size(),0.));
        vals[k][k]=IVals[k];
    }
}

void BasisFunctionIndex::valDer(const std::complex<double> &X, std::vector<std::complex<double> > &Val, std::vector<std::complex<double> > &Der, bool ZeroOutside) const
{
    int idx=int(round(X.real()));
    if(abs(imag(X))>1.e-12 or abs(X-double(idx))>1.e-12)ABORT("only integer index numbers allowed, is: "+tools::str(X,12));
    if(idx>=vals.size())ABORT("cannot evaluate index outside index range: "+tools::str(X)+" idx: "+tools::str(idx)+" size: "+tools::str(vals.size()));
    Val=vals[idx];
    Der.assign(Val.size(),0.);
}
