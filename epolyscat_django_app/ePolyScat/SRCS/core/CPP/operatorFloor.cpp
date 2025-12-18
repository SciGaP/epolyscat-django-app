// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "operatorFloor.h"

#include "qtEigenDense.h"
#include <Core>

#include "mpiWrapper.h"

#include "tools.h"
#include "str.h"
#include "readInput.h"
#include "algebra.h"
#include "useMatrix.h"
#include "index.h"
#include "discretization.h"
#include "operatorSingle.h"
#ifdef _USE_HACC_
#include "basisfunctionciion.h"
#endif
#include "operatorTensor.h"
#include "operatorDefinition.h"

#include "operatorZG.h"
#include "operatorZGSblock.h"
#include "operatorZD.h"
#include "operatorRG.h"
#include "operatorZGxZG.h"
#include "operatorZDxZG.h"
#include "operatorZGxZD.h"
#include "operatorZDxZD.h"
#include "operatorFloorEE.h"
#include "operatorFloorGP.h"
#include "operatorHF.h"
#include "operatorMeanEE.h"
#include "operatorFloor3d.h"
#include "operatorFloorInverse.h"
#include "operatorHartree.h"
#include "operatorFloorXC.h"
#include "printOutput.h"
#include "basisProd.h"
#include "basisNdim.h"
#include "basisMat1D.h"
#ifdef _USE_HACC_
#include "basisMatCI.h"
#endif
#include "basisMat2D.h"
#include "basisMatMulti.h"
#include "basisBesselCoulomb.h"
#include "operatorBesselCoulomb.h"
#include "basisMatDependent.h"

#include "parallel.h"
#include "constrainedView.h"

using namespace std;
#include "eigenNames.h"

double OperatorFloor::UNDEFINED=-DBL_MAX;

double OperatorFloor::bandedRatio=0.3; // only if at most 0.3 non-zeros, a matrix is considered banded
std::map <std::string,std::vector<std::complex<double> > > OperatorFloor::complexData;
std::map <std::string,std::vector<double> > OperatorFloor::realData;
std::vector<std::complex<double> > OperatorFloor::tempComplex;


std::map<std::string,int> OperatorFloor::_packCode;

int OperatorFloor::diagnoseDataSize(){

    int siz=0;
    for (auto d: complexData)siz+=d.second.size()*sizeof(std::complex<double>);
    return siz;
}

void OperatorFloor::axpy(const std::complex<double> &Alfa, const std::complex<double> *X, unsigned int SizX,
                         const std::complex<double> &Beta, std::complex<double> *Y, unsigned int SizY) const{
    DEVABORT("for some reasons, shared library needs this to be there, but should not be executed");
};

// generate a unique integer code for easy MPI'ing
int OperatorFloor::packCode(string Kind) {
    if(_packCode.size()==0){
        _packCode["DUM"]=_packCode.size();
        _packCode["Zero"]=_packCode.size();
        _packCode["ZG"]=_packCode.size();
        _packCode["ZD"]=_packCode.size();
        _packCode["ZGxZG"]=_packCode.size();
        _packCode["ZDxZG"]=_packCode.size();
        _packCode["ZGxZD"]=_packCode.size();
        _packCode["ZGSblock"]=_packCode.size();
        _packCode["FloorEE"]=_packCode.size();
        //        _packCode["TensorProduct"]=_packCode.size();
        _packCode["constrained"]=_packCode.size();
    }
    if(_packCode.count(Kind)==0)
        ABORT("no pack code defined for \""+Kind+"\", add to list");
    return _packCode[Kind];
}

UseMatrix OperatorFloor::UseMatrixTensor(std::vector<const UseMatrix*> Dat){
    if(Dat.size()!=2)ABORT("only for 2 matrices");
    UseMatrix prod(Dat[0]->rows()*Dat[1]->rows(),Dat[0]->cols()*Dat[1]->cols());
    for(unsigned int m=0,mb=0;m<Dat[0]->rows();m++,mb+=Dat[1]->rows())
        for(unsigned int n=0,nb=0;n<Dat[0]->cols();n++,nb+=Dat[1]->cols()){
            prod.block(mb,nb,Dat[1]->rows(),Dat[1]->cols())=*Dat[1];
            prod.block(mb,nb,Dat[1]->rows(),Dat[1]->cols())*=(*Dat[0])(m,n).complex();
        }
    return prod;
}

string OperatorFloor::failAbsorb;


static unsigned int cnt1=0;
void OperatorFloor::packBasic(std::vector<int> &Info, vector<complex<double> > & Buf) const{
    if(Info!=vector<int>(5,0))DEVABORT("must enter pack with vector of 5 zeros");
    Info[0]=packCode(kind());
    Info[1]=_rows;
    Info[2]=_cols;
    Info[3]=Buf.size();
    if(_rows*_cols!=0 and Buf.size()==0)ABORT("cannot pack empty floor");
    // append iWeights to end of buffer (if any)
    if(iWeights!=0)Buf.insert(Buf.end(),iWeights->begin(),iWeights->end());
    Info[4]=Buf.size();
}

void OperatorFloor::unpackBasic(const std::vector<int> &Info, const std::vector<std::complex<double> > &Buf){
    if(Info[0]!=packCode(kind()))ABORT(Str("mismatch Info[0] and packCode():")+Info[0]+" vs. "+packCode(kind())+"("+kind()+")");
    _rows=Info[1];
    _cols=Info[2];
    iWeights=0;
    if(Info[3]<Info[4])iWeights.reset(new std::vector<complex<double>>(Buf.begin()+Info[3],Buf.end()));
}

string OperatorFloor::tensorType(string Kinds){
    vector<string> facs(tools::splitString(Kinds,' '));
    if(facs.size()!=2)return "UNKNOWN";
    if((facs[0].length()==2 and facs[1].length()==2))return "ZGxZG";
    if((facs[0].length()==2 and facs[1].length()==3) and facs[1][2]=='d')return "ZGxZD";
    if((facs[1].length()==2 and facs[0].length()==3) and facs[0][2]=='d')return "ZDxZG";
    if((facs[0].length()==3 and facs[0][2]=='d' and facs[1].length()==3) and facs[1][2]=='d')return "ZDxZD";
    return "ZGxZG";
}

void OperatorFloor::apply(std::complex<double> Alfa, const std::complex<double>* X, unsigned int SizX,
                          const std::complex<double> & Beta, std::complex<double>* Y, unsigned int SizY) const{
    if(X==0){
        // input may be sparse
        if(Beta!=1.)for(size_t k=0;k<SizY;k++)*Y*=Beta;
        return;
    }
    if(Y==0)DEVABORT("lhs not assigned");

    if(timeDepFac!=0)Alfa*=*timeDepFac;
    axpy(Alfa,X,SizX,Beta,Y,SizY);
}

TIMER(setupFloor,)
//OperatorFloor* OperatorFloor::factory(const OperatorTensor &OT){
//    if(OT.mats.size()==0)
//        return new OperatorZero(OT.iIndex->sizeStored(),OT.jIndex->sizeStored());
//    return factory(OT.mats,OT.definition+tools::str(OT.iIndex->label())+"|"+tools::str(OT.jIndex->label()));
//}
OperatorFloor* OperatorFloor::factory(const std::string& TermOper,
                                      const Index * IIndex, const Index *JIndex, std::complex<double> Multiplier)
{
#ifdef _USE_HACC_
    BasisMatCI matCI;
#endif
    BasisMat1D opMat;
    BasisMat2D o2Mat;
    BasisMatDependent depMat;
    BasisMatMulti multMat;
    vector<const Eigen::MatrixXcd*>m;

    if(TermOper.find("HartreeFock")!=string::npos)return new OperatorFloorHF(tools::stringInBetween(TermOper,"[","]"),IIndex,JIndex,Multiplier);
    if(TermOper.find("Hartree")!=string::npos)return new OperatorHartree(tools::stringInBetween(TermOper,"[","]"),IIndex,JIndex,Multiplier);
    if(TermOper.find("XC")!=string::npos)     return new OperatorFloorXC(tools::stringInBetween(TermOper,"[","]"),IIndex,JIndex,Multiplier);
    if(TermOper.find("MeanEE")!=string::npos)return new OperatorMeanEE(tools::stringInBetween(TermOper,"[","]"),IIndex,JIndex,Multiplier);
    if(TermOper.find("GrossPitaevskii")!=string::npos)return new OperatorFloorGP(tools::stringInBetween(TermOper,"[","]"),IIndex,JIndex,Multiplier);


    // this should go into a more abstract factory
    if(IIndex->basis()->name().find("besselCoulomb")==0){
        OperatorBesselCoulomb *opBC = new OperatorBesselCoulomb(TermOper,dynamic_cast<const BasisBesselCoulomb*>(IIndex->basis()->integrable()),dynamic_cast<const BasisBesselCoulomb*>(JIndex->basis()->integrable()),Multiplier);
        return dynamic_cast<OperatorFloor*>(opBC);
    }
    const BasisMatMatrix* bmm;
    if(not (opMat=BasisMat1D(TermOper,IIndex,JIndex)).isEmpty())m=opMat.mats();
    if(not (depMat=BasisMatDependent(TermOper,IIndex,JIndex)).isEmpty())m=depMat.mats();
    else if(not (o2Mat=BasisMat2D(TermOper,IIndex,JIndex)).isEmpty())m=o2Mat.mats();
    else if(not (multMat=BasisMatMulti(TermOper,IIndex,JIndex)).isEmpty())m=multMat.mats();
    else if((bmm=BasisMatMatrix::factory(TermOper,IIndex->basis(),JIndex->basis())))m=bmm->mats();
#ifdef _USE_HACC_
    else if(not (matCI=BasisMatCI(TermOper,IIndex,JIndex)).isEmpty())m=matCI.mats();
#endif
    if(m.size()>0){
        Eigen::MatrixXcd m0;
        if(Multiplier!=1.){
            m0=*m[0]*Multiplier;
            m[0]=&m0;
        }
        string hash=TermOper+IIndex->hash()+JIndex->hash()+"*"+tools::str(Multiplier);
        return OperatorFloor::factory(m,hash);
    }
    else {
        if(IIndex->root()->isHybrid() || IIndex->root()->axisName() == "Ion" || IIndex->root()->axisName() == "Neutral" ||
                JIndex->root()->isHybrid() || JIndex->root()->axisName() == "Ion" || JIndex->root()->axisName() == "Neutral" ||
                IIndex->axisName() == "Ion" || JIndex->axisName() == "Ion" || IIndex->axisName() == "Neutral" || JIndex->axisName() == "Neutral"){
//            return factory(OperatorTensor(TermOper,0,0,IIndex,JIndex,Multiplier));
            DEVABORT("disabled");
        }
        ABORT(Str("failed to construct operator (part)"," ")+TermOper+"for"+IIndex->strNode()+"|"+JIndex->strNode()+", check syntax"
              +"\nlist of standard functions in Algebra:\n"+Algebra::listStandard());
    }
    return 0;
}
OperatorFloor* OperatorFloor::copyFactory(const OperatorFloor *Other){
    OperatorFloor*oFloor=const_cast<OperatorFloor*>(Other);
    if     (dynamic_cast<const OperatorZG*>(Other)!=0)oFloor=new OperatorZG(*dynamic_cast<const OperatorZG*>(Other));
    else if(dynamic_cast<const OperatorZD*>(Other)!=0)oFloor=new OperatorZD(*dynamic_cast<const OperatorZD*>(Other));
    else if(dynamic_cast<const OperatorZGSblock*>(Other)!=0)oFloor=new OperatorZGSblock(*dynamic_cast<const OperatorZGSblock*>(Other));
    else if(dynamic_cast<const OperatorZGxZG*>(Other)!=0)oFloor=new OperatorZGxZG(*dynamic_cast<const OperatorZGxZG*>(Other));
    else if(dynamic_cast<const OperatorZDxZG*>(Other)!=0)oFloor=new OperatorZDxZG(*dynamic_cast<const OperatorZDxZG*>(Other));
    else if(dynamic_cast<const OperatorZGxZD*>(Other)!=0)oFloor=new OperatorZGxZD(*dynamic_cast<const OperatorZGxZD*>(Other));
    else if(dynamic_cast<const OperatorZero*>(Other)!=0)oFloor=new OperatorZero(*dynamic_cast<const OperatorZero*>(Other));
    else if(dynamic_cast<const OperatorDUM*>(Other)!=0)oFloor=new OperatorDUM(*dynamic_cast<const OperatorDUM*>(Other));
    else
        PrintOutput::DEVwarning("in OperatorTree::nodeCopy -- "
                                +Other->kind()
                                +" not copied, only pointed to -- may break, if one of the terms is deleted",1);
    return oFloor;
}

//TIMER(factory1,);
//TIMER(factory2,);
//TIMER(factory3,);
OperatorFloor* OperatorFloor::factory(const std::vector<const Eigen::MatrixXcd*> &PMats, string Hash){
    // not very elegant
    //    STARTDEBUG(factory1);
    vector<const UseMatrix*>pmats;
    for(size_t k=0;k<PMats.size();k++)pmats.push_back(new UseMatrix(*PMats[k]));
    //    STOPDEBUG(factory1);
    //    STARTDEBUG(factory2);
    OperatorFloor* of=factory(pmats,Hash);
    //    STOPDEBUG(factory2);
    //    STARTDEBUG(factory3);
    for(size_t k=0;k<PMats.size();k++)delete pmats[k];
    //    STOPDEBUG(factory3);
    return of;
}

OperatorFloor* OperatorFloor::factoryInverse(const Index* Idx, unsigned int SubD, unsigned int SuperD, bool BandOvr){
    return new OperatorFloorInverse(Idx, SubD, SuperD, BandOvr);
}

//TIMER(factory4,)
//TIMER(factory5,)
//TIMER(factory6,)
//TIMER(floor1,)
//TIMER(floor2,)
//TIMER(floor3,)
//TIMER(floor4,)
//TIMER(floor5,)
//TIMER(floor6,)
//TIMER(floor7,)
OperatorFloor* OperatorFloor::factory(const std::vector<const UseMatrix *> &PMats, string Hash){
    // analyze the matrices
    //    STARTDEBUG(factory4);
    if(PMats.size()==0)ABORT("factory requires at least one matrix");

    string kinds;
    for(unsigned int k=0;k<PMats.size();k++){
        kinds+=Eigen::Map<Eigen::MatrixXcd>(PMats[k]->data(),PMats[k]->rows(),PMats[k]->cols()).isDiagonal()?"dia ":"oth ";
        //        kinds+=PMats[k]->type(bandedRatio)+" ";
    }

    // production form of the operator
    OperatorFloor* ret;
    //    STOPDEBUG(factory4);
    //    STARTDEBUG(factory5);
    if(        kinds=="oth "
               or kinds== "rs "
               or kinds=="rgb "
               or kinds=="rsb "
               or kinds=="rg "
               or kinds=="ig "
               or kinds=="cs "
               or kinds=="ch "
               or kinds=="ih "
               or kinds=="cg "
               or kinds=="is "
               or kinds=="rh "
               or kinds=="csb " // Added temporarily
               or kinds=="cgb "
               or kinds=="igb "
               ){
        // check for shape cphase*real
        cnt1++;

        if(not PMats[0]->isFull()){
            //            STARTDEBUG(floor5);
            UseMatrix mat(*PMats[0]);
            ret=new OperatorZG(&mat.expand(),Hash);
            //            STOPDEBUG(floor5);
        } else {
            //            STARTDEBUG(floor6);
            ret=new OperatorZG(PMats[0],Hash);
            //            STOPDEBUG(floor6);
        }
    }

    else if(   kinds=="dia "
               or kinds=="csd "
               or kinds=="rsd "
               or kinds=="isd "
               ){
        //        STARTDEBUG(floor7);
        if(PMats[0]->isFull())
            ret=new OperatorZD(Eigen::Map<Eigen::MatrixXcd>(PMats[0]->data(),PMats[0]->rows(),PMats[0]->cols()),Hash);
        else
            ret=new OperatorZD(PMats[0],Hash);
        //        STOPDEBUG(floor7);
    }

    else if(kinds=="dia dia " or tensorType(kinds)=="ZDxZD"){
        //        STARTDEBUG(floor1);
        ret=new OperatorZDxZD(PMats,Hash);
        //        STOPDEBUG(floor1);
    }

    else if(kinds=="dia oth " or tensorType(kinds)=="ZDxZG"){
        //        STARTDEBUG(floor2);
        ret=new OperatorZDxZG(PMats,Hash);
        //        STOPDEBUG(floor2);
    }

    else if(kinds=="oth dia " or tensorType(kinds)=="ZGxZD"){
        //        STARTDEBUG(floor3);
        ret=new OperatorZGxZD(PMats,Hash);
        //        STOPDEBUG(floor3);
    }

    else if("oth oth " or tensorType(kinds)=="ZGxZG"){
        //        STARTDEBUG(floor4);
        ret=new OperatorZGxZG(PMats,Hash);
        //        STOPDEBUG(floor4);
    }
    else {
        //        STARTDEBUG(factory6);
        for(unsigned int k=0;k<PMats.size();k++)
            kinds+=PMats[k]->type(bandedRatio)+" ";
        if(kinds.length()<5 and kinds.substr(2,1)=="s"){
            // for now compress sparse to block to the extend possible
            ret=new OperatorZGSblock(PMats[0],Hash);
        } else {
            PMats[0]->print("floor matrix: "+Hash,1);
            PMats[0]->print("floor matrix: "+Hash,0);
            ABORT("unknown type \""+kinds+"\", diagnosed as "+tensorType(kinds));
        }
        //        STOPDEBUG(factory6);
    }

    //    STOPDEBUG(factory5);
    if(ret->oNorm==DBL_MAX)ret->setNorm(); // !!! may be slow... !!!
    return ret;
}

OperatorFloor* OperatorFloor::specialFactory(const std::string Name, const std::string Def,
                                             const Index *IIndex, const Index *JIndex, std::complex<double> Multiplier)
{
    // either directly compute  OperatorFloor
    //     or set up matrices for standard OperatorFloor::factory
    OperatorFloor* of=0;
    vector<UseMatrix> mat;

    if( IIndex->basis()->ndim()!=0 or JIndex->basis()->ndim()!=0 ){
        mat.resize(1);
        BasisNdim::matrix(Def,IIndex,JIndex,mat[0],1.);
    }
    else if(Def=="[[allOnes]]"){
        mat.push_back(UseMatrix::Constant(IIndex->sizeCompute(),JIndex->sizeCompute(),1.));
    }
    else if(Def=="[[Id]]"){
        mat.push_back(UseMatrix::Identity(IIndex->sizeCompute(),JIndex->sizeCompute()));
    }
#ifdef _USE_HACC_
    else if(Def.find("haCC")!=string::npos){
        DEVABORT("disabled");
//        for(const Index* s=IIndex; not s->isRoot(); s=s->parent()){
//            if(s->basis()->name().find("CIion")!=string::npos) {
//                mat.push_back(const_cast<BasisFunctionCIion* >(dynamic_cast<const BasisFunctionCIion*>(s->basisSet()->PointerToFunctions()))->matrixWithNeutrals(Def,IIndex,JIndex));
//                break;
//            }
//        }
//        if(mat.size()==0){
//            for(const Index* s=JIndex; not s->isRoot(); s=s->parent()){
//                if(s->basis()->name().find("CIion")!=string::npos) {
//                    mat.push_back(const_cast<BasisFunctionCIion* >(dynamic_cast<const BasisFunctionCIion*>(s->basisSet()->PointerToFunctions()))->matrixWithNeutrals(Def,JIndex,IIndex).adjoint());
//                    break;
//                }
//            }
//        }
//        if(mat[0].isZero(1.e-12))mat.clear();
    }
#endif

    else if(Def=="[[eeInt6DHelium]]"){
        if(Multiplier!=1.)ABORT("cannot have any factors for "+Def);
        of=new OperatorFloorEE(Name,Def,IIndex,JIndex);
    }
    else if(Def.find("[[HartreeFock")==0){
        if(Multiplier!=1.)ABORT("cannot have any factors for "+Def);
        of=new OperatorFloorHF("",IIndex,JIndex,Multiplier);
    }
    else if(Def=="[[Pot3d]]"){
        mat.clear(); // no matrices needed
        if(Multiplier!=1.)ABORT("cannot have any factors for "+Def);
        of=new OperatorFloor3d(IIndex,JIndex);
    }
#ifdef _USE_HACC_
    else if(IIndex->basis()->ci() or JIndex->basis()->ci()){
        mat.push_back(UseMatrix::Constant(IIndex->sizeCompute(),JIndex->sizeCompute(),0.));
        PrintOutput::DEVwarning("use zero for matrix elements <"+IIndex->axisName()+"|"+Def+"|"+JIndex->axisName()+">",1);
    }
#endif
    else
        ABORT("special operator "+Def+" undefined on "+IIndex->strNode()+"|"+JIndex->strNode());

    if(mat.size()>0){
        mat[0]*=Multiplier; // multiply first factor matrix by multiplier
        vector<const UseMatrix*> pMat;
        for(size_t k=0;k<mat.size();k++)pMat.push_back(&mat[k]);
        of=OperatorFloor::factory(pMat,Def+IIndex->hash()+JIndex->hash());
    }
    //    if(of!=0)of->iWeights=new vector<complex<double> >(0);
    return of;
}

void OperatorFloor::replace(OperatorFloor *&Floor, const std::vector<int> &Info, const std::vector<std::complex<double> > &Buf){
    OperatorFloor * rep;
    if(Info.size()==0)rep=new OperatorDUM();
    else rep=OperatorFloor::unpackFactory(Info,Buf);
    rep->oNorm=Floor->oNorm;
    rep->_cost=Floor->_cost;
    rep->timeDepFac=Floor->timeDepFac;
    delete Floor;
    Floor=rep;
}

// cannot construct through macro, as it would be switched off by -D_TIMER_OFF_
static Timer TIMERcost("cost","INTERNAL","operatorFloor.cpp");
double OperatorFloor::applicationCost(bool Bcast) const {
    if(_cost==UNDEFINED){
        vector<complex<double> > x(_cols),y(_rows);
        for(unsigned int k=0;k<x.size();k++)x[k]=complex<double>(exp(complex<double>(0,k)));
        TIMERcost.start();
        TIMERcost.stopTimer();
        double startSecs=TIMERcost.secs();
        TIMERcost.start();
        // do a substantial number of applications
        for(unsigned int k=0;k<1;k++){
            axpy(1.,x.data(),x.size(),0.,y.data(),y.size());
            axpy(1.,x.data(),x.size(),0.,y.data(),y.size());
            axpy(1.,x.data(),x.size(),0.,y.data(),y.size());
            axpy(1.,x.data(),x.size(),0.,y.data(),y.size());
            axpy(1.,x.data(),x.size(),0.,y.data(),y.size());
            axpy(1.,x.data(),x.size(),0.,y.data(),y.size());
            axpy(1.,x.data(),x.size(),0.,y.data(),y.size());
            axpy(1.,x.data(),x.size(),0.,y.data(),y.size());
            axpy(1.,x.data(),x.size(),0.,y.data(),y.size());
            axpy(1.,x.data(),x.size(),0.,y.data(),y.size());
        }
        TIMERcost.stopTimer();
        double cst=TIMERcost.secs()-startSecs;
        if(Bcast)MPIwrapper::Bcast(&cst,1,MPIwrapper::master()); // make sure all processes have the same _cost
        const_cast<OperatorFloor*>(this)->_cost=1.e5*cst;
    }
    return _cost;
}

OperatorFloor* OperatorFloor::unpackFactory(const std::vector<int> &Info, const std::vector<std::complex<double> > &Buf){
    OperatorFloor* of(0);
    if     (Info[0]==packCode("DUM"))of=         new OperatorDUM();
    else if(Info[0]==packCode("ZG"))of=          new OperatorZG(Info,Buf);
    else if(Info[0]==packCode("ZD"))of=          new OperatorZD(Info,Buf);
    else if(Info[0]==packCode("ZGxZG"))of=       new OperatorZGxZG(Info,Buf);
    else if(Info[0]==packCode("ZDxZG"))of=       new OperatorZDxZG(Info,Buf);
    else if(Info[0]==packCode("ZGxZD"))of=       new OperatorZGxZD(Info,Buf);
    else if(Info[0]==packCode("ZGSblock"))of=    new OperatorZGSblock(Info,Buf);
    else if(Info[0]==packCode("FloorEE"))of=     new OperatorFloorEE(Info,Buf);
    else if(Info[0]==packCode("constrained"))of= new ConstrainedZero::Floor(Info,Buf);
    else
        ABORT("not defined for packCode()="+tools::str(Info[0]));
    of->setNorm();
    return of;
}


string OperatorFloor::hashString(unsigned int Rows, unsigned int Cols){
    return tools::str(Rows)+"x"+tools::str(Cols);
}

string OperatorFloor::str(int Digits) const{
    string s=strInfo();
    if(Digits>=0){
        UseMatrix mat;
        matrix(mat);
        return s=mat.str(s,Digits);
    }
    return s;
}

void OperatorFloor::write(ofstream &File) const{
    std::vector<int>info(5,0);
    std::vector<std::complex<double>> buf;
    pack(info,buf);
    tools::write(File,info);
    for(auto a: buf)if(std::isnan(a.real()) or std::isnan(a.imag()))ABORT("corrupted write");
    tools::write(File,buf.size());
    tools::write(File,buf);
}
OperatorFloor* OperatorFloor::readFactory(std::ifstream &File){
    std::vector<int> info(5);
    tools::read(File,info);
    size_t siz;
    tools::read(File,siz);
    std::vector<std::complex<double>> buf(siz);
    tools::read(File,buf);
    for(auto a: buf)if(std::isnan(a.real()) or std::isnan(a.imag()))ABORT("corrupted read");
    return unpackFactory(info,buf);
}


string OperatorFloor::strInfo() const {
    std::string s(kind());
    s+=" "+tools::str(_rows,3)+"x"+tools::str(_cols,3);
    if(_cost!=UNDEFINED)s+=" "+tools::str(_cost,3);
    else s+="  * ";
    return s;
}

bool OperatorFloor::isDiagonal() const {
    return     dynamic_cast<const OperatorDUM*>(this)
            or dynamic_cast<const OperatorZD*>(this);
}

bool OperatorFloor::isAbsorbable() const {
    if(dynamic_cast<const OperatorZG*>(this)!=0)return true;
    if(dynamic_cast<const OperatorZD*>(this)!=0)return true;
    if(dynamic_cast<const OperatorZGSblock*>(this)!=0)return true;
    if(dynamic_cast<const OperatorDUM*>(this)!=0)return true;
    if(dynamic_cast<const OperatorZero*>(this)!=0)return true;
    if(dynamic_cast<const OperatorTensorProduct*>(this)!=0)return false;
    if(dynamic_cast<const OperatorFloorEE*>(this)!=0)return false;
    if(dynamic_cast<const OperatorFloorGP*>(this)!=0)return false;
    if(dynamic_cast<const OperatorMeanEE*>(this)!=0)return false;
    if(dynamic_cast<const OperatorBesselCoulomb*>(this)!=0)return true;
    ABORT("please, define isAbsorbable() for OperatorFloor kind="+kind());
    return true; // shut up the compiler
}

//TIMER(absorb1,)
//TIMER(absorb2,)
//TIMER(absorb3,)
//TIMER(absorb4,)
//TIMER(absorb5,)
//TIMER(absorb,)
bool OperatorFloor::absorb(OperatorFloor *&Into, OperatorFloor *&Other, string Def){

    //    STARTDEBUG(absorb);
    // to maintain operator tree structure, this must be done in parallel
    // admits for the cases: a single proces or all processes hold the data

    int hinder=1;
    {
        if(Into==Other)goto AllHinder;
        if(Into->kind()=="DUM" and Other->kind()=="DUM")goto AllHinder;
        if(Into->kind()=="DUM" or Other->kind()=="DUM")goto AllHinder;//ABORT("must not have floor and dummy on same iIndex,jIndex block");
        if(Into->_rows!=Other->_rows)goto AllHinder;
        if(Into->_cols!=Other->_cols)goto AllHinder;
        if(Into->timeDepFac!=Other->timeDepFac)goto AllHinder;
        if(((Into->iWeights==0) != (Other->iWeights==0)))goto AllHinder;
        if(Into->iWeights!=0 and *Into->iWeights!=*Other->iWeights)goto AllHinder;
        if(not Into->isAbsorbable())goto AllHinder;
        if(not Other->isAbsorbable())goto AllHinder;

        // form matrices and get new floor
        //        STARTDEBUG(absorb2);
        Eigen::MatrixXcd mat,other;
        mat=Into->matrix();
        other=Other->matrix();
        mat+=other;
        //        STOPDEBUG(absorb2);
        //        STARTDEBUG(absorb3);
        EigenTools::purge(mat,1.e-13);
        //        STOPDEBUG(absorb3);
        //        STARTDEBUG(absorb4);
        //        mat.purge(1.e-13,1.e-13);
        OperatorFloor* oAbs=factory(std::vector<const Eigen::MatrixXcd*>(1,&mat),Def);
        //        STOPDEBUG(absorb4);

        // test for successful new floor
        if(oAbs==0)goto AllHinder;

        // make new host known to floor host list
        delete Into;Into=oAbs;
        Into->timeDepFac=Other->timeDepFac;

        // all passed
        hinder=0;
    }

AllHinder:

    //    STARTDEBUG(absorb1);
    MPIwrapper::AllreduceSUM(&hinder,1);
    //    STOPDEBUG(absorb1);

    // all or a single process did absorb
    if(hinder==0 or hinder==MPIwrapper::Size()-1){
        delete Other;Other=0;
        Into->forceNorm(); // reset norm
        Into->forceCost(); // reset cost
    }
    else if(hinder!=MPIwrapper::Size())
        ABORT(Str("non-absorbable floors on")+hinder+"processes out of"+MPIwrapper::Size()+"(illegal)");

    //    STOPDEBUG(absorb);
    return Other==0;
}

void OperatorFloor::applyLeftOverlap(std::complex<double> *Data){
    PrintOutput::DEVwarning("applyLeftOverlap: not meant to be used");
    if(iWeights!=0)
        for(unsigned int k=0;k<iWeights->size();k++,Data++)*Data*=iWeights->data()[k];
}


void OperatorFloor::addWeights(const Index* IIndex, const std::vector<std::complex<double> > &IWeights){
    PrintOutput::DEVwarning("addWeights: not meant to be used");
    if(IWeights.size()==0)iWeights=0;
    else
        iWeights=addComplex(tools::str(IIndex),IWeights);
}

void OperatorFloor::scale(std::complex<double> Beta, std::vector<std::complex<double> > &Y){
    if(Beta!=1.){
        if(Beta==0.)for(unsigned int k=0;k<Y.size();k++)Y[k]=0.;
        else        for(unsigned int k=0;k<Y.size();k++)Y[k]*=Beta;
    }
}

void OperatorFloor::scale(std::complex<double> Beta,std::complex<double>*Y, unsigned int SizY){
    if(Beta!=1.){
        if(Beta==0.)for(unsigned int k=0;k<SizY;k++)Y[k]=0.;
        else        for(unsigned int k=0;k<SizY;k++)Y[k]*=Beta;
    }
}

OperatorFloor::~OperatorFloor(){}

std::shared_ptr<std::vector<std::complex<double>>> OperatorFloor::addComplex(const std::string &Hash, const std::vector<std::complex<double> > &Dat)
{
    std::shared_ptr<std::vector<std::complex<double>>>res(new std::vector<std::complex<double>>(Dat.begin(),Dat.end()));
    return res;
}

vector<double> * OperatorFloor::addReal(const string &Hash, const std::vector<double> &Dat){
    vector<double> * added;
    string hash=Hash+"|"+tools::str(Dat.size())+":";
    added=&realData[hash];
    unsigned int c=0;
    while(added->size()!=0 and (added->size()!=Dat.size() or *added!=Dat)){
        c++;
        added=&realData[Hash+tools::str(c)];
    }
    if(added->size()==0)*added=Dat;
    return added;
}

Eigen::MatrixXcd OperatorFloor::matrixFactor(int D) const {
    if(dynamic_cast<const OperatorTensorProduct*>(this)!=0)
        ABORT("matrixFactors not implemented for "+kind());
    return matrix();
}

Eigen::MatrixXcd OperatorFloor::matrix() const {
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(_rows,_cols);
    vector<complex<double> > x(_cols,0.),y(_rows);
    for(unsigned int j=0;j<_cols;j++){
        x[j]=1.;
        axpy(1.,x.data(),x.size(),0.,y.data(),y.size());
        for (unsigned int i=0;i<_rows;i++) {mat(i,j)=y[i];}
        x[j]=0.;
    }
    if(iWeights!=0){
        if(iWeights->size())
            for(unsigned int k=0;k<iWeights->size();k++){
#ifdef _DEVELOP_
                PrintOutput::DEVwarning("iWeights not meant to be used");
#endif
                mat.row(k)*=iWeights->at(k);
//#else
//                if(iWeights->at(k)!=1.)DEVABORT("not meant to be used");
//#endif
            }
    }
    return mat;
}

void OperatorFloor::matrix(UseMatrix & mat) const{
    // default matrix setup through axpy operation
    mat=UseMatrix(_rows,_cols);
    vector<complex<double> > x(_cols,0.),y(_rows);
    for(unsigned int j=0;j<_cols;j++){
        x[j]=1.;
        axpy(1.,x.data(),x.size(),0.,y.data(),y.size());
        for (unsigned int i=0;i<_rows;i++) {mat(i,j)=y[i];}
        x[j]=0.;
    }
}

void OperatorFloor::setNorm() const {
    UseMatrix mat;
    matrix(mat);
    const_cast<OperatorFloor*>(this)->oNorm=mat.maxAbsVal();
}
void OperatorDUM::axpy(const std::complex<double> & Alfa, const std::complex<double>*X, unsigned int SizX,
                       const std::complex<double> & Beta,       std::complex<double>*Y, unsigned int SizY) const
{DEVABORT("cannot apply dummy");}

void OperatorZero::axpy(const std::complex<double> & Alfa, const std::complex<double>*X, unsigned int SizX,
                        const std::complex<double> & Beta,       std::complex<double>*Y, unsigned int SizY) const
{for(size_t k=0;k<SizY;k++)Y[k]*=Beta;}
