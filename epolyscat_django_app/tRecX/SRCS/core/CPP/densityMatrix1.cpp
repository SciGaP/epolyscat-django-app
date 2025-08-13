// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "densityMatrix1.h"

#include "index.h"
#include "basisIntegrable.h"

#include "coefficients.h"
#include "basisOrbitalNumerical.h"
#include "eigenTools.h"


#include "timer.h"

static void extractValues(const BasisOrbitalNumerical & OrbJ, const BasisIntegrable & BasA, std::map<int,std::map<int,Eigen::MatrixXcd> > & MlJa){

    // sort into matrices
    MlJa.clear();
    std::string firstAxis=OrbJ.orbital(0)->idx()->axisName();
    if(firstAxis=="Phi"){
        for(size_t j=0;j<OrbJ.size();j++){
            const Coefficients * oj=OrbJ.orbital(j);
            for(size_t km=0;km<oj->childSize();km++){
                int m=oj->idx()->basis()->physical(km);
                const Coefficients* orbjm=oj->child(km);
                for(size_t kl=0;kl<orbjm->childSize();kl++){
                    int l=orbjm->idx()->basis()->physical(kl);
                    for (size_t n=0;n<orbjm->child(kl)->childSize();n++)
                    {
                        if(*orbjm->child(kl)->child(n)->idx()->basis()->integrable()==BasA){
                            if(MlJa[m][l].size()==0)MlJa[m][l]=Eigen::MatrixXcd::Zero(OrbJ.size(),BasA.size());
                            for(size_t a=0;a<BasA.size();a++)
                                MlJa[m][l](j,a)=orbjm->child(kl)->child(n)->data()[a];
                        }
                    }
                }
            }
        }
    }
    else if(firstAxis=="Rn"){
        size_t nRad=0;
        for(;nRad<OrbJ.orbital(0)->childSize() ;nRad++)
            if(*OrbJ.orbital(0)->child(nRad)->firstLeaf()->idx()->basis()->integrable()==BasA)break;
        if(nRad==OrbJ.orbital(0)->childSize())DEVABORT("radial basis not found");

        for(size_t j=0;j<OrbJ.size();j++){
            const Coefficients * oj=OrbJ.orbital(j)->child(nRad);
            for(size_t km=0;km<oj->childSize();km++){
                int m=oj->idx()->basis()->physical(km);
                const Coefficients* orbjm=oj->child(km);
                for(size_t kl=0;kl<orbjm->childSize();kl++){
                    int l=orbjm->idx()->basis()->physical(kl);
                    if(MlJa[m][l].size()==0)MlJa[m][l]=Eigen::MatrixXcd::Zero(OrbJ.size(),BasA.size());
                    for(size_t a=0;a<BasA.size();a++)
                        MlJa[m][l](j,a)=orbjm->child(kl)->data()[a];
                }
            }
        }
    }
    else
        DEVABORT("only for hiearchies Rn.Phi.Eta or Phi.Eta.Rn, is: "+OrbJ.orbital(0)->idx()->hierarchy());

}

DensityMatrix1::DensityMatrix1(const std::vector<std::vector<Eigen::MatrixXcd> > &RhoIJ)
    :_s(-1),_r(-1){
    // SVD RhoIJ
    _cacb=Eigen::MatrixXi::Constant(RhoIJ.size(),RhoIJ[0].size(),-1);
    if(_cacb.rows()!=_cacb.cols())DEVABORT(Sstr+"RhoIJ blocks must be square, found"+_cacb.rows()+"x"+_cacb.cols());

    for(size_t i=0;i<RhoIJ.size();i++)
        for(size_t j=0;j<RhoIJ[i].size();j++){
            if(RhoIJ[i][j].size()>0 and not RhoIJ[i][j].isZero(1.e-12)){
                Eigen::JacobiSVD<Eigen::MatrixXcd> svd(RhoIJ[i][j],Eigen::ComputeThinU | Eigen::ComputeThinV);
                _cacb(i,j)=_iSvd.size();
                _iSvd.push_back(svd.matrixU()*svd.singularValues().asDiagonal());
                _jSvd.push_back(svd.matrixV());
            }
        }
}

void DensityMatrix1::patch(const BasisIntegrable &BasA, const BasisIntegrable &BasB,
                           const BasisOrbitalNumerical & OrbI, const BasisOrbitalNumerical & OrbJ){
    // store all orbitals on IdxA/IdxB for fast density computation
    _s=-1;
    _r=-1;
    reset();
    extractValues(OrbI,BasB,_mlIs);
    extractValues(OrbJ,BasA,_mlJr);

    // m-ranges
    _imMin= INT_MAX;
    _imMax=-INT_MAX;
    for (auto p: _mlIs){
        _imMin=std::min(_imMin,p.first);
        _imMax=std::max(_imMax,p.first);
    }
    _jmMin= INT_MAX;
    _jmMax=-INT_MAX;
    for (auto p: _mlJr){
        _jmMin=std::min(_jmMin,p.first);
        _jmMax=std::max(_imMax,p.first);
    }
}
DensityMatrix1::M::M(DensityMatrix1 & Dens,int Mi, int Mj):_dens(&Dens),_lIs(&Dens._mlIs[Mi]),_lJr(&Dens._mlJr[Mj]){

    _ilMin= INT_MAX;
    _ilMax=-INT_MAX;
    for (auto p: *_lIs){
        _listLi.push_back(p.first);
        _ilMin=std::min(_ilMin,p.first);
        _ilMax=std::max(_ilMax,p.first);
    }
    _jlMin= INT_MAX;
    _jlMax=-INT_MAX;
    for (auto p: *_lJr){
        _listLj.push_back(p.first);
        _jlMin=std::min(_jlMin,p.first);
        _jlMax=std::max(_jlMax,p.first);
    }
}

DensityMatrix1::M* DensityMatrix1::operator ()(int Mi,int Mj){
    if(_s<0 or _r<0)DEVABORT("Must set radial points DensitMatrix1::set(a,b) before setting Mi,Mj");
    if(not _m.count(Mi) or not _m[Mi].count(Mj))_m[Mi][Mj].reset(new DensityMatrix1::M(*this,Mi,Mj));
    return _m[Mi][Mj].get();
}

TIMER(densLGet,)
TIMER(densLRet,)
DensityMatrix1::L* DensityMatrix1::M::operator ()(int Li){
    // somewhat time-critical - optimize
    auto _pi=_l.find(Li);
    if(_l.find(Li)==_l.end()){
        DensityMatrix1::L * d=new DensityMatrix1::L(*this,Li);
        _l[Li].reset(d);
        return d;
    }
    else
        return _pi->second.get();
}
DensityMatrix1::LL* DensityMatrix1::L::operator ()(int Lj){
    if(not _ll[Lj])_ll[Lj].reset(new DensityMatrix1::LL(*this,Lj));
    return _ll[Lj].get();
}

TIMER(density,)
DensityMatrix1::LL::LL(const L &DensL, int Lj)
    :_densM(DensL._densM),_Is(DensL._Is)
{
    auto p=DensL._densM->_lJr->find(Lj);
    if(p==DensL._densM->_lJr->end())return; // combination with l1 does not occur
    _Jr=&(p->second);
    if(densC.size()==0){
        START(density);
        for(size_t cc=0;cc<_densM->_dens->_iSvd.size();cc++){
            Eigen::MatrixXcd mm
                    =(_Jr->col(_densM->_dens->_r).transpose()*(_densM->_dens->_jSvd[cc]))*(_Is->col(_densM->_dens->_s).transpose()*(_densM->_dens->_iSvd[cc])).adjoint();
                    
            densC.push_back(mm(0,0));
        }
        STOP(density);
    }
}
DensityMatrix1::L::L(M &DensM, int Li)
    :_densM(&DensM),_Is(&DensM._lIs->at(Li)),_lJr(DensM._lJr){
    for(auto p: *_lJr)_listLj.push_back(p.first);
    _ll.resize(40);

}

void DensityMatrix1::reset(){
    for (auto &ldi: _m)
        for (auto &ldij: ldi.second)
            ldij.second->_l.clear();
    _m.clear();
}

