// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "basisNdim.h"
#include "readInput.h"

#include "operatorDefinition.h"
#include "useMatrix.h"
#include "index.h"
#include "tools.h"
#include "coordinateTrans.h"
#include "operatorNdim.h"
#include "basisProd.h"
#include "basisIntegrable.h"
#include "basisPolarOff.h"
#include "basisPolar2D.h"
#include "printOutput.h"
#include "str.h"
#include "operatorData.h"
#include "inverseDvr.h"
#include "discretizationHybrid.h"

using namespace std;

std::map<const std::string,std::map<const Index*,std::map<const Index*,UseMatrix > > > BasisNdim::opers;
std::map<std::string,const BasisAbstract*> BasisNdim::bases;

std::map<std::string,const OperatorTree*> BasisNdim::basicOp;

// evaluate basis-tree at a single coordinate (will go into class Index)
std::vector<std::complex<double>> BasisNdim::IndexAt(const Index* Idx, const std::vector<double> Coors)
{
    if(Coors.size()==0)return {1.};

    std::vector<std::complex<double>>res(Idx->size(),0.);
    std::vector<std::complex<double>>bval;
    size_t pos=0;
    for(size_t k=0;k<Idx->childSize();k++){
        const BasisIntegrable* bas=Idx->basis()->integrable();
        auto resK=IndexAt(Idx->child(k),{Coors.begin()+bool(bas),Coors.end()});
        if(bas and std::find_if(resK.begin(),resK.end(),[](std::complex<double> v){return v!=0.;})!=resK.end()){
            if(bas->lowBound()<=Coors[0] and Coors[0]<=bas->upBound()){
                if(bval.size()==0)bval=bas->val(Coors[0]);
                for(size_t p=0;p<resK.size();p++)res[pos+p]=resK[p]*bval[k];
            }
        }
        pos+=Idx->child(k)->size();
    }
    return res;
}

const BasisAbstract* BasisNdim::factory(std::string Name){
    if(bases.count(Name)==0){
        if(allNdim(Name).size())return allNdim(Name)[0];
        return 0;
    }
    return bases[Name];
}

std::vector<const BasisNdim *> BasisNdim::allNdim(std::string Name){
    std::vector<const BasisNdim*> res;
    for(auto b: bases){
        if(b.second->name().find(Name)==0 and b.second->size()>0){
            res.push_back(dynamic_cast<const BasisNdim*>(b.second));
//            for(size_t k=1;k<res.size();k++){
//                if(res[0]->hasOverlap(*res[k])){
//                    PrintOutput::warning("Ndim function "+res[0]->str()+" has overlap with "+res[k]->str(),2);
//                }
//            }
            PrintOutput::DEVwarning("overlaps not checked");
        }
    }
    return res;
}

void BasisNdim::read(ReadInput &Inp){

    for(int line=1;line<100;line++){
        BasisPolarOff* b=new BasisPolarOff(Inp,"PolarOffCenter",line);
        bases[b->name()]=b;
        if(b->size()==0)break;
    }

    bases["Polar2D"]=new BasisPolar2D(Inp);

}

static void basisNodes(const Index *Idx,vector<vector<double> > & Node){
    if(Idx->isLeaf())return;

    if(Idx->basis()->integrable()!=0){
        if(Node.size()==0)Node.resize(1,vector<double>());

        for(size_t k=Node.size()-1;k>=0;k--){
            Node.push_back(Node[k]);
            Node[k].push_back(Idx->basis()->integrable()->lowBound()); // attach new lower boundary
            Node.back().push_back(Idx->basis()->integrable()->upBound());
        }
        basisNodes(Idx->child(0),Node);
    }
    else
        for(size_t k=0;k<Idx->childSize();k++)basisNodes(Idx,Node);
}

void BasisNdim::matrix(const string &Op, const Index* IIndex, const Index* JIndex, UseMatrix &Mat, complex<double> Multiplier){

    if(opers.count(Op)==0 or
            opers[Op].count(IIndex)==0 or
            opers[Op][IIndex].count(JIndex)==0)
    {

        const Index* iRoot=IIndex;
        const Index* jRoot=JIndex;


        // try cast to multi-dimensional
        const BasisNdim * iNdim=iRoot->basis()->ndim();
        const BasisNdim * jNdim=jRoot->basis()->ndim();
        if(iNdim==0 and jNdim==0)PrintOutput::warning("at least one basis should be BasisNdim - else use standard product");

        if(iNdim==0){
            vector<vector<double> > node;
            basisNodes(IIndex,node);
        }

        if(iNdim==0){iRoot=jNdim->rootIndex(IIndex);iNdim=BasisProd::factory(iRoot,jRoot);}
        if(jNdim==0){jRoot=iNdim->rootIndex(JIndex);jNdim=BasisProd::factory(jRoot,iRoot);}



        OperatorNdim opNd(Op,iNdim->_quadCoor,jNdim->_quadCoor);

        // auxiliary indexing for fast use in loop
        vector<const Index*> iFloor,jFloor;
        vector<int> iPos,jPos;
        vector<vector<complex<double>*> >mats;
        for(const Index * jF=jRoot->firstFloor();jF!=0;jF=jF->nodeRight(jRoot)){
            jFloor.push_back(jF);
            jPos.push_back(jF->posIndex(jRoot));
        }
        for(const Index * iF=iRoot->firstFloor();iF!=0;iF=iF->nodeRight(iRoot)){
            iFloor.push_back(iF);
            iPos.push_back(iF->posIndex(iRoot));
        }
        for(size_t iF=0;iF<iFloor.size();iF++){
            mats.push_back(vector<complex<double>*>());
            for(size_t jF=0;jF<jFloor.size();jF++){
                opers[Op][iFloor[iF]][jFloor[jF]]=UseMatrix::Zero(iFloor[iF]->size(),jFloor[jF]->size());
                mats.back().push_back(opers[Op][iFloor[iF]][jFloor[jF]].data());
            }
        }

        // loop through all points
        for(size_t pt=0;pt<iNdim->quadSize();pt++){
            // loop through all terms
            for(size_t trm=0;trm<opNd.terms();trm++){

                // evaluate function for term
                vector<complex<double> > cgrid;
                complex<double> fij=1.;
                for(size_t k=0;k<iNdim->quadGrid(pt).size();k++)
                    cgrid.push_back(iNdim->_comSca[k].xScaled(iNdim->quadGrid(pt)[k]));
                fij=jNdim->quadWeig(pt)*opNd.term(trm).kernel(cgrid);
                int ivd=opNd.term(trm).ivd();
                int jvd=opNd.term(trm).jvd();

                // add contribution to all matrix elements
                for(size_t jF=0;jF<jFloor.size();jF++){
                    for(size_t j=0;j<jFloor[jF]->size();j++){
                        complex<double> fijBas=fij*jNdim->valPartial(pt,jPos[jF]+j,jvd);
                        for(size_t iF=0;iF<iFloor.size();iF++){
                            int ij=j*iFloor[iF]->size();
                            for(size_t i=0;i<iFloor[iF]->size();i++,ij++){
                                mats[iF][jF][ij]+=fijBas*iNdim->valPartial(pt,iPos[iF]+i,ivd);
                            }
                        }
                    }
                }
            }
        }

        for(size_t iF=0;iF<iFloor.size();iF++)
            for(size_t jF=0;jF<jFloor.size();jF++){
                opers[Op][iFloor[iF]][jFloor[jF]].purge();
            }
    }
    if(opers[Op].count(JIndex)!=0 and opers[Op][JIndex].count(IIndex)!=0){
        if(not
                (opers[Op][JIndex][IIndex]-opers[Op][IIndex][JIndex].transpose()).isZero(1.e-12*max(1.,opers[Op][JIndex][IIndex].maxAbsVal()))
                and not
                (opers[Op][JIndex][IIndex]+opers[Op][IIndex][JIndex].transpose()).isZero(1.e-12*max(1.,opers[Op][JIndex][IIndex].maxAbsVal()))){
            PrintOutput::DEVmessage("Operator is neither symmetric nor anti-symmetric "+Op);
        }
    }

    // get operator from storage and multiply
    Mat=opers[Op][IIndex][JIndex]*Multiplier;
}

void BasisNdim::productGridWeig(const Index* Idx, std::vector<std::vector<double> > & Grid, std::vector<double> & Weig,
                                int MinQuad,std::vector<double>GridK,double WeigK)
{
    if(GridK.size()==0){
        Grid.clear();
        Weig.clear();
    }
    if(Idx->isLeaf()){
        Grid.push_back(GridK);
        Weig.push_back(WeigK);
        return;
    }

    UseMatrix pts,wgs;
    //NOTE: orders are somewhat sloppy...
    //    int addN=max(int(Idx->basis()->order()+2),MinQuad-int(Idx->basis()->order()));
    const BasisIntegrable* bas=Idx->basis()->integrable();
    int addN=max(int(bas->order()),MinQuad);
    bas->quadRule(bas->order()+addN,pts,wgs);
    bas->jacobian()->operator()(pts,wgs);

    GridK.resize(GridK.size()+1);
    for(size_t k=0;k<pts.size();k++){
        GridK.back()=pts(k).real();
        productGridWeig(Idx->descend(),Grid,Weig,MinQuad,GridK,WeigK*wgs(k).real());
    }
}

const Index* BasisNdim::rootIndex(const Index *Idx) const{

    // find coordinates in sequence, starting from last
    vector<string> coor=tools::splitString(_quadCoor,'.');
    Idx=Idx->firstFloor();
    for(size_t k=coor.size();k>0;k--)
        while(Idx!=0 and Idx->axisName()!=coor[k-1])Idx=Idx->parent();
    return Idx;
}

void BasisNdim::mainQuadValDer() {
    // transform _quadGrid, _quadWeig, _valDer to main coordinates
    // defining transformation is in derived class's "toCartesian"
    _comSca.resize(dim()); // cannot complex scale off-center

    CoordinateMap toMain=CoordinateTrans::fromCartesian(_quadCoor);
    JacobianMap jacToMain=CoordinateTrans::jacCartesian(_quadCoor);
    IntegrationFactor facMain=CoordinateTrans::integrationFactor(_quadCoor);
    NablaFactor nabMain=CoordinateTrans::nablaFactor(_quadCoor);

    for(size_t pt=0;pt<_quadGrid.size();pt++){
        // map from ndim to main coordinates (going through shared reference, usually cartesian)
        vector<double>ndimCoor(_quadGrid[pt]);
        _quadGrid[pt]=(toMain(toCartesian(ndimCoor)));

        // adjust quadrature weights to refer to main coordinates (origin at 0 0 0)
        double factorRatio=facMain(_quadGrid[pt])/absFactor(ndimCoor);
        _quadWeig[pt]/=pow(factorRatio,2);

        // create the transformation at point
        Eigen::MatrixXd trans=jacobianToNdim(ndimCoor);
        trans=Eigen::Map<Eigen::MatrixXd>(jacToMain(_quadGrid[pt]).data(),dim(),dim())*trans.inverse().eval();
        trans*=factorRatio; // for integrations wrt main coordinate system

        // transform to value and gradient wrt main coordinate system
        for(size_t ibas=0;ibas<size();ibas++){
            for(size_t i=0;i<dim();i++)valPartial(pt,ibas,1+i)-=
                    value(pt,ibas)*nablaFactor(ndimCoor,i)/absFactor(ndimCoor); // d/dr' f - f/r'
            Eigen::Map<Eigen::VectorXcd>(_valDer[pt][ibas].data()+1,dim())=
                    trans*Eigen::Map<Eigen::VectorXcd>(_valDer[pt][ibas].data()+1,dim());
            _valDer[pt][ibas][0]*=factorRatio; // for integrations wrt main coordinate system
            for(size_t i=0;i<dim();i++)_valDer[pt][ibas][1+i]+=
                    _valDer[pt][ibas][0]*nabMain(_quadGrid[pt],i)/facMain(_quadGrid[pt]); //  + f*r/(r')^2
        }
    }
}
std::string BasisNdim::mainCoor(ReadInput & Inp){
    // get the reference coordinates from input
    std::vector<std::string> s;
    for(int l=1;not Inp.endCategory("Axis",l);l++){
        if(std::find(s.begin(),s.end(),Axis::readName(Inp,l))==s.end())
            s.push_back(Axis::readName(Inp,l));
    }
    std::string ndims="Ndim Orbital";
    for(std::string n: tools::splitString(ndims,' ')){
        if(s[0]==n)return tools::joinString(std::vector<std::string>({s.begin()+1,s.end()}),".");
    }
    ABORT("need any of"+ndims+" as first axis, is: "+s[0]);
    return "";
}


std::string BasisNdim::str(int Level) const{
    std::string s="ndim("+_quadCoor+")";

    if(Level==0)return s;

    for(size_t pt=0;pt<_quadGrid.size();pt++){
        s+="\n"+tools::str(pt)+": "+tools::str(_quadGrid[pt]);
        for(size_t ibas=0;ibas<size();ibas++){
            s+="\n"+tools::str(ibas,5)+": "+tools::str(_valDer[pt][ibas]);
        }
    }
    return s;
}

void BasisNdim::test(){
    return;
    vector<double> cart;
    for (size_t k=0;k<dim();k++)cart.push_back(rand()%2048/(2048.));

    CoordinateMap toMain=CoordinateTrans::fromCartesian(quadCoor());
    CoordinateMap fromMain=CoordinateTrans::toCartesian(quadCoor());

    vector<double> orig;
    orig=fromMain(toMain(cart));
    for(size_t k=0;k<dim();k++)
        if(abs(orig[k]-cart[k])>1.e-12){
            ABORT(Str("transformation to coordinates ")+quadCoor()+"failed, trans="+orig+"init="+cart);
        }
    PrintOutput::DEVmessage("OK coordinate trans to quad="+quadCoor());

    orig=toCartesian(CoordinateTrans::fromCartesian(_ndimCoor)(cart));
    for(size_t k=0;k<dim();k++)
        if(abs(orig[k]-cart[k])>1.e-12){
            ABORT(Str("transformation to coordinates ")+_ndimCoor+"failed, trans="+orig+"init="+cart);
        }
    PrintOutput::DEVmessage("OK coordinate trans to ndim="+_ndimCoor);
}







