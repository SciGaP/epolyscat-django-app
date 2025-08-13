// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "discretizationSurface.h"

#include "readInput.h"

#include "threads.h"
#include "folder.h"
#include "algebra.h"
#include "index.h"
#include "indexSurface.h"
#include "operatorTensor.h"
#include "pulse.h"
#include "mapGauge.h"
#include "discretizationHybrid.h"
#include "discretizationSurfaceHybrid.h"
#include "discretizationSurface2D.h"
#include "discretizationtsurffspectra.h"
#include "axisTree.h"
#include "operatorFloor.h"
#include "basisSub.h"
#include "basisVector.h"
#include "basisIntegrable.h"
#include "basisExpIm.h"
#include "basisAssocLeg.h"
#include "basisGrid.h"
#include "eigenTools.h"
#include "parallelOperator.h"
#include "tsurffSource.h"

using namespace std;

string DiscretizationSurface::prefix="surface_";
string DiscretizationSurface::prefixAx="surf"; // NOTE: this is hard-coded in many places, keep "surf" for now

std::string DiscretizationSurface::readDef(ReadInput& Inp){
    std::string def;
    Inp.read("Surface","definition",def,"not given","convert to: Spherical[R,Lmax{,Mmax}]...spherical basis with maximal angular momenta")
            .texdocu(R"tex(
                     if Mmax is not specified, it will be inferred from input discretization
                     )tex");
    return def;
}

std::string DiscretizationSurface::read(ReadInput &Inp, std::vector<double> & Points, bool & WriteAscii){

    Inp.exclude("Surface","points","Surface","definition");
    std::string def=readDef(Inp);

    if(def=="not given"){
        Points.clear();
        Inp.read("Surface","points",Points,"","save values and derivatives at points (blank-separated list), default 2 au befor absorption");
    }
    ReadInput::main.read("Surface","ascii",WriteAscii,ReadInput::flagOnly,"produce ascii file with surface values",0,"surfaceAscii");
    if(Points.size()==0 and ReadInput::main.found("Spectrum"))Points.assign(1,DBL_MAX);
    return def;
}

std::vector<DiscretizationSurface*> DiscretizationSurface::all(ReadInput& Inp, std::shared_ptr<Discretization> D, std::vector<double> & surf, std::string Region){
    // messy, sorry...
    std::vector<DiscretizationSurface*> discSurf;
    if(Region=="" and (TsurffSource::readSymmetry12(Inp,D->idx()->hierarchy())
                       or DiscretizationSurface::readDef(Inp)!="not given")){
        discSurf.push_back(DiscretizationSurface::factory(D.get(),surf,0));
    }
    else
        for(unsigned int c=0;surf.size()>0 and D->idx()->continuity(c)!=Index::npos;c++)
            discSurf.push_back(DiscretizationSurface::factory(D.get(),surf,c)); // Surface disc corresponding to each continuity level
    return discSurf;
}


string DiscretizationSurface::surfacePath(std::string RunDir, const Index * Idx, std::string Region){

    std::string src=RunDir;
    if(Region.find('.')==std::string::npos)
        src+=DiscretizationSurface::prefix+Region;
    else
        src+="S_"+Region.substr(0,Region.find('.'))+"/"+DiscretizationSurface::prefix+Region.substr(Region.rfind('.')+1);
    if(not folder::exists(src))
        ABORT("Could not find surface file "+src+" needed as source for subregion "+Region
              +", first run with -region="+Region.substr(0,Region.find('.')));

    return src;
}

DiscretizationSurface::DiscretizationSurface(std::string SrcFile)
{
    std::ifstream surfStream;
    surfStream.open( SrcFile.c_str(), ios::in|ios::binary);
    if(not surfStream.is_open())ABORT("coul not find surface file" +SrcFile);
    idx()=new Index(surfStream);
}


DiscretizationSurface * DiscretizationSurface::factory(const Discretization *Parent, std::vector<double> &Rad, unsigned int Nsurf,
                                                       std::string SrcFile){
    if(SrcFile!="")
        return new DiscretizationSurface(SrcFile);

    bool bdum;
    std::string def;
    if((def=read(ReadInput::main,Rad,bdum))!="not given")
        return new DiscretizationSurface2D(Parent->idx(),def);       // surface was defined explicitly
    else if(dynamic_cast<const DiscretizationHybrid*>(Parent)!=0)
        return new DiscretizationSurfaceHybrid(dynamic_cast<const DiscretizationHybrid*>(Parent),Rad,Nsurf);
    else
        return new DiscretizationSurface(Parent,Rad,Nsurf);
}

// return the NSurf'th continuity index
const Index* getContinuityLevel(const Index* Idx, int NSurf){

    const Index* idx=Idx;
    while(not idx->isLeaf()){
        //HACK for hybrid: second part contains continuity
        if(idx->axisName().find("&")!=string::npos)idx=idx->child(1);

        const Index* idxFloor=idx->descend()->axisIndex(idx->axisName()); // find matching axis name
        if(idx->basis()->isIndex() and idxFloor!=0 and idxFloor->basis()->integrable())
        {
            if(NSurf==0)return idx;
            NSurf--;
        }
        idx=idx->descend();
    }
    PrintOutput::DEVwarning(Sstr+"no continuity level"+NSurf+"in index\n"+Idx->str());
    return 0; // to shut up some over-ambitious compilers (CLANG)
}

DiscretizationSurface::DiscretizationSurface(const Discretization *Parent, const std::vector<double> &Rad, int NSurf, string DerivativeSide)
{
    if(Rad.size()!=1)
        ABORT(Str("specify exactly 1 surface radius, found")+Rad.size()+"at"+Rad);

    // create new hierarchy and continuitLevel: add value/derivative index just above floor
    parent = dynamic_cast<const Discretization *>(Parent);
    name = Parent->name;
    // get the name of the NSurf'th continuity axis

    if(Rad.size()==0)ABORT("specify at least one surface radius")
            wIdx = 0;

    // construct global index
    idx() = new IndexS(Parent->idx(),Rad,NSurf,DerivativeSide);

    _mapFromParent.reset(new MapSurface(idx(),Parent->idx()));
    ParallelOperator::setDistribution(mapFromParent());

    _construct();

}

DiscretizationSurface::DiscretizationSurface(const Discretization *Parent, const std::vector<double> &Rad, const std::string SurfaceFile, string DerivativeSide)
{
    if(Rad.size()!=1)
        ABORT(Str("specify exactly 1 surface radius, found")+Rad.size()+"at"+Rad);

    // create new hierarchy and continuitLevel: add value/derivative index just above floor
    parent = dynamic_cast<const Discretization *>(Parent);
    name = Parent->name;
    std::string surfCoor=SurfaceFile.substr(SurfaceFile.rfind("_")+1);
    //    name.insert(name.find(surfCoor)+surfCoor.length(),"_surf");

    if(Rad.size()==0)ABORT("specify at least one surface radius");
    ifstream surfStream(SurfaceFile.c_str(), ios::in|ios::binary);
    idx() = new Index(surfStream);

    _mapFromParent.reset(new MapSurface(idx(),Parent->idx()));
    _construct();
}

void DiscretizationSurface::_construct(){
    ParallelOperator::setDistribution(mapFromParent());

    // mixed gauge: apply back-gauge transform before write
    // (do not back-transform if length gauge)
    if(Algebra::isAlgebra("Rg") and Algebra("Rg").val(0.).real()!=0 and Algebra("Rg").val(0.).real()<DBL_MAX/2){
        std::shared_ptr<OperatorAbstract> m(new MapGauge(this));
        _mapFromParent=m;
    }

}

DiscretizationSurface::MapSurface::MapSurface(const Index *SurfI, const Index *FromI)
    :OperatorTree("mapTo"+SurfI->hierarchy(),SurfI,FromI)
{
    // set top level name that is used as surface file name
    if(SurfI->parent()==0)
        name=prefix+SurfI->findAxisStarts(prefixAx)->axisName().substr(prefixAx.length());

    if(FromI->hasFloor()){
        // one of the coordinates of Floor is converted to surface
        std::vector<const Eigen::MatrixXcd*> pMats;
        Eigen::MatrixXcd * pmat;
        for(const Index*from=FromI,*surf=SurfI;;from=from->descend(),surf=surf->descend()){
            if(not from->subEquivalent())
                DEVABORT("surfaces only for floors with tensor product structure, is:\n"+FromI->str());

            // tensor product of identity operators and conversion to value and derivative
            if(surf->axisName().substr(0,6)=="ValDer"){
                std::vector<std::complex<double> > val,der;
                if(surf->basis()->grid()->mesh().size()!=2)DEVABORT(Sstr+"not 2 grid points at ValDer:"+surf->basis()->grid()->mesh());
                from->basis()->integrable()->valDerD({surf->basis()->grid()->mesh()[0]},val,der);

                if(from->basis()->integrable()->eta()!=1.)
                    ABORT(Sstr+"Illegal: surface"+surf->basis()->grid()->mesh()[0]+"in complex scaled region, basis"+from->basis()->str());

                pmat=new Eigen::MatrixXcd(2,val.size());
                for(size_t k=0;k<from->basis()->size();k++){
                    //HACK true radial derivative here - this should be changed for consistency, ugly HACK?? for XYZ
                    if(FromI->axisName().find("Rn")!=string::npos or FromI->axisName().find_first_of("XYZ")!=string::npos)
                        der[k]-=val[k]/surf->basis()->grid()->mesh()[0];
                    pmat->operator()(0,k)=val[k];
                    pmat->operator()(1,k)=der[k];
                }
            } else {
                pmat=new Eigen::MatrixXcd(surf->basis()->size(),surf->basis()->size());
                *pmat=Eigen::MatrixXcd::Identity(surf->basis()->size(),surf->basis()->size());
            }
            pMats.push_back(pmat);

            if(from->isBottom())break;
        }
        floor()=OperatorFloor::factory(pMats,"toSurf"+FromI->hash());
        for(auto m: pMats)delete m;
    }
    else if(SurfI->axisName().find(prefixAx)==0){
        // FE level of surface and matching branch of SurfI
        const BasisSub* sub=dynamic_cast<const BasisSub*>(SurfI->basis());
        if(!sub and SurfI->basis()->strDefinition()!="Vector:1")
            DEVABORT("expected BasisSub or Vector:1, got: "+SurfI->basis()->str()+"\n"+SurfI->root()->str());
        for(size_t k=0;sub and k<SurfI->childSize();k++)
            childAdd(new MapSurface(SurfI->child(k),FromI->child(sub->subset()[k])));
    } else {
        // identity on all other axes, except for hybrid where terms will have been omitted
        for(unsigned int k=0,kOrig=0;k<SurfI->childSize();k++,kOrig++){
            // for hybrid axes
            while(FromI->axisName().find("&")!=string::npos and FromI->child(kOrig)->axisName()!=SurfI->child(k)->axisName())kOrig++;
            childAdd(new MapSurface(SurfI->child(k),FromI->child(kOrig)));
        }
    }
}

static void purgeIdx(Index* Idx){
    // get maximal height
    size_t maxH=0;
    for(Index* idx=Idx;idx;idx=idx->nodeNext())maxH=std::max(maxH,size_t(idx->depth()));

    // standard Index purge
    Idx->purge(maxH);

    // replace subsets by full
    for(;Idx;Idx=Idx->nodeNext()){
        if(Idx->basis()->sub()){
//            if(Idx->axisName()=="Eta" or Idx->axisName()=="Phi"){
//                std::vector<int> sub=Idx->basis()->sub()->subset();
//                for(size_t k=0;k<sub.size();k++)
//                    if(sub[k]!=int(k))DEVABORT(Sstr+"only for lower subsets, got:"+Idx->basis()->str());
//                const BasisAbstract* bas=BasisSub::superBas(Idx->basis());
//                std::string def=bas->strDefinition();
//                if(dynamic_cast<const BasisAssocLeg*>(bas)){
//                    def=def.substr(0,def.find(",")+1)+tools::str(sub.size());
//                }
//                else if(dynamic_cast<const BasisTrigon*>(Idx->basis())){
//                    std::vector<std::string> ms=tools::splitString(def.substr(0,def.find(":")+1),',');
//                    def=def.substr(0,def.find(":")+1);
//                    for(auto m: ms)def+=m+",";
//                    def.pop_back();
//                }
//                else
//                    DEVABORT("only for AssocLeg or BasisTrigon, got:"+def);
//                Idx->setBasis(BasisAbstract::factory(def));
//            }
        }
    }

}

// Algorithm:
// duplicate structure in From index until the NSurf'th continuity level is met
// attach a copy of From branch that contains the surface
// replace the node at the continuity level for coordinate C with surfC (grid, at present limited to single point)
// replace the basis covering the surface coordinate C with ValeDerC consiting of value and derivative at the surface point
DiscretizationSurface::IndexS::IndexS(const Index *From, std::vector<double> Radius, unsigned int NSurf, string DerivativeSide)
{
    //    if(From->isLeaf())ABORT(From->root()->str()+"\n\nNo unbounded axis found in Index");

    Index::nodeCopy(From,false); // true copy, not view
    fromIndex.push_back(From);

    const Index* fIndex=getContinuityLevel(From,0);
    if(fIndex and (fIndex!=From or NSurf>0)){
        // not the continuity level for desired surface - descend
        if(fIndex==From)NSurf--; // found continuity - reduce surface count

        // hybrid axis
        int kMin=From->axisName().find("&")!=string::npos ? 1 : 0;
        if(kMin){
            // must remove stub of first axis
            std::vector<int> subs(1,kMin);
            for(size_t k=kMin+1;k<From->basis()->size();k++)subs.push_back(k);
            setBasis(BasisAbstract::factory(BasisSub::strDefinition(From->basis(),subs)));
        }

        for(unsigned int k=kMin;k<From->childSize();k++)
            childAdd(new IndexS(From->child(k),Radius,NSurf));
    }

    else if(fIndex)
    {
        // indicate surface level
        setAxisName(prefixAx+From->axisName());
        fromIndex.push_back(From);
        vector<int> elemN;
        if(Radius.size()!=1)ABORT("specify only single surface radius");
        double surf=Radius[0];

        // loop through sections below From (to locate Radius)
        const Index* elem=From->descend();
        for(;elem!=0;elem=elem->rightSibling()){
            Index* elemFunc=elem->axisIndex(From->axisName());
            double lB=elemFunc->basis()->integrable()->lowBound(),uB=elemFunc->basis()->integrable()->upBound();
            double eps=min((uB-lB)*1.e-12,1.e-10); // eps is tricky, as boundaries may be at +-infty;

            // maximal surface specfied, make automatic selection
            const Index* next=elemFunc->upperNeighbor();
            if(abs(surf)>DBL_MAX/2 and next){
                if(elemFunc->basis()->integrable()->eta()==1. and next->basis()->integrable()->eta()!=1.)
                    surf=uB;
            }

            // locate inside element; if on boundary
            // derivatives are taken from the side closer to zero, unless specified explicitly at boundary
            bool selectThis=false;
            if(DerivativeSide=="boundaryFromBelow")
                selectThis=abs(uB-surf)<eps;
            else if (DerivativeSide=="boundaryFromAbove")
                selectThis=abs(lB-surf)<eps;
            else if (DerivativeSide=="fromInside")
                selectThis=
                        (lB>=0 and lB< surf     and surf<=uB+eps) or
                        (lB< 0 and lB<=surf+eps and surf< uB);
            else
                DEVABORT("illegal value of DerivativeSide="+DerivativeSide
                         +", allowed are boundaryFromBelow,boundaryFromAbove,fromInside(default)");

            if(selectThis){
                // attach copy of tree that contains surface
                for(unsigned int l=childSize();l>0;l--)childErase(l-1);
                elemN.push_back(elem->nSibling());
                childAdd(new Index(*From->child(elemN.back())));
                // insert v/d at matching level
                for(Index *vd=childBack()->axisIndex(From->axisName());vd!=0;vd=vd->nodeRight(),elemFunc=elemFunc->nodeRight()){
                    // erase function level
                    for(int l=vd->childSize();l>0;l--)vd->childErase(l-1);
                    // insert v/d level
                    // NOTE: we need all possible levels below to be equivalent
                    if(elemFunc->descend())for(int l=0;l<2;l++)vd->childAdd(new Index(*elemFunc->descend()));
                    else                   {vd->leafAdd();vd->leafAdd();} // for now, surface index DOES carry dummies
                    vd->setBasis(BasisGrid::factory(std::vector<double>(2,surf)));
                    vd->setAxisName("ValDer"+From->axisName());
                }
                break;
            }
        }
        if(elem==0){
            if(DerivativeSide.substr(0,8)=="boundary"){
                ABORT(Str("surface does not coincide with any element boundary:")+surf);
            } else if(abs(surf)<DBL_MAX/2){
                ABORT(Str("surface outside all elements:")+surf);
            }
        } else {
            setBasis(BasisAbstract::factory(BasisSub::strDefinition(From->basis(),elemN)));
        }
    }
    purgeIdx(this); // some branches may not end in surfaces
    sizeCompute();
}
