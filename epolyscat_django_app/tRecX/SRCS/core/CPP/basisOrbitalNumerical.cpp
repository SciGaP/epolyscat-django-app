// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "basisOrbitalNumerical.h"
#include "basisNdim.h"
#include "basisIntegrable.h"
#include "basisGrid.h"
#include "basisSub.h"

#include "indexGrid.h"
#include "operatorMap.h"
#include "indexNew.h"
//#include "indexGrid.h"
#include "inverse.h"

#include "vectorValuedFunction.h"
#include "coordinateTrans.h"
#include "coorSystem.h"

#include "readInput.h"
#include "printOutput.h"
#include "plot.h"

#include "eigenTools.h"
#include "integrate.h"

static double integrationAccuracy=1e-7;
static std::string integrationMethod="singleQuadrature";
static int points=4;

/// return order at all integrable axes or at specified Axes
static void volumeOrders(std::vector<std::vector<double>> &Vol, std::vector<size_t> &Ord, const Index* Idx, const std::vector<std::string> Axes={} /** include only these */){

    if(Idx->descend())volumeOrders(Vol,Ord,Idx->descend(),Axes);
    if(Idx->basis()->integrable() and (Axes.size()==0 or std::find(Axes.begin(),Axes.end(),Idx->axisName())!=Axes.end())){
        size_t ord(0);
        double lb=+DBL_MAX,ub=-DBL_MAX;
        for(const Index* idx=Idx;idx!=0;idx=idx->nodeRight()){
            ord=std::max(ord,size_t(idx->basis()->integrable()->order()));
            lb=std::min(lb,idx->basis()->integrable()->lowBound());
            ub=std::max(ub,idx->basis()->integrable()->upBound());
        }
        Ord.insert(Ord.begin(),ord);
        Vol.insert(Vol.begin(),{lb,ub});
    }
}

/// return volume containing all integrable axes or at specified Axes
static std::vector<std::vector<double>> volume(const Index* Idx, const std::vector<std::string> Axes={} /** include only these */){
    std::vector<size_t>ord;
    std::vector<std::vector<double>> vol;
    volumeOrders(vol,ord,Idx,Axes);
    return vol;
}

/// return order at all integrable axes or at specified Axes
static std::vector<size_t> orders(const Index* Idx, size_t AddOrd, const std::vector<std::string> Axes={} /** include only these */){
    std::vector<size_t>ord;
    std::vector<std::vector<double>> vol;
    volumeOrders(vol,ord,Idx,Axes);
    for(auto &o: ord)o+=AddOrd;
    return ord;
}

void BasisOrbitalNumerical::setup(){
    std::vector<const BasisAbstract*> bases
            =BasisAbstract::select([&](const BasisAbstract* bas){return dynamic_cast<const BasisOrbitalNumerical*>(bas);});

    bool hasCO2=false;
    const BasisOrbitalNumerical* b;
    for(auto bb: bases){
        b=dynamic_cast<const BasisOrbitalNumerical*>(bb);
        b->orbitals();
        hasCO2=hasCO2 or b->name().find("CO2")!=std::string::npos;
    }
}

static std::string _printList;
static const std::vector<std::pair<std::string,std::string>> _printProperties={
    {"T","1/2<<Laplacian>>"},
    {"L(L+1)","<<AngularMomentumSquared>>"},
    {"Lz^2","<<AngularMomentumZsquared>>"},
    //    {"Lz","<1_d><1><1>-<d_1><1><1>"},
    {"cos^8Phi","<pow[8](cos(Q))><1><1>"},
    {"cos^2Phi","<cos(Q)*cos(Q)><1><1>"},
    {"sin^2Phi","<sin(Q)*sin(Q)><1><1>"},
    //    {"Parity","<<Parity>>"}, // nice to have, but not implemented
    {"X","<<X>>"},
    {"Y","<<Y>>"}
};

void BasisOrbitalNumerical::read(ReadInput &Inp){
    std::string mess;
    for(auto p: _printProperties)mess+=p.first+" ";
    Inp.read("Orbital","properties",_printList,"","print properties (blank-separted list), choices are "+mess);
    for(int line=1;line<10;line++){
        std::string def;
        Inp.read("Orbital","definition",def,"","defintion of the orbitals"+mess,line);
        if(def=="")break;
        BasisAbstract::factory("Orbital:"+def);
    }
    Inp.read("Orbital","integrationAccuracy",integrationAccuracy,"1e-7","accuracy target for numerical integration, values<1e-7 may lead to slow convergence");
    Inp.read("Orbital","method",integrationMethod,"singleQuadrature",
             "methods: singleQuadrature...evaluate on overall quadrature grid,"
             "volumeSplit...quadratures on subvolumes");
    Inp.read("Orbital","integrationPoints",points,"5","number of points beyond standard quadrature");
}

static std::string strFromDef(const BasisSetDef & Def){
    std::string s(tools::stringInBetween(Def.funcs,"[","]"));
    s+=":"+tools::str(int(Def.lowBound()));
    s+=":"+tools::str(int(Def.lowBound())+Def.order);
    return s;
}

BasisOrbitalNumerical::BasisOrbitalNumerical(const BasisSetDef &Def)
    :BasisOrbitalNumerical(strFromDef(Def)){}

std::string BasisOrbitalNumerical::strDefinition() const{
    return "Orbital:"+_funcDef+":"+_refName+":"+tools::str(_select,3,",");
}

BasisOrbitalNumerical::BasisOrbitalNumerical(const std::string Def)
    :BasisOrbital("NumOrb")
{
    std::string def(Def);
    if(def.find("Orbital:")==0)def=def.substr(8);
    else if(def.find("Orbital[")==0)def=tools::stringInBetween(Def,"[","]");
    std::vector<std::string> part=tools::splitString(def,':',"[(<","])>");

    // a subset of orbitals may be defined with the orbital name
    std::vector<int> subset;
    std::string ssel=tools::stringInBetween(part[0],"{","}");
    if(ssel!=part[0])subset=tools::string_to_intVec("{"+ssel+"}");

    // determine which orbitals to use in the basis
    if(part.size()==4){
        if(_select.size())DEVABORT("_select has been set earlier");
        _select=tools::string_to_intVec(part[2]+":"+part[3]);
    }
    else if(part.size()==3){
        _select=tools::string_to_intVec("{"+part[2]+"}");
    }
    else if (part.size()>3)
        ABORT("need formats Func, Func:Refidx, or Func:Refidx:ibeg:iend, got: "+def);

    if(subset.size()>0){
        // select from subset

        if(_select.size()==0)
            _select=subset;
        else{
            if(int(subset.size())<1+*std::max_element(_select.begin(),_select.end()))
                ABORT(Sstr+"Subset size smaller than maximal selection"+_select);
            for(int &s: _select)s=subset[s];
        }
    }

    if(_select.size()==0)
        ABORT("need to specify selection of orbitals as Name:Ref:ibeg:iend, Name:Ref:1,2,6,4 or Name{1,2,3..}:Ref, got: "+Def);


    _funcDef=part[0].substr(0,part[0].find("{")); // selection is in _select now
    if(part.size()>1)_refName=part[1];
    else _refName="auto";
    _orb.resize(_select.size());
    _name="NumOrb["+_funcDef+(_refName!="auto"?"_"+_refName:"")+"]";
}

static std::vector<double> gridPoint(const Index* Idx){
    std::vector<double> pts;
    for(const Index * idx=Idx;idx->parent()!=0;idx=idx->parent()){
        const BasisGrid* b=idx->parent()->basis()->grid();
        if(b!=0)pts.push_back(b->mesh()[idx->nSibling()]);
    }
    std::reverse(pts.begin(),pts.end());
    return pts;
}

std::string BasisOrbitalNumerical::str(int Level) const{
    Str s("NumOrb[","");
    if(Level)s+=""; // just to shut up the compiler about unused parameter
    s+=_funcDef+"] ";
    if(size()==1)     s=s+"no."+_select[0];
    else if(size()<10)s=s+"{"+_select+"}";
    else              s=s+_select.front()+",..,"+_select.back()+"} ["+size()+"]";
    if(orbitalFunctions()==0)return "function not found: "+_funcDef;
    return s+" out of "+orbitalFunctions()->info();
}

/// set CoefficientFloors with f.innerProduct(f)<Eps*C.innerProduct(C) to =0
static void purge(Coefficients & C, double Eps=1.e-20){
    double eps=Eps*abs(C.innerProduct(&C));
    for(Coefficients* f=C.firstLeaf();f!=0;f=f->nextLeaf())
        if(abs(f->innerProduct(f))<eps)f->setToZero();
}

static void mapNdim(const Index* Idx, const VectorValuedFunction * Orbs, std::vector<int> Select,std::vector<Coefficients*> Orb){
    CoordinateMap toCartesian=CoordinateTrans::toCartesian(Idx->basis()->ndim()->quadCoor());
    CoordinateMap toFunctionArg=CoordinateTrans::fromCartesian(Orbs->coordinates());
    IntegrationFactor factorIdx=CoordinateTrans::integrationFactor(Idx->basis()->ndim()->quadCoor());
    // for coordinates X.Y.Z factorFunc = 1
    IntegrationFactor factorFunc=CoordinateTrans::integrationFactor(Orbs->coordinates());

    std::vector<Coefficients> gVal;

    if(Idx->basis()->ndim()->quadSize()==0)DEVABORT("no quadGrid");
    for(size_t k=0;k<Idx->basis()->ndim()->quadSize();k++){
        std::vector<double> point=Idx->basis()->ndim()->quadGrid()[k];
        std::vector<double>arg=toFunctionArg(toCartesian(point));
        if(not Orbs->validArguments(arg))
            Orbs->abortInvalidArguments(tools::str(arg)+" converted from "+tools::str(point));
        double ratio=factorFunc(arg);
        ratio=ratio==0.? 1. : factorIdx(point)/ratio; //quadrature point=0 may be included, but is not needed here...

        std::vector<std::complex<double> > orbVal=(*Orbs)(arg);

        double ratioWeig=ratio*Idx->basis()->ndim()->quadWeig()[k];
        for(size_t l=0;l<Orb.size();l++){
            std::complex<double> orbValWeig=orbVal[Select[l]]*ratioWeig;
            for(size_t b=0;b<Idx->basis()->ndim()->size();b++){
                Orb[l]->orderedData()[b]+=std::conj(Idx->basis()->ndim()->valDers()[k][b][0])*orbValWeig;
            }
        }
    }
    for(auto o: Orb)purge(*o,1.e-20);
    Eigen::MatrixXcd vals(Orb[0]->size(),Orb.size());
    for(size_t j=0;j<Orb.size();j++)
        vals.col(j)=Eigen::Map<Eigen::MatrixXcd>(Orb[j]->values().data(),vals.rows(),1);
    PrintOutput::matrix(vals.real(),5);
}

void BasisOrbitalNumerical::mapTreeProduct(const Index* Idx, const VectorValuedFunction * Orbs, std::vector<int> Select,std::vector<Coefficients*> Orb,
                                           std::vector<unsigned int> Points){

    CoordinateMap fromIdx=CoordinateTrans::toCartesian(Idx->coordinates());
    IntegrationFactor factorIdx=CoordinateTrans::integrationFactor(Idx->coordinates());
    CoordinateMap toFunctionArg=CoordinateTrans::fromCartesian(Orbs->coordinates());
    IntegrationFactor factorFunc=CoordinateTrans::integrationFactor(Orbs->coordinates());

    std::vector<std::complex<double> > val;
    std::vector<Coefficients> gVal;

    std::vector<std::string> ax=tools::splitString(Idx->coordinates(),'.');
    std::vector<unsigned int> pts(ax.size(),0);
    if(Points.size()==0){
        // get excess size quadrature points
        auto pt(orders(Idx,points,ax));
        for(size_t k=0;k<pt.size();k++)pts[k]=pt[k];
    }
    else {
        pts=Points;
    }

    // get discretization on a quadrature grid
    std::vector<std::vector<double>>gg(ax.size()),ww(ax.size());
    IndexGrid gdx(Idx,ax,gg,ww);
    const Index * idxLeaf=gdx.firstLeaf();
    // set up the vectors
    if(Select.size()==0)
        for(size_t k=0;k<Orbs->length();k++)Select.push_back(k);
    for(size_t sel: Select)if(sel>=Orbs->length())ABORT(Str("cannot select"," ")+Select+"from only"+Orbs->length()+"orbitals");
    Orb.resize(Select.size()); // no size was given, use full function
    gVal.resize(Orb.size());
    //CAUTION: the Coefficients copy constructor broken: does not copy the contiguous storage, assignemt does
    for(size_t k=0;k<Orb.size();k++)gVal[k]=Coefficients(&gdx);

    for (int k=0;idxLeaf!=0;idxLeaf=idxLeaf->nodeNext()){
        idxLeaf->bottomExpand();
        for(const Index *idxBottom=idxLeaf->firstLeaf();idxBottom!=0;idxBottom=idxBottom->rightSibling(),k++){
            std::vector<double>arg,point;
            point=gridPoint(idxBottom);
            arg=toFunctionArg(fromIdx(point));
            if(not Orbs->validArguments(arg))
                Orbs->abortInvalidArguments(tools::str(arg)+" converted from "
                                            +tools::str(gridPoint(idxBottom)));
            double ratio=factorFunc(arg);
            ratio=ratio==0.? 1. : factorIdx(point)/ratio; //quadratur point=0 may be included, but is not needed here...
            val=(*Orbs)(arg);
            for(size_t l=0;l<gVal.size();l++)gVal[l].orderedData()[k]=val[Select[l]]*ratio;
        }
        idxLeaf->bottomUnexpand();
    }

    Coefficients tmp(Idx);
    Plot plt(Idx,ReadInput::main);
    OperatorMap map(Idx,&gdx); // map from grid to original (w/o inverse overlap)
    for(size_t k=0;k<Orb.size();k++){
        map.apply(1.,gVal[k],0.,tmp);
        tmp.makeContinuous();
        *Orb[k]=tmp;
        purge(*Orb[k],1.e-20);
    }
}

// get Index reduced to element next after Isub on a FE axis
// return first if Isub is =0
// return 0 if no further element
// return copy, if not FE axis is found
// abort if multiple finite element axes
static const Index* nextElement(const Index* Idx, const Index* Isub){
    // remove elements except for first above Isub (use first, if Isub==0)
    // return 0 if none above

    if(Isub==Idx)return 0; // already all

    // full copy
    IndexNew* res(new IndexNew(*Idx));

    // find continuity level
    Index*feLev=res;
    while(feLev and feLev->depthOfDuplicate()==Index::npos){
        feLev=feLev->descend();
    }
    if(feLev==0)return res;

    const Index* nxt=feLev->descend();
    while(nxt and nxt->depthOfDuplicate()==Index::npos)nxt=nxt->descend();
    if(nxt)DEVABORT("can only be used on single FE axes, found: "+feLev->axisName()+" and "+nxt->axisName());

    // depth of element below child of feLev
    int childToElement=feLev->depthOfDuplicate()-feLev->depth()-1;
    // upper boundary of previous Isub
    double upBound=Isub==0?-DBL_MAX:Isub->descend(feLev->depth()+childToElement+1)->basis()->integrable()->upBound();

    // go through all nodes on this level
    for(Index* f=feLev;f!=0;f=f->nodeRight()){
        // select first branch above upBound
        int sub(0);
        while(Isub!=0 and sub<int(f->childSize()) and
              f->child(sub)->descend(childToElement)->basis()->integrable()->lowBound()<upBound)sub++;
        if(sub==int(f->childSize())){
            delete res;
            return 0;
        }

        for(size_t k=f->childSize();k>0;k--)
            if(int(k-1)!=sub)f->childErase(k-1);

        f->setBasis(BasisAbstract::factory(BasisSub::strDefinition(f->basis(),{sub})));
    }
    delete Isub;
    res->sizeCompute();
    return res;
}


static void addSubbranch(std::vector<Coefficients*> C, std::vector<Coefficients*> B){
    // add subtrees B into C where bases match
    if(C.front()->idx()->treeEquivalent(B.front()->idx(),{})){
        // equivalent trees - add whole
        for(size_t k=0;k<B.size();k++)*C[k]+=*B[k];
    }
    else{
        if(B.front()->isLeaf())DEVABORT("B does not seem to be subbrach of C");
        const Index*cI=C.front()->idx(),*bI=B.front()->idx();
        if(not (*BasisSub::superBas(cI->basis())==*BasisSub::superBas(bI->basis())))
            DEVABORT("B is not sub-tree of C");
        std::vector<int> cSub=BasisSub::subset(cI->basis());
        std::vector<int> bSub=BasisSub::subset(bI->basis());
        for(size_t l=0;l<bI->childSize();l++){
            size_t kMatch=std::find(cSub.begin(),cSub.end(),bSub[l])-cSub.begin();
            if(kMatch==cSub.size())DEVABORT("no branch in C matching B");
            // vectors of matching children

            std::vector<Coefficients*>cK,bL;
            for(auto e: C)cK.push_back(e->child(kMatch));
            for(auto e: B)bL.push_back(e->child(l));
            addSubbranch(cK,bL);
        }
    }
}

// return orbitals wrt basis with quadrature on fixed quadratur rule
// NOTE: can probably be accelerated quite a bit by avoiding going throug OperatorMap
Eigen::MatrixXcd BasisOrbitalNumerical::intTreeProduct(const Index* Idx, const VectorValuedFunction * Orbs,std::vector<int> Select,
                                                       std::vector<std::vector<double>> Grids,std::vector<std::vector<double>> Weigs){

    CoordinateMap fromIdx=CoordinateTrans::toCartesian(Idx->coordinates());
    IntegrationFactor factorIdx=CoordinateTrans::integrationFactor(Idx->coordinates());
    CoordinateMap toFunctionArg=CoordinateTrans::fromCartesian(Orbs->coordinates());
    IntegrationFactor factorFunc=CoordinateTrans::integrationFactor(Orbs->coordinates());

    // get grid Index for given quadrature rule
    std::vector<std::string> axes=tools::splitString(Idx->coordinates(),'.');
    IndexGrid gdx(Idx,axes,Grids,Weigs);

    // set up the vectors
    std::vector<std::complex<double> > val(Orbs->length());
    if(Select.size()==0)for(size_t k=0;k<val.size();k++)Select.push_back(k);
    for(size_t sel: Select)if(sel>=val.size())ABORT(Str("cannot select"," ")+Select+"from only"+val.size()+"orbitals");

    //CAUTION: the Coefficients copy constructor broken: does not copy the contiguous storage, assignemt does
    std::vector<Coefficients> gVal;
    gVal.resize(Select.size());
    for(size_t k=0;k<Select.size();k++)gVal[k]=Coefficients(&gdx);

    int k=0;
    for (const Index * idxLeaf=gdx.firstLeaf();idxLeaf!=0;idxLeaf=idxLeaf->nodeNext()){
        idxLeaf->bottomExpand();
        for(const Index *idxBottom=idxLeaf->firstLeaf();idxBottom!=0;idxBottom=idxBottom->rightSibling(),k++){
            std::vector<double>arg,point;
            point=gridPoint(idxBottom);
            arg=toFunctionArg(fromIdx(point));
            if(not Orbs->validArguments(arg))
                Orbs->abortInvalidArguments(tools::str(arg)+" converted from "
                                            +tools::str(gridPoint(idxBottom)));
            double ratio=factorFunc(arg);
            ratio=ratio==0.? 1. : factorIdx(point)/ratio; //quadratur point=0 may be included, but is not needed here...

            val=(*Orbs)(arg);

            for(size_t l=0;l<gVal.size();l++)gVal[l].orderedData()[k]=val[Select[l]]*ratio;
        }
        idxLeaf->bottomUnexpand();
    }

    // this we may move into the error criterion
    OperatorMap map(Idx,&gdx);
    Coefficients tmp(Idx);
    Eigen::MatrixXcd res(tmp.size(),Select.size());
    for(size_t k=0;k<Select.size();k++){
        map.apply(1.,gVal[k],0,tmp);
        res.col(k)=Eigen::Map<Eigen::VectorXcd>(tmp.data(),tmp.size());
    }
    return res;
}

class IntParameters
{
    double _epsAbs,_epsRel;
    size_t _maxLevel;
    std::vector<std::vector<double>> _grid,_weig;
public:
    /// \brief Integrate parameters for integration
    /// \param NQuad     number of quadrature points
    /// \param AccRel    desired relative accuracy
    /// \param AccAbs    desired absolute accuracy
    IntParameters(std::vector<size_t> NQuad, double AccRel=1.e-12, double AccAbs=1.e-12,size_t MaxLevel=20)
        :_epsAbs(AccAbs),_epsRel(AccRel),_maxLevel(MaxLevel){
        // n-point Gauss-Legendre Quadratures on [0,1]
        for(auto n: NQuad){
            OrthogonalLegendre l;
            _grid.push_back({});
            _weig.push_back({});
            l.quadratureGauss(n,_grid.back(),_weig.back());
            for(auto &g: _grid.back())g=0.5*(g+1.);
            for(auto &w: _weig.back())w*=0.5;
        }
    }
    std::vector<std::vector<double>> grid() const {return _grid;};
    std::vector<std::vector<double>> weig() const {return _weig;};
    size_t maxLevel() const {return _maxLevel;};
    double epsAbs() const {return _epsAbs;}
    double epsRel() const {return _epsRel;}
};

// scale quadature Grids/Weigs to current volume
static void scaledGridWeig(const BasisOrbitalNumerical::IntParameters &IntPars,const std::vector<std::vector<double>> &Vol,std::vector<std::vector<double>> &Grids,std::vector<std::vector<double>>& Weigs){
    std::vector<std::vector<double> > sub;
    std::vector<Eigen::MatrixXcd> subInts;

    // scale grids to subvolume
    Weigs=IntPars.weig();
    Grids=IntPars.grid();
    for(size_t l=0;l<Vol.size();l++){
        for(auto &w: Weigs[l])w=w*(Vol[l][1]-Vol[l][0]);
        for(auto &g: Grids[l])g=g*(Vol[l][1]-Vol[l][0])+Vol[l][0];
    }
}

// this, unfortunately, somewhat duplicates Integrate::Recursive
// recursively compute integral to pre-defined accuracy by binary splitting of volume
Eigen::MatrixXcd BasisOrbitalNumerical::intTreeRecursive(const Index *Idx, const VectorValuedFunction *Orbs, std::vector<int> Select,
                                                         const IntParameters &IntPars, std::vector<std::vector<double>> Vol,size_t Level, Eigen::MatrixXcd Previous)
{

    if(IntPars.epsAbs()<1e-7 or IntPars.epsRel()<1e-7){
        PrintOutput::DEVwarning(Sstr+"requesting very high integration relative and absolut accuracies "
                                +IntPars.epsAbs()+IntPars.epsRel()+"- may not converge",1);
    }

    std::vector<std::vector<double>> grids,weigs;
    if(Level==0){
        // initialize upon entry
        scaledGridWeig(IntPars,Vol,grids,weigs);
        Previous=BasisOrbitalNumerical::intTreeProduct(Idx,Orbs,Select,grids,weigs);
    }
    Eigen::MatrixXcd integral(Previous);

    // systematically run through sub-volumes
    integral*=0.;
    std::vector<std::vector<double> > sub;
    std::vector<Eigen::MatrixXcd> subInts;
    while (Integrate::NextSubvolume(Vol,sub).size()>0){
        scaledGridWeig(IntPars,sub,grids,weigs);
        subInts.push_back(BasisOrbitalNumerical::intTreeProduct(Idx,Orbs,Select,grids,weigs));
        integral+=subInts.back();
    }
    // return if absolute and relative error estimates are met
    if((integral-Previous).lpNorm<2>()<std::max(IntPars.epsAbs(),IntPars.epsRel()*integral.lpNorm<2>())){
        return integral;
    }

    if(Level>IntPars.maxLevel()){
        Integrate::showVol(Vol);
        PrintOutput::message(Sstr+integral.lpNorm<2>()+" - "+Previous.lpNorm<2>()+" >? "+IntPars.epsAbs());
        if(Integrate::measureVol(Vol)*1.e5<integral.lpNorm<2>())PrintOutput::message("!!! integrand may be singular !!!");
        ABORT("exceeded maximal recursion depth="+tools::str(IntPars.maxLevel())
              +"\ncheck whether integrand has singularities or non-analyticities"
              +"\nrelax accuracy or increase order");
    }

    // get accurate integrals on each subvolume
    integral=Eigen::MatrixXcd::Zero(integral.rows(),integral.cols());
    for(unsigned int k=0;k<subInts.size();k++){
        // adjust tolerance according to importance of sub-volume
        integral+=intTreeRecursive(Idx,Orbs,Select,IntPars,Integrate::NextSubvolume(Vol,sub),Level+1,subInts[k]);
    }
    return integral;
}

// insert <orb[i]|VectorValuedFunction> into Coefficients Orb
// orb[i] is assumed to be expanded into coordinates that map to VectorValuedFunction
void BasisOrbitalNumerical::mapOrbitals(const Index* Idx, const VectorValuedFunction * Orbs, std::vector<int> Select,std::vector<Coefficients*> Orb)
{
    if(not Idx->basis()->orbital())ABORT("not an orbital basis");
    const Index* idx=Idx->basis()->orbital()->orbital(0)->idx();
    std::vector<Coefficients*>orbx;
    for(size_t k=0;k<Orb.size();k++)orbx.push_back(new Coefficients(idx));
    mapTreeProductExact(orbx[0]->idx(),Orbs,Select,orbx);
    for(size_t k=0;k<Orb.size();k++){
        Orb[k]->data()[k]=orbx[k]->innerProduct(Idx->basis()->orbital()->orbital(k));
    }
    for(auto &o: orbx)delete o;
}

// return orbatils at sub-branch where idx() matches Idx
std::vector<Coefficients*> orbsAtIdx(std::vector<Coefficients>& InOrbs, const Index* Idx){
    std::vector<Coefficients*> out;
    for(Coefficients &i: InOrbs){
        out.push_back(&i);
        for(;out.back()!=0 and out.back()->idx()!=Idx;out.back()=out.back()->nodeNext());
        if(out.back()==0)DEVABORT("index not found in orbital");
    }
    return out;
}

// create coefficients for orthonormal numerical orbitals
// driver covers various index structures
void BasisOrbitalNumerical::generateOrbitals(const Index *Idx){
    const Index* idx=Idx;
    if(idx==0){
        if(_refName=="auto"){
            idx=referenceIndex["main"];
            while(idx->hierarchy().find("&")!=std::string::npos)idx=idx->child(1);
            _refName="auto("+idx->hierarchy()+")";
            referenceIndex[_refName]=idx;
        }
        if(referenceIndex.find(_refName)==referenceIndex.end())
            ABORT("reference index not defined: "+_refName+"\navailable: "+tools::listMapKeys(referenceIndex));
        idx=referenceIndex[_refName];
    }

    // get the function
    const VectorValuedFunction * orbs=orbitalFunctions();
    if(orbs==0)DEVABORT("not defined "+_funcDef+", do VectorValuedFunction::add(...) first\nAvailable:\n"
                        +VectorValuedFunction::list());

    if(not idx->inverseOverlap()){
        const_cast<Index*>(idx)->localOverlapAndInverse(0,0);
        Inverse::factory(idx);
    }

    // get selection size and create empty orbitals
    (*orbs)({0.,0.,0.});
    if(_select.size()==0)
        for(size_t k=0;k<orbs->length();k++)_select.push_back(k);
    for(size_t sel: _select)if(sel>=orbs->length())ABORT(Str("cannot select"," ")+_select+"from only"+orbs->length()+"orbitals");
    _orb.assign(_select.size(),Coefficients(idx,0.));

    // descend index transformable coordinates are found
    std::vector<std::vector<Coefficients*>> minusC;
    for(const Index* tdx=idx;tdx!=0;tdx=tdx->nodeNext()){
        if(tdx->basis()->ndim()){
            DEVABORT("cannot use Ndim at present");
            minusC.push_back(std::vector<Coefficients*>());
            Coefficients tmp(tdx);
            std::vector<Coefficients*> tdxC(orbsAtIdx(_orb,tdx));
            mapNdim(tdx,orbs,_select,tdxC);
            for(size_t k=0;k<tdxC.size();k++){
                if(tdxC[k]->norm()>1.e-20){
                    minusC.back().push_back(tdxC[k]);
                    tmp=*tdxC[k];
                    tdx->inverseOverlap()->apply(1.,tmp,0,*minusC.back().back());
                }
                else tdxC[k]->setToZero();
            }
        }
    }

    for(const Index* tdx=idx;tdx!=0;tdx=tdx->nodeNext()){
        if(tdx->basis()->orbital()){
            std::vector<Coefficients*> tdxC(orbsAtIdx(_orb,tdx));
            mapOrbitals(tdx,orbs,_select,tdxC);
        }
    }

    for(const Index* tdx=idx;tdx!=0;tdx=tdx->nodeNext()){
        if(CoordinateTrans::toCartesian(tdx->coordinates())!=CoordinateTrans::undefined){
            Coefficients tmp(tdx);
            std::vector<Coefficients*>tdxC(orbsAtIdx(_orb,tdx));
            if(integrationMethod=="singleQuadrature")
                mapTreeProduct(tdx,orbs,_select,tdxC);
            else if(integrationMethod=="volumeSplit")
                mapTreeProductExact(tdx,orbs,_select,tdxC);
            else
                DEVABORT("unknown integration method: "+integrationMethod);
            break;
        }
    }
    for(auto &o: _orb)purge(o);

    Coefficients tmp(idx);
    for(size_t k=0;k<_orb.size();k++){
        tmp=_orb[k];
        idx->inverseOverlap()->apply(1.,tmp,0.,_orb[k]);
    }
    _overOnBasis=overlap();

    std::vector<double>nrm=norms();
    for(int k=nrm.size()-1;k>=0;k--){
        if(nrm[k]<0.99){
            if(nrm[k]<1.e-12){
                nrm.erase(nrm.begin()+k);
                _orb.erase(_orb.begin()+k);
                ABORT(Str("empty orbital no."," ")+k+"- reference basis may be too small, e.g. angular momenta");
            }
        }
    }
    orthonormalize(true);

    // warn about complex scaled region
    for(size_t k=0;k<_orb.size();k++){
        double scaledExp=abs( idx->overlap()->matrixElement(_orb[k],_orb[k])
                              -idx->overlap()->matrixElementUnscaled(_orb[k],_orb[k]));
        if(scaledExp>1.e-9){
            PrintOutput::warning(Str("OrbitalNumerical"," ")+k+"has expectation value"+scaledExp+"in scaled region - may not be meaningful");
        }
    }
    plot();
}

bool BasisOrbitalNumerical::operator==(const BasisAbstract & Other) const {
    const BasisOrbitalNumerical* other=dynamic_cast<const BasisOrbitalNumerical*>(&Other);
    if(not other)return false;
    return size()==other->size()
            and _refName==other->_refName
            and _funcDef==other->_funcDef
            and _select==other->_select;
}

void BasisOrbitalNumerical::print(std::string Title) const{

    // NOTE: use operator setup sparingly (memory in absence of symmetry)
    std::vector<std::vector<double>> prop;
    std::vector<std::string> title;
    for(auto p: _printProperties)
        if(_printList.find(p.first)!=std::string::npos){
            prop.push_back(expectationValues(p.second));
            title.push_back(p.first);
        }

    if(Title!="")PrintOutput::title(Title);

    PrintOutput::newRow();
    PrintOutput::rowItem("No.");
    PrintOutput::rowItem("{orig}");
    for(auto t: title)PrintOutput::rowItem("<"+t+">");
    for(size_t k=0;k<_orb.size();k++){
        PrintOutput::newRow();
        PrintOutput::rowItem(int(k));
        PrintOutput::rowItem(_select[k]);
        for(size_t l=0;l<prop.size();l++)PrintOutput::rowItem(prop[l][k]);
    }
    PrintOutput::paragraph();
    PrintOutput::subTitle("Overlaps after expansion (before orthonormalization)");
    PrintOutput::matrix(_overOnBasis);
    PrintOutput::paragraph();

}

const VectorValuedFunction* BasisOrbitalNumerical::orbitalFunctions() const {
    return VectorValuedFunction::get(_funcDef.substr(0,_funcDef.find("{")));
}

// true if volume is finite, i.e. none of the boundaries is beyond +-DBL_MAX/2
static bool finiteVolume(std::vector<std::vector<double>> Vol){
    for(auto v: Vol)
        if(v[0]<-DBL_MAX/2 or v[1]>DBL_MAX/2)return false;
    return true;
}

// recursively compute integrals to preset accuracy
void BasisOrbitalNumerical::mapTreeProductExact(const Index *Idx, const VectorValuedFunction *Orbs, std::vector<int> Select, std::vector<Coefficients *> Orb){
    std::vector<std::string> ax=tools::splitString(Idx->coordinates(),'.');


    // loop through all subtrees for single elements
    const Index* isub;
    for(isub=nextElement(Idx,0);isub!=0;isub=nextElement(Idx,isub)){
        auto v(volume(isub));
        if(finiteVolume(v)){
            Eigen::MatrixXcd mat=intTreeRecursive(isub,Orbs,Select,IntParameters(orders(isub,points),integrationAccuracy,integrationAccuracy),v);
            Coefficients csub(isub);
            for(int j=0;j<mat.cols();j++){
                for(size_t k=0;k<csub.size();k++)csub.data()[k]=mat(k,j);
                addSubbranch({Orb[j]},{&csub});
            }
        } else{
            std::vector<Coefficients*> vsub;
            for(size_t k=0;k<Select.size();k++)vsub.push_back(new Coefficients(isub));
            mapTreeProduct(isub,Orbs,Select,vsub);
            addSubbranch(Orb,vsub);
            for(auto &o: vsub)delete o;
        }
    }
    delete isub;
}
























