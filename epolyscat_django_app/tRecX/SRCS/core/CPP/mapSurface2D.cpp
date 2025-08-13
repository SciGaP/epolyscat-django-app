#include "../mapSurface2D.h"


#include "index.h"
#include "indexGrid.h"
#include "indexNew.h"
#include "algebra.h"

#include "basisGrid.h"
#include "basisExpIm.h"
#include "operatorFloor.h"
#include "coordinateTrans.h"
#include "discretizationSurface.h"
#include "discretizationGrid.h"
#include "operatorMap.h"
#include "inverse.h"
#include "coorSystem.h"
#include "basisSub.h"

#include "eigenTools.h"

// should go to stringTools
static Str toStr(const std::map<std::string,double> Map){
    Str res("","");
    for(const auto & m: Map)res=res+m.first+": "+m.second+", ";
    res.resize(res.size()-2);
    return res;
}

Index* indexSphericalSurface(const Index* Idx,int LMax, double R){
    // construct a spherical surface index from (parabolic) index
    // - copy structure up to and including Phi-level
    // - insert BasisSub-level, this is expected by the tSurff routines, as it arises as
    //   leftover from the FE-DVR discretization of the Rn-axis
    // - replace branches with desired Eta branches (AssociatedLegendre basis up to LMax)
    // - attach a ValDer floor level to each Eta branch

    Index* res=new IndexNew();
    res->nodeCopy(Idx,false);
    // descend to Phi-level
    if(res->axisName()!="Phi")
        for(size_t k=0;k<Idx->childSize();k++)
            res->childAdd(indexSphericalSurface(Idx->child(k),LMax,R));

    // attach spherical basis
    for(size_t k=0;k<Idx->childSize();k++){

        int m=dynamic_cast<const BasisTrigon*>(Idx->basis())->mValueAbs(k);
        if(m>LMax)ABORT("lMax too small, smaller than m-values in original basis ");

        // attach a surface level
        res->childAdd(new Index());
        Index* surf=res->childBack();
        surf->setAxisName("surfRn");

        // surface expects "BasisSub" here
        surf->setBasis(BasisAbstract::factory("BasisSub: 1 of "+BasisAbstract::factory("Vector: 2")->strDefinition()));

        Index* iEta=new Index();
        surf->childAdd(iEta);

        iEta->setAxisName("Eta");
        iEta->setBasis(BasisAbstract::factory(Sstr+"AssociatedLegendre: "+m+","+(LMax-m)));
        for(size_t l=0;l<iEta->basis()->size();l++){
            iEta->childAdd(new Index());
            iEta->childBack()->setAxisName("ValDerRn");
            iEta->childBack()->setBasis(BasisAbstract::factory(Sstr+"Grid: "+R+","+R));
            iEta->childBack()->setFloor();
        }
    }
    res->sizeCompute();
    return res;
}

MapSurface2D::~MapSurface2D(){delete _gIndex;}

MapSurface2D::MapSurface2D(const Index* Idx,std::string SurfDef)
    :OperatorAbstract("mapSurface"+SurfDef,0,Idx){

    std::vector<std::shared_ptr<Algebra>> algs;
    std::vector<std::string> parts=tools::splitString(tools::stringInBetween(SurfDef,"Spherical[","]",true),',');
    for(auto s: parts){
        algs.push_back(std::shared_ptr<Algebra>(new Algebra(s)));
        if(not algs.back()->isAlgebra())ABORT("strings between [...] in "+SurfDef+" are no proper algebras");
    }

    std::vector<std::string> axGrid,surfCoor;
    if(SurfDef.find("Spherical[")==0){
        surfCoor={"Phi","Eta","Rn"};
        axGrid={"Eta"};

        // construct surface index: first algs[0] is Radius, algs[1] is Lmax
        iIndex=indexSphericalSurface(Idx,int(algs[1]->val(0.).real()),algs[0]->val(0.).real());
        name="surface_Rn";
    }
    else
        ABORT("for now, need Spherical[lmax,r], got:"+SurfDef);

    std::string inSys,refSys;
    std::vector<std::string> inCoor=tools::splitString(Idx->coordinates(),'.');
    for(auto c: inCoor  )if(std::find(surfCoor.begin(),surfCoor.end(),c)==surfCoor.end()){ inSys+=c+".";inSys.pop_back();}
    for(auto c: surfCoor)if(std::find(inCoor.begin(),    inCoor.end(),c)==  inCoor.end()){refSys+=c+".";refSys.pop_back();}

    //Note: CoorSystem pointers are handled internally, do NOT delete (not nice...)
    CoorSystem* coorSys(CoorSystem::factory(refSys,inSys));
    _gIndex=new IndexGrid(iIndex,axGrid,std::vector<std::vector<double>>(1),{});
    _toSurf.reset(new OperatorMap(iIndex,_gIndex));
    _toGrid.reset(new MapSurfaceGrid(_gIndex,Idx,*coorSys));
}

MapSurface2D::MapSurfaceGrid::MapSurfaceGrid(const Index* ISurf, const Index* Idx, const CoorSystem& CoorSys)
    :OperatorTree("Map "+ISurf->coordinates()+"<-"+Idx->coordinates(),ISurf,Idx)
{
    // - descend while both index structures match (i.e. identity on wrt to these factors)
    // - descend through possible FE levels on each side
    // - once the above is exhausted, descend in surface index until ValDer level is met
    // - on ValDer level
    //   + determine (Eta,R) of present surface  point
    //   + compute PXi(Eta,R), PEta(Eta,R) coordinates, and dPXi/dR, dPEta/dR
    //   + temporarily move floor level into Index leafs (single coefficients)  in ISurf and Idx
    //   + compute map beteween current branches (different constructor)
    //   + reset the floors to ValDer and FE levels

    if(idx()->nodeEquivalent(jdx())){
        // identical on this level - no transform
        for(size_t k=0;k<idx()->childSize();k++){
            childAdd(new MapSurfaceGrid(idx()->child(k),jdx()->child(k),CoorSys));
        }
    }
    else if(jdx()->isFem()){
        // descend into coordinate section
        for(size_t k=0;k<jdx()->childSize();k++)
            childAdd(new MapSurfaceGrid(idx(),jdx()->child(k),CoorSys));
    }
    else if(idx()->isFem() or idx()->axisName().find("surf")==0){
        // descend into coordinate section
        for(size_t k=0;k<idx()->childSize();k++)
            childAdd(new MapSurfaceGrid(idx()->child(k),jdx(),CoorSys));
    }

    else if (ISurf->axisName().find("ValDer")==0){
        // arrived

        // get target coordinate values
        std::vector<std::string> AxisOut=tools::splitString(CoorSys.name(),'.');
        std::vector<double> tarCoors;
        for(const Index* idx=ISurf;idx->parent()!=0;idx=idx->parent()){
            if(std::find(AxisOut.begin(),AxisOut.end(),idx->parent()->axisName())!=AxisOut.end()
                    and idx->parent()->basis()->grid())
                tarCoors.push_back(idx->parent()->basis()->grid()->mesh()[idx->nSibling()]);
        }
        tarCoors.push_back(ISurf->basis()->grid()->mesh()[0]);
        if(tarCoors.size()!=AxisOut.size())DEVABORT("did not find coordinate values");

        // compute original coordinate values
        std::vector<std::string> AxisIn=tools::splitString(CoorSys.refSystem(),'.');
        std::vector<double> coors=CoorSys.toRef(tarCoors);
        std::map<std::string,double> inCoors;
        for(size_t k=0;k<AxisIn.size();k++)inCoors[AxisIn[k]]=coors[k];

        // compute transform for gradient
        Eigen::MatrixXd jac=CoorSys.jacRefdCoor(tarCoors); // d[inCoor]/d[tarCoors](tarCoors), column-index is tarCoor
        std::vector<std::complex<double>> jacD;
        double fac(1.);
        if(CoorSys.refSystem()=="PXi.PEta" and CoorSys.name()=="Eta.Rn"){
            // tRecX uses r dPsi/dr for surface derivative values from Rn
            fac=0.5*(inCoors["PXi"]+inCoors["PEta"]);
            jacD.push_back(fac*jac(0,1));
            jacD.push_back(fac*jac(1,1));
        }
        else
            DEVABORT("not implemented for "+CoorSys.name()+" <- "+CoorSys.refSystem());

        // expand and lower floors
        ISurf->bottomExpandAll();
        Idx->bottomExpandAll();
        size_t iDepth=idx()->depth();
        size_t jDepth=jdx()->depth();
        size_t iFloor=iDepth+ISurf->heightAboveFloor();
        size_t jFloor=jDepth+Idx->heightAboveFloor();
        idx()->resetFloor(iDepth+idx()->height());
        jdx()->resetFloor(jDepth+jdx()->height());

        childAdd(new MapSurfaceGrid(iIndex->child(0),jdx(),inCoors,fac,{},0));
        childAdd(new MapSurfaceGrid(iIndex->child(1),jdx(),inCoors,fac,jacD,0));

        // remove empty branches (point may not fall into current axis sections)
        for(size_t k=childSize();k>0;k--)
            if(child(k-1)->isLeaf() and not child(k-1)->floor())childErase(k-1);

        // re-set floor level to input
        idx()->sizeCompute();
        jdx()->sizeCompute();
        reFloor(iFloor,jFloor);

        // restore entry indices
        const_cast<Index*>(idx())->bottomUnexpandAll();
        const_cast<Index*>(jdx())->bottomUnexpandAll();
    }

    else {
        // if not on ValDer level yet, descend in ISurf
        for(size_t k=0;k<ISurf->childSize();k++){
            childAdd(new MapSurfaceGrid(ISurf->child(k),Idx,CoorSys));
        }
    }
}

MapSurface2D::MapSurfaceGrid::MapSurfaceGrid(const Index* ISurf, const Index* Idx, std::map<std::string, double> Coors, std::complex<double> Val,
                                             std::vector<std::complex<double>>Ders, size_t Pos)
    :OperatorTree("Map "+ISurf->coordinates()+"<-"+Idx->coordinates(),ISurf,Idx)
{
    // descend through rhs index
    // - descent into all FEs
    // - if is BasisIntegrable (must belong to coordinates),
    //   and point in support of basis, for all rhs basis functions k:
    //   + multiply value by current factor value                      Val <- Val*val[k]
    //   + multiply matching derivatives by current factor derivative  Der[m] <- Der[m]*der[m]
    //   + multiply non-matching derivatives by current factor value:  Der[k] <- Der[k]*val[k]
    // - on 1x1 floor level, place value or derivative into trivial floor Operator1X1

    if(jdx()->isFem()){
        // descend into coordinate section
        for(size_t k=0;k<jdx()->size();k++)
            childAdd(new MapSurfaceGrid(idx(),jdx()->child(k),Coors,Val,Ders,Pos));
    }
    else if(Idx->basis()->integrable()) {
        if(not Coors.count(Idx->axisName()))
            DEVABORT("axis"+Idx->axisName()+"is not among coordinates: "+tools::listMapKeys(Coors));

        //        if(tools::doubleInside(Coors[Idx->axisName()],Idx->basis()->integrable()->lowBound(),Idx->basis()->integrable()->upBound())){
        if(Idx->basis()->integrable()->lowBound()<Coors[Idx->axisName()] and Coors[Idx->axisName()]<=Idx->basis()->integrable()->upBound()){
            if(Idx->basis()->integrable()->eta()!=1.){
                ABORT("surface is in complex scaled region: "+toStr(Coors));
            }
            std::vector<std::complex<double>> val,der;
            jdx()->basis()->integrable()->valDer({Coors[Idx->axisName()]},val,der);
            for(size_t k=0;k<jdx()->childSize();k++){
                std::vector<std::complex<double>> ders(Ders);
                for(size_t p=0;p<ders.size();p++)ders[p]*=p==Pos?der[k]:val[k];
                childAdd(new MapSurfaceGrid(idx(),jdx()->child(k),Coors,Val*val[k],ders,Pos+1));
            }
        }
    }
    else if(jdx()->isLeaf()){
        // compute value or derivative
        if(Ders.size()>0){
            Val=0;
            for(const auto &d: Ders)Val+=d;
        }
        // add as trivial floor
        oFloor=new Operator1X1(Val);
    }
    else {
        Sstr+jdx()->strNode()+jdx()->isBottom()+jdx()->isLeaf()+Sendl;
        DEVABORT("this is not supposed to happen");
    }

    // remove empty branches
    for(size_t k=childSize();k>0;k--)
        if(child(k-1)->isLeaf() and not child(k-1)->floor())childErase(k-1);
}

void MapSurface2D::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    // map to grid
    _toGrid->apply(A,Vec,0,*_toGrid->tempLHS());

    // further expand to full angular grid (or we may do this from the beginning)
    // apply gauge factors (use MapGauge)

    // map to original basis
    Coefficients* y=_inverse?_toSurf->tempLHS():&Y;
    _toSurf->apply(1,*_toGrid->tempLHS(),B,*y);

    // apply invers (if any)
    if(_inverse)_inverse->apply(1.,*y,0,Y);
}
