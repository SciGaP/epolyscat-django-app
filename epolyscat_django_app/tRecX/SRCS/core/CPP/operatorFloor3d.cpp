// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorFloor3d.h"

#include "labelled.h"
#include "readInput.h"
#include "index.h"
#include "qtEigenDense.h"
#include "printOutput.h"
#include "operatorNdim.h"
#include "coordinateTrans.h"
#include "inverseDvr.h"
#include "basisDvr.h"
#include "algebra.h"
#include "eigenTools.h"

using namespace std;

// all static variables need to re-appear here
unsigned int OperatorFloor3d::mExpansionSize=INT_MAX;
unsigned int OperatorFloor3d::lExpansionSize=INT_MAX;
unsigned int OperatorFloor3d::mLev;
unsigned int OperatorFloor3d::lLev;
unsigned int OperatorFloor3d::femLev;

vector<double> OperatorFloor3d::_potShift;
vector<vector<complex<double> > > OperatorFloor3d::potIJ; // potential for IJ-block at quadrature grid as well as radius
vector<vector<complex<double> > > OperatorFloor3d::phiBas; // phi-basis at quadrature grid
vector<vector<vector<complex<double> > > > OperatorFloor3d::etaBas; // eta-basis at quadrature grid
vector<vector<complex<double> > > OperatorFloor3d::invOvr; // Inverse overlap for radial grid


vector<int> OperatorFloor3d::quad;
potential3d OperatorFloor3d::potOrigin0=0;
string OperatorFloor3d::potDef;
double OperatorFloor3d::chargeC=0.;
double OperatorFloor3d::screenHydrogen=100.;
double OperatorFloor3d::screenCarbon=0.265;
double OperatorFloor3d::methaneSize=2.042794;
std::string OperatorFloor3d::methaneAlign="z";

complex<double> OperatorFloor3d::pot3d(const std::vector<std::complex<double> > & Coor){
    if(_potShift.size()==0)return potOrigin0(Coor);
    else return potOrigin0(CoordinateTrans::shiftPolar(Coor,_potShift));
}

static std::unique_ptr<Algebra> _potRadial;
complex<double> OperatorFloor3d::radial(const std::vector<std::complex<double> > & Coor){
    return _potRadial->val(Coor[2]);
}
complex<double> OperatorFloor3d::methane(const std::vector<std::complex<double> > & Coor){
    double Phi=Coor[0].real();
    double Eta=Coor[1].real();
    complex<double> R=Coor[2];

    //0. Conversion to cartesian coordinates
    double sintheta=sqrt(1-Eta*Eta);
    Eigen::Vector3cd X(R*sintheta*cos(Phi),R*sintheta*sin(Phi),R*Eta);

    //1. Settings
    double R0=methaneSize; //radial distance in a.u.
    double phi1=acos(-1)*2./3.; //the angle

    //2. Define the radial coordinates
    vector<Eigen::Vector3cd> RR;

    if(methaneAlign=="x"){
        RR.push_back(Eigen::Vector3cd(R0         ,0.                        ,0.                        ));
        RR.push_back(Eigen::Vector3cd(R0*(-1./3.),R0*sqrt(8./9.)*cos(  phi1),R0*sqrt(8./9.)*sin(  phi1)));
        RR.push_back(Eigen::Vector3cd(R0*(-1./3.),R0*sqrt(8./9.)*cos(2*phi1),R0*sqrt(8./9.)*sin(2*phi1)));
        RR.push_back(Eigen::Vector3cd(R0*(-1./3.),R0*sqrt(8./9.)*cos(3*phi1),R0*sqrt(8./9.)*sin(3*phi1))); // Algined in xy plane
    }else if(methaneAlign=="y"){
        RR.push_back(Eigen::Vector3cd(0.                        ,R0         ,0.                        ));
        RR.push_back(Eigen::Vector3cd(R0*sqrt(8./9.)*sin(  phi1),R0*(-1./3.),R0*sqrt(8./9.)*cos(  phi1)));
        RR.push_back(Eigen::Vector3cd(R0*sqrt(8./9.)*sin(2*phi1),R0*(-1./3.),R0*sqrt(8./9.)*cos(2*phi1)));
        RR.push_back(Eigen::Vector3cd(R0*sqrt(8./9.)*sin(3*phi1),R0*(-1./3.),R0*sqrt(8./9.)*cos(3*phi1))); // Aligned in yz plane
    }else if(methaneAlign=="z"){
        RR.push_back(Eigen::Vector3cd(0.                        ,0.                        ,R0         ));
        RR.push_back(Eigen::Vector3cd(R0*sqrt(8./9.)*cos(  phi1),R0*sqrt(8./9.)*sin(  phi1),R0*(-1./3.)));
        RR.push_back(Eigen::Vector3cd(R0*sqrt(8./9.)*cos(2*phi1),R0*sqrt(8./9.)*sin(2*phi1),R0*(-1./3.)));
        RR.push_back(Eigen::Vector3cd(R0*sqrt(8./9.)*cos(3*phi1),R0*sqrt(8./9.)*sin(3*phi1),R0*(-1./3.))); // Aligned in zx plane
    }else{
        ABORT("Unknown methaneAlign: "+methaneAlign);
    }

    //3. Calculate the potential
    complex<double> val=-chargeC*(1.-exp(-R/screenCarbon))/R;
    for(int k=0;k<RR.size();k++){
        complex<double> r=(X-RR[k]).transpose()*(X-RR[k]);
        r=sqrt(r);
        val-=(1.-chargeC)*(0.25/r)+chargeC*0.75*exp(-r/screenHydrogen);
    }
    return val;

}

void OperatorFloor3d::print(){
    if(potOrigin0==0)return;

    PrintOutput::title("3d potential");
    PrintOutput::lineItem("name",potDef);
    PrintOutput::newLine();
    PrintOutput::lineItem("quadOrders",tools::str(quad));
    PrintOutput::paragraph();
}


void OperatorFloor3d::read(ReadInput &Inp){
    quad.clear();
    Inp.read("Pot3d","quadOrder",quad,"40 40 0","quadrature order for interpolating the potential, last=0 indicates DVR)");
    if(quad.size()!=3)ABORT("specify single for all or 3 quadOrder's for phi,eta,r; found: "+tools::str(quad))
    if(quad[2]!=0)ABORT("non-dvr radial integration does not work")

    Inp.read("Pot3d","chargeC",chargeC,"1/5","fraction or Coulombic charge on C");
    Inp.read("Pot3d","methaneSize",methaneSize,"2.042794","Length scale of Methane molecule (au)");
    Inp.read("Pot3d","methaneAlign",methaneAlign,"z","Align one symmetry axis along either 'x', 'y' or 'z'");
    Inp.read("Pot3d","screenHydrogen",screenHydrogen,"3.100","screen length for Hydrogens");
    Inp.read("Pot3d","screenCarbon",screenCarbon,"2.","screen length for Carbon");

    Inp.read("Pot3d","potential",potDef,"undefined",": methane,hydrogen,harmOsc,identity,radial[x,y,z,V(r)]...radial around x,y,z");
    Inp.read("Pot3d","origin",_potShift,"","postion of potential origin (give in cartesian coordinates)");

    if(not Inp.found("Pot3d"))return;
    if(potDef.find("radial[")==0){
        potOrigin0=radial;
        if(Inp.found("Pot3d","origin"))ABORT("when using Pot3d:potential=radial[x,y,z,V], do not specify Pot3d:origin")
        std::vector<std::string> arg=tools::splitString(tools::stringInBetween(potDef,"[","]"),',');
        if(arg.size()!=4)ABORT("need radial[x,y,z,V], got: "+potDef)
        _potShift={Algebra::constantValue(arg[0]).real(),Algebra::constantValue(arg[1]).real(),Algebra::constantValue(arg[2]).real()};
        _potRadial.reset(new Algebra(arg[3]));
        if(not _potRadial->isAlgebra())
            ABORT("potential (4th argument) of"+potDef+" is not a well-formed algebra\nerror "+Algebra::failures);
    }
    else if(potDef=="undefined")potOrigin0=undefined;
    else if(potDef=="methane") potOrigin0=methane;
    else if(potDef=="identity")ABORT("Pot3d:potential=identity absolete, use radial[0,0,0,1] instead")
    else if(potDef=="harmOsc")ABORT("Pot3d:potential=harmOsc absolete, use radial[0,0,0,Q*Q/2] instead")
    else if(potDef=="hydrogen")ABORT("Pot3d:potential=hydrogen absolete, use radial[0,0,0,-1/Q] instead")
    else if(potDef=="charB1")potOrigin0=experimental;
    else ABORT("undefined Pot3d: "+potDef)

    if(_potShift==vector<double>(3,0.))_potShift.clear();
    else for(size_t k=0;k<_potShift.size();k++)_potShift[k]*=-1.;

    // make potential available to OperatorNdim
    OperatorNdim::addPotNdim("Pot3d",pot3d);
}

void OperatorFloor3d::setup(const Index *Idxfloor){

    if(potOrigin0==0)ABORT("call OperatorFloor3d::read before setup");

    //1. Get the top level index to make sense of precalculating values for speed
    const Index *iRoot=Idxfloor; // find top of polar
    for(;iRoot and iRoot->coordinates()!="Phi.Eta.Rn";iRoot=iRoot->parent());
    if(not iRoot)
        ABORT("for now only for hierarchy ending in Phi.Eta.Rn but is \""+Idxfloor->root()->hierarchy()+"\"");

    // locate the levels
    mLev=iRoot->axisIndex("Phi")->depth()-iRoot->depth();
    lLev=iRoot->axisIndex("Eta")->depth()-iRoot->axisIndex("Phi")->depth();
    femLev=iRoot->axisIndex("Rn")->depth()-iRoot->depth();
    size_t rLev=iRoot->firstFloor()->depth()-iRoot->depth();
    if(femLev==rLev)femLev=Index::npos;

    int mExpansionSize2=iRoot->childSize();
    int lExpansionSize2=mExpansionSize2-2;

    UseMatrix phiGrid,phiWeig;
    iRoot->descend(mLev)->basis()->integrable()->quadRule(quad[0],phiGrid,phiWeig);
    phiBas=basVal(phiGrid,phiWeig,*iRoot->basis()->integrable(),mExpansionSize2);

    //HACK: computation of quadrature order needs checking
    UseMatrix etaGrid,etaWeig;
    iRoot->descend(lLev)->basis()->integrable()->quadRule(quad[1],etaGrid,etaWeig);

    //4. the eta basis at the quadrature points times sqrt(weight) for different m's: size(nm,lgrid,netaq)
    for(int mm=0;mm<iRoot->descend(mLev)->childSize();mm++){
        Index* lIdx=iRoot->descend(mLev)->child(mm);
        etaBas.push_back(basVal(etaGrid,etaWeig,*lIdx->basis()->integrable(),lExpansionSize2));
    }

    //5. the r basis at the quadrature points for the different FEM elements size(nFEM,nfunc,rgrid)
    potIJ.clear();

    Index* rIdx=iRoot->descend(rLev);
    int nEle=1;
    if(femLev!=Index::npos)nEle=iRoot->descend(femLev)->childSize();
    vector<complex<double> > coor(3);
    for(int n=0; n<nEle; n++)
    {
        if(femLev!=Index::npos)rIdx=iRoot->descend(femLev)->child(n);
        // get an accurate(?) quadrature (input or default=2*order)
        std::vector<double> grid,weig;
        const BasisDVR* dvr=dynamic_cast<const BasisDVR*>(rIdx->basis());
        if(not dvr)DEVABORT("must use DVR basis, is:"+rIdx->basis()->str());
        if(quad[2])DEVABORT("must not specify quadrature for r-coordinate");
        dvr->dvrRule(grid,weig);

        // transform from quadrature grid s to dvr basis
        // trans(r,s) = (<m|n>)^-1 val[m](s)w_s
        UseMatrix rVal,lVal,ovr;

        rVal=dvr->val(grid);
        lVal=rVal.transpose();
        for(int s=0;s<lVal.cols();s++)lVal.col(s)*=weig[s];

        ovr=lVal*rVal;
        lVal=ovr.solve(lVal);

        UseMatrix dvrGrid,dvrWeig;
        if(not dvr)DEVABORT("need dvr basis here");
        dvr->dvrRule(dvrGrid,dvrWeig);

        // store transpose for easier summation in loop
        UseMatrix trans;
        trans=lVal.transpose();
        trans*=dvr->val(dvr->nodes());


        // complex scaled grid points
        // rIdx->basisSet()->def.comSca.coordinates(grid);
        //6. precompute the potFunc at the phi,eta,n,rn grid
        potIJ.push_back(vector<complex<double> >());
        for(size_t p=0; p<phiGrid.size();p++){
            coor[0]=phiGrid(p).complex();
            for(size_t q=0; q<etaGrid.size();q++){
                coor[1]=etaGrid(q).complex();
                // pot a given phi,eta and all scaled radial points
                vector<complex<double> > potr;
                for(size_t r=0;r<dvr->size();r++){
                    coor[2]=grid[dvr->nBeg()+r];
                    potIJ.back().push_back(pot3d(coor));
                }
            }
        }
        for(auto &a: potIJ.back())
            if(std::isnan(a.real()) or std::isnan(a.imag()))
                ABORT("non-number in potential, check for singularites, especially at element boundaries");
    }
}

std::vector<std::vector<std::complex<double> > > OperatorFloor3d::basVal(const UseMatrix & Grid, const UseMatrix & Weig, const BasisIntegrable & Bas, unsigned int PotPoints){

    //1. values of the different basis functions at the quadrature points
    std::vector<std::vector<std::complex<double> > > vals;

    UseMatrix val=Bas.val(Grid);//get the values of the basis at the grid
    vals.resize(Bas.size());//so vals(nbasis,nquadpoints)
    for(size_t m=0;m<vals.size();m++){
        vals[m].resize(Grid.size());
        for(size_t k=0;k<Grid.size();k++)
            vals[m][k]=val(k,m).complex()*sqrt(Weig(k).real());// absorb square-root of weight
    }
    return vals;
}


OperatorFloor3d::OperatorFloor3d(const Index *IIndex, const Index *JIndex)
{
    /* For given IIndex (assuming DVR in k-direction)
     * f(ml,pq) = aI_mp bI[m]_lq
     *
     * matrix element at DVR point r_k:
     *     Pot[IIndex,JIndex,k]=sum[p,q]  aI^*_mp bI^*_lq V_qpk aJ_m'p bJ_l'q
     */
    const Index* iRoot=JIndex;
    for(;iRoot and iRoot->coordinates()!="Phi.Eta.Rn";iRoot=iRoot->parent());
    while(iRoot and iRoot->axisName().find("&")!=std::string::npos)iRoot=iRoot->descend();
    if(not iRoot)DEVABORT("only for floor in hierarchy Phi.Eta.Rn, found"+JIndex->root()->coordinates());

    PrintOutput::DEVwarning("re-setup floor, needs fixing",2);
    setup(IIndex);
    invOvr.clear();

    vector<unsigned int> iIdx(IIndex->index());
    vector<unsigned int> jIdx(JIndex->index());
    iIdx=std::vector<unsigned int>(iIdx.begin()+iRoot->depth(),iIdx.end());
    jIdx=std::vector<unsigned int>(jIdx.begin()+iRoot->depth(),jIdx.end());

    _rows=IIndex->sizeCompute();
    _cols=JIndex->sizeCompute();
    oNorm=0;

    size_t iFem=0,jFem=0;
    if(femLev!=Index::npos){
        iFem=iIdx[femLev];
        jFem=jIdx[femLev];
    }

    if(iFem==jFem){
        // Setup DVR rule
        while(invOvr.size()<=iFem) invOvr.push_back(std::vector<std::complex<double> >());
        if(invOvr[iFem].size()==0){
            const InverseDVR* invDVR = dynamic_cast<const InverseDVR*>(iRoot->inverseOverlap());
            int pos=IIndex->posIndex(iRoot);
            for(unsigned int i=0; i<IIndex->size(); i++)invOvr[iFem].push_back(invDVR->diagonal()[pos+i]);
        }


        // get aIaJ[p] := aI^*_mp aJ_m'p
        vector<complex<double> > aIaJ(phiBas[jIdx[mLev]]);
        for(size_t k=0;k<aIaJ.size();k++)aIaJ[k]*=std::conj(phiBas[iIdx[mLev]][k]);
        vector<complex<double> > bIbJ(etaBas[jIdx[mLev]][jIdx[lLev]]);
        for(size_t k=0;k<bIbJ.size();k++)bIbJ[k]*=std::conj(etaBas[iIdx[mLev]][iIdx[lLev]][k]);

        //actually add the data to the operator
        const BasisDVR* dvr=dynamic_cast<const BasisDVR*>(IIndex->basis());
        if(!dvr)DEVABORT("need dvr basis");
        vector<complex<double> > diag(IIndex->size(),0.);
        for (size_t p=0,pqk=0;p<aIaJ.size();p++)
            for (size_t q=0;q<bIbJ.size();q++){
                complex<double> ab=aIaJ[p]*bIbJ[q];
                for (size_t k=0;k<IIndex->size();k++,pqk++){
                    diag[k]+=ab*potIJ[iFem][pqk];
                }
            }
        string hash="Pot3d"+IIndex->hash()+JIndex->hash();
        for(size_t k=0;k<diag.size();k++){
            oNorm=max(oNorm,abs(diag[k]));
            diag[k]/=invOvr[iFem][k];
        }
        dat=addComplex(hash,diag);
    }
}
