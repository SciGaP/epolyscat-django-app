// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatormapchannelssurface.h"
#include "operatorGradient.h"
#include "discretization.h"
#include "discretizationFactor.h"
#include "discretizationSpectral.h"
#include "discretizationSurface.h"
#include "discretizationGrid.h"
#include "discretizationHybrid.h"
#include "wavefunction.h"
#include "coefficients.h"
#include "readInput.h"
#include "index.h"
//#include "operator.h"
#include "operatorGradient.h"
#include "timePropagator.h"
#include "timePropagatorOutput.h"
#include "pulse.h"
#include "derivativeFlat.h"
#include "derivativeLocal.h"
#include "operatorDiagonal.h"
#include "timer.h"
#include "printOutput.h"
#include "operatorDefinition.h"
#include "parameters.h"
#include "parallelOperator.h"

using namespace std;

void OperatorMapChannelsSurface::initializeIonWfs()
{
    if(spec==0 or ionCh==0) ABORT("Cannot initilize wavefunctions");

    if(ionWf==0) ionWf = new Wavefunction(ionCh);
    ionWf->setToZero();

    // Compute initial states
    Coefficients specCoeff(spec->idx());
    Coefficients ionTemp(ion->idx());
    vector<complex<double>* > specCoeffPointers;
    specCoeff.pointerToC(specCoeffPointers);
    for(unsigned int k=0;k<specCoeffPointers.size();k++){
        specCoeff.setToZero();
        *specCoeffPointers[k] = 1.0;
        ionTemp.setToZero();
        spec->mapToParent()->apply(1.,specCoeff,1.,ionTemp);
        *ionWf->coefs->child(k) = ionTemp;
    }

}

void OperatorMapChannelsSurface::printInfo() const {
    PrintOutput::title(Str("IONIC STATES AND CHANNELS"));
    if(ionCh->name+freeCh->name!="")PrintOutput::subTitle("     "+ionCh->name+" and "+freeCh->name);
    PrintOutput::paragraph();
    PrintOutput::newRow();
    PrintOutput::rowItem("Channel");
    PrintOutput::rowItem("Real(Energy)");
    PrintOutput::rowItem("Imag(Energy)");

    UseMatrix mat;
    chanOp->matrix(mat);
    for(unsigned int k=0;k<chanOp->iIndex->childSize();k++){
        PrintOutput::newRow();
        PrintOutput::rowItem(k);
        PrintOutput::rowItem(mat(k,k).real(),6);
        PrintOutput::rowItem(mat(k,k).imag(),6);
    }
    PrintOutput::newLine();
    PrintOutput::flush();
}

void OperatorMapChannelsSurface::setup(ReadInput & In){

    string ChanOp,ChanInt,IonAxes,shift;
    double Emin,Emax;
    read(In,ChanOp,ChanInt,IonAxes,shift,Emin,Emax);

    //    const Discretization* _mainDisc;
    hyb=dynamic_cast<const DiscretizationHybrid*>(_mainDisc);
    if(hyb!=0){
        //HACK temporarily detach index such that setup works
        const_cast<BasicDisc*>(hyb->comp[1])->idx()->parentRef()=0;
        _mainDisc=hyb->comp[1];
    }

    ion = new DiscretizationFactor(_mainDisc,IonAxes);

    OperatorTree ionOp("ChannelOperator",OperatorDefinition(ChanOp,ion->idx()->hierarchy()),ion->idx(),ion->idx());
    UseMatrix mat;
    ionOp.matrix(mat);

    if(shift == "max") { // shift by the ionic channel with highest energy
        chanEmax = mat(0, 0).real();
        for(unsigned i = 1; i < ionOp.iIndex->childSize(); ++i)
            if(std::abs(mat(i, i).real()) < std::abs(chanEmax)) chanEmax = mat(i, i).real();
    }
    else if(shift == "none")
        chanEmax = 0.;
    else
        ABORT("unknown global energy shift option \"" + shift + "\"");

    if(MPIwrapper::Size() > 1) {
        std::complex<double> dat{ chanEmax };
        MPIwrapper::AllreduceSUM(&dat, 1);
        chanEmax = dat.real();
    }


    spec=new DiscretizationSpectral(ion,&ionOp,Emin,Emax);
    if(spec==0 or spec->idx()==0 or spec->idx()->sizeStored()==0)
        ABORT("no channels found in energy range ["+tools::str(Emin)+","+tools::str(Emax)+"]");

    ionCh = new DiscretizationFactor(_mainDisc,IonAxes,false,spec->idx()->sizeStored());
    freeCh = new DiscretizationFactor(_mainDisc,IonAxes,true,spec->idx()->sizeStored());

    //HACK - should get consistent handling of ChanOp definition
    chanOp = new OperatorTree("ChannelOperator Total",
                              OperatorDefinition("<1>"+ChanOp+"+<1>"+ChanInt,ionCh->idx()->hierarchy()),ionCh->idx(),ionCh->idx()); //Ionic channel Hamiltonian
    Parameters::update(Pulse::gettBegin()); // set time-dependent parameters

    if(Algebra::isAlgebra("Rg")){
        size_t k=0;
        for(;k<_mainDisc->getAxis().size();k++)
            if(_mainDisc->getAxis()[k].name=="Rn")break;
        double box=abs(_mainDisc->getAxis()[k].boxsize());
        if(box-Parameters::pointer("Rg")->real()<-1.e-10*box)
            ABORT(Str("gauge radius Rg=")+Parameters::pointer("Rg")->real()+"larger than box-size of"+box);
    }
    else
        PrintOutput::warning("no gauge radius found - use of haCC requires mixed gauge");

    ionWf=0;
    initializeIonWfs();
    iniState = new Wavefunction(ionCh);
    *iniState = *ionWf;                      //save a copy for untwisting
    freeWf = new Wavefunction(freeCh);

    printInfo();

    // create channel surface
    vector<double> surf;
    In.read("Surface","points",surf,"","save values and derivatives at surface points (blank-separated list)");
//    if(freeCh->continuityLevel.size()!=1)ABORT("Which continuity?");

    // write either surfaces deriatives for spectra or full gradient
    if(0!=(grad=OperatorGradient::read(freeCh,In))){
        freeChSurf=0;
        name="grad_"+freeCh->name;
        iIndex = grad->iIndex;
        jIndex = grad->jIndex;
    }
    else {
        // name of current operator, decides surface ouput file name
        freeChSurf = new DiscretizationSurface(freeCh,surf,0);
        vector<string> pars;
        pars = tools::splitString(freeChSurf->mapFromParent()->name,'<');
        name=pars[0];     // Informative surface file name !
        // Its a map from jIndex to iIndex
        iIndex = freeChSurf->idx();        // Map To disc
        jIndex = _mainDisc->idx();              // Map From disc
    }

    // Case: Spectrum code HACK
    if(OperatorAbstract::useTensor==false and OperatorAbstract::fuseOp==false and OperatorAbstract::useOperatorFloor==false){  // These flags are set in main_spectrum
        prop=0;
        odeLoc=0;
    }
    else {

        // Defining the propagator
        // Replicating the code form main file, may be avoided
        double accuracy,cutE,fixStep;
        double tBeg=0.,tEnd=0.,tPrint,tStore;
        string method;
        TimePropagator::read(In,tBeg,tEnd,tPrint,tStore,accuracy,cutE,fixStep,applyThreshold,method);

        string preconDef;
        ReadInput::main.read("Operator","projection",preconDef,"","spectrally project in time propagation wrt to this operator");

        //set all wfs time to tBeg
        iniState->time = tBeg;
        ionWf->time = tBeg;
        freeWf->time = tBeg;

        prop=0;
    }

    // re-attach index (if needed)
    if(hyb!=0)const_cast<BasicDisc*>(hyb->comp[1])->idx()->parentRef()=hyb->idx();
}

double OperatorMapChannelsSurface::energyShift() const {
    return chanEmax;
}

static DerivativeFlat* df;
std::string OperatorMapChannelsSurface::axis(ReadInput &In){
    string ChanOp,ChanInt,IonAxes="",shift;
    double Emin,Emax;
    // force all reads for completing input
    if(In.found("Channel"))read(In,ChanOp,ChanInt,IonAxes,shift,Emin,Emax);
    return IonAxes;
}

void OperatorMapChannelsSurface::read(ReadInput & In, string & ChanOp, string & ChanInt, string & IonAxes, string &Shift, double& Emin, double& Emax){
    if(In.found("Channel")){
        In.read("Channel","Axis",IonAxes,ReadInput::noDefault,"Ionic channel Axis");
        In.read("Channel","ChanOp",ChanOp,ReadInput::noDefault,"Channel hamiltonian");
        In.read("Channel","ChanInteraction",ChanInt,ReadInput::noDefault,"Channel field interaction");
        In.read("Channel","Emin",Emin,ReadInput::noDefault,"Lower energy");
        In.read("Channel","Emax",Emax,ReadInput::noDefault,"Upper energy");
        In.read("Channel","Shift",Shift,"none","Constant global energy shift");
    }
}


OperatorMapChannelsSurface::OperatorMapChannelsSurface(ReadInput& In, const Discretization* Disc)
    :_mainDisc(Disc){setup(In);}

void OperatorMapChannelsSurface::setProp() const {

    const OperatorAbstract * ovr=ionCh->idx()->localOverlap(),*chan=chanOp;
    DiscretizationSpectral* proj=0;
    OperatorTree* asTree = new OperatorTree(chanOp);
    if(energyShift() != 0.) {
        std::string shiftE=tools::str(-0.5*chanEmax);
        OperatorTree* shift = new OperatorTree("Ionic Shift", OperatorDefinition(shiftE+"(<<Id>>+<<Id>>)",asTree->iIndex->hierarchy()), asTree->iIndex, asTree->iIndex); // odd HACK in definition string, as otherwise the parsing leads to a hierarchy structure incompatible with all other haCC operators
        ParallelOperator par(shift);
        par.bcast();
        asTree->add(shift);
    }
    DerivativeFlat* der=new DerivativeFlat(asTree,applyThreshold,proj);
    df = der;
    DerivativeLocal* derLoc=new DerivativeLocal(der);
    const_cast<OperatorMapChannelsSurface*>(this)->odeLoc=new OdeRK4<LinSpaceMap<CoefficientsLocal>,CoefficientsLocal>(derLoc);

    TimePropagatorOutput* out = new TimePropagatorOutput(DBL_MAX);
    out->addExpec(const_cast<OperatorAbstract*>(ovr));
    out->addExpec(const_cast<OperatorAbstract*>(chan));

    const_cast<OperatorMapChannelsSurface*>(this)->prop=new TimePropagator(odeLoc,out,1.e-10,0);

}

OperatorMapChannelsSurface::~OperatorMapChannelsSurface()
{
    if(ionWf!=0)     {delete ionWf; ionWf=0;}
    if(freeWf!=0)    {delete freeWf; freeWf=0;}
    if(iniState!=0)  {delete iniState; iniState=0;}
    if(prop!=0)      {delete prop; prop=0;}
    if(chanOp!=0)    {delete chanOp; chanOp=0;}
    //HACK (segfaults)    if(odeLoc!=0)    {delete odeLoc; odeLoc=0;}
}

void OperatorMapChannelsSurface::apply(std::complex<double> A, const Coefficients &X, std::complex<double> B, Coefficients &Y) const{
    // set up propagator
    if(prop==0)setProp();
    if(A!=1.)DEVABORT("must call with A=1");
    if(B!=0.)DEVABORT("must call with B=0");

    //logically, the next few lines should go into "update(Time)"
    CoefficientsGlobal * globalC=CoefficientsGlobal::view(ionWf->coefs);
    CoefficientsLocal * localC  =CoefficientsLocal::view(ionWf->coefs);

    MPIwrapper::Bcast(globalC->storageData(),globalC->size(),MPIwrapper::master());
    odeLoc->step(*localC,ionWf->time,_time-ionWf->time);
    ionWf->time = _time;
    
    // Contract with incoming coefficient X with ionWfs
    //HACK for hybrid case
    if(hyb==0 and X.idx()->axisName().find("&")==string::npos)
        freeCh->contract(X,*ionWf->coefs,*freeWf->coefs);
    else
        freeCh->contract(*X.child(1),*ionWf->coefs,*freeWf->coefs);

    freeWf->time = ionWf->time;

    // Map to surface disc or to gradient
    if(freeChSurf!=0)
        const_cast<OperatorAbstract*>(freeChSurf->mapFromParent())->axpy(1., *freeWf->coefs, 0., Y, _time);
    else
        grad->axpy(1., *freeWf->coefs,0., Y, _time);

}

void OperatorMapChannelsSurface::currentUntwistingMatrix(UseMatrix &res)
{
    res = UseMatrix::Zero(ionCh->idx()->childSize(),ionCh->idx()->childSize());

    for(unsigned int i=0;i<ionCh->idx()->childSize();i++)
        for(unsigned int k=0;k<ionCh->idx()->childSize();k++)
            res(i,k).complex() = ionCh->idx()->child(0)->localOverlap()->matrixElement(*iniState->coefs->child(i),
                                                                                     *ionWf->coefs->child(k));

    res = res.inverse();
}

DiscretizationSurface *OperatorMapChannelsSurface::pointerToChannelSurface()
{
    return freeChSurf;
}
