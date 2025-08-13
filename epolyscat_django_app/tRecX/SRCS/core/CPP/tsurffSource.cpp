// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "tsurffSource.h"

#include "mpiWrapper.h"
#include "threads.h"

#include "surfaceFlux.h"
#include "discretization.h"
#include "discretizationSurface.h"
#include "discretizationtsurffspectra.h"
#include "operatorDefinition.h"
#include "parameters.h"
#include "index.h"
#include "timer.h"
#include "printOutput.h"

#include "gaunt.h"

#include "operatorTreeVolkov.h"
#include "volkovGrid.h"
#include "spectrumISurff.h"

#include "log.h"
#include "operatorTreeExplorer.h"
#include "pulse.h"
#include "readInputList.h"
#include "basisGrid.h"
#include "coefficientsPermute.h"
#include "fileTools.h"
#include "timeCritical.h"

#include "debugInfo.h"
#include "basisNdim.h"
#include "asciiFile.h"

static CoefficientsLocal* _viewSource=0;
static CoefficientsLocal* _viewSourceDamped=0;


static int surfaceNumber(std::string Ax, const Index* Idx){
    const Index* idx=Idx;
    int surfNum=0;
    while(idx!=0 and idx->axisName()!=Ax){
        if(idx->continuity()!=Index::npos){
            surfNum++;
            idx=idx->descend();
        } else {
            idx=idx->nodeNext();
        }
    }
    if(idx==0)ABORT("axis "+Ax+" is not unbound axis in "+Idx->hierarchy());
    return surfNum;
}

const OperatorTree* TsurffSource::commutator() const {return dynamic_cast<OperatorTree*>(_commutator);}
const OperatorMap* TsurffSource::surfaceToK() const {
    return dynamic_cast<const OperatorMap*>(smoothSpecDisc->mapFromParent());}
const OperatorAbstract* TsurffSource::volkovPhase(double Time) const {
    _volkovPhase->update(Time,0);
    return _volkovPhase;
}


bool TsurffSource::recompute(){
    return ReadInput::main.flag("recompute","recompute all amplitudes and spectra");
}

// wrapper form momenutum grid (from DiscretizationTsurffSpectra)
std::vector<double> momenta(ReadInput & Inp, int radialPoints=DiscretizationTsurffSpectra::_defaultRadialPoints){
    // get the input flags
    // default value for radialPoints
    double maxEnergy=std::max(3*Pulse::current.omegaMax(),12.*Pulse::current.uPonderomotive()),minEnergy=0.;
    bool kGrid=true;
    DiscretizationTsurffSpectra::readFlags(Inp,radialPoints,minEnergy,maxEnergy,kGrid,true);

    std::vector<double> _momenta(radialPoints);
    DiscretizationTsurffSpectra::computeMomenta(minEnergy,maxEnergy,_momenta,kGrid);
    if (_momenta.size()==0) { ABORT("specify at least one momentum"); }
    return _momenta;
}

double TsurffSource::readTurnoff(ReadInput& Inp){
    double res;
    std::string deflt=SpectrumISurff::readUse(Inp)?"0":tools::str(4.*math::pi/Pulse::current.omegaMin());
    Inp.read("Source","turnOff",res,deflt,"linearly turn off source over this time-interval, default = 1 OptCyc",1,"tAver")
            .matchesLinp(Inp.output()+"linp");
    if(SpectrumISurff::readUse(Inp) and res!=0.)ABORT("cannot use Source:turnOff != 0 with iSurff");
    return res;
}
double readEndIntegration(ReadInput& Inp, bool & Matches){
    double res;
    bool match=Inp.read("Source","end",res,Algebra::Infty,"stop time-integration here",1,"tMax")
            .matchesLinp(Inp.outputTopDir()+"linp");
    Matches=Matches and match;
    if(not Matches and not TsurffSource::recompute())
        ABORT("-tMax does not match previous, specify -recompute");
    return res;
}

void TsurffSource::allReads(ReadInput & Inp){
    readTurnoff(Inp);
    bool match(true);
    readEndIntegration(Inp,match);
    nextRegion(Inp);
    readSymmetry12(Inp);
    SurfaceFlux::read(Inp);
}

void TsurffSource::setTurnOff(ReadInput &Inp, double &EndProp){

    _turnoffRange=readTurnoff(Inp);
    bool match(true);
    _endIntegration=readEndIntegration(Inp,match);
    if(_endIntegration<DBL_MAX/2 and _endIntegration>EndProp)
        PrintOutput::warning(Str("spectral integration stops before time-propagtion:")+_endIntegration+"<"+EndProp);
    _endIntegration=std::min(EndProp,_endIntegration);
    EndProp=_endIntegration;

    if(_turnoffRange>0.){
        _currentSourceDamped.reset(new Coefficients(idx()));
        //HACK for parallel, need CoefficientsLocal view (to be replaced)
        _viewSourceDamped=CoefficientsLocal::view(_currentSourceDamped.get());
    }
}

void TsurffSource::generate(const std::vector<std::string> AxisPath, ReadInput &Inp, const Index *Parent, Index *&Surface, Index *&Source){
    // Surface discretization
    std::vector<double> radii;
    int nMomenta=DiscretizationTsurffSpectra::_defaultRadialPoints;
    Inp.texdocuCategoryAdd("Surface","points","","points");
    if(AxisPath.size()>1){
        // get surface and spectral data from previous calculations
        ReadInputList linp(Inp.outputTopDir()+"S_"+AxisPath[AxisPath.size()-2]+"/linp");
        std::string linpStr;

        linpStr=linp.readValue("Surface","points",ReadInput::noDefault,"retrieve for TsurffSource",1,"dum1","");
        for(std::string rStr: tools::splitString(linpStr,' '))radii.push_back(tools::string_to_double(rStr));
        linpStr=linp.readValue("Spectrum","radialPoints",ReadInput::noDefault,"retrieve for TsurffSource",1,"dum2","");
        nMomenta=tools::string_to_int(linpStr);
    } else {
        Inp.read("Surface","points",radii,"","save values and derivatives at surface points (blank-separated list)")
                .texdocu(R"tex(
                         {\bf Note:} at present, use exactly one value (multiple surfaces are disabled).
                         )tex");
    }
    Surface=new DiscretizationSurface::IndexS(Parent,radii,surfaceNumber(AxisPath[0],Parent));
    Source= new DiscretizationTsurffSpectra::IndexTsurffSpectra(Surface,momenta(Inp,nMomenta));
    if(AxisPath.size()==1)return;

    delete Surface;
    Index* oldSource=Source;
    generate(std::vector<std::string>(AxisPath.begin()+1,AxisPath.end()),Inp,oldSource,Surface,Source);
    delete oldSource;
}

std::string TsurffSource::str() const {

    Str s("","");
    for(const Index* ix=idx();ix;ix=ix->descend()){
        if(ix->axisName().find("k")==0){
            s=s+"\n"+ix->axisName()+SEP(", ");
            s=s+ix->basis()->size()+tools::str(pow(ix->basis()->grid()->mesh()[0],2)*0.5,3);
            s=s+tools::str(pow(ix->basis()->grid()->mesh().back(),2)*0.5,3);
            s=s+tools::str(_endIntegration,3);
            s=s+SEP(" ");
            if(_turnoffRange>0.)s=s+"[turn off: "+_turnoffRange+"]";
            std::string space="";
            for(size_t k=1;space=="" and k<ix->basis()->grid()->size()-1;k++)
                if(     std::abs(ix->basis()->grid()->mesh()[k-1]+ix->basis()->grid()->mesh()[k+1]-2.*ix->basis()->grid()->mesh()[k])>
                        std::abs(ix->basis()->grid()->mesh()[k-1]-ix->basis()->grid()->mesh()[k+1])*1.e-12)space="NOT k-equidistant";
            s=s+space;
        }
    }
    if(s!="")s=s.substr(1);
    return std::move(s);
}

std::vector<std::string> TsurffSource::unboundAxes(ReadInput & Inp){
    std::vector<Axis> ax;
    Axis::fromFile(Inp,ax);
    std::vector<std::string> unbound;
    for(const Axis & a: ax){
        if(Coordinate::isDiscrete(a.name))continue;
        //HACK: never consider axes named "Ndim" as unbounded (should test for basisNdim)
        if(a.upperEnd()>DBL_MAX/2. or a.lowerEnd()<-DBL_MAX/2)
            unbound.push_back(a.name);
        else if(a.comsca.r0up()<a.upperEnd())
            unbound.push_back(a.name);
    }
    return unbound;

}

bool TsurffSource::readSymmetry12(ReadInput &Inp, const std::string Hierarchy){
    bool symm12;
    Inp.read("Spectrum","symmetry12",symm12,"false","assume 1<->2 exchange symmetry and avoid unnecessary regions")
            .texdocu(R"tex(
                     Specific for 6D two-electron systems. Assumes two polar coordinate systems numbered as 1 and 2, as in
                     \nameref{docu:tutorial:23}
                     )tex");
    if(symm12 and Hierarchy!=""){
        std::vector<std::string> hier=tools::splitString(Hierarchy,'.');
        for(auto &c: hier)
            if(c.back()!='1' and c.back()!='2')
                ABORT("for using Spectrum:symmetry12 all need all coordinates numbered 1 or 2, have: "+Hierarchy);
    }
    return symm12;
}

void TsurffSource::moveSubdirs(ReadInput & Inp){
    std::vector<Axis> ax;
    Axis::fromFile(Inp,ax);
    std::vector<std::string> unbound(TsurffSource::unboundAxes(Inp));
    std::string backup("back");
    for(int i=0;i<100;i++){
        backup=Inp.outputTopDir()+"backup"+tools::str(i,2,'0')+"/";
        if(not folder::exists(backup)){
            folder::create(backup);
            break;
        }
    }
    for(std::string reg1: unbound){
        std::string outf=Inp.outputTopDir()+"S_"+reg1+"/"+PrintOutput::outExtension;
        if(folder::exists(outf) and tools::findFirstLine(outf,{"accept/reject","TimePropagator"})!=""){
            PrintOutput::warning("move S_"+reg1+" --> "+backup);
            std::rename((Inp.outputTopDir()+"S_"+reg1).c_str(),(backup+"S_"+reg1).c_str());
        }
    }
    for(std::string reg1: unbound){
        for(std::string reg2: unbound){
            std::string outf=Inp.outputTopDir()+"S_"+reg1+"."+reg2+"/"+PrintOutput::outExtension;
            if(folder::exists(outf) and tools::findFirstLine(outf,{"accept/reject","TimePropagator"})!=""){
                PrintOutput::warning("move S_"+reg1+"."+reg2+" --> "+backup);
                std::rename((Inp.outputTopDir()+"S_"+reg1+"."+reg2).c_str(),(backup+"S_"+reg1+"."+reg2).c_str());
            }
        }
    }

}

static std::string _nextRegion(ReadInput &Inp){

    if(MPIwrapper::isMaster()){
        bool symm12=TsurffSource::readSymmetry12(Inp);


        std::string outRegion="";

        // if not set yet, try reading:
        std::string inRegion;
        Inp.read("Source","region",inRegion,"","axis perpendicular to surface of subregion",1,"region");
        if(inRegion!=""){
            PrintOutput::message("calculation forced in region "+inRegion);
            return inRegion;
        }


        // check for completed initial run
        std::string outf=Inp.outputTopDir()+PrintOutput::outExtension;
        if(not  folder::exists(outf) or
                tools::findFirstLine(outf,{"accept/reject","TimePropagator"})==""){
            return "";
        }

        std::vector<std::string> unbound;
        //HACK
        std::vector<double> vdum;
        if(DiscretizationSurface::readDef(Inp)!="not given"){
            unbound.push_back("Rn");
        }
        else {
            std::vector<Axis> ax;
            Axis::fromFile(Inp,ax);
            unbound=TsurffSource::unboundAxes(Inp);
        }


        // single coordinate regions
        std::string skip="";
        for(std::string reg: unbound){
            std::string outf=Inp.outputTopDir()+"S_"+reg+"/"+PrintOutput::outExtension;
            if(not  folder::exists(outf) or
                    tools::findFirstLine(outf,{"accept/reject","TimePropagator"})==""){
                if(symm12 and reg.back()!='1')skip+=reg+",";
                else                  return reg;
            }
        }

        // double coordinate regions
        for(std::string reg1: unbound)
            for(std::string reg2: unbound){
                if(reg1==reg2)continue;
                std::string outf=Inp.outputTopDir()+"S_"+reg1+"."+reg2+"/"+PrintOutput::outExtension;
                if(not folder::exists(outf) or
                        tools::findFirstLine(outf,{"accept/reject","TimePropagator"})==""){
                    if(symm12 and reg2.back()!='2')skip+=reg1+"."+reg2;
                    else                   return reg1+"."+reg2;
                }
            }

        // all subregions have been done
        return Sstr+"unbound axes:"+unbound;
    }
    else{
        return "";
    }
}
std::string TsurffSource::nextRegion(ReadInput &Inp){
    std::string region;
    if(MPIwrapper::isMaster())region=_nextRegion(Inp);
    MPIwrapper::Bcast(region,MPIwrapper::master());
    return region;
}

void TsurffSource::print(std::string SourceFile) const {

    //    PrintOutput::title("Integration for spectral amplitude");
    PrintOutput::lineItem("source(s)",SourceFile);
    PrintOutput::newLine();
    PrintOutput::lineItem("end",_endIntegration);
    if(_turnoffRange>0.)
        PrintOutput::lineItem("turnOff",_turnoffRange);
    else
        PrintOutput::lineItem("turnOff","ABRUPT (Caution)");

    PrintOutput::paragraph();
    PrintOutput::paragraph();
    bool head=true;
    for(const Index* ix=idx();ix;ix=ix->descend()){
        if(ix->axisName().find("k")==0){
            if(head){
                head=false;
                PrintOutput::newRow();
                PrintOutput::rowItem("axis");
                PrintOutput::rowItem("Points");
                PrintOutput::rowItem("Emin");
                PrintOutput::rowItem("Emax");
                PrintOutput::rowItem("spacing");
            }
            PrintOutput::newRow();
            PrintOutput::rowItem(ix->axisName());
            PrintOutput::rowItem(ix->basis()->grid()->size());
            PrintOutput::rowItem(pow(ix->basis()->grid()->mesh()[0],2)*0.5);
            PrintOutput::rowItem(pow(ix->basis()->grid()->mesh().back(),2)*0.5);
            std::string space="k-equidistant";
            for(size_t k=1;space=="" and k<ix->basis()->grid()->size()-1;k++)
                if(     std::abs(ix->basis()->grid()->mesh()[k-1]+ix->basis()->grid()->mesh()[k+1]-2.*ix->basis()->grid()->mesh()[k])>
                        std::abs(ix->basis()->grid()->mesh()[k-1]-ix->basis()->grid()->mesh()[k+1])*1.e-12)space="NOT k-equidistant";
            PrintOutput::rowItem(space);
        }
    }

}

TsurffSource::TsurffSource(const std::vector<std::string> & AxisPath, const Index *Parent, ReadInput &Inp){
    //    readTurnOff(Inp,DBL_MAX);
    //    if(_endIntegration>_endTurnOff)ABORT("source turn off time after end time");
    DEVABORT("not activated");
    Index *surfI, *sourceI;
    //NOTE: surfI should be recoverd from surface_* file rather than reconstructed from input
    generate(AxisPath,Inp,Parent,surfI,sourceI);

    //    flux=new SurfaceFlux(AxisPath,Inp,surfI);
    _commutator=new OperatorTree("Commutator",OperatorDefinition("<<Commutator>>", flux->idx()->hierarchy()),flux->idx(), flux->idx());

    wfSurface = new Coefficients(flux->idx());
    wfSurface->treeOrderStorage();

    wfSmoothSpecDisc = new Coefficients(sourceI);
    wfSmoothSpecDisc->treeOrderStorage();

    //    tempSmoothSpecDisc = new Wavefunction(smoothSpecDisc);
    tempSmoothSpecDisc = new Coefficients(*wfSmoothSpecDisc);
    tempSmoothSpecDisc->treeOrderStorage();

    _volkovPhase = new VolkovGrid(flux.get(),idx());

}

TsurffSource::TsurffSource(const DiscretizationTsurffSpectra* SmoothSpecDisc, ReadInput& Inp, DiscretizationSurface *SurfD, std::shared_ptr<SurfaceFlux> &Flux):
    flux(Flux),_turnoffRange(0.),_endIntegration(DBL_MAX),smoothSpecDisc(SmoothSpecDisc)
{
    timeCritical::suspend();
    const Discretization* D = dynamic_cast<Discretization*>(SurfD);
    while(D->parent)D=D->parent;
    //AS this works, but is logically not every clean, cast should be to discretization derived, and stop when fails

    // Initialize Surface Flux
    LOG_PUSH("SurfaceFlux");
    // extract file name from surface discretization
    std::string file;
    std::vector<std::string> axes=tools::splitString(SurfD->idx()->hierarchy(),'.');
    for(std::string ax: axes)
        if(ax.find("k")==0)file+="_"+ax.substr(1);
    if(file!="")file="S"+file+"/";
    for(std::string ax: axes)
        if(ax.find("ValDer")==0)file+=DiscretizationSurface::prefix+ax.substr(6);
    file=Inp.outputTopDir()+file;
    if(flux==0)flux.reset(new SurfaceFlux(file,Inp));

    LOG_POP();

    LOG_PUSH("smoothSpecDisc");
    // determine number of k-grid points for present source
    // from Index's of remaining surfD (i.e. surfD[k], k!=SurfNumber)
    if(smoothSpecDisc==0){
        smoothSpecDisc = dynamic_cast<DiscretizationTsurffSpectra*>(
                    DiscretizationTsurffSpectra::factoryTsurff(flux->idx(),Inp));
    }

    LOG_POP();

    _commutator=new OperatorTree("Commutator",OperatorDefinition("<<Commutator>>", flux->idx()->hierarchy()),flux->idx(), flux->idx());


    LOG_PUSH("wavefunctions");
    wfSurface = new Coefficients(flux->idx());
    wfSurface->treeOrderStorage();

    //    wfSmoothSpecDisc = new Wavefunction(smoothSpecDisc);
    wfSmoothSpecDisc = new Coefficients(smoothSpecDisc->idx());
    wfSmoothSpecDisc->treeOrderStorage();

    //    tempSmoothSpecDisc = new Wavefunction(smoothSpecDisc);
    tempSmoothSpecDisc = new Coefficients(*wfSmoothSpecDisc);
    tempSmoothSpecDisc->treeOrderStorage();
    LOG_POP();

    Parameters::update(flux->FluxBufferBeginTime());
    _currentTime=flux->FluxBufferBeginTime();

    LOG_PUSH("Volkov");

    bool useGrid=true;

    // Not pretty, but should detect Helium 6D
    if(idx()->hierarchy().find("Eta1") != std::string::npos and idx()->hierarchy().find("Eta2") != std::string::npos){
        useGrid = false;
    }
    if(MPIwrapper::Size()>1 and not useGrid){
        useGrid=true;
        PrintOutput::DEVwarning("forcing useGrid in parallel code");
    }
    if(ReadInput::main.found("_EXPERT_","usegrid","DEBUGusegrid")){
        PrintOutput::warning("forcing grid usage");
        ReadInput::main.read("_EXPERT_","usegrid",useGrid,ReadInput::noDefault,"force grid true/false",0,"DEBUGusegrid");
    }

    //    useGrid=true;

    if(!useGrid)
        PrintOutput::message("Applying Volkov Phase without grid");

    if(useGrid){
        _volkovPhase = new VolkovGrid(D, flux.get(), smoothSpecDisc);
    }else{
        std::string kLevel;
        for(const Index *s=flux->idx(),*t=idx();;s=s->child(0), t=t->child(0)){

            if(s->axisName().size()>6 and t->axisName().size()>1){
                if(s->axisName().substr(0,6)=="ValDer" and t->axisName().substr(0,1)=="k"){
                    kLevel = t->axisName();
                    break;
                }
                else if(s->axisName().substr(0,6) == "ValDer"){
                    if(t->child(0)->axisName().substr(0,1) == "k"){
                        kLevel = t->child(0)->axisName();
                        break;
                    }
                }
            }
            if(t->isBottom() or s->isBottom())break;
        }

        /*
                 * Not exactly pretty, but should do
                 */
        OperatorTreeVolkov::AxisDescriptor ax = {"NONE", "NONE", kLevel};
        for(const Index* s=idx(); not s->isLeaf(); s=s->descend()){
            if(kLevel.find("kRn")==0 and s->axisName().find("Phi")==0){
                if(kLevel.size()>3 and s->axisName().size()>3){
                    if(kLevel.substr(3,1) == s->axisName().substr(3,1)) ax.phi = s->axisName();
                }else if(kLevel.size() == 3 and s->axisName().size() == 3){
                    ax.phi = s->axisName();
                }
            }
            if(kLevel.find("kRn")==0 and s->axisName().find("Eta")==0){
                if(kLevel.size()>3 and s->axisName().size()>3){
                    if(kLevel.substr(3,1) == s->axisName().substr(3,1)) ax.eta = s->axisName();
                }else if(kLevel.size() == 3 and s->axisName().size() == 3){
                    ax.eta = s->axisName();
                }
            }
        }

        PrintOutput::DEVmessage("Applying Volkov Phase on "+ax.phi+"."+ax.eta+"."+ax.k);
        // Don't use idx() due to parallelization
        _volkovPhase = new OperatorTreeVolkov(wfSmoothSpecDisc->idx(), wfSmoothSpecDisc->idx(), ax);
    }
    delete _viewSource; delete _viewSourceDamped; //HACK get rid of this soon
    _viewSource=CoefficientsLocal::view(wfSmoothSpecDisc);
    if(_currentSourceDamped)_viewSourceDamped=CoefficientsLocal::view(_currentSourceDamped.get());
    LOG_POP();
    timeCritical::resume();
}


TIMER(fluxUpdate,)
TIMER(applyCommutator,)
TIMER(smoothSpec,)
TIMER(volkov,)
CoefficientsLocal* TsurffSource::UpdateSource(double Time){
    //HACK for use in parallel

    STARTDEBUG(fluxUpdate);
    _currentTime=Time;
    if(not flux->update(Time)) return 0;
    STOPDEBUG(fluxUpdate);


    Parameters::update(Time);

    STARTDEBUG(applyCommutator);

    _commutator->apply(1., *flux->coefs, 0., *wfSurface);
    STOPDEBUG(applyCommutator);
    STARTDEBUG(smoothSpec);
    //HACK until index is consistent
    *smoothSpecDisc->mapFromParent()->tempRHS()=*wfSurface;
    smoothSpecDisc->mapFromParent()->apply(1.,*smoothSpecDisc->mapFromParent()->tempRHS(), 0., *tempSmoothSpecDisc);
    STOPDEBUG(smoothSpec);

    STARTDEBUG(volkov);
    _volkovPhase->update(Time);
    _volkovPhase->apply(1., *tempSmoothSpecDisc, 0., *wfSmoothSpecDisc);
    STOPDEBUG(volkov);

    return CurrentSource();
}

int x12count=0;
CoefficientsLocal *TsurffSource::CurrentSource(){
    if(_x12Source){
        if(++x12count>1000){
            ABORT("x21");
        }
        *wfSmoothSpecDisc+=_x12Source->fromOrig(*wfSmoothSpecDisc);
        wfSmoothSpecDisc->scale(0.5);
    }
    if(_turnoffRange==0. or _currentTime<_endIntegration-_turnoffRange){
        return _viewSource;
    }
    if(_currentTime>_endIntegration){
        PrintOutput::warning("time beyond  integration time",1);
        wfSmoothSpecDisc->setToZero();
        return _viewSource;
    }
    // linearly turn off source
    *_currentSourceDamped=*wfSmoothSpecDisc;
    _currentSourceDamped->scale((_endIntegration-_currentTime)/(_turnoffRange));
    return _viewSourceDamped;
}

double TsurffSource::SourceBufferBeginTime()
{
    return flux->FluxBufferBeginTime();     // Should be used only in the beginning
}

const Index* TsurffSource::idx() const{return wfSmoothSpecDisc->idx();}

std::vector<double> TsurffSource::surfaceRadius(const std::string Axis) const {
    return flux->idx()->grid("ValDer"+Axis);
}

