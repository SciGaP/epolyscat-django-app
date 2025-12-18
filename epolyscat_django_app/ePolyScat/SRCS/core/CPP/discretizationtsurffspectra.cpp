// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "discretizationtsurffspectra.h"
#include "pulse.h"
#include "inverseDvr.h"
#include "overlapDVR.h"
//#include "operator.h"
#include "operatorMap.h"
#include "readInputList.h"
#include "discretizationSurface.h"
#include "discretizationCoulXNew.h"
#include "surfaceFlux.h"
//#include "basisMat.h"
#include "operatorDefinition.h"
#include "readInput.h"
#include "log.h"
#include "basisGrid.h"
#include "basisGridQuad.h"
#include "threads.h"
#include "indexOverlap.h"
#include "tsurffSource.h"

// for std::find()
#include <algorithm>
#include "eigenTools.h"

using namespace std;

int DiscretizationTsurffSpectra::_defaultRadialPoints=50;


void DiscretizationTsurffSpectra::read(ReadInput &Inp, const Discretization *D, int unboundDOF, std::string &Region, int &NSurf){
    if(Region=="")return;

    std::vector<std::string> infCoorNames;
    for(unsigned int i=0; i<D->continuityLevel.size(); i++){
        infCoorNames.push_back(D->getAxis()[D->continuityLevel[i]].name);
    }
    if(Region==""){ // no flag
        DiscretizationTsurffSpectra::getNextRegion(unboundDOF, Region, NSurf, infCoorNames, Inp);
    }
    else{
        // what is this doing hsurfaceFilelike it does nothing?
        DiscretizationSurface::surfacePath(Inp.outputTopDir(),D->idx(),Region);
    }
}

// could be moved to "stringTools.cpp"
std::vector<std::string> DiscretizationTsurffSpectra::splitRegion(const std::string s){
    std::string tempString = s;
    std::vector<std::string> elements;
    while(tempString.size()>0){
        bool found=false;
        // reverse iterator (start from end of Coordinate::list)
        for(auto it = Coordinate::list.rbegin(); it!=Coordinate::list.rend(); ++it){
            if(tempString.substr(0,it->first.length())==it->first){
                found=true;
                tempString=tempString.substr(it->first.length(), tempString.length());
                std::string number="";
                while(tempString.substr(0,1).find_first_of("0123456789")!=string::npos){
                    number+=tempString.substr(0,1);
                    tempString=tempString.substr(1,tempString.length());
                }
                elements.push_back(it->first+number);
                break;
            }
        }
        if(found==false){
            ABORT("Cannot split region into its components: "+s);
        }
    }
    return elements;
}

void DiscretizationTsurffSpectra::checkRegion(int &unbound, std::string &Region, int &NSurf, std::vector<std::string> &infCoorNames, ReadInput &inp){
    // List of all infinite coordinates appearing in Region
    std::vector<std::string> infCoorList = DiscretizationTsurffSpectra::splitRegion(Region);

    unbound = infCoorList.size();

    if(unbound>infCoorNames.size()) ABORT("Too many axis names specified in input flag");

    // check whether the elements of infCoorList are in infCoorNames
    // Maybe, this should be done using "Coordinate::list"
    for(unsigned int i=0; i<infCoorList.size(); i++){
        std::string tempName = infCoorList[i];
        if(std::find(infCoorNames.begin(), infCoorNames.end(), tempName)==infCoorNames.end()){
            ABORT("axis "+infCoorList[i]+" is not valid");
        }
    }

    // Determine whether all necessary subregions exist
    if(unbound>1){
        std::vector<std::string> subregions, tempNames(unbound-1);
        // Determine the names of all subregions going into "region"
        DiscretizationTsurffSpectra::getRegions(0, 0, unbound-1, infCoorList, subregions, tempNames);
        for(unsigned int i=0; i<subregions.size();i++){
            if(folder::exists(inp.outputTopDir()+"S_"+subregions[i])==false) ABORT("Input region "+Region+" is invalid as subdirectory S_"+subregions[i]+" is missing");
        }
    }

    if(unbound==infCoorNames.size()) ABORT("All subregions computed - use Spectrum to generate final spectra");

    // determine NSurf
    std::vector<std::string> regions, names(unbound);
    DiscretizationTsurffSpectra::getRegions(0, 0, unbound, infCoorNames, regions, names);
    for(unsigned int j=0; j<regions.size(); j++){
        if(regions[j]==Region){
            NSurf=j;
            break;
        }
    }
}

void DiscretizationTsurffSpectra::getNextRegion(int &unbound, std::string &Region, int &NSurf, std::vector<std::string> &infCoorNames, ReadInput &inp){
    for(unsigned int i=1; i<infCoorNames.size(); i++){
        unbound = i;
        std::vector<std::string> regions, names(unbound);
        DiscretizationTsurffSpectra::getRegions(0, 0, unbound, infCoorNames, regions, names);
        // check, whether the directories exist
        for(unsigned int j=0; j<regions.size(); j++){
            if(folder::exists(inp.outputTopDir()+"S_"+regions[j])==false){
                Region=regions[j];
                NSurf = j;
                return;
            }
        }
    }
    ABORT("All subregions computed - use Spectrum to generate final spectra");
}

void DiscretizationTsurffSpectra::getRegions(int pos, unsigned int ibeg, int unbound, std::vector<std::string> &infCoorNames, std::vector<std::string> &regions, std::vector<std::string> &names){
    for(unsigned int i=ibeg; i<infCoorNames.size()-unbound+pos+1; i++){
        names[pos] = infCoorNames[i];
        if(pos+1<unbound){
            DiscretizationTsurffSpectra::getRegions(pos+1, i+1, unbound, infCoorNames, regions, names);
        }
        else{
            std::string region="";
            for(unsigned int j=0; j<names.size(); j++){
                region+=names[j];
            }
            regions.push_back(region);
        }
    }
}

void DiscretizationTsurffSpectra::readFlags(ReadInput &Inp, int &radialPoints, double &MinEnergy, double &maxEnergy, bool &kGrid,bool AllowFlags){
    //function is called >=twice for parabolic coords
    
    
    // this input logics is hacked badly and very confusing
    if(radialPoints==0)radialPoints=_defaultRadialPoints;
    maxEnergy=max(3*Pulse::current.omegaMax(),12.*Pulse::current.uPonderomotive());
    bool eGrid;
    ReadInput::main.texdocuCategoryAdd("Spectrum","radialPoints,maxEnergy,minEnergy,plot,symmetry12,energyGrid,spectrum,amplitude,benchmark",
                                       R"tex(
                                       Controls grid on which spectral amplitudes are computed and how spectra are extracted.
                                       \\ Get different plots by re-running >tRecX runDir/00xx -\nameref{docu:Spectrum:plot}=XXX. This is cheap.
                                       \\ Changing grids requires re-computing amplitudes from surface files
                                       with the -\nameref{docu:Flag:recompute} flag. This can be very expensive for 2-particle spectra.
                                       )tex","09,20,23,220");
    for(int k=0;k<2;k++){
        ReadInput *inp;
        if(k==0)inp=&ReadInput::main;
        if(k==1)inp=&Inp;
        std::string linpFile=inp->output()+"linp";
        bool matches(true);
        matches=matches and
                inp->read("Spectrum", "maxEnergy", maxEnergy,tools::str(maxEnergy), "spectral energy range, default = max(12Up,3 Omega)",1,"Emax")
                .matchesLinp(linpFile);
        matches=matches and
                inp->read("Spectrum", "minEnergy", MinEnergy,"0", "lowest energy",1,"Emin")
                .texdocu(R"tex(
                         Note: energy =0 will NOT be in the grid.
                         )tex")
                .matchesLinp(linpFile);
        // if Spectrum is not specified, there is also no default number of radial points, but can input by flag
        if(not inp->found("Spectrum"))radialPoints=0;
        matches=matches and
                inp->read("Spectrum", "radialPoints", radialPoints,tools::str(radialPoints), "number or momentum points in energy range",1,"nR")
                .texdocu(R"tex(
                         Number of grid points on the continous spectrum axis. Refers to any  continous axis.
                         At present all continuous axes will have the same number of grid points.
                         For different points axes contact developers.
                         ("radial" is legacy from the initial 3d spherical problem)
                         )tex")
                .matchesLinp(linpFile);
        matches=matches and
                inp->read("Spectrum", "energyGrid",eGrid,tools::str(false),"grid equidistant in energy (slower, not recommended)",1,"eGrid")
                .texdocu(R"tex(
                         Default is a k-grid (modulus of momentum), which allows faster evaluation of phases.
                         eGrid=true may lead to further slow-down in time-critical calculations.
                         )tex")
                .matchesLinp(linpFile);
        if(radialPoints>256)
            PrintOutput::warning(Sstr+"Large number of"+radialPoints
                                 +"radial spectral points, may lead to large memory, recommended value <=256");
        if(not TsurffSource::recompute() and not matches){
            ABORT("spectrum grid inputs do not match originals: use matching flags or -recompute");
        }
    }

    kGrid=not eGrid;


    if( not AllowFlags and (
                ReadInput::main.found("Spectrum","benchmark","benchmark") or
                ReadInput::main.found("Spectrum","minEnergy","Emin") or
                ReadInput::main.found("Spectrum","maxEnergy","Emax") or
                ReadInput::main.found("Spectrum","radialPoints","nR") or
                ReadInput::main.found("Spectrum","energyGrid","eGrid"))){
        ABORT("must not specify flags with incoming current from subregion - already fixed by sub-region calculation");
    }

    // zero energy is not admissible
    if(MinEnergy<=0.){
        if(kGrid)MinEnergy=pow(sqrt(2*maxEnergy)/radialPoints,2)*0.5;
        else     MinEnergy=maxEnergy/radialPoints;
    }

    // special case for Esry benchmark
    bool dum;
    Inp.read("Spectrum", "benchmark",dum,"false","standard spectra for Bret Esry's benchmark",1,"benchmark")
            .texdocu("[DEVELOPER,TEST]");
    if(ReadInput::main.found("Spectrum","benchmark","benchmark")){
        maxEnergy=10*Pulse::current.uPonderomotive();
        MinEnergy=maxEnergy*0.01;
        radialPoints=2;
    }
}

DiscretizationTsurffSpectra * DiscretizationTsurffSpectra::factoryTsurff(const Index *ParentIdx, ReadInput &Inp, int radialpoints){
    if(radialpoints == 0)
        return new DiscretizationTsurffSpectra(ParentIdx, Inp);
    else
        return new DiscretizationTsurffSpectra(ParentIdx, Inp, radialpoints);
}

//should be moved into AxisTree
AxisTree* fromVec(const std::vector<Axis> Ax){
    if(Ax.size()==0)return 0;
    AxisTree * tree=new AxisTree(Ax[0]);
    for(int k=1;k<Ax.size();k++)tree->childAdd(new AxisTree(Ax[k]));
    return tree;
}



//DiscretizationTsurffSpectra::DiscretizationTsurffSpectra(const DiscretizationSurface *Parent, ReadInput &Inp, int RadialPoints)
//{
//    parent=Parent;
//    name=Parent->name;

//    LOG_PUSH("DiscretizationTsurffSpectra");

//    // get the input flags
//    // default value for radialPoints
//    int radialPoints=_defaultRadialPoints;
//    radialPoints=RadialPoints;
//    double maxEnergy=max(3*Pulse::current.omegaMax(),12.*Pulse::current.uPonderomotive()),minEnergy=0.;
//    bool kGrid=true;
//    //reached
//    readFlags(Inp,radialPoints,minEnergy,maxEnergy,kGrid,true);
//    DEVABORT("debug");
//    // not reached
//    // replace default value for radialPoints by input value
//    if(RadialPoints!=_defaultRadialPoints and radialPoints==_defaultRadialPoints)radialPoints=RadialPoints;
    
//    vector<double> momenta(radialPoints);
//    computeMomenta(minEnergy,maxEnergy,momenta,kGrid);

//    _construct(Parent->idx(),momenta);
//    LOG_POP();
//}

//DiscretizationTsurffSpectra::DiscretizationTsurffSpectra(const SurfaceFlux *Flux, ReadInput &Inp, int RadialPoints)
//{
//    parent=0;
//    name=Flux->idx()->hierarchy();
//    LOG_PUSH("DiscretizationTsurffSpectra");

//    // get the input flags
//    // default value for radialPoints
//    int radialPoints=_defaultRadialPoints;
//    radialPoints=RadialPoints;
//    double maxEnergy=max(3*Pulse::current.omegaMax(),12.*Pulse::current.uPonderomotive()),minEnergy=0.;
//    bool kGrid=true;
//    readFlags(Inp,radialPoints,minEnergy,maxEnergy,kGrid,true);
//    // replace default value for radialPoints by input value
//    if(RadialPoints!=_defaultRadialPoints and radialPoints==_defaultRadialPoints)radialPoints=RadialPoints;

//    vector<double> momenta(radialPoints);
//    computeMomenta(minEnergy,maxEnergy,momenta,kGrid);

//    _construct(Flux->idx(),momenta);
//    LOG_POP();
//}

DiscretizationTsurffSpectra::DiscretizationTsurffSpectra(const Index *FluxIdx, ReadInput &Inp, int RadialPoints)
{
    parent=0;
    name=FluxIdx->hierarchy();
    LOG_PUSH("DiscretizationTsurffSpectra");

    // get the input flags
    // default value for radialPoints
    int radialPoints=_defaultRadialPoints;
    radialPoints=RadialPoints;
    double maxEnergy=max(3*Pulse::current.omegaMax(),12.*Pulse::current.uPonderomotive()),minEnergy=0.;
    bool kGrid=true;
    readFlags(Inp,radialPoints,minEnergy,maxEnergy,kGrid,true);
    // replace default value for radialPoints by input value
    if(RadialPoints!=_defaultRadialPoints and radialPoints==_defaultRadialPoints)radialPoints=RadialPoints;

    vector<double> momenta(radialPoints);
    computeMomenta(minEnergy,maxEnergy,momenta,kGrid);

    _construct(FluxIdx,momenta);
    LOG_POP();
}

//DiscretizationTsurffSpectra::DiscretizationTsurffSpectra(const DiscretizationSurface *Parent, const std::vector<double> Momenta)
//{
//    parent=Parent;
//    name=Parent->name;
//    _construct(parent->idx(),Momenta);
//}

void DiscretizationTsurffSpectra::_construct(const Index *IParent, const std::vector<double> Momenta){
    if (Momenta.size()==0) {ABORT("specify at least one momentum"); }

    LOG_PUSH("IndexTsurffSpectra");

    const Index* id=IParent;
    while(id and id->axisName().find("ValDer")!=0)id=id->descend();
    if(not id)DEVABORT("no ValDer-level in hiearchy "+IParent->hierarchy());

    // find the new k-level (=surface level)
    idx() = new IndexTsurffSpectra(IParent,Momenta);
    const Index * six=IParent;
    while(six and six->axisName().find("ValDer")!=0)six=six->descend();

    idx()=Threads::fork(idx(),"k"+six->axisName().substr(6));
    idx()=Threads::detach(idx()); // keep only present

    LOG_POP();
    _mapFromParent.reset(new OperatorMap(idx(),IParent));

    if(idx()->hierarchy().find("kRn1.kRn2")+9!=idx()->hierarchy().size()){
        IndexOverlap::set(idx());
        idx()->setInverseOverlap(Inverse::factory(idx()));
    }
}

void DiscretizationTsurffSpectra::computeMomenta(double minEnergy, double maxEnergy, vector<double>& momenta, bool KGrid)
{
    if(not KGrid)
        for (int k=0;k<momenta.size();k++)
            momenta[k]=sqrt(2.*(minEnergy+k*(maxEnergy-minEnergy)/(momenta.size()-1)));
    else{
        double kmax=sqrt(2*maxEnergy);
        double kmin=sqrt(2*minEnergy);
        for (size_t k=0;k<momenta.size();k++)momenta[k]=(kmin+k*(kmax-kmin)/(momenta.size()-1));
    }
}

std::vector<double> quadWeights(std::string AxName,const std::vector<double> Momenta){
    std::vector<double> w(Momenta.size(),1.);

    for(size_t k=1;k<w.size()-1;k++)
        w[k]=0.5*(Momenta[k+1]-Momenta[k-1]);
    w[0]=w[1];
    w.back()=w[w.size()-2];

    if(AxName.find("kX")==0 or AxName.find("kY")==0 or AxName.find("kZ")==0)
        return w;
    else if(AxName.find("kRn")==0){
        for(int k=0;k<w.size();k++)w[k]*=std::pow(Momenta[k],2);
    }
    else
        ABORT("cannot handle quadrature weight for Axis= "+AxName);
    return w;
}

DiscretizationTsurffSpectra::IndexTsurffSpectra::IndexTsurffSpectra(const Index *From, const std::vector<double> momenta)
{
    Index::nodeCopy(From,false); // true copy, not view
    fromIndex.push_back(From);

    string s=From->axisName();
    if(From->axisName().substr(0,4)=="surf")setAxisName(s.replace(0,4,"spec"));
    if(From->axisName().substr(0,6)!="ValDer")
    {   // not desired valder level - descend, unless truncated
        if(From->basis()->strDefinition()!="Vector:1"){
        if(From->isBottom())ABORT(From->root()->str()+"\n\nNo unbounded axis found in Index");
        for(unsigned int k=0;k<From->childSize();k++)
            childAdd(new IndexTsurffSpectra(From->child(k),momenta));
        }
    }
    else{
        setAxisName("k"+From->axisName().substr(6,From->axisName().size()));
        // get proper quadrature weights
        setBasis(BasisGridQuad::factory(momenta,quadWeights(axisName(),momenta)));

        if(From->basis()->size()!=2) ABORT("More than value and derivative??");
        if(not From->isBottom())
            for(unsigned int l=0;l<momenta.size();l++)
                childAdd(new Index(*From->child(0)));    // Assumes the branches below are not dependent on v/d index.
    }
    sizeCompute();
}
