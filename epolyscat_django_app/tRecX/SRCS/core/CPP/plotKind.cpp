// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "plotKind.h"
#include "abort.h"
#include "tools.h"

#include "str.h"
#include "algebra.h"
#include "constants.h"
#include "index.h"
#include "indexGrid.h"
#include "discretizationGrid.h"
#include "spectrumLinCom.h"
#include "readInput.h"

using namespace std;

map<std::string,std::map<std::string,std::shared_ptr<PlotKind>> > PlotKind::plotKinds;

std::string PlotKind::defaultKind(const std::string Hierarchy){
    std::string coors=Index::coordinates(Hierarchy);
    if(coors.find("kRn1")!=string::npos and coors.find("kRn2")!=string::npos)return "total";
    if(coors.find("kX")!=string::npos and coors.find("kY")!=string::npos)return "total";
    if(coors.find("kX1")!=string::npos and coors.find("kX2")!=string::npos)return "total";
    if(coors==("kX"))return "total";
    if(coors==("kY"))return "total";
    if(coors==("kZ"))return "total";
    return "partial";
}

PlotKind::PlotKind(const string Kind, ReadInput& Inp)
{
    std::string cat=Kind==""? "Plot":"Plot_"+Kind;
    int p,line=0;
    string axStr="";
    double lb,ub;
    Inp.texdocuCategoryAdd("Plot","axis,lowerBound,upperBound,points,digits,usage,append,currents,interval,tagMax,tagMin",
                           R"tex(
                           Defines how to generate plots of spectra
                           \[
                           \si_{i_s\ldots i_{m-1}}(q_m,\ldots,q_{N})=\sum_{i_{s+1},\ldots,i_{m-1}}|\Psi_{i_0,\ldots,i_s,i_{m-1}}(q_m,\ldots,q_{N})|^2
                           \]
                           )tex","03,04");
    while (axStr!="BLANK") {
        line++;
        Inp.read(cat,"axis",axStr,"BLANK","(mutable) axis to plot (must match an axis name in discretization)",line)
                .texdocu(R"tex(
                         Name of plot axis or plot index.\\
                         Sequence of intput lines determinint the sorting or the plot axes.
                         The first axis is always the first plot axis, with \nameref{docu:Plot:usage} \lcode{grid}.
                         If a second axis with \lcode{usage=grid} is defined, the plot will be two-dimensional.
                         Further axes must be either summed or separate 2d plots, see \nameref{docu:Plot:usage}.
                         If an axis is not specified, it will be summed/integrated over the complete axis.
                         )tex");
        Inp.read(cat,"points",p,"0","(mutable) number of points, 0...automatic",line)
                .texdocu(R"tex(
                         Number of equidistant grid points on the axis including \lcode{lowerBound} and \lcode{upperBound}.
                         Not applicable (ingored?) for discrete axes.
                         )tex");
        Inp.read(cat,"lowerBound",lb,"0","(mutable) upper grid boundary",line)
                .texdocu(R"tex(
                         Lower boundary of the plot axis. Depending on \lcode{usage} plot, grid, separate columns, or summation will start from here.
                         )tex");
        Inp.read(cat,"upperBound",ub,tools::str(lb),"(mutable) lower grid boundary; defaults to upper=lower - single point or no grid conversion",line)
                .texdocu(R"tex(
                         Upper boundary of the plot axis, see \lcode{lowerBound}
                         )tex");

        string usage,defUse="grid";
        if(line>2)defUse="separate";
        if(p==0)defUse="separate";
        Inp.read(cat,"usage",usage,defUse,"(mutable) axis usage: grid...plot axis, separate...separate plots, sum...sum over coordinates",line,"",
                 "grid,separate,sum")
                .texdocu(R"tex(
                         Determines how the (generalized) coordinate $q_i$ in $\Psi(q_0,q_1,\ldots)$ is used.
                         Note that $q_i$ can be a discrete or continuous coordinate.
                         \begin{itemize}
                         \item[\lcode{grid}] plot $|\Psi|^2$ on an equidistant grid in \lcode{axis}, see \lcode{lowerEnd,upperEnd,points}.
                         \item[\lcode{separate}] make a separate plot for each point of the equidistant grid or for each discrete indes
                         \item[\lcode{sum}] sum or integrate over given section of axis
                         \end{itemize}
                         )tex");

        if(axStr=="BLANK")break;

        if(usage=="grid")_use.push_back("g");
        else if(usage=="separate")_use.push_back("p");
        else if(usage=="sum")_use.push_back("s");
        else ABORT("illegal Plot:usage="+usage+", allowed: grid,separate,sum");

        _axes.push_back(axStr);
        _points.push_back(p);
        _bounds.push_back(vector<double>(2,lb));
        _bounds.back()[1]=ub;
    }
    plotKinds["input"][Kind].reset(new PlotKind(*this));
}


void PlotKind::resize(const string Axis, unsigned int NewSize){
    size_t k=std::find(_axes.begin(),_axes.end(),Axis)-_axes.begin();
    if(k==_axes.size())ABORT(Sstr+"cannot resize PlotKind"+Axis+"not in"+_axes);
    _points[k]=NewSize;
}

// compact listing of available plots
const PlotKind *PlotKind::definePlot(const string Kind, const string Hierarchy){
    string pars;
    string kind=Kind=="" ? defaultKind(Hierarchy) : Kind;
    if(kind.find("[")!=string::npos){
        kind=Kind.substr(0,Kind.find("["));
        pars=tools::stringInBetween(Kind,"[","]");
    }

    if(plotKinds["input"].count(Kind)){
        const PlotKind* p=plotKinds["input"][Kind].get();
        std::vector<std::string> ax(tools::splitString(Hierarchy,'.'));
        for(std::string a: p->axes())
            if(std::find(ax.begin(),ax.end(),a)==ax.end())
                ABORT("Plot_"+kind+": cannot find axis "+a+" in "+Hierarchy);
        return p;
    }

    // dummies for 1d
    plotKinds["total"]["kX"].reset(new PlotKind({"kX"},{"g"},{0},vector<vector<double> >(1,{0.,0.})));
    plotKinds["total"]["kY"].reset(new PlotKind({"kY"},{"g"},{0},vector<vector<double> >(1,{0.,0.})));
    plotKinds["total"]["kZ"].reset(new PlotKind({"kZ"},{"g"},{0},vector<vector<double> >(1,{0.,0.})));

    //--------------------------------------------------------------------------------------------------
    // cartesian 2d
    plotKinds["sumX"]["specY.kX.kY"].reset(new PlotKind({"kY","kX"},{"g","s"},{0,0},vector<vector<double> >(2,{0.,0.})));
    plotKinds["sumY"]["specY.kX.kY"].reset(new PlotKind({"kX","kY"},{"g","s"},{0,0},vector<vector<double> >(2,{0.,0.})));
    plotKinds["sumX"]["specZ.kX.kZ"].reset(new PlotKind({"kZ","kX"},{"g","s"},{0,0},vector<vector<double> >(2,{0.,0.})));
    plotKinds["sumZ"]["specZ.kX.kZ"].reset(new PlotKind({"kX","kZ"},{"g","s"},{0,0},vector<vector<double> >(2,{0.,0.})));
    plotKinds["sumY"]["specZ.kY.kZ"].reset(new PlotKind({"kZ","kY"},{"g","s"},{0,0},vector<vector<double> >(2,{0.,0.})));
    plotKinds["sumZ"]["specZ.kY.kZ"].reset(new PlotKind({"kY","kZ"},{"g","s"},{0,0},vector<vector<double> >(2,{0.,0.})));

    //--------------------------------------------------------------------------------------------------
    // polar coordinates
    plotKinds["sumRn1"]["Phi1.Eta1.kRn1.Phi2.Eta2.kRn2"]
            .reset(new PlotKind({"kRn2"},{"g"},{0},vector<vector<double> >(1,{0.,0.})));
    plotKinds["partialRn1"]["Phi1.Eta1.kRn1.Phi2.Eta2.kRn2"]
            .reset(new PlotKind({"kRn1","Eta1"},{"g","p"},{0,0},vector<vector<double> >(2,{0.,0.})));

    plotKinds["sumRn2"]["Phi1.Eta1.kRn1.Phi2.Eta2.kRn2"]
            .reset(new PlotKind({"kRn1"},{"g"},{0},vector<vector<double> >(1,{0.,0.})));

    // Joint angular distributions
    if(Kind.substr(0,3)=="JAD")
        plotKinds["JAD"]["Phi1.Eta1.kRn1.Phi2.Eta2.kRn2"].reset(new PlotJAD(Kind));

    // total
    plotKinds["total"]["kRn"].reset(new PlotKind({"kRn"},{"g"},{0},vector<vector<double> >(1,{0.,0.})));
    plotKinds["total"]["Phi.Eta.kRn"]=plotKinds["total"]["kRn"];

    // total single-ionization spectra
    plotKinds["total"]["Phi1.Eta1.Phi2.Eta2.kRn1.Rn2"].reset(new PlotKind({"kRn1"},{"g"},{0},vector<vector<double> >(1,{0.,0.})));
    plotKinds["total"]["Phi1.Eta1.Phi2.Eta2.Rn1.kRn2"].reset(new PlotKind({"kRn2"},{"g"},{0},vector<vector<double> >(1,{0.,0.})));


    // multichannel - separate by channels
    plotKinds["total"]["Vec.Phi.Eta.kRn"]
            .reset(new PlotKind({"kRn","Vec"},{"g","p"},{0,0},vector<vector<double> >(2,{0.,0.})));

    // 2d
    plotKinds["total"]["Phi1.Phi2.Eta1.Eta2.kRn1.kRn2"]
            .reset(new PlotKind({"kRn1","kRn2"},{"g","g"},{0,0},vector<vector<double> >(2,{0.,0.})));
    plotKinds["total"]["kX.kY"]
            .reset(new PlotKind({"kX","kY"},{"g","g"},{0,0},vector<vector<double> >(2,{0.,0.})));
    plotKinds["total"]["kX1.kX2"]
            .reset(new PlotKind({"kX1","kX2"},{"g","g"},{0,0},vector<vector<double> >(2,{0.,0.})));

    // lm-partial waves
    plotKinds["partial"]["Phi.Eta.kRn"].reset(new PlotKind({"kRn","Eta","Phi"},{"g","p","p"},{0,0,0},vector<vector<double> >(3,{0.,0.})));

    // multichannel
    plotKinds["partial"]["Vec.Phi.Eta.kRn"].reset(new PlotKind({"kRn","Eta","Phi","Vec"},{"g","p","p","p"},{0,0,0,0},vector<vector<double> >(4,{0.,0.})));

    // l-partial waves, phi integrated
    plotKinds["lPartial"]["Phi.Eta.kRn"].reset(new PlotKind({"kRn","Eta"},{"g","p"},{0,0},vector<vector<double> >(2,{0.,0.})));
    plotKinds["lPartial"]["Phi1.Eta1.kRn1.Phi2.Eta2.kRn2"].reset(new PlotKind({"kRn1","Eta1","kRn2","Eta2"},{"g","p","p","p"},{0,0,0,0},
                                                                              vector<vector<double> >(4,{0.,0.})));
    // cut at Eta
    unsigned int nCut=tools::string_to_int(pars);
    plotKinds["cutEta"]["Phi.Eta.kRn"].reset(new PlotKind({"kRn","Phi","Eta"},{"g","g","p"},{0,91,nCut},{{0.,0.},{0.,2*math::pi},{-1.,1.}}));
    plotKinds["cutPhi"]["Phi.Eta.kRn"].reset(new PlotKind({"kRn","Eta","Phi"},{"g","g","p"},{0,91,nCut},{{0.,0.},{-1.,1.},{0.,2*math::pi}}));

    // integral over Phi
    plotKinds["intPhi"]["Phi.Eta.kRn"].reset(new PlotKind({"kRn","Eta"},{"g","g"},{0,101},{{0.,0.},{-1.,1.}}));
    plotKinds["intPhi"]["Phi1.Phi2.Eta1.Eta2.kRn1.kRn2"]
            .reset(new PlotKind({"kRn1","Eta1","kRn2","Eta2"},{"g","g","g","g"},{0,0,0,0},{{0.,0.},{-1.,1.},{-1.,1.},{-1.,1.}}));

    // cartesian kZ - input must have been converted already
    plotKinds["kZ"]["kZ"].reset(new PlotKind({"kZ"},{"g"},{0},{{0.,0.}}));
    plotKinds["kZ"]["kZ1.kZ2"].reset(new PlotKind({"kZ1","kZ2"},{"g","g"},{0,0},{{0.,0.},{0.,0.}}));

    double eta0(-1.),eta1(1.);
    unsigned int pts=181;
    if(kind.find("zone")==0){
        auto parts=tools::splitString(pars,',');
        if(parts.size()!=2)ABORT("need zone[th0,th1], got: "+kind);
        double th0=Algebra(parts[0]).constVal().real();
        double th1=Algebra(parts[1]).constVal().real();
        pts=std::max(21,std::min(int(pts),int(abs(th1-th0)/math::pi*181*4)));
        eta0=cos(th0);
        eta1=cos(th1);
    }
    plotKinds["zone"]["Phi.Eta.kRn"].reset(new PlotKind({"kRn","Eta"},{"g","s"},{0,pts},{{0.,0.},{eta0,eta1}}));

    if(plotKinds.count(kind)==0)ABORT("no such spectrum kind "+kind
                                      +"\navailable kinds are\n"+tools::listMapKeys(plotKinds));

    //------------------------------------------------------------------------------------------------------------

    // remove spec-levels
    string coors=Index::coordinates(Hierarchy);
    //    size_t specPos;
    //    while((specPos=coors.find("spec"))!=string::npos)
    //        coors.erase(coors.begin()+specPos,coors.begin()+coors.find(".",specPos)+1);

    // remove hybrid coordintates
    std::vector<std::string> cc=tools::splitString(coors,'.');
    coors="";
    for(auto c: cc)
        if(c.find("&")==std::string::npos)coors+=c+".";
    coors.pop_back();


    // permute until coordinates match
    for(auto c=plotKinds[kind].begin();c!=plotKinds[kind].end();c++){
        if(tools::permutation(tools::splitString(c->first,'.'),
                              tools::splitString(coors,'.')).size()>0){
            return c->second.get();
        }
    }

    ABORT("no spectral plot \""+kind+"\" defined for hierarchy "+Hierarchy+" (coordinates "+coors+")"
          +"\navailable coordinates are "+tools::listMapKeys(plotKinds[kind])
          +"\nsupplement desired plot in PlotKind"
          +"\nor define a Plot_user ");
}

std::vector<std::vector<double>> PlotKind::grid() const{
    std::vector<std::vector<double>> g,w;
    IndexGrid::gridWeight(g,w,_axes,_points,_bounds);
    return g;
}
std::vector<std::vector<double>> PlotKind::weig() const{
    std::vector<std::vector<double>> g,w;
    IndexGrid::gridWeight(g,w,_axes,_points,_bounds);
    return w;
}

std::string PlotJAD::format(){return "JAD[thMin{:thMax:n},eMin{:eMax:n}{~eV}{,phiMin{:phiMax:n}}]";}
PlotJAD::PlotJAD(std::string Def){
    string messageFormat="need format "+format()+", {...} is optional, got: "+Def;
    if(Def.substr(0,4)!="JAD[")ABORT(messageFormat);

    std::vector<std::string> defs=tools::splitString(tools::stringInBetween(Def,"[","]"),',');
    if(defs.size()<2)ABORT(messageFormat);

    // there is not need to have the grids here...
    eta1Grid=tools::rangeToGrid(defs[0],0,{"",":",""});
    for(double & theta: eta1Grid)theta=cos(theta);

    for(int k=0;k<181;k++)eta2Grid.push_back(cos(math::pi*k/180.));

    std::string range=defs[1].substr(0,defs[1].find_first_of(" ~"));
    kGrid=tools::rangeToGrid(range,0,{"",":",""});

    std::string uni=defs[1].find_first_of(" ~")!=std::string::npos?defs[1].substr(defs[1].find_first_of(" ~")+1):"DEFAULT_SYSTEM";
    for(double & ener: kGrid)ener=sqrt(2.*Units::convert(ener,uni));

    std::vector<double>phiGrid;
    phiGrid={0.,math::pi};
    if(defs.size()>2)phiGrid=tools::rangeToGrid(defs[2],0,{"",":",""});

    _axes={"Eta2","Phi2","Eta1","Phi1","kRn1","kRn2"};
    _use={"g","p","p","p","p","p"};
    _points={181,(unsigned int)phiGrid.size(),(unsigned int)(eta1Grid.size()),1,(unsigned int)(kGrid.size()),(unsigned int)(kGrid.size())};
    _bounds={{-1.,1.},{phiGrid[0],phiGrid.back()},{eta1Grid[0],eta1Grid.back()},{0.,0.},{kGrid[0],kGrid.back()},{kGrid[0],kGrid.back()}};

}
std::vector<std::vector<double>> PlotJAD::grid() const{
    std::vector<std::vector<double>> u(PlotKind::grid());
    for(size_t k=0;k<_axes.size();k++){
        if(_axes[k]=="Eta1")       u[k]=eta1Grid;
        if(_axes[k]=="Eta2")       u[k]=eta2Grid;
        if(_axes[k].find("kRn")==0)u[k]=kGrid;
    }
    return u;
}
std::string PlotJAD::strDef() const {
    std::string s("eV");
    s+=                      tools::str(Units::convert(0.5*std::pow(kGrid.front(),2),"au","eV"),2);
    if(kGrid.size()>1)s+=":"+tools::str(Units::convert(0.5*std::pow(kGrid.back(),2),"au","eV"),2)+":"+tools::str(kGrid.size());
    s+="_th";
    if(eta1Grid.front()!=-1. or eta1Grid.back()!=1.){
        s+=                        tools::string(SpectrumLinCom::strFractionOfPi(acos(eta1Grid.front())));
        if(eta1Grid.size()>1)s+=":"+tools::string(SpectrumLinCom::strFractionOfPi(acos(eta1Grid.back())))+":"+tools::str(eta1Grid.size());
    }
    return s;
}
