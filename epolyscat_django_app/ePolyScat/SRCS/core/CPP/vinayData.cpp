// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "vinayData.h"

#include <vector>
#include "readInput.h"
#include "index.h"
#include "coefficients.h"
#include "algebra.h"
#include "interpolate.h"
#include "basisAbstract.h"

using namespace std;

vector<const Index*> VinayData::_allIdx;

void VinayData::readVec(ifstream &Inp, int Size, vector<double> & Vec){
    Vec.resize(Size);
    for(int k=0;k<Size;k++)Inp>>Vec[k];
}
VinayData::~VinayData(){delete _c;}

/// converts NpiM to N*pi/M
double strToAngle(std::string NumPiDenom){
    if(NumPiDenom.find("pi")==string::npos)return Algebra::realConstant(NumPiDenom);
    double res=math::pi;
    if(NumPiDenom.find("pi")!=0)res*=tools::string_to_double(NumPiDenom.substr(0,NumPiDenom.find("pi")));
    if(NumPiDenom.find("pi")+2!=NumPiDenom.length())res/=tools::string_to_double(NumPiDenom.substr(NumPiDenom.find("pi")+2));
    return res;
}


VinayData::VinayData(std::string & FileName)
{
    // default ranges for channels
    // one XUV peak must be included for phase fixing
    std::map<std::string,std::string> ranges;
    ranges["0"]="0.53,0.78,100";
    ranges["1"]="0.53,0.78,100";
    ranges["2"]="0.30,0.55,100";
    ranges["3"]="0.30,0.55,100";
    ranges["4"]="0.25,0.50,100";
    ranges["5"]="-infty,infty,100";

    std::map<std::string,std::vector<double>> sideBand;
    sideBand["0"]={0.73,0.78};
    sideBand["1"]=sideBand["0"];
    sideBand["2"]={0.49,0.55};
    sideBand["3"]=sideBand["2"];
    sideBand["4"]={0.44,0.50};

    string range,chan;
    bool interpolate;
    ReadInput::main.read("VinayData","interpolate",interpolate,"false","interpolate in grid",1,"interpolate");
    ReadInput::main.read("VinayData","channel",chan,ReadInput::noDefault,"channel (file name will be ampl_channel",1,"channel");
    ReadInput::main.read("VinayData","kRange",range,ranges[chan],"sub-range for plot -kRange=[kMin,kMax,pts], pts is optional",1,"kRange");
    FileName+="/ampl_"+chan;
    _sideBand=sideBand[chan];

    if(range=="-infty,infty,100" and interpolate)ABORT("specify -kRange=[kMin,kMax,pts] for interpolation");

    vector<string> vRange=tools::splitString(tools::stringInBetween(range,"[","]"),',');
    vector<double> kRange(2,0.);
    kRange[0]=Algebra::realConstant(vRange[0]);
    kRange[1]=Algebra::realConstant(vRange[1]);
    int maxSize(100);
    if(vRange.size()==3)maxSize=int(Algebra::realConstant(vRange[2]));

    if(ReadInput::main.flag("Emax","(forbidden with VinayData)") or
            ReadInput::main.flag("nR","(forbidden with VinayData)") )
        ABORT("cannot change parameters of VinayData");

    ifstream inp(FileName.c_str(),std::ios_base::in);
    if(not inp.is_open())ABORT("failed to open Vinay file "+FileName);


    string s=FileName.substr(0,FileName.rfind("/"));
    while(s.back()=='/')s=s.substr(0,s.length()-1);
    s=s.substr(s.rfind("/")+1);
    vector<string> sPars=tools::splitString(s,'_');
    _polarAngle=math::pi*strToAngle(sPars[0])/180.;
    _phiCEO=strToAngle(sPars[1]);
    _channel=FileName.substr(FileName.rfind("ampl_")+5);

    size_t code,sizePhi,sizeEta,sizeK;
    inp>>code;
    if(code!=1)ABORT(Str("was expecting code=1, is =")+code);
    inp>>sizePhi;
    inp>>sizeEta;
    inp>>sizeK;

    vector<double> gPhi,gEta,wEta,gK;
    readVec(inp,sizePhi,gPhi);
    readVec(inp,sizeEta,gEta);
    readVec(inp,sizeEta,wEta);
    readVec(inp,sizeK,gK);


    size_t kMin,kMax,sampl;
    vector<double> gKsampl,gKfull;
    if(not interpolate){
        kMin=std::lower_bound(gK.begin(),gK.end(),kRange[0])-gK.begin();
        kMax=std::upper_bound(gK.begin(),gK.end(),kRange[1])-gK.begin();
        sampl=(kMax-kMin-1)/maxSize+1;
        for(size_t k=0;k<gK.size();k++)
            if(k%sampl==0 and kMin<=k and k<kMax)gKsampl.push_back(gK[k]);
    }
    else {
        kMin=std::lower_bound(gK.begin(),gK.end(),kRange[0]-0.1*(kRange[1]-kRange[0]))-gK.begin();
        kMax=std::upper_bound(gK.begin(),gK.end(),kRange[1]+0.1*(kRange[1]-kRange[0]))-gK.begin();
        sampl=1;
        gKsampl=tools::rangeToGrid(range);
        gKfull=std::vector<double>(gK.begin()+kMin,gK.begin()+kMax);
    }

    _fullGrid=Str("")+sizePhi+sizeEta+sizeK+gK.front()+gK.back();
    _sampleGrid=Str("")+sizePhi+sizeEta+gKsampl.size()+gKsampl.front()+gKsampl.back();
    if(interpolate)_sampleGrid+=" (interpolated)";

    _idx=getIndex({{gPhi,{}},{gEta,wEta},{gKsampl,{}}},vector<const BasisAbstract*>(),0);
    _allIdx.push_back(_idx); // keep track of pointer (may clean up if desired)

    _c=new Coefficients(_idx);
    if(not _c->anyData())DEVABORT("need this continguous");

    complex<double> a;

    if(_c->childSize()!=sizePhi or _c->descend(2)->levelSize()!=sizePhi*sizeEta)
        DEVABORT(Str("incorrect index structure")+_c->childSize()+sizePhi+_c->descend(2)->levelSize()+sizeEta);

    for(Coefficients* cc=_c->descend(2);cc!=0;cc=cc->nodeNext()){
        if(not interpolate){
            for(size_t k=0,ks=0;k<sizeK;k++){
                inp>>a;
                if(k%sampl==0 and kMin<=k and k<kMax)cc->anyData()[ks++]=a;
            }
        }
        else {
            vector<complex<double> >cFull;
            for(size_t k=0;k<sizeK;k++){
                inp>>a;
                if(kMin<=k and k<kMax)cFull.push_back(a);
            }
            InterpolateNewton<double,complex<double> >intPol(gKfull,cFull,4);
            for(size_t k=0;k<gKsampl.size();k++)
                intPol.val(gKsampl[k],cc->anyData()[k]);
        }
    }
}

Index* VinayData::getIndex(std::vector<std::vector<std::vector<double> > > Grids, std::vector<const BasisAbstract*> Bas, int Level){
    vector<string> axNam={"Phi","Eta","kRn","NONE"};

    if(Bas.size()==0){
        string basDef;
        for(int l=0;l<3;l++){
            basDef="Grid:";
            if(l<2)basDef="GridQuad:";
            for(double g: Grids[l][0])basDef+=tools::str(g)+",";
            for(double w: Grids[l][1])basDef+=tools::str(w)+",";
            if(axNam[l]=="Phi"){
                string phiWeig=tools::str(2*math::pi/double(Grids[l][0].size()))+",";
                for(size_t k=0;k<Grids[l][0].size();k++)basDef+=phiWeig;

            }
            basDef.resize(basDef.length()-1);
            Bas.push_back(BasisAbstract::factory(basDef));
        }
    };

    Index * iLevel(new Index());
    iLevel->setAxisName(axNam[Level]);
    iLevel->setBasis(Bas[Level]);
    if(iLevel->axisName()=="kRn")iLevel->setFloor();
    for(size_t k=0;k<iLevel->basis()->size();k++)
        if(Level<2)iLevel->childAdd(getIndex(Grids,Bas,Level+1));
        else       iLevel->leafAdd();

    iLevel->sizeCompute();
    return iLevel;
}
