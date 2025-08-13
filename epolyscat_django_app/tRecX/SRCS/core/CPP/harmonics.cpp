// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "harmonics.h"

#include <vector>

#include "abort.h"
#include "asciiFile.h"
#include "algebra.h"
#include "interpolate.h"
#include "readInput.h"
#include "vectorReal.h"
#include "interpolate.h"
#include "pulse.h"
#include "units.h"
#include "constants.h"
#include "operatorData.h"
#include "plot.h"
//#include "printOutput.h"

#ifdef _USE_FFTW_
#include <fftw3.h>
#endif

typedef double fftw_complex[2];

using namespace std;

string Harmonics::dipoleDefinitions(ReadInput &Inp, const string & Coor, const string &Hamiltonian,vector<string>&Names){

    Names.clear();
    vector<string> ax;
    if(Hamiltonian.find("<<MixedGauge")!=string::npos)
        ax={"Z","Y","X"};
    else {
        if     (Hamiltonian.find("iLaserA0[t]")!=string::npos)ax.push_back("Z"); // old component naming
        else if(Hamiltonian.find("iLaserAz[t]")!=string::npos)ax.push_back("Z");
        if(Hamiltonian.find("iLaserAy[t]")!=string::npos)ax.push_back("Y");
        if(Hamiltonian.find("iLaserAx[t]")!=string::npos)ax.push_back("X");
    }
    if(Inp.found("Dipole") and ax.size()==0 and Hamiltonian.find("<<Laser")==std::string::npos)
        PrintOutput::warning("no dipoles constructed, as no iLaserAx,y,z in Hamiltionan definition = "+Hamiltonian);
    string accel,velo,leng;
    Inp.texdocuCategoryAdd("Dipole","velocity,acceleration,length",
                           R"tex(
                           Time-dependent dipole expectation values will be computed in any of the specfied forms.
                           This works only for selected coordinate systems.
                           If the coordinate system is unknown, you can specify the explicit form of the operator.
                           The time-dependent expectation values will be used to compute harmonices.
                           )tex","07,17,18");
    Inp.read("Dipole","acceleration",accel,"dum","compute dipole expectation")
            .texdocu(R"tex(
                     most numerically robust where it can be compueted.
                     \lcode{Dipole: acceleration} will compute acceleration form for selected coordinate systems and potentials
                     )tex");
    Inp.read("Dipole","velocity",velo,"dum","compute dipole expectation")
            .texdocu(R"tex(
                     numerically rather robust and independent of potential.
                     )tex");
    Inp.read("Dipole","length",leng,"dum","compute dipole expectation")
            .texdocu(R"tex(
                     performs generally poor, cheap to compute and good to see just HOW poorly it works)tex");
    if(accel!="dum" and accel.find("DipAcc")!=0 and accel.find(":")!=6)ABORT("need format DipAccZ:<..definition>");

    // try to automatically determine dipole operator
    if(Inp.found("Dipole","acceleration")){
        if(accel=="dum"){
            if(Hamiltonian.find("<<Coulomb>>")==string::npos and Hamiltonian.find("<<PotSolid>>")==string::npos)
                ABORT(string("dipole acceleration only predefined for <<Coulomb>> and <<PotSolid>>")
                      +",            manually define for potential in\nHamiltonian = "+Hamiltonian
                      +",            use format \"Harmonics: acceleration=STRING\"");
            for(unsigned int k=0;k<ax.size();k++)Names.push_back("DipAcc"+ax[k]);
        }
        else Names.push_back(accel);
    }
    if(Inp.found("Dipole","velocity"))
        for(unsigned int k=0;k<ax.size();k++)Names.push_back("DipVel"+ax[k]);

    if(Inp.found("Dipole","length"))
        for(unsigned int k=0;k<ax.size();k++)Names.push_back("DipLen"+ax[k]);
    string s;
    if(Names.size()>0){
        s=Names[0];
        for(unsigned int k=1;k<Names.size();k++)s+=","+Names[k];
    }
    return s;
}

void Harmonics::addDipolesToPlot(std::string DipoleDefinitions, Plot & PlotDef, const std::vector<OperatorAbstract*> & OpList ){
    if(PlotDef.isEmpty() or DipoleDefinitions=="")return;
    vector<string> names,seps;
    tools::splitString(DipoleDefinitions,",",names,seps,"<[(",">])");
    for(size_t k=0;k<names.size();k++){
        vector<OperatorAbstract*>::const_iterator op=OpList.begin();
        for(;op!=OpList.end();op++){
            if((*op)->name.find(names[k].substr(0,names[k].find(':')))==0)break;
        }
        if(op==OpList.end())ABORT("did not find "+names[k]+" in list of operators");
        PlotDef.addOperator(*op);
    }
    PrintOutput::message("added dipoles to density plots: "+DipoleDefinitions);
}

// returns  the first dipole type found in Head
std::string dipoleColumn(std::string Head){
    vector<string> dName={
        "DipVelX","Op[Map_D(X)]",
        "DipVelY","Op[Map_D(Y)]",
        "DipVelZ","Op[Map_D(Z)]",
        "DipAccX",
        "DipAccY",
        "DipAccZ",
        "DipLenX",
        "DipLenY",
        "DipLenZ"
    };
    string dip;
    for(size_t k=0;k<dName.size();k++){
        if(dName[k].find("Dip")!=string::npos)dip=dName[k];
        if(Head.find(dName[k])!=string::npos)return dip;
    }
    return "";
}

static double zeroFunc(double T){return 0.;}
static double Ax(double T){return Pulse::iAx(T).imag();}
static double Ay(double T){return Pulse::iAy(T).imag();}
static double Az(double T){return Pulse::iAz(T).imag();}
static double Fz(double T){return Pulse::F0(T).real();}
static double Fx(double T){return Pulse::F1(T).real();}
static double Fy(double T){return Pulse::F2(T).real();}

static int cnt=0;
Harmonics::Harmonics(ReadInput &Inp, string Exten, string Sep, bool RowWise, double UnitT)
{


    typedef double(*pulseFunc)(double T);
    vector<pulseFunc>addFunc;
    vector<int>powOmega;

    if(0==cnt++)PrintOutput::title("COMPUTE HIGH HARMONIC SPECTRA");

    // try list of posible extensions
    vector<string>headLines,outHeader(1);
    vector<unsigned int> col;

    // read multiple colums
    AsciiFile dataFile(Inp.output()+Exten);
    headLines.clear();
    dataFile.readComments(headLines);

    if(headLines.size()>0 and dipoleColumn(headLines.back())!=""){
        vector<string>nam,sep;
        tools::splitString(headLines.back(),Sep,nam,sep);
        for(unsigned int l=0;l<nam.size();l++){
            string cnam=tools::cropString(nam[l]);
            string onam=dipoleColumn(cnam);
            if(onam!=""){
                col.push_back(l); // note: first "column" will be #
                outHeader[0]+=" <"+onam+">";
                if(     onam=="DipVelX")addFunc.push_back(Ax);
                else if(onam=="DipVelY")addFunc.push_back(Ay);
                else if(onam=="DipVelZ")addFunc.push_back(Az);
                else if(onam=="DipAccX")addFunc.push_back(Fx);
                else if(onam=="DipAccY")addFunc.push_back(Fy);
                else if(onam=="DipAccZ")addFunc.push_back(Fz);
                else if(onam=="DipLenX")addFunc.push_back(zeroFunc);
                else if(onam=="DipLenY")addFunc.push_back(zeroFunc);
                else if(onam=="DipLenZ")addFunc.push_back(zeroFunc);
                else ABORT("no function defined for dipole type "+onam);

                powOmega.push_back(0);
                if(onam.find("Vel")!=string::npos)powOmega.back()=2;
                if(onam.find("Len")!=string::npos)powOmega.back()=4;
            }
        }
    }

    if(col.size()==0){
        if(folder::exists(dataFile.name()))
            PrintOutput::warning(Sstr+"no dipole found on file root"+dataFile.name()
                                 +"--- cannot compute harmonics from that file");
        return;
    }

    vector<double> time;
    vector<vector<double> > cols;
    vector<Eigen::MatrixXd> mTime;
    vector<Eigen::MatrixXd> mCols;
    vector<unsigned int> ax(1,0);
    if(Exten=="wf")ax.push_back(1);
    dataFile.readCols(mTime,ax,Sep,RowWise);

    mTime[0]*=UnitT;
    for(int k=0;k<mTime[0].rows();k++)time.push_back(mTime[0](k,0));
    dataFile.readCols(mCols,col,Sep,RowWise);


    // get some higher power of 2
    // ==============================================================
    // NOTE: there appears to be some aliasing effect
    // that creates equidistant large spikes in the spectra
    // seems to depend on the combination of
    // sample points, interpolation polynomial, and FFT points
    // keeping with FFT points well above sample points, it disappeared
    // ================================================================

    ReadInput::main.read("Harmonics","log2(Points)",nPoints,"16","use 2^log2points points",1,"log2points");
    nPoints=pow(2,nPoints);
    if(nPoints/time.size()>4)
        PrintOutput::warning(Str("oversampling by"," ")+nPoints+"points, with only"+time.size()+"dipoles points on file");

    double xmin=-DBL_MAX;
    double xmax= DBL_MAX;
    string range;
    if(Exten=="wf"){
        ReadInput::main.read("Harmonics","sumRange",range,"-infty,infty","sum local harmonics over this range",1,"sumRange");
        vector<string>sx=tools::splitString(tools::stringInBetween(range,"[","]"),',');
        if(sx.size()!=2)ABORT(Str("need input as -sumRange=[-3,4], got: -sumRange=")+range);
        xmin=tools::string_to_double(sx[0]);
        xmax=tools::string_to_double(sx[1]);
        PrintOutput::message(Str("summing local responses over [","")+xmin+","+xmax+"]");
    }

    // true resolution is limited by propagation time
    double dOmega=2*math::pi/(time.back()-time.front());

    // maximal frequency (absolute limit by time-spacing)
    double maxOmega=2.+10.*Pulse::current.uPonderomotive();
    maxOmega=min(maxOmega,dOmega*time.size());

    // rearrange for interpolation, smoothly turn off after end of pulse
    Algebra*read;
    if(time.back()>Pulse::current.gettEnd()){
        read=new Algebra ("trunc["+tools::str(Pulse::current.gettEnd())+","+tools::str(time.back())+"]");
    } else {
        if(time.back()-Pulse::current.gettEnd()>1.e-3*Units::convert(1.,"OptCyc"))
            PrintOutput::warning("dipole data end at "+tools::str(time.back())
                                 +",  before end of pulse at "+tools::str(Pulse::current.gettEnd()));
        read=new Algebra("1.");
    }

    // select times as needed
    vector<VectorReal> rows;
    std::vector<double> tSub;


    // required time-spacing
    double dTime=(time.back()-time[0])/nPoints;
    double tNext=time[0]-dTime;
    // select values at the coarsest meaningful grid
    for(unsigned int l=0;l<mCols[0].rows();l++){
        if(time[l]>tNext){
            tSub.push_back(time[l]);
            tNext=tSub.back()+dTime/2;
            rows.push_back(VectorReal());
            for(unsigned int k=0;k<mCols.size();k++){
                double addTerm=addFunc[k](tSub.back());
                for(unsigned int m=0;m<mCols[0].cols();m++){
                    rows.back().push_back(mCols[k](l,m)*read->val(tSub.back()).real()+addTerm);
                }
            }
        }
    }
    PrintOutput::message(Str("selected"," ")+tSub.size()+"out of"+time.size()+"time points (<dt> ~"+(tSub.back()-tSub[0])/tSub.size()+")");

    size_t k=Interpolate<double,VectorReal>::notIncreasing(tSub,rows,true);
    if(k!=0)ABORT("non-increasing support points "+tools::str(time[k-1],12)+" > "+tools::str(time[k],12));


    // interpolate
    VectorReal vals(rows[0].size());
    cols.assign(vals.size(),vector<double>());
    InterpolateNewton<double,VectorReal> intTab(tSub,rows,2);
    double t=tSub.front();
    dTime=(tSub.back()-tSub.front())/nPoints;
    vector<double> tGrid;
    for(int k=0;k<nPoints;k++,t+=dTime){
        tGrid.push_back(t);
        vals=intTab.val(t,vals);
        for(unsigned int l=0;l<cols.size();l++)cols[l].push_back(vals[l]);
    }


#ifdef _USE_FFTW_
    // fundamental frequency
    double omega=2.*math::pi/Units::convert(1.,"OptCyc","au");

    // limit printout to reasonable range of frequencies
    vector<double>oGrid;
    for(unsigned int k=0;k*dOmega<maxOmega;k++)oGrid.push_back(dOmega/omega*k);

    // fourier transform
    fftw_complex *in, *out;
    fftw_plan p(0);
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * cols[0].size() );
    out =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * cols[0].size() );

    vector<vector<double> > colsX(mCols.size());
    vector<vector<double> > colsSum;
    //    double qnrm=pow(1./double(nPoints),2);
    double qnrm=pow(dTime,2);
    for(int k=0,km=0;k<mCols.size();k++){
        colsX[k].resize(mCols[0].cols()*oGrid.size());
        vector<complex<double> > sumK(oGrid.size(),0.);
        for(int m=0;m<mCols[k].cols();m++,km++){
            if(cols[km].size()!=cols[0].size())ABORT(Str("unequal length columns")+cols[0].size()+cols[km].size()+"at"+km);

            p = fftw_plan_dft_1d(cols[0].size(), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

            for(fftw_complex* c=in;c<in+cols[0].size();c++){(*c)[0]=0.;(*c)[1]=0.;}
            for (unsigned int l=0;l<cols[km].size();l++)
            {
                in[l][0]=cols[km][l];
                in[l][1]=0;
            }
            fftw_execute(p); /* repeat as needed */

            // local responses
            for (unsigned int l=0;l<oGrid.size();l++){
                colsX[k][m*oGrid.size()+l]=std::norm(complex<double>(out[l][0],out[l][1]))*qnrm*pow(dOmega*l,powOmega[k]);
            }
            // coherently sum up the contributions
            if(Exten=="wf" and
                    tools::doubleInside(mTime[1](0,m),xmin,xmax)){
                for (unsigned int l=0;l<oGrid.size();l++){
                    sumK[l]+=complex<double>(out[l][0],out[l][1]);
                }
            }
        }
        if(Exten=="wf"){
            colsSum.push_back(vector<double>(oGrid.size(),0.));
            for(int l=0;l<oGrid.size();l++){
                colsSum.back()[l]=std::norm(sumK[l])*qnrm*pow(dOmega*l,powOmega[k]);
            }
        }
    }
    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);

    outHeader.insert(outHeader.begin(),2,string());
    outHeader[0]="       Omega0 = "+tools::str(omega)+" (au)";
    outHeader[1]="Upondermotive = "+tools::str(Pulse::current.uPonderomotive())+" (au)";
    outHeader[2]="HarmonicOrder "+outHeader[2];

    // write to file
    AsciiFile harmonicFile(Inp.output()+"harmonics_"+Exten);
    PrintOutput::title("harmonics on file "+harmonicFile.name());

    if(mTime.size()==1){
        colsX.insert(colsX.begin(),oGrid);
        harmonicFile.writeComments(outHeader);
        harmonicFile.writeCols(colsX);
    }
    else  if(mTime.size()==2){
        // sum of 2d data
        colsSum.insert(colsSum.begin(),1,vector<double>(oGrid.begin(),oGrid.end()));
        AsciiFile sumFile(Inp.output()+"harmonics_sumX");
        vector<string> sumHeader(outHeader);
        sumHeader.insert(sumHeader.begin(),"local dipoles summed over "+range);
        sumFile.writeComments(sumHeader);
        sumFile.writeCols(colsSum);
        PrintOutput::message("sum of local harmonics on file "+sumFile.name());

        // 2d data
        outHeader[2]="Coor "+outHeader[2];
        harmonicFile.blankIfChange(1);
        colsX.insert(colsX.begin(),2,vector<double>());
        for(int k=0;k<mTime[0].cols();k++)colsX[0].insert(colsX[0].end(),oGrid.begin(),oGrid.end());
        for(int k=0;k<mTime[1].cols();k++)colsX[1].insert(colsX[1].end(),oGrid.size(),mTime[1](0,k));
        harmonicFile.writeComments(outHeader);
        harmonicFile.writeCols(colsX);


    }
    else
        ABORT(Str("data can only be only 1d or 2d, found")+mTime.size());

#else
    PrintOutput::warning("only dipoles will be computete - for harmonics compile with -D_USE_FFTW_");
#endif


    if(Exten=="expec"){
        AsciiFile dipoleFile(Inp.output()+"dipoles");
        PrintOutput::title("dipoles at equidistant time grid on file "+dipoleFile.name());
        vector<string> dipHeader(outHeader);
        dipHeader.back()="Time "+outHeader.back();
        cols.insert(cols.begin(),1,tGrid);
        dipoleFile.writeComments(dipHeader);
        dipoleFile.writeCols(cols,12);
    } else if(Exten=="wf"){
        AsciiFile dipoleFile(Inp.output()+"sumDip");
        PrintOutput::message("sum of local dipoles on "+dipoleFile.name());
        vector<string> dipHeader(outHeader);
        dipHeader.back()="Time "+outHeader.back();
        vector<vector<double> > sumDip;
        sumDip.push_back(time);
        for(size_t k=0;k<mCols.size();k++){
            sumDip.push_back(vector<double>(sumDip[0].size(),0.));
            for(int c=0;c<mCols[k].cols();c++){
                for (int l=0;l<mCols[k].rows();l++){
                    sumDip.back()[l]+=mCols[k](l,c);
                }
            }
        }
        dipoleFile.writeComments(dipHeader);
        dipoleFile.writeCols(sumDip,12);

    }


}

void Harmonics::compute(int argc,char* argv[]){



    if(argc<2 or argv[1][0]=='-'){
        AsciiFile("Harmonics.doc").copy(cout);
        PrintOutput::paragraph();
        PrintOutput::message("usage: >Harmonics runDir/0123 [flags]");
        exit(0);
    }

    ReadInput::openMain(ReadInput::flagOnly,argc,argv);

    Units::setDefault("au");

    ReadInput dataInp(string(argv[1])+"/inpc",argc,argv,false);
    dataInp.setUnits("au");

    Pulse::read(dataInp,true);

    string file;
    ReadInput::main.read("Harmonics","file",file,"all","which file to use: expec,wf,all...(default)",1,"file");

    // data format on various files
    vector<string> fileExt={"wf","expec"};
    vector<string> fileSep={ ","," "};
    vector<bool> rowWise={true,false};
    vector<double> unitTime={0.,1.};
    if(Units::isDefined("OptCyc"))unitTime[0]=Units::convert(1.,"OptCyc");


    for(size_t k=0;k<fileExt.size();k++){
        if(file!="all" and fileExt[k]!=file)continue;
        Harmonics(dataInp,fileExt[k],fileSep[k],rowWise[k],unitTime[k]);
    }
    ReadInput::main.writeDoc();
}
