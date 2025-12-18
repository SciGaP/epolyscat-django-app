// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "pulse.h"

#include <memory>
#include "readInput.h"
#include "readInputList.h"
#include "abort.h"
#include "mpiWrapper.h"
#include "constants.h"
#include "parameters.h"
#include "printOutput.h"
#include "algebra.h"
#include "multiParam.h"

using namespace std;
using namespace tools;

Pulse Pulse::current=Pulse(); // declare variable for current pulse
bool Pulse::checkNew=false; // for now, check new pulse


void Pulse::setCurrent(ReadInput & Inp, bool Print, string Kind)
{
    std::string File="";
    if(Print)File=Inp.output()+Kind;
    setCurrent(Pulse(Inp,Kind),File);
}

static std::complex<double> _Asqu(std::complex<double> Time){return complex<double>(Pulse::current.Apot2(Time.real()),0.);}

static const Algebra* LaserAlgebraFactory(std::string FunctionName){
    if(FunctionName=="LaserAsq[t]")return new AlgebraExternalFunction(FunctionName,_Asqu);

    return 0;
}

void Pulse::setCurrent(const Pulse & Pulse, string PrintFile){
    current=Pulse;
    Parameters::add("LaserF0[t]",1.,Pulse::F0);
    Parameters::add("LaserFz[t]",1.,Pulse::F0);
    Parameters::add("LaserFx[t]",1.,Pulse::F1);
    Parameters::add("LaserFy[t]",1.,Pulse::F2);

    // alternate, more consistent naming
    Parameters::add("iLaserA0[t]",1.,Pulse::iAz);
    Parameters::add("iLaserAz[t]",1.,Pulse::iAz);
    Parameters::add("iLaserAx[t]",1.,Pulse::iAx);
    Parameters::add("iLaserAy[t]",1.,Pulse::iAy);

    Parameters::add("LaserAxx[t]",1.,Pulse::Axx);
    Parameters::add("LaserAyy[t]",1.,Pulse::Ayy);
    Parameters::add("LaserAzz[t]",1.,Pulse::Azz);
    Parameters::add("LaserAxy[t]",1.,Pulse::Axy);
    Parameters::add("LaserAxz[t]",1.,Pulse::Axz);
    Parameters::add("LaserAyz[t]",1.,Pulse::Ayz);

//    Parameters::add("LaserAsq[t]",1.,Pulse::Asqu);
    Parameters::add("LaserAsqHalf[t]",1.,[](double d) { return Pulse::Asqu(d) / 2.; });
    Parameters::add("iLaserAabs[t]",1.,Pulse::iAabs);
    Parameters::add("iLaserAdTheta[t]",std::complex<double>(0.,1.),Pulse::idThetaA);
    Parameters::add("iLaserFdTheta[t]",std::complex<double>(0.,1.),Pulse::idThetaF);

    Parameters::add("LaserFabs[t]",1.,Pulse::Fabs);
    // only printout laser file when on master process, or it will be very slow
    if(PrintFile!="" and MPIwrapper::isMaster())current.print(PrintFile, 0., 0., 2000);

    Algebra::addExternalFactory(LaserAlgebraFactory,"LaserAlgebraFactory");
}


/// electric field
std::complex<double> Pulse::F0(double Time){return current.Field(Time,0);}
std::complex<double> Pulse::F1(double Time){return current.Field(Time,1);}
std::complex<double> Pulse::F2(double Time){return current.Field(Time,2);}

/// vector potential (time-integral of electric field
std::complex<double> Pulse::iAz(double Time){return complex<double>(0,current.Apot(Time,0));}
std::complex<double> Pulse::iAx(double Time){return complex<double>(0,current.Apot(Time,1));}
std::complex<double> Pulse::iAy(double Time){return complex<double>(0,current.Apot(Time,2));}

std::complex<double> Pulse::iAabs(double Time){return complex<double>(0,std::sqrt(current.Apot2(Time)));}
std::complex<double> Pulse::idThetaA(double Time){return complex<double>(0,current.dThetaA(Time));}
std::complex<double> Pulse::idThetaF(double Time){return complex<double>(0,current.dThetaF(Time));}
//std::complex<double> Pulse::idTheta(double Time){return complex<double>(0,1.);}

std::complex<double> Pulse::Axx(double Time){return complex<double>(std::pow(current.Apot(Time,1),2),0.);}
std::complex<double> Pulse::Ayy(double Time){return complex<double>(std::pow(current.Apot(Time,2),2),0.);}
std::complex<double> Pulse::Azz(double Time){return complex<double>(std::pow(current.Apot(Time,0),2),0.);}

std::complex<double> Pulse::Axy(double Time){return complex<double>(current.Apot(Time,1)*current.Apot(Time,2),0.);}
std::complex<double> Pulse::Axz(double Time){return complex<double>(current.Apot(Time,1)*current.Apot(Time,0),0.);}
std::complex<double> Pulse::Ayz(double Time){return complex<double>(current.Apot(Time,2)*current.Apot(Time,0),0.);}

std::complex<double> Pulse::Asqu(double Time){return complex<double>(current.Apot2(Time),0.);}

std::complex<double> Pulse::Fabs(double Time){
    return std::sqrt(std::pow(current.Field(Time,0),2)
                     +std::pow(current.Field(Time,1),2)
                     +std::pow(current.Field(Time,2),2));
}

static bool allowJump=false;
static double angleLastA=0.;
static double angleLastF=0.;
double Pulse::thetaA(double Time){
    double vX=current.Apot(Time,1),vZ=current.Apot(Time,0);
    if(std::abs(vX)>1e-7 or std::abs(vZ)>1e-7)angleLastA=std::atan2(vX,vZ);
    return angleLastA;
}
double Pulse::thetaF(double Time){
    double vX=current.Field(Time,1),vZ=current.Field(Time,0);
    if(std::abs(vX)>1e-7 or std::abs(vZ)>1e-7)angleLastF=std::atan2(vX,vZ);
    return angleLastF;
}

Pulse::Pulse( ReadInput &Inp,string PulseName):t0F(0),t0A(0),t0D(0),pulseName(PulseName)
{
    //    comp.assign(3,vector<PulseSingle>());

    string shape;
    double peakIntensity=0.;
    double fwhm=0.;
    double lambda=0.;
    double phiCeo=0.;
    double t0=0.;
    double polAng;
    double aziAng=0.;
    double ellip=0.,ellipAng=0.;
    double photonEnergy=0.;
    double energyChirp=0.;
    double amplitude=0.;
    print0=DBL_MAX;
    print1=-DBL_MAX;

    Inp.read("Pulse","check",checkNew,"false","check new against old pulse checking code");

    // a series of pulses is read in
    int line=0;
    Inp.texdocuCategoryAdd(PulseName,"shape,lambda(nm),lambda@I/2(nm),I(W/cm2),FWHM,phiCEO,polarAngle,azimuthAngle",
                           R"tex(
                           Pulse (mostly Laser) is composed of single pulses specified in subesequent lines.
                           The wavelength (if any) in the first line determines the value of \lcode{OptCyc}
                           that can be used as a unit of time anywhere in the input. \ref{docu:Laser}
                           The vector potential $\vA$ of a single pulse has the form
                           \[
                           \vA(t) = \vec{\ep} \frac{A_0(t)}{\om} \sin(\om t - \phi_{CEO})
                           \]
                           $A_0$ is determined by \nameref{docu:Laser:shape},\nameref{docu:Laser:I_W_cm2_)}, and \nameref{docu:Laser:FWHM}
                           )tex","06,10,13,15");
    while (true) {
        line++;
        Inp.read(PulseName,"shape",shape,"BLANK","pulse shapes: gauss,cos8,cos4,cos2,sin2,flatTop,train,BLANK",line)
                .texdocu(R"tex(
                         Functional forms for $A_0$ (see \nameref{docu:category:Laser}). $T_{FWHM}$ is dertermined from \nameref{docu:Laser:FWHM}
                         such that the desired pulse duration FWHM in intensity is reached (Caution: really clean only for cosN, gauss).
                         \begin{itemize}
                         \item[gauss] $\exp(-\frac{(t-t_0)^2}{T_{FWHM}})$
                         \item[cos8] cos4, cos2: n=8,4,2 $\cos^n(-\frac{(t-t_0)^2}{T_{FWHM}})$
                         \\{\bf Note:} use of cos8 is recommended. Good compromise for a finite-duration, still similar to gauss.
                         Smaller powers tend to create artefacts for short pulses.
                         \item[flatTop][4 OptCyc] linear rise/constant/linear decay, rise and decay over the [4 OptCyc] (in this example).
                         Total pulse duration is controlled by \nameref{docu:Laser:FWHM}.
                         \item[train] Train of pulses (not currently maintained)
                         \end{itemize}
                         )tex");
        if(shape=="BLANK" or shape==ReadInput::notFound) break;
        if (PulseName=="Pulse") {
            Inp.read(PulseName,"peakAmpl",amplitude,ReadInput::noDefault,"peak amplitude",line);
            Inp.read(PulseName,"circFreq",photonEnergy,"0","circular carrier frequency (omega)",line);
        }
        else if(PulseName=="Laser"){
            Inp.read(PulseName,"I(W/cm2)",peakIntensity,ReadInput::noDefault,"peak intensity",line)
                    .texdocu(R"tex(
                             Peak intensity is defined as $[\max_t A_0(t)/\om]^2/2$, which is exact laser intensity in a cw pulse.
                             For $A_0$ see \nameref{docu:Laser:shape}.
                             )tex");
            Inp.read(PulseName,"lambda(nm)",lambda,"-1","wave length",line)
                    .texdocu(R"tex(
                             central wave length (defaul input in nm). See \nameref{docu:category:Laser}.
                             )tex");
            Inp.read(PulseName,"lambda@I/2(nm)",energyChirp,tools::str(lambda),"linear chirp - wave length at half intensity after peak",line)
                    .texdocu(R"tex(
                             linear chirp such that the frequency at half-intensity = value given here
                             )tex");
            photonEnergy=Units::convert(physics::planck_constant*physics::speed_of_light/(lambda*1.e-9),"SI|energy","au");
            energyChirp=Units::convert(physics::planck_constant*physics::speed_of_light/((energyChirp)*1.e-9),"SI|energy","au");
            Inp.read(PulseName,"ePhoton",photonEnergy,tools::str(photonEnergy),"photon energy",line)
                    .texdocu(R"tex(
                             alternative input of central wave length \nameref{docu:Laser:lambda_nm_} in terms of photon energy.
                             lambda(nm) and ePhoton are mutually exclusive.
                             )tex");
            Inp.read(PulseName,"ePhoton@I/2",energyChirp,tools::str(energyChirp),"linear chirp - photon energy at half intensity after peak",line)
                    .texdocu(R"tex(
                             linear chirp in terms of photon-energy instead of \nameref{docu:Laser:lambda_I_2_nm_}
                             )tex");
            energyChirp-=photonEnergy;
            if(abs(energyChirp)<photonEnergy*1.e-7)energyChirp=0.;

            Inp.exclude("Laser","Laser","ePhoton","lambda(nm)");
            Inp.exclude("Laser","Laser","ePhoton@I/2","lambda@I/2(nm)");
            amplitude=sqrt(Units::convert(peakIntensity,"W/cm2","au"))/photonEnergy;
        }
        else
            ABORT("no specifics defined for pulse type: "+PulseName);
        if(photonEnergy<0.)ABORT("must specify ePhoton or lambda(nm)");



        if(photonEnergy!=0. and line==1)
            Units::addUnit("OptCyc|time",2*math::pi*Units::convert(1/photonEnergy,"au|time","SI"));

        Inp.read(PulseName,"FWHM",fwhm,ReadInput::noDefault,"intensity FWHM of pulse",line);
        Inp.read(PulseName,"peak",t0,"0","time of peak envelope",line)
                .texdocu(R"tex(
                         shift laser peak time to given value (default is t=0)
                         )tex");
        Inp.read(PulseName,"phiCEO",phiCeo,"0","(rad) carrier-envelope offset",line,"","[-6.2832,6.2832]")
                .texdocu(R"tex(
                         $\phi_{CEO}$, see \nameref{docu:category:Laser}
                         )tex");
        Inp.read(PulseName,"polarAngle",polAng,"0","(deg) polar angle (measured from z-axis)",line)
                .texdocu(R"tex(
                         polar angle $\th$ (relative to z-axis) of the polarzation vector $\vec{\ep}$, see \nameref{docu:category:Laser}.
                         )tex");
        Inp.read(PulseName,"azimuthAngle",aziAng,"0","(deg) azimuthal angle (measured from x-axis)",line)
                .texdocu(R"tex(
                         polar angle $\phi$ (around z-axis) of the polarzation vector $\vec{\ep}$, see \nameref{docu:category:Laser}.
                         )tex");
        if((polAng!=0. and abs(polAng)<5.) or (aziAng!=0. and abs(aziAng)<5.))
            PrintOutput::warning(Str("ARE YOU SURE? polar and azimuthAngle are in degrees, found: ")+polAng+aziAng);
        Inp.read(PulseName,"ellip",ellip,"0","ellip%2 gives -1=left circular, 0=linear, +1=right cirular ",line)
                .texdocu(R"tex(
                         ellipticity, inquire for details, if needed.
                         )tex");
        if(checkNew and ellip!=0.){
            checkNew=false;
            PrintOutput::warning("checking cannot be done for elliptic polarization, switched off");
        }
        Inp.read(PulseName,"ellipAng",ellipAng,"0","angle of ellipticity plane arount main polarization axis",line);
        Inp.read(PulseName,"allowAngleJump",allowJump,"false","allow jump of rotation angle in rotating frem",line)
                .texdocu(R"tex(
                         When the field becomes 0, its direction becomes undefined. This may cause jump in the
                         rotation angle in a rotating coordinate system and produce artefacst. tRecX will stop.
                         If you want to admit it anyway, set to true.
                         )tex");

        Inp.obsolete(PulseName,"LaserCycles","enter FWHM in OptCyc instead",line);

        energyChirp/=(fwhm/2);
        if(abs(energyChirp)<1.e-10)energyChirp=0.;
        if(energyChirp!=0)ABORT("pulse chirp is broken - fix or do not use");

        sing.push_back(std::shared_ptr<PulseSingle>(
                           PulseSingle::factory(shape,amplitude,t0,fwhm,photonEnergy,phiCeo,
                                                ellip,ellipAng,polAng,aziAng,energyChirp)));
        print0=min(print0,sing.back()->print0);
        print1=max(print1,sing.back()->print1);
    }
    // possibly read parameter range
    //parameterRange(Inp);


    Inp.exclude("Pulse","Laser");
    if(sing.size()==0)PrintOutput::warning("no pulse given - dummy zero pulse");

    t0int=tBegin();
    if(t0int<-DBL_MAX/2){
        t0int=print0;
        PrintOutput::DEVwarning("infinite pulse, start integration at beginning of default print range");
    }

    intAsqT0=0.;
    vecA(t0A+1.e-10);
    vecF(t0F+1.e-10);
    intAsq(t0int+1.e-10);
}

double Pulse::tBegin() const {
    double t=DBL_MAX;
    for (size_t k=0;k<sing.size();k++)t=min(t,sing[k]->tBegin);
    return t;
}

double Pulse::tEnd() const {
    double t=-DBL_MAX;
    for (size_t k=0;k<sing.size();k++)t=max(t,sing[k]->tEnd);
    return t;
}

std::vector<double> Pulse::wavelengths(){
    std::vector<double> res;
    for(auto p: sing)
        res.push_back(p->omega==0?0:math::pi/p->omega);
    return res;
}
std::vector<double> Pulse::intensities(){
    std::vector<double> res;
    for(auto p: sing)
        res.push_back(std::pow(p->aPeak*p->omega,2));
    return res;
}
std::string Pulse::polarization(){
    std::string res="z-linear";
    for (auto p: sing){
        if(p->polarMain[1]!=0 or p->polarMain[2]!=0 or p->polarPerp[1]!=0 or p->polarPerp[2]!=0){
            res="general";
            break;
        }
    }
    return res;
}
double Pulse::duration(){
    double tmin=+DBL_MAX,tmax=-DBL_MAX;
    for (auto p: sing){
        tmin=std::min(tmin,p->t0-p->fwhm()/2.);
        tmax=std::max(tmax,p->t0+p->fwhm()/2.);
    }
    return tmax-tmin;
}

void Pulse::alignFieldZ(double Time){
    // align polarization such that pulse Field(Time) points in positive Z-direction
    VectorReal newZ(vecF(Time));
    double cosA=newZ[0]/sqrt(newZ.normSqu());
    double sinA=abs(sqrt(1.-cosA*cosA));
    if(newZ[1]<0.)sinA=-sinA;
    for(size_t k=0;k<sing.size();k++){
        VectorReal oldMain(sing[k]->polarMain),oldPerp(sing[k]->polarPerp);
        if((abs(oldMain[2])+abs(oldPerp[2]))>1.e-10)ABORT("only for polarization in 0-1 plane");
        sing[k]->polarMain =oldMain*cosA;
        sing[k]->polarMain-=oldPerp*sinA;
        sing[k]->polarPerp =oldMain*sinA;
        sing[k]->polarPerp+=oldPerp*cosA;
    }
}


void Pulse::output(string Title,string Operator) const{

    if(Title=="")PrintOutput::title(current.pulseName);
    else         PrintOutput::title(Title);

    PrintOutput::paragraph();
    if(apotMax()==0.){
        PrintOutput::message("no non-zero component");
        return;
    }

    PrintOutput::newRow();
    PrintOutput::rowItem("polAxis");
    PrintOutput::rowItem("Shape");
    if(current.pulseName=="Laser"){
        PrintOutput::rowItem("Intens(W/cm2)");
        PrintOutput::rowItem("lambda(nm)");
    } else {
        PrintOutput::rowItem("peakAmpl");
        PrintOutput::rowItem("photonEnergy");

    }
    PrintOutput::rowItem("FWHM");
    PrintOutput::rowItem("begin");
    PrintOutput::rowItem("end");
    PrintOutput::rowItem("phiCEO");
    vector<string> ax(3,"z");ax[1]="x";ax[2]="y";
    vector<string> ax_alt(ax);ax_alt[0]="0";
    double omega0;
    bool multicolor=false;
    for(unsigned int i=0;i<sing.size();i++){
        if(i==0)omega0=sing[i]->omega;
        multicolor=multicolor or omega0!=sing[i]->omega;
        if (abs(sing[i]->aPeak)<1.e-12)continue;
        PrintOutput::newRow();
        PrintOutput::rowItem(Str("","")+"["+sing[i]->polarMain.purge()+"]");
        PrintOutput::rowItem(sing[i]->shape);
        if(current.pulseName=="Laser"){
            PrintOutput::rowItem(Units::convert((sing[i]->aPeak*sing[i]->omega)*(sing[i]->aPeak*sing[i]->omega),"au|intensity","W/cm2"));
            PrintOutput::rowItem(physics::planck_constant/Units::convert(sing[i]->omega,"au|energy","SI")*physics::speed_of_light*1.e9);
        } else {
            PrintOutput::rowItem(sing[i]->aPeak);
            PrintOutput::rowItem(sing[i]->omega);
        }
        PrintOutput::rowItem((sing[i]->fwhm()));
        PrintOutput::rowItem((sing[i]->tBegin));
        PrintOutput::rowItem((sing[i]->tEnd));
        if(sing[i]->phiCeo!=0)PrintOutput::rowItem((sing[i]->phiCeo));
        else                  PrintOutput::rowItem("");
        if(abs(sing[i]->deltaOmega)>sing[i]->omega*1.e-10)
            PrintOutput::rowItem(Sstr+"chirp"+(sing[i]->omega+sing[i]->deltaOmega*sing[i]->fwhm()/2.));

    }
    // check consistency of operator with pulse
    if(Operator!="NONE"){
        for(size_t i=0;i<3;i++){
            if(vecAmax(i)<1.e-10)continue;
            if(Operator.find("LaserA"+ax[i])==string::npos
                    and Operator.find("LaserF"+ax[i])==string::npos
                    and Operator.find("LaserA"+ax_alt[i])==string::npos
                    and Operator.find("LaserF"+ax_alt[i])==string::npos
                    and Operator.find("<<MixedGaugeDipole")==string::npos
                    and Operator.find("<<Coriolis")==string::npos
                    and Operator.find("<<LaserLength")==string::npos
                    and Operator.find("<<LaserVelocity")==string::npos
                    ){
                PrintOutput::warning("Operator "+Operator+" does not couple to "+ax[i]+" component of pulse");
                if(ReadInput::main.found("Spectrum") and
                        not ReadInput::main.flag("DEBUGforceInconsistentSpectra",
                                                 "procede with computation although interation is incompatible with spectra"))
                    ABORT("all non-zero components must dipole-couple, else spectra will be inconsistent\n"
                          "if you want to run anyway, specify -DEBUGforceInconsistentSpectra");
            }
        }
    }
    PrintOutput::paragraph();
    PrintOutput::lineItem("ponderomotive potential (eV)",Units::convert(current.uPonderomotive(),"au","eV"));
    PrintOutput::newLine();
    if(multicolor){
        PrintOutput::paragraph();
        PrintOutput::message("Multi-color pulse, OptCyc is for 1st color, value="
                             +tools::str(Units::convert(1.,"OptCyc|time","au"),7)+" (au)");
    }



    if(sing.size()==0)PrintOutput::message("NO PULSE SPECIFIED");

}

string Pulse::str(unsigned int Brief) const {

    if(apotMax()==0.)return "noPulse";

    string s;
    if(Brief==0){
        for(unsigned int i=0;i<sing.size();i++){
            if (sing[i]->aPeak==0.)continue;
            if(i>0)s+=".";
            s+=sing[i]->shape+"["+tools::str(i)+"]";
            if(current.pulseName=="Laser"){
                s+=tools::str(int(sing[i]->fwhm()*sing[i]->omega/(2.*math::pi)*(1.+1.e-10)))+"c ";
                s+=tools::str(Units::convert((sing[i]->aPeak*sing[i]->omega)*(sing[i]->aPeak*sing[i]->omega),"au|intensity","W/cm2"),2);
                s+="@"+tools::str(physics::planck_constant/Units::convert(sing[i]->omega,"au|energy","SI")*physics::speed_of_light*1.e9,4);
            } else {
                s+=tools::str(sing[i]->fwhm())+"d ";
                s+=tools::str(sing[i]->aPeak,2);
                s+="@"+tools::str(sing[i]->omega,3);
            }
        }

    }
    else if(Brief>0){
        for(size_t l=0;l<sing.size();l++)
            if(sing[l]->aPeak>1.e-12)s+="\n"+sing[l]->str();
        if(s.length()>0)s.erase(s.begin(),s.begin()+1);
        else s=" --- empty pulse ----";
    }
    return s;
}

double Pulse::gettBegin()
{

    if(current.sing.size()==0)return -DBL_MAX; // empty pulse lasts forever
    double t = current.sing[0]->tBegin;
    for (size_t j=0; j<current.sing.size(); ++j)
        t=min(t,current.sing[j]->tBegin);
    return t;
}
double Pulse::gettEnd()
{
    if(current.sing.size()==0)return DBL_MAX; // empty pulse lasts forever
    double t = current.sing[0]->tEnd;
    for (size_t j=0; j<current.sing.size(); ++j)
        t=max(t,current.sing[j]->tEnd);
    return t;
}

double Pulse::uPonderomotive() const{
    double uP=0.;
    for(unsigned int l=0;l<sing.size();l++){
        if(sing[l]->omega==0.){ABORT("no meaningful ponderomotive potential for frequency=0");}
        uP+=0.25*sing[l]->aPeak*sing[l]->aPeak;
    }
    return uP;
}

std::complex<double> Pulse::aIntensityPowN::val(const std::complex<double> Q) const{
    return std::pow(aFieldPowN(2).integral(Q.real()-math::pi/_pulse->omegaMin(),Q.real()+math::pi/_pulse->omegaMin()),n);
}

double Pulse::intensityMoment(int N) const{
    return aIntensityPowN(this,N).integral(gettBegin(),gettEnd()).real();
}

double Pulse::nPhotonCrossSection(double Yield, int N, std::string Units) const{
    if(sing.size()!=1 or sing[0]->shape!="cos2" or N!=2){
        PrintOutput::warning(Sstr+"only for single cos2 pulse - dummy returned");
        return -1.;
    }
    if(omegaMax()!=omegaMin()){
        PrintOutput::warning(Sstr+"multiple frequencies"+omegaMax()+"..."+omegaMax()+" - no meaningful n-photon crossection");
        return -1.;
    }

    double tEff=(gettEnd()-gettBegin())*35./128.; // for a cos2 pulse only!

    //alternative conversion
    double ePhotB=Units::convert(omegaMax(),"au|energy","SI");
    double iPeakB=Units::convert(pow(sing[0]->aPeak*omegaMax(),2),"au|intensity","W/cm2");
    double tEffB=Units::convert(tEff,"au|time","SI");
    double sigB=Yield*std::pow(ePhotB/iPeakB,N)/tEffB;

    return sigB;


}

double Pulse::omegaMax() const {
    double res=0.;
    for(unsigned int l=0;l<sing.size();l++)
        res=max(res,sing[l]->omega);
    return res;
}
double Pulse::apotMax() const {
    double res=0.;
    for(unsigned int l=0;l<sing.size();l++)
        res=max(res,sing[l]->aPeak);
    return res;
}
double Pulse::omegaMin() const {
    double res=DBL_MAX;
    for(unsigned int l=0;l<sing.size();l++)
        res=min(res,sing[l]->omega);
    return res;
}

static double dTheta (double Ax, double Az, double Fx, double Fz){
    double eps=1.e-12;
    if(abs(Fx)+abs(Fz)<2*eps)return 0.;

    double Asq=Ax*Ax+Az*Az;
    if(Asq>eps*eps){
        double d=Az*Fx-Ax*Fz;
        if(abs(d)>abs(Az*Fx)*1.e-8)
            return d/Asq;
    }

    if(abs(Az*Fz)<eps*eps)return 0.;

    if(std::abs(std::abs(Az*Fx)-std::abs(Ax*Fz))<eps)return 0.;

    return DBL_MAX;
}

double Pulse::dThetaA(double Time) const{
    double res=DBL_MAX;
    // Time may hit an accidental 0, in that case use nearby value
    for(double dt:{0.,1.e-5,-1.e-5}){
        res=dTheta(Apot(Time+dt,1),Apot(Time+dt,0),Field(Time+dt,1),Field(Time+dt,0));
        if(res!=DBL_MAX)return res;
    }
    return 0.; // if not defined, return 0
    ABORT("failed to uniquely compute dTheta: t="+tools::str(Time));
}

double Pulse::dThetaF(double Time) const{
    double res=DBL_MAX;
    // Time may hit an accidental 0, in that case use nearby value
    for(double dt:{0.,1.e-5,-1.e-5}){
        res=dTheta(Field(Time+dt,1),Field(Time+dt,0),dField(Time+dt,1),dField(Time+dt,0));
        if(res!=DBL_MAX)return res;
    }
    return 0.; // if not defined, return 0
    ABORT("failed to uniquely compute dTheta: t="+tools::str(Time));
}

double Pulse::Apot2(double Time) const {
    double res=0;
    for(unsigned int i=0;i<3;i++)res+=std::pow(Apot(Time,i),2);
    return res;
}

VectorReal Pulse::vecA(double Time, bool Check) const {
    if(Time!=t0A or vecAt0.size()!=3){
        t0A=Time;
        // new time-recalculate
        vecAt0.assign(3,0.);
        for(size_t k=0;k<sing.size();k++){
            double env=sing[k]->aPeak*sing[k]->valEnv(Time-sing[k]->t0);
            double car=sing[k]->carrier(Time,0);
            if(sing[k]->omega!=0.)
                vecAt0+=sing[k]->polarMain*env*car*sing[k]->compMain;
            else
                vecAt0+=sing[k]->polarMain*env*sing[k]->compMain;

            if(sing[k]->polarPerp.size()!=0)
                vecAt0+=sing[k]->polarPerp*env*sqrt(1-car*car)*sing[k]->compPerp;

        }
    }
    return vecAt0;
}

double Pulse::vecAmax(int Component) const {
    VectorReal vecAmax;
    vecAmax.assign(3,0.);
    for(size_t k=0;k<sing.size();k++){
        double env=sing[k]->aPeak;
        vecAmax+=sing[k]->polarMain*env*sing[k]->compMain;

        if(sing[k]->polarPerp.size()!=0)
            vecAmax+=sing[k]->polarPerp*env*sing[k]->compPerp;

    }
    return std::abs(vecAmax[Component]);
}

VectorReal Pulse::vecF(double Time,bool Check) const {
    if(Time!=t0F or vecFt0.size()!=3){
        t0F=Time;
        vecFt0.assign(3,0.);
        for(size_t k=0;k<sing.size();k++){
            double der=sing[k]->aPeak*sing[k]->derEnv(Time-sing[k]->t0);
            double env=sing[k]->aPeak*sing[k]->valEnv(Time-sing[k]->t0);
            double car=sing[k]->carrier(Time,0),dCar=sing[k]->carrier(Time,1);
            vecFt0+=sing[k]->polarMain*sing[k]->compMain*(dCar*env+car*der); // derivative of carrier

            if(sing[k]->polarPerp.size()>0 and abs(sing[k]->compPerp)>1.e-10){
                if(sing[k]->deltaOmega!=0)ABORT("chirp not implemented for elliptic");
                vecFt0+=sing[k]->polarPerp*sing[k]->compPerp*(-car*env+dCar*der);
            }
        }
    }
    return vecFt0;
}

VectorReal Pulse::vecD(double Time,bool Check) const {
    if(Time!=t0D or vecDt0.size()!=3){
        t0D=Time;
        vecDt0.assign(3,0.);
        for(size_t k=0;k<sing.size();k++){
            double dde=sing[k]->aPeak*sing[k]->ddEnv(Time-sing[k]->t0);
            double der=sing[k]->aPeak*sing[k]->derEnv(Time-sing[k]->t0);
            double env=sing[k]->aPeak*sing[k]->valEnv(Time-sing[k]->t0);
            double car=sing[k]->carrier(Time,0),dCar=sing[k]->carrier(Time,1),ddCar=sing[k]->carrier(Time,2);
            vecDt0+=sing[k]->polarMain*sing[k]->compMain*(ddCar*env+2*dCar*der+car*dde); // derivative of carrier

            if(sing[k]->polarPerp.size()>0 and abs(sing[k]->compPerp)>1.e-10){
                ABORT("elliptic not implemented");
                if(sing[k]->deltaOmega!=0)ABORT("chirp not implemented for elliptic");
                vecDt0+=sing[k]->polarPerp*sing[k]->compPerp*(-car*env+dCar*der);
            }
        }
    }
    return vecDt0;
}

double Pulse::intAsq(double Time) const {
    if(Time!=t0int){
        intAsqT0+=aSq().integral(complex<double>(t0int),complex<double>(Time)).real();
        t0int=Time;
    }
    return intAsqT0;
}

void Pulse::print(string File,double T0, double T1, int Points){
    ofstream laser;
    vector<double>intF(3,0.);
    vector<double>intD(3,0.);
    laser.open((File).c_str());
    laser<<"# laser dipole field in 0,1,2 = z,x,y directions, and integrals"<<endl;
    laser<<"# dThetaA/dt angular velocity of A around y-coordinate"<<std::endl;
    laser<<"# dThetaF/dt angular velocity of F around y-coordinate"<<std::endl;
    laser<<"# Int(Field) and Int(dField/dt) can be used for consistency checks ?= A,F"<<std::endl;
    std::string head="#\n# time,|A|, dThetaA/dt, dThetaF/dt";
    for(int i: {0,1,2})head+=Str("","")+", Apot["+i+"], Field["+i+"], Int(Field["+i+"]), Int(dField["+i+"]/dt)";
    laser<<head<<endl;
    if(T0>=T1){T0=print0;T1=print1;};
    double dt=(T1-T0)/std::max(Points,1);
    for (double t=T0;t<=T1;t+=dt){
        laser<<setw(8)<<setprecision(10)<<t<<setw(18)<<sqrt(Apot2(t))<<setw(18)<<dThetaA(t)<<setw(18)<<dThetaF(t);
        for (unsigned int i=0;i<3;i++) {
            laser<<setw(13)<<setprecision(10)<<setw(18)<<Apot(t,i)<<" "<<setw(18)<<Field(t,i)
                <<" "<<setw(18)<<intF[i]<<setw(18)<<intD[i];
            intF[i]+=aField(i).integral(complex<double>(t),complex<double>(t+dt)).real();
            intD[i]+=aDfield(i).integral(complex<double>(t),complex<double>(t+dt)).real();
        }
        laser<<endl;
    }
    laser.close();
}
