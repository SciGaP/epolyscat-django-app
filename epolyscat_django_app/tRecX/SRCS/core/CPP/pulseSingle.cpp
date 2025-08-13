// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "pulseSingle.h"

#include "math.h"
#include "constants.h"
#include "pulse.h"


#ifdef _USE_FFTW_
#include "fft.h"
#endif

using namespace std;

PulseSingle* PulseSingle::factory(string Enve, double APeak, double T0, double Dur, double Ome, double Phi,
                                  double Ellip, double AngleEllip, double PolAng, double AziAng, double DeltaOmega){

    if(Enve.find("flatTop")==0)  return new PulseFlatTop(Enve,APeak,T0,Dur,Ome,Phi,Ellip,AngleEllip,PolAng,AziAng,DeltaOmega);
    else if (Enve.find("cos")==0)return    new PulseCosN(Enve,APeak,T0,Dur,Ome,Phi,Ellip,AngleEllip,PolAng,AziAng,DeltaOmega);
    else if (Enve=="gauss")      return   new PulseGauss(Enve,APeak,T0,Dur,Ome,Phi,Ellip,AngleEllip,PolAng,AziAng,DeltaOmega);
    else ABORT("undefined pulse: "+Enve);
}

PulseSingle::PulseSingle(string Enve, double APeak, double T0, double Dur, double Ome, double Phi,
                         double Ellip, double AngleEllip, double PolAng, double AziAng,double DeltaOmega)
    :shape(Enve),aPeak(APeak),t0(T0),omega(Ome),phiCeo(Phi),print0(t0),print1(t0),tBegin(t0),tEnd(t0),deltaOmega(DeltaOmega){

    // polarization and ellipticity, new style
    polarMain.clear();
    polarMain.push_back(cos(math::pi*PolAng/180.)                          ); //z-axis
    polarMain.push_back(sin(math::pi*PolAng/180.)*cos(math::pi*AziAng/180.)); //x-axis
    polarMain.push_back(sin(math::pi*PolAng/180.)*sin(math::pi*AziAng/180.)); //y-axis

    // get perpendicular unit vector in polar-to-z plane: |z>-|p><p|z>
    polarPerp.clear();
    polarPerp=polarMain;
    polarPerp*=-polarMain[0];
    polarPerp[0]+=1.;
    if(polarPerp.normSqu()<1.e-10)polarPerp[1]=1.; // main in z-direction

    // renormalize
    polarMain*=1./sqrt(polarMain.normSqu());
    polarPerp*=1./sqrt(polarPerp.normSqu());

    if(AngleEllip!=0)ABORT("ellipticity implemented only in plane to z-axis");
    resetParameter("ellip",Ellip);
}

double PulseSingle::carrier(double T, int Derivative) const {
    double res;
    if(Derivative==0){
        res=sin((T-t0)*(omega+deltaOmega*(T-t0))+phiCeo);
    }else if(Derivative==1){
        res=cos((T-t0)*(omega+deltaOmega*(T-t0))+phiCeo)*(omega+2*deltaOmega*(T-t0));
    }else if(Derivative==2){
        if(deltaOmega!=0.)ABORT("cannot use chirp with derivative of field");
        res=-sin((T-t0)*omega+phiCeo)*std::pow(omega,2);
    } else {
        res=0.;
        ABORT("can have at most 2nd derivative of carrier");
    }
    return res;
}

double PulseSingle::fwhmTransformLimit() const {
    ABORT("temporarily disabled");
    //    deltaOmega=0.;
    //    return fwhm()*p0.spectralWidth()/spectralWidth();
}

double PulseSingle::spectralWidth() const {
#ifndef _USE_FFTW_
    PrintOutput::DEVwarning("compiled w/o FFTW - cannot compute spectral width");
    return 1.;
#else
    // evaluate on good time grid
    int npow=14;
    double t=2*print0-print1,dt=2*(print1-print0)/pow(2,npow);
    vector<complex<double> > puls;
    for(int k=0;k<pow(2,npow);k++,t+=dt)
        puls.push_back(carrier(t,false)*valEnv(t-t0));

    // Fourier transform
    puls=Fft(puls.size(),true).transform(puls);

    // compute variance
    double integ=0,omInt=0.,sqInt=0.;
    double om=0,dom=2*math::pi/(dt*puls.size());
    for(int k=0;k<puls.size();k++,om+=dom){
        integ+=norm(puls[k]);
        omInt+=norm(puls[k])*om;
        sqInt+=norm(puls[k])*om*om;
    }
    return sqrt(sqInt/integ-pow(omInt/integ,2));
#endif
}

void PulseSingle::resetParameter(string Name, double Value){
    if(Name=="phiCEO")phiCeo=Value;
    else if(Name=="peak")t0=Value;
    else if(Name=="ellip"){
        compMain=cos(Value*math::pi/4.);
        compPerp=sin(Value*math::pi/4.);
    }
    else
        ABORT("parameter \""+Name+"\" cannot be reset");
}

string PulseSingle::str() const{
    string s=shape;
    double inten=Units::convert(pow(aPeak*omega,2),"DEFAULT_SYSTEM","W/cm2");
    s+="\t"+tools::str(inten,3);
    s+="\t"+tools::str(t0,3);
    s+="\t"+tools::str(omega,3);
    s+="\t"+tools::str(phiCeo,3);
    if(deltaOmega!=0){
        s+="\t"+tools::str(deltaOmega,3);
        s+="\t"+tools::str(fwhm(),3);
        //        s+="\t"+tools::str(fwhmTransformLimit(),3);
    }
    return s;
}

PulseFlatTop::PulseFlatTop(std::string Enve, double APeak, double T0, double Dur, double Ome, double Phi,
                           double Ellip, double AngleEllip, double PolAng, double AziAng, double DeltaOmega)
    :PulseSingle(Enve,APeak,T0,Dur,Ome,Phi,Ellip,AngleEllip,PolAng,AziAng,DeltaOmega)
{
    string ramp=tools::stringInBetween(Enve,"[","]");
    // flat top needs an extra parameter
    if(ramp==shape)ABORT("define ramp in the format \"flatTop[2 OptCyc]\", is: "+shape);
    vector<string> parts=tools::splitString(ramp,' ');
    if(parts[0].find_first_not_of("0123456789e.+-")!=string::npos)
        ABORT("illegal ramp format: "+Enve+", did you forget the blank between value and units?");
    par1=tools::string_to_double(parts[0]);
    if(parts.size()==2)par1=Units::convert(par1,parts[1],"au");
    par0=Dur-par1;

    der=new AlgebraTrunc("trunc["+tools::str(par0/2.)+","+tools::str(par0/2.+par1)+",1]");
    env=new AlgebraTrunc("trunc["+tools::str(par0/2.)+","+tools::str(par0/2.+par1)+"]");

    tBegin-=par0/2.+par1;
    tEnd  +=par0/2.+par1;
    print0=tBegin;
    print1=tEnd;
}

PulseGauss::PulseGauss(string Enve, double APeak, double T0, double Dur, double Ome, double Phi,
                       double Ellip, double AngleEllip, double PolAng, double AziAng, double DeltaOmega)
    :PulseSingle(Enve,APeak,T0,Dur,Ome,Phi,Ellip,AngleEllip,PolAng,AziAng,DeltaOmega)
{
    par0=2*log(2.)/(Dur*Dur); // FWHM(intensiy)^2= 2log2/par0
    print0-=3*Dur;
    print1+=3*Dur;
    tBegin=-DBL_MAX;
    tEnd  =+DBL_MAX;
}

PulseCosN::PulseCosN(string Enve, double APeak, double T0, double Dur, double Ome, double Phi,
                     double Ellip, double AngleEllip, double PolAng, double AziAng, double DeltaOmega)
    :PulseSingle(Enve,APeak,T0,Dur,Ome,Phi,Ellip,AngleEllip,PolAng,AziAng,DeltaOmega)
{
    cosPow=tools::string_to_int(Enve.substr(Enve.find("cos")+3));
    double nCyc=Dur*Ome/(2*math::pi);
    if(cosPow<8 and Dur*Ome<2*math::pi*8)
        PrintOutput::warning("Pulse shape "+Enve+" can produce artefacts for short pulses (here: "+tools::str(nCyc,2)
                             +" cycles), recommended shape is cos8",1);
    switch(cosPow){
    default: ABORT("admissible cosN, N=2,4,8; is: "+Enve);
    case 2:par0=2*acos(sqrt(0.5))/Dur;break;                   // FWHM(intensity)=2 acos((0.5)^(1/2))/par0
    case 4:par0=2*acos(sqrt(sqrt(0.5)))/Dur;break;             // FWHM(intensity)=2 acos((0.5)^(1/4))/par0
    case 8:par0=2*acos(sqrt(sqrt(sqrt(sqrt(0.5)))))/Dur;break; // FWHM(intensity)=2 acos((0.5)^(1/16))/par0
    }
    par1=math::pi/(2*par0);
    print0-=par1;
    print1+=par1;
    tBegin-=par1;
    tEnd  +=par1;
}

double PulseCosN::fwhm() const{
    switch(cosPow){
    case 8: return (2*acos(sqrt(sqrt(sqrt(sqrt(0.5))))))/par0;
    case 4: return (2*acos(sqrt(sqrt(0.5))))/par0;
    case 2: return (2*acos(sqrt(0.5)))/par0;
    default: ABORT(Str("need to define fwhm for cos")+cosPow);
    }
}
