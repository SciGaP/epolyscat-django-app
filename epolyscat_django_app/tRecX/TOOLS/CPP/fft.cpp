// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "fft.h"
#include "abort.h"

using namespace std;

#ifdef _USE_FFTW_

Fft::Fft(unsigned int Size,bool Forward)
{
    // [FFTW] reserve storage of the internal C type fftw_complex
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Size);
    out =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Size);

    // forward or backward, optimize (FFTW_MEASURE) for repeated applications
    if(Forward)plan = fftw_plan_dft_1d(Size, in, out, FFTW_FORWARD, FFTW_MEASURE);
    else       plan = fftw_plan_dft_1d(Size, in, out, FFTW_BACKWARD,FFTW_MEASURE);
}

vector<complex<double> > & Fft::transform(vector<complex<double> > & Vec){

    for (unsigned int k=0;k<Vec.size();k++)
    {
        // [FFTW] fftw_complex stores complex as two subsequent doubles
        in[k][0]=Vec[k].real();
        in[k][1]=Vec[k].imag();
    }

    // do the FFT
    fftw_execute(plan);

    // include a normalization for conservation of sum_k |F[f](k)|^2=sum_l |f(l)|^2
    double qnrm=1./sqrt(double(Vec.size()));
    for (unsigned int k=0;k<Vec.size();k++)
    {
        Vec[k]=complex<double>(out[k][0],out[k][1])*qnrm;
    }
    return Vec;
}

Fft::~Fft(){
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
}

void Fft::Test(unsigned int Size){

    // set up forward and backward transformations
    Fft forw(Size,true),back(Size,false);

    // put a random vector
    vector<complex<double> >ran;
    for(int k=0;k<Size;k++){
        double rRand=rand()/double(RAND_MAX);
        double iRand=rand()/double(RAND_MAX);
        ran.push_back(complex<double>(rRand,iRand));
    }

    vector<complex<double> >work(ran);
    work=forw.transform(work);
    work=back.transform(work);
    for(int k=0;k<ran.size();k++)work[k]/=ran[k];

    for(int k=0;k<work.size();k++){
        if(abs(work[k]-1.)>1.e-10){
            cout<<"error: "<<work[k]-1.<<",  ratios (sample):\n";
            for(int l=0;l<work.size();l+=work.size()/13)cout<<" "<<work[l];
            cout<<endl;
            ABORT("fft pair does not match");
        }
    }
    cout<<"OK Fft with "<<Size<<endl;

}
#endif
