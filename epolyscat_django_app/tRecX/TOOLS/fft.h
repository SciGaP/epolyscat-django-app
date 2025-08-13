// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef FFT_H
#define FFT_H

#include <vector>
#include <complex>
#ifdef _USE_FFTW_
#include <fftw3.h>
#endif

class Fft
{
#ifdef _USE_FFTW_
    fftw_complex *in, *out;
    fftw_plan plan;
public:
    ~Fft();
    Fft(unsigned int Size, bool Forward);

    /// normalized DFT (normalizations 1/sqrt(N) included)
    std::vector<std::complex<double> > &transform(std::vector<std::complex<double> > & Vec);

    static void Test(unsigned int Size);
#endif
};

#endif // FFT_H
