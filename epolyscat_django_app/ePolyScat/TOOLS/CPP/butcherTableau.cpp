// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "butcherTableau.h"

#include "abort.h"

ButcherTableau::ButcherTableau(std::string Name)
    :_name(Name)
{
    if(Name=="RK4"){
        _consistency=4;
        _a.resize(4);
        _a[0]={0.,0.,0.,0.};
        _a[1]={0.5,0.,0.,0.,};
        _a[2]={0.,0.5,0.,0.};
        _a[3]={0.,0.,1.0,0.};
        _b={1./6.,1./3.,1./3.,1./6};
        _c={0.,0.5,0.5,1.};
    }
    else if (Name=="Butcher67")
    {// source: legacy
        _consistency=6;
        _a.resize(7);
        _a[0].assign(7,0.);
        _a[1]={1./3.,      0.,     0.,     0.,   0.,  0.,    0.};
        _a[2]={0.,      2./3.,     0.,     0.,   0.,  0.,    0.};
        _a[3]={1./12.,  1./3.,-1./12.,     0.,   0.,  0.,    0.};
        _a[4]={-1./16., 9./8.,-3./16., -3./8.,   0.,  0.,    0.};
        _a[5]={    0.,  9./8., -3./8., -3./4.,1./2.,  0.,    0.};
        _a[6]={9./44.,-9./11.,63./44.,18./11.,   0.,-16./11.,0.};
        _c={1./3.,     2./3.,  1./3.,  1./2., 1./2.,       1.};
        _b={11./120.,27./40.,27./40.,-4./15.,-4./15.,11./120.};

//        _c[1]=1./3.;
//        _c[2]=2./3.;
//        _c[3]=1./3.;
//        _c[4]=1./2.;
//        _c[5]=1./2.;
//        _c[6]=1.;


//        _b[0]=11./120.;
//        _b[2]=27./40.;
//        _b[3]=27./40.;
//        _b[4]=-4./15.;
//        _b[5]=-4./15.;
//        _b[6]=11./120.;

    }
    else if(Name=="CrankNicolson"){
        _consistency=2;
        _a.resize(2);
        _a[0]={0., 0.};
        _a[1]={0.5,0.5};
        _b={0.5,0.5};
        _c={0.,1.};
    }
    else if(Name=="GaussLegendre2"){
        // from Wikepedia: Gauss-Legendre method
        _consistency=4;
        _a.resize(2);
        _a[0]={0.25,             0.25-sqrt(3.)/6.};
        _a[1]={0.25+sqrt(3.)/6., 0.25            };
        _b={0.5,             0.5};
        _c={0.5-sqrt(3.)/6., 0.5+sqrt(3.)/6.};
    }

    else
        ABORT("no default Butcher tableau with name \""+Name+"\"");
}

