// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#ifndef MAPGAUGE_H
#define MAPGAUGE_H

#include "operatorAbstract.h"
#include <deque>
#include <memory>
#include "tree.h"

class DiscretizationGrid;
class DiscretizationSurface;

/// \ingroup OperatorData
/// \brief map surface to velocity gauge
class MapGauge:public OperatorAbstract
{
    // use deque to ensure data can be pointed to
    mutable double _time;
    std::vector<std::vector<double> > rgrid; // 0,1,2...unit vector, 3...min(length,GaugeRadius)
    std::deque<std::complex<double> > phase;
    std::deque<double> aDotUnit;

    // tree with pointers to the Phases
    class GaugePhase: public Tree<GaugePhase>{
    public:
        std::complex<double>* pPhase;
        double * pADotUnit;
        GaugePhase(const Index* IdxAng,double Phi,double Eta,
                   std::vector<std::vector<double> > &Unit,
                   std::deque<double> &ADotUnit,
                   std::deque<std::complex<double> > & Phase
                   );
    };
    GaugePhase* gPhas;
    std::shared_ptr<const OperatorAbstract> mapToSurface;
    std::shared_ptr<Coefficients> _cSurf,_cGrid;
    DiscretizationGrid* gridAngular;
    double setNorm(){return 1.;}
    void applyPhase(Coefficients* Vec, const GaugePhase *Grid) const;
    void _construct();

    // these are private, apply only with time set
    void apply(std::complex<double> Alfa, const Coefficients &X, std::complex<double> Beta, Coefficients &Y) const ;
    void update(double Time, const Coefficients* CurrentVec=0);
public:
    virtual ~MapGauge();
    MapGauge(const Index* SurfI);
    MapGauge(const DiscretizationSurface *S);
    void axpy(std::complex<double> Alfa, const Coefficients &X, std::complex<double> Beta, Coefficients &Y,double Time){
        if(_time!=Time)update(Time);
        _time=Time;
        apply(Alfa,X,Beta,Y);
    }
};


#endif // MAPGAUGE_H
