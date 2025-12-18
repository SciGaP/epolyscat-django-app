// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#ifndef TSURFF_SOURCE_NEW_H
#define TSURFF_SOURCE_NEW_H

#include <vector>

#include "operatorTree.h"
#include "operatorFloor.h"
#include <memory>

class SurfaceFlux;
class ReadInput;
class DiscretizationSurface;
class Wavefunction;
class DiscretizationTsurffSpectra;
class OperatorAbstract;
class OperatorMap;
class VolkovGrid;
class OperatorTreeVolkov;
class Coefficients;
class CoefficientsPermute;

class TsurffSource{
    std::shared_ptr<SurfaceFlux> flux;

    Coefficients* tempSmoothSpecDisc;
    Coefficients* wfSmoothSpecDisc;
    std::shared_ptr<Coefficients> _currentSourceDamped;
    std::shared_ptr<CoefficientsPermute> _x12Source
    ;
    Coefficients* wfSurface;

    OperatorAbstract* _commutator;
    mutable OperatorAbstract* _volkovPhase;

    double _currentTime,_endIntegration;
    mutable double _turnoffRange; ///< linearly turn off source across this time interval

public:
    TsurffSource(const DiscretizationTsurffSpectra* SmoothSpecDisc,
                 ReadInput& Inp, DiscretizationSurface* surfD, std::shared_ptr<SurfaceFlux> &Flux);
    TsurffSource(const std::vector<std::string> &AxisPath, const Index* Parent, ReadInput& Inp);

    ///\brief generates Index's of Surface and Source according to AxisPath
    ///
    /// source will be from surface defined by all except last AxisPath
    static void generate(const std::vector<std::string> AxisPath, ReadInput& Inp, const Index* Parent, Index* & Surface, Index* & Source);

    /// moves amplitude calculasion to backup
    static bool recompute();

    CoefficientsLocal *UpdateSource(double Time);
    CoefficientsLocal *CurrentSource();

    const DiscretizationTsurffSpectra* smoothSpecDisc; //!< discretization for radial k-value and spherical harmonics
    double SourceBufferBeginTime();
    const Index* idx() const; ///< source Index (i.e. surface converted to spectrum)
    std::vector<double> surfaceRadius(const std::string Axis) const; ///< surface radii (at present only 1)
    static std::string nextRegion(ReadInput & Inp); ///< find next spectral region to compute
    static std::vector<std::string> unboundAxes(ReadInput & Inp);

    double currentTime() const {return _currentTime;}
    double &turnOffTime() const {return _turnoffRange;}
    std::string str() const; ///< energy range and points, turn-off interval
    void print(std::string SourceFile) const; ///< pretty print details of the source

    static void allReads(ReadInput & Inp);
    static double readTurnoff(ReadInput& Inp);
    static bool readSymmetry12(ReadInput & Inp, std::string Hierarchy="");
    static void moveSubdirs(ReadInput& Inp);
    void setTurnOff(ReadInput & Inp, double &EndProp); ///< set time for lineaer source-turnoff
    const OperatorTree* commutator() const;
    const OperatorMap *surfaceToK() const;
    const OperatorAbstract* volkovPhase(double Time) const;

};
#endif
