// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef QUANTUMCHEMICALDATA_H
#define QUANTUMCHEMICALDATA_H

#include <vector>
#include "qtEigenDense.h"
#include "my4darray_shared.h"
#include "sod.h"

struct columbus_input_file{
  std::string path;
  double det_tol,sce_cutoff,multipole_cutoff;
  int numI,multF;
  columbus_input_file():det_tol(1.e-12),sce_cutoff(1.e-10),multipole_cutoff(1.e-12),numI(0),multF(1){}
};

class mo;

/// \ingroup ChemStruc
/// \brief data shared by new QuantumChemicalInput and old columbus_data (TEMPORARY)
class QuantumChemicalData
{
public:
    QuantumChemicalData():charge(0),no_atoms(1),no_electrons(0),det_cutoff(0),det_threshold(0.){}

    columbus_input_file cif;

    std::vector<double > charge;
    std::vector<Eigen::RowVector3d > cord;
    int no_atoms;
    int no_electrons;
    unsigned int det_cutoff;
    double det_threshold;

    Eigen::MatrixXd Overlap, KineticEnergy, Potential;
    My4dArray_shared* vee;
    std::string source(){return cif.path;}
};

#endif // QUANTUMCHEMICALDATA_H
