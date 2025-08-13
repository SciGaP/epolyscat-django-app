// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef DISCRETIZATIONSPECTRALPRODUCT_H
#define DISCRETIZATIONSPECTRALPRODUCT_H

#include <memory>

#include "discretizationSpectral.h"
#include "operatorAbstract.h"

/**
 * Extension of DiscretizationSpectral to allow for situations where the Hilbert space is of
 * pseudo-tensor-product structure and the operator used for creating the spectral projection (the Laplacian)
 * is separable. In the tensor-product case (with obvious generalisations to pseudo-tensor-product),
 * a projection \f$P_{\mathrm{low}}\f$ removing energies higher than cutE can be written as
 * \f[
 *
 * P_{\mathrm{low}} = P_{\mathrm{low},1}\otimes P_{\mathrm{low},2} =
 *      (P_{\mathrm{low},1}\otimes 1)(1\otimes P_{\mathrm{low},2})
 *
 * \f]
 *
 * with \f$P_{\mathrm{low},i}\f$ given by \f$1 - \mathrm{getMapTo()}\circ \mathrm{mapFromParent}\f$ of the factor
 * discretizations (DiscretizationSpectral).
 */
class DiscretizationSpectralProduct : public DiscretizationSpectral{
    class Projector: public OperatorAbstract{
        const DiscretizationSpectralProduct* parent;

    public:
        Projector(const DiscretizationSpectralProduct* Parent);
        void apply(std::complex<double> A, const Coefficients& Vec, std::complex<double> B, Coefficients& Y) const;
    };

    mutable std::unique_ptr<Projector> proj;

public:
    std::vector<std::unique_ptr<DiscretizationSpectral>> factors;

    DiscretizationSpectralProduct(const Discretization *D,
                                  const std::string OpSeparable,
                                  double Emin=-DBL_MAX,
                                  double Emax=DBL_MAX,
                                  bool excludeEnergyRange = false
            );


    const OperatorAbstract* projector() const;

    /**
     * Perform various consistency checks on the discretization (see implementation).
     *
     * This requires the full projector (f. e. set up Hamiltonian == Laplacian)
     */
    void checkFull(const OperatorAbstract* Projector, double cutE) const;
};

#endif // DISCRETIZATIONSPECTRALPRODUCT_H
