#ifndef INDEXORTHO_H
#define INDEXORTHO_H

#include "index.h"
#include "basisMO.h"
#include "basisOrbitalNumerical.h"

/// Index that is orthogonalized to a set of molecular orbitals
class IndexOrtho: public Index
{
    std::unique_ptr<Index> _idx;
    std::unique_ptr<BasisMO> _mo;
    std::unique_ptr<BasisOrbitalNumerical> _oNum; // _mo expressed on a sufficiently dens expansion
public:
    IndexOrtho();
};

#endif // INDEXORTHO_H
