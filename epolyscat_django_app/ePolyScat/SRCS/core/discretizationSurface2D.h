#ifndef DISCRETIZATIONSURFACE2D_H
#define DISCRETIZATIONSURFACE2D_H

#include "discretizationSurface.h"

class ReadInput;
class DiscretizationSurface2D : public DiscretizationSurface
{
public:
    /// interprete spherical surface defintion string
    static void sphericalParameters(std::string Def, double &R, int &LMax, int &MMax);
    DiscretizationSurface2D(const Index* Idx, std::string SurfaceDef);
};

#endif // DISCRETIZATIONSURFACE2D_H
