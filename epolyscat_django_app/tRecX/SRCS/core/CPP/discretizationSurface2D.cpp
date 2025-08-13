#include "../discretizationSurface2D.h"

#include "mapSurface2D.h"
#include "readInput.h"
#include "algebra.h"

DiscretizationSurface2D::DiscretizationSurface2D(const Index* Idx, std::string SurfaceDef)
{
   _mapFromParent.reset(new MapSurface2D(Idx,SurfaceDef));
   idx()=const_cast<Index*>(mapFromParent()->jdx());
}

void DiscretizationSurface2D::sphericalParameters(std::string Def, double &R, int &LMax, int &MMax){
    if(Def.find("Spherical[")==std::string::npos)ABORT("need format Spherical[R,Lmax{,Mmax}], got: "+Def);
    std::vector<std::string> parts=tools::splitString(tools::stringInBetween(Def,"Spherical[","]",true),',');
    Algebra surfR(parts[0]);
    if(not surfR.isAlgebraOfConsts())ABORT("specify surface radius as (algebra of) const, is: "+parts[0]+" in "+Def);
    R=surfR.val(0.).real();
    LMax=tools::string_to_int(parts[1]);
    MMax=parts.size()>2?tools::string_to_int(parts[2]):LMax;
}
