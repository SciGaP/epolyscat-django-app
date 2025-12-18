#ifndef MAPSURFACE2D_H
#define MAPSURFACE2D_H

#include "operatorTree.h"

class Index;
class Inverse;
class OperatorMap;
class CoorSystem;
class MapSurface2D : public OperatorAbstract
{
    class MapSurfaceGrid: public OperatorTree{
        MapSurfaceGrid(const Index* ISurf, const Index* Idx, std::map<std::string,double> Coors,
                       std::complex<double> Val, std::vector<std::complex<double>>Ders, size_t Pos);
    public:
        MapSurfaceGrid(){};
        MapSurfaceGrid(const Index* ISurf, const Index* Idx, const CoorSystem & CoorSys);
    };
    std::unique_ptr<MapSurfaceGrid> _toGrid;
    std::unique_ptr<OperatorMap>_toSurf;
    std::unique_ptr<Inverse>_inverse;
    Index* _gIndex;

public:
    ~MapSurface2D();
    /// map Idx to a surface defined by SurfDef
    ///
    /// at present, only SurfDef="Spherical[R,Lmax]" and Idx in parabolic coordinates is recognized
    MapSurface2D(const Index* Idx,std::string SurfDef);
    void apply(std::complex<double> A, const Coefficients& Vec, std::complex<double> B, Coefficients& Y) const;
};

#endif // MAPSURFACE2D_H
