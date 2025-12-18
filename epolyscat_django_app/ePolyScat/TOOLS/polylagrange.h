// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef POLYLAGRANGE_H
#define POLYLAGRANGE_H

#include <vector>

/// Lagrange polynomials
template<class Arg>
class PolyLagrange
{
    std::vector<Arg> _valInit;
    std::vector<Arg> _mesh;
public:
    PolyLagrange(){}
    PolyLagrange(std::vector<Arg> Mesh)
        :_mesh(Mesh)
    {
        _valInit.assign(Mesh.size(),1.);
        for (int i=0;i<_valInit.size();i++)
            for(int j=0;j<_mesh.size();j++)
                if(i!=j)_valInit[i]=_valInit[i]*(1./(_mesh[i]-_mesh[j]));
    }

    const std::vector<Arg>& mesh(){return _mesh;}

    /// values of all polynomials at Q
    std::vector<Arg> values(Arg Q) const
    {
        std::vector<Arg> v(_valInit);
        for (int i=0;i<_valInit.size();i++)
            for(int j=0;j<_mesh.size();j++)
                if(i!=j)v[i]*=(Q-_mesh[j]);
        return v;
    }

    std::vector<Arg> derivatives(Arg Q) const
    {
        std::vector<Arg> v(_valInit);
        std::vector<Arg> d(v.size(),0.);
        for (int i=0;i<_valInit.size();i++)
            for(int j=0;j<_mesh.size();j++)
                if(i!=j){
                    d[i]=d[i]*(Q-_mesh[j])+v[i];
                    v[i]=v[i]*(Q-_mesh[j]);
                }
        return d;
    }

};

#endif // POLYLAGRANGE_H
