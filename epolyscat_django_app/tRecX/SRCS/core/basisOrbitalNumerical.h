// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#ifndef BASISORBITALNUMERICAL_H
#define BASISORBITALNUMERICAL_H

#include "basisOrbital.h"

class Index;
class VectorValuedFunction;

/** \ingroup Basissets */
///@brief Generate BasisOrbital from function
class BasisOrbitalNumerical : public BasisOrbital
{
public:
    class IntParameters
    {
        double _epsAbs,_epsRel;
        size_t _maxLevel;
        std::vector<std::vector<double>> _grid,_weig;
    public:
        /// \brief Integrate parameters for integration
        /// \param NQuad     number of quadrature points
        /// \param AccRel    desired relative accuracy
        /// \param AccAbs    desired absolute accuracy
        IntParameters(std::vector<size_t> NQuad, double AccRel=1.e-12, double AccAbs=1.e-12,size_t MaxLevel=20)
            :_epsAbs(AccAbs),_epsRel(AccRel),_maxLevel(MaxLevel){
            // n-point Gauss-Legendre Quadratures on [0,1]
            for(auto n: NQuad){
                OrthogonalLegendre l;
                _grid.push_back({});
                _weig.push_back({});
                l.quadratureGauss(n,_grid.back(),_weig.back());
                for(auto &g: _grid.back())g=0.5*(g+1.);
                for(auto &w: _weig.back())w*=0.5;
            }
        }
        std::vector<std::vector<double>> grid() const {return _grid;};
        std::vector<std::vector<double>> weig() const {return _weig;};
        size_t maxLevel() const {return _maxLevel;};
        double epsAbs() const {return _epsAbs;}
        double epsRel() const {return _epsRel;}
    };
protected:
    std::string _refName;
    std::string _funcDef;
    std::vector<int> _select;
    Eigen::MatrixXd _overOnBasis;
    static void mapTreeProduct(const Index* Idx, const VectorValuedFunction* Orbs, std::vector<int> Select,std::vector<Coefficients*> Orb, std::vector<unsigned int> Points={});
    static void mapOrbitals(const Index* Idx, const VectorValuedFunction * Orbs, std::vector<int> Select,std::vector<Coefficients*> Orb);
    void convergeAxis(std::string Ax, double LowB, double UpB, const VectorValuedFunction *Orbs, std::vector<int> Select,
                      Eigen::MatrixXcd Morb,
                      std::vector<std::string> Axes, std::vector<std::vector<double>> Grids, std::vector<std::vector<double>> Weigs);
    static void convergeAxis(std::string Axis, const VectorValuedFunction *Orbs, std::vector<int> Select, std::vector<Coefficients *> Orb,
                             std::vector<std::vector<double> > &Grids, std::vector<std::vector<double> > &Weigs);
    static Eigen::MatrixXcd intTreeProduct(const Index* Idx, const VectorValuedFunction * Orbs, std::vector<int> Select,
                                           std::vector<std::vector<double> > Grids, std::vector<std::vector<double> > Weigs);

    static Eigen::MatrixXcd intTreeRecursive(const Index *Idx, const VectorValuedFunction *Orbs, std::vector<int> Select,
                                             const IntParameters &IntPars, std::vector<std::vector<double>> Vol, size_t Level=0, Eigen::MatrixXcd Previous=Eigen::MatrixXcd());
    static void mapTreeProductExact(const Index *Idx, const VectorValuedFunction *Orbs, std::vector<int> Select, std::vector<Coefficients *> Orb);
public:
    static void read(ReadInput & Inp);
    static void  setup(); ///< setup all defined BasisOrbitalNumerical
    BasisOrbitalNumerical(const BasisSetDef &Def);
    unsigned int size() const {return _select.size();}
    unsigned int order() const {return size();}
    const std::vector<int> & select() const {return _select;}

    ///@brief Constructor
    ///
    /// Def is FuncDef:RefIndex:ibeg:iend, :ibeg:iend is optional,
    /// <br>e.g. Phi*Eta*exp(-Rn):main (main reference Index is set up in main_trecx)
    /// <br>MO:complement...MO must have been defined by VectorValuedFunctions::add
    /// <br>             ..."complement" index must have been set previously, e.g. for a hybrid Index
    /// Coordinates must match the Index axes
    BasisOrbitalNumerical(std::string Def);
    void generateOrbitals(const Index* Idx=0);
    void print(std::string Title="") const;
    //    void plot() const; ///< plot, if Plot is defined in ReadInput::main
    std::string str(int Level=0) const;
    bool operator==(const BasisAbstract & Other) const;
    const VectorValuedFunction* orbitalFunctions() const;
    std::string strDefinition() const;
};

#endif // BASISORBITALNUMERICAL_H
