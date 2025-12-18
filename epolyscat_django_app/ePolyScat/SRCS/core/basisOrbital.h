// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISORBITAL_H
#define BASISORBITAL_H

#include "basisAbstract.h"
#include "coefficients.h"
#include "operatorAbstract.h"

class OperatorTree;
class OperatorAbstract;
class OperatorDefinition;
class ReadInput;

/** \ingroup Basissets */

///@brief Several orbitals w.r.t. a reference Index
class BasisOrbital: public BasisAbstract
{
    std::vector<std::string> _plotDef;
    void generate() const;
    bool _isOrthonormal;
protected:
    mutable std::vector<Coefficients> _orb;
    void orthonormalize(bool Warn=true); ///< (pseudo)-Schmidt-orthonormalize, warn if not orthogonal
    BasisOrbital(std::string Name, ReadInput* Inp=0);
    std::vector<double> expectationValues(std::string OpDef) const;
    Eigen::MatrixXd overlap() const;

    class Map: public OperatorAbstract{
        std::vector<Coefficients> _mapV;
        std::unique_ptr<const Index> _idx;
        const BasisOrbital* _orb;
    public:
        Map(const Index* Grid, const Index *OrbJdx, const BasisOrbital* Orb);
        void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const;
    };
    mutable std::shared_ptr<const OperatorAbstract> _mapOrbs;
public:
    void setMap(const Index* Grid,const Index* OrbJdx) const;
    std::shared_ptr<const OperatorAbstract> map() const {return _mapOrbs;}

    virtual ~BasisOrbital(){}
    /// table of names, as Index may not exist yet upon creation of BasisOrbital (should go protected)
    static std::map<std::string,const Index*> referenceIndex;
    static void addIndex(std::string Name, const Index* Idx);
    static void setReference(const Index* Idx); // adds all suitable indices from Idx to reference
    static const Index* getReference(std::string Name){return referenceIndex.count(Name)?referenceIndex[Name]:0;}; // adds all suitable indices from Idx to reference

    virtual void generateOrbitals(const Index* Idx)=0;
    unsigned int size() const {return _orb.size();}
    const Coefficients * orbital(int N) const;
    std::vector<const Coefficients*> orbitals() const;

    ///\brief attach operator tree to Node, where at least IIndex, JIndex or both have BasisOrbital
    ///
    /// i.e. construct <orbBasis| OperDef | referenceBasis> etc.
    static void attachOperator(OperatorTree * Node, std::string Name, const OperatorDefinition & Def, const Index* IIndex, const Index* JIndex,
                               std::complex<double> Multiplier, std::vector<std::complex<double> *> TFac);
    static void clear();
    void plot() const; ///< plot, if plot is defined in input
    std::vector<double> kineticEnergy() const {return expectationValues("0.5<<Laplacian>>");}
    std::vector<double> norms() const  {return expectationValues("<<Overlap>>");}
    bool isIndex() const {return true;}
    std::vector<size_t> inScaledRegion() const;
    bool isOrthonormal() const {return _isOrthonormal;}

};

#endif // BASISORBITAL_H
