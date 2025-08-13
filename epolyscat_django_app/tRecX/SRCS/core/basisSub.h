// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef BASISSUB_H
#define BASISSUB_H

#include "str.h"
#include "qtEigenDense.h"
#include "basisAbstract.h"

/** \ingroup Basissets */
/// \brief subset of the functions from a BasisAbstract
class BasisSub: public BasisAbstract
{
    friend class BasisAbstract;
    friend class BasisThread;
    std::vector<int> _subset;
    const BasisAbstract* _bas;

    /// create subset, default use all, for empty susbset specify Subset={}
    BasisSub(const BasisAbstract *Bas, const std::vector<int> & Subset={});
    BasisSub(std::string StrDefinition);
public:
    static const BasisAbstract* superBas(const BasisAbstract* Sub); ///< first super-basis that is not class BasisSub
    static std::vector<int> subset(const BasisAbstract* Sub); ///< subset wrt. superBas
    std::vector<int> subset() const { return _subset;}
    static Eigen::MatrixXcd subMatrix(const Eigen::MatrixXcd & Mat, std::vector<int> ISub, std::vector<int> JSub); ///< (to toolsEigen)

    void valDer(const UseMatrix &X, UseMatrix &Val, UseMatrix &Der, bool ZeroOutside) const;

    unsigned int size() const override {return _subset.size();}
    unsigned int order() const;
    std::string str(int Level) const override;//{return Str("subset{","")+_subset.size()+"} of "+_bas->str(Level);}
    std::string strDefinition()  const override;
    static std::string strDefinition(const BasisAbstract* Bas, std::vector<int> Subset);

    static Eigen::MatrixXcd map(const BasisAbstract * Ibas, const BasisAbstract * Jbas);
    double physical(int Index) const override;
    bool operator==(const BasisAbstract &other) const override;
    bool isOrthonormal() const override {return _bas->BasisAbstract::isOrthonormal();}

};

#endif // BASISSUB_H
