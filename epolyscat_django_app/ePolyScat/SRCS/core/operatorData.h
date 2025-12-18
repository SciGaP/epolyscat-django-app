// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORDATA_H
#define OPERATORDATA_H

#include "basisAbstract.h"
#include "useMatrix.h"
#include "discretization.h"
#include "index.h"
#include "operatorDefinition.h"


/** \ingroup OperatorData */

/// @brief (OBSOLESCENT, use OperatorDefinition instead) parse strings and generate operator (sub-)matrices
///
/// String syntax example: <br>
/// <Op0><Op1>parameter1<Op2><Op3> - <Op4><Op5>parameter2<Op6><Op7> for the tensor product
///
/// \f$
/// \hat{O}^{(0)}\otimes\hat{O}^{(1)}\otimes \left[par_1\times\hat{O}^{(2)}\right]\otimes\hat{O}^{(3)}
/// -\hat{O}^{(4)}\otimes\hat{O}^{(5)}\otimes \left[par_2\times\hat{O}^{(6)}\right]\otimes\hat{O}^{(7)}
///\f$
///
/// "parameter1" and 2 are strings that are interpreted by  <br>
/// class Parameters or derived classes of class FunctionOneArg
///
/// "Op0", "Op1" describe a sum of terms of one-dimensional operators of the form
/// <br>OpK = Term0+Term1-Term2+Term3 etc.
///
/// each term can consist of a parameter and a (product of) factors
/// <br>TermL = par.Fac0[Pars0]*Fac1*Fac2[Pars2]
///
/// par is interpreted by  class Parameters or derived classes of class FunctionOneArg
/// <br>FacM is either an (old style) standard BasisMat::kind: only a single factor allowed
/// <br> or a basisMatFunc: a product of factors are allowed
class OperatorData
{
public:
    /// a list of standard operators for a range of coordinate systems
    static std::map<std::string,std::map<std::string,OperatorDefinition> >standardOperators;
    static OperatorDefinition standardDefs;
    static void setStandard(); ///< set up pre-defined standard operators
    static std::string expandStandard(const std::string Def,const Discretization * Disc);
    static std::string expandStandard(const std::string Def,std::string Hierarchy);
    static std::string expandStandard(const std::string Def, const std::string &IHierarchy, const std::string &JHierarchy);

    static bool isStandard(const std::string & Def); ///< true if standard operator defintion string
    static bool isMultip(const std::string & Def); ///< true if multiplication operator (to be applied by transformation to quadratur grid)
    static bool isZero(const std::string & Def, const Index* IIndex, const Index* JIndex); ///< true if operator is identical to zero
    //    static std::string opIJ(const std::string & Def, int I, int J); ///< extract matrix entries
    static void modifyDefinition(std::string &Def, const Index *IIndex, const Index *JIndex); ///< supplement definition
    static void operatorData(std::string Def, const Index * const IFloor, const Index * const JFloor, std::vector<UseMatrix> & Mats); ///< HACK from Discretization

    /// keeps only those parts of Definition that start with <I,J> (or <1>), terminates if structure does not match
    static std::string extractBlock(std::string Definition, unsigned int I,unsigned int J);

    /// split into factors and remainders: <factor_i><remainder_i>
    static void factorize(std::string Def,Discretization*IDisc, Discretization*JDisc,Index*IIndex,Index*JIndex,
                          std::vector<std::string>&Factors, std::vector<std::string>&Remainders);
    static std::vector<std::string> terms(std::string Def); ///< return vector of terms extracted from definition
    static std::vector<std::string> singleTerms(std::string Def, const std::string &IHierarchy, const std::string &JHierarchy); ///< expand such that first term is factor<def> or factor[[special]]
    static void dependence(const std::string & Term, std::vector<std::string> & Dep); ///< extract dependencies from a single term
    static std::string parameter(const std::string & Term); ///< extract the parameter from a single term
    static std::string sign(const std::string & Term); ///< return sign (+, -, or blank) of term
    static std::string first(const std::string & Def,bool unsign=false); ///< return the first <..> factor of Def (if any)
    static std::string remainder(const std::string & Def); ///< return Def without the first <..> factor (if any)
    static void constrain(UseMatrix & Mult,std::string Constraint,unsigned int Level); ///< return Def without the first <..> factor (if any)

    enum substitutionSpec {
        index,
        surfaceOnFloor
    };
    //! \brief find indices to substitute in the definition
    //! save position and length of occurence of name of left and right axes in Def in leftSub and rightSub;
    //! save if information lies on the same level or the floor in spec.
    static void substitutionIndices(const std::string &Def, substitutionSpec& spec, const std::string &leftName, const std::string &rightName,
                                    std::vector<size_t> &leftSub, std::vector<size_t> &rightSub);
};

#endif // OPERATORDATA_H
