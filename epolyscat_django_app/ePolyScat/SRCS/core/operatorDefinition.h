// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#ifndef OPERATORDEFINITION_H
#define OPERATORDEFINITION_H

#include <string>
#include <map>
#include <vector>

class Index;
class UseMatrix;

/** @defgroup OperatorData Definitions
 *  \ingroup Operators
 *  \brief definition strings, syntax, manipulation, and matrix evaluation
 *  @{
*/

/// define and manipulate standard operator strings (Laplacian, derivatives etc) for a range of coordinate systems
class OperatorDefinition: public std::string
{
    // place constraint string C after the Pos'th "." in after $ConstType
    static std::string addConstraint(std::string Def, std::string ConstType, int Pos, std::string C);
    ///HACK - do not use
    static std::string opPermuted(const std::string &Op, const std::vector<unsigned int> Perm);
    OperatorDefinition axisToCoor(const std::string Hierarchy);///< if top-level is continuity-level, move factor to end
    std::string parentHierarchy(std::string Hierarchy);///< replace certain coordinates by what they were derived from

    static OperatorDefinition remainder(const OperatorDefinition & Def);

    /// return vector of terms, single blank term if empty
    static std::vector<std::string> terms(const std::string Def);
    /// extract parameter para from Term, if Param!="" return product with proper sign
    std::string parameter(std::string Term, std::string Param);
    /// return sign of parameter times OtherSign (+ or -)
    std::string sign(std::string Par, const std::string OtherSign);
    /// remainder after first factor was removed
    std::string remainder(std::string Term, std::string First);
    /// return Term without parameter
    std::string noParameter(std::string Term);
    /// return Term without sign
    std::string noSign(std::string Term);
    /// single factor of the form <string>, empty if not first factor
    std::string firstSingleFactor(const std::string Term) const;
    /// single predefined operator <<string>>, empty if not first factor
    std::string firstPredefined(std::string Term) const;
    /// single predefined operator <<string>>, empty if not first factor
    std::string firstSpecial() const;
    /// single first (...) expression that contains any <...>, <<...>> or [[...]]
    std::string firstBracketed() const;

    /// any first factor, <single>, <<predefined>> or [[special]]
    OperatorDefinition firstFactor() const;

    /// coordinates as given in input
    static std::string inputCoordinates(std::string Hierarchy);
    /// compose factor and Term, include update of constraint if there is any
    std::string compose(std::string Factor, std::string Term);
    /// permute single term
    std::string permute(std::string Term, std::string InCoor, std::string OutCoor);
    std::string permuteFEM(std::string Hierarchy);

    /// return operator where an identity factor that appears in all terms is dropped, attach coordinates:
    /// coordinates=X.Y.Z, <A><1><1>+<1><1><B> ->  (<A><1>+<1><B>):X.Z
    std::string removeId(const std::string Hierarchy) const;

    OperatorDefinition & construct(std::string InputDef, std::string Hierarchy);
    bool _abort;
public:
    static const OperatorDefinition identity;
    ~OperatorDefinition(){}
    OperatorDefinition():_abort(true){}

    /// @brief transparent constructor for operator definition
    ///
    /// legal operator strings are of the form
    ///
    /// OpStr = (+/-Term1 +/- Term2 +/-Term3 ... )
    ///
    /// (no blanks are needed, round brackets optional)\n
    ///
    /// TermX has one of the following forms\n
    ///  -  param<singleFactor>OpStrX
    ///  -  param<<PredefinedOp>>OpStrX
    ///  -  OpStrX
    ///
    /// where:
    ///  +  OpStrX....... a legal operator string (including blank)
    ///  +  param ....... string that can be interpreted by class Algebra
    ///  +  singleDef.... string that can be interpreted by class BasisMat1D
    ///  +  PredefinedOp. name of an operator that is predefined for Coordinates
    ///  +  Coordinates.. string of coordinate axis names, separated by '.' (operator factors will be properly permuted)
    ///
    /// there is limited support for handling constraints of the type $BAND.0.*.2  etc.
    /// transitional signature
    OperatorDefinition(const std::string String, std::string Hierarchy, bool Abort=true);
    OperatorDefinition(const OperatorDefinition Def, std::string Hierarchy, bool Abort=true)
        :OperatorDefinition(Def.str(),Hierarchy,Abort){}

    OperatorDefinition(std::string Def):OperatorDefinition(){assign(Def);} ///< use operator "as is" assuming it will match hierarchy

    /// false if axis name is not to be counted into operator defintion (e.g. FE axis)
    static bool skipAxis(const Index* Idx);

    /// get overalp definition for Idx
    static std::string defaultOverlap(const Index* Idx);

    std::string str() const {return *this;}

    /// return string where certain terms are droped (for deriving operators)
    OperatorDefinition dropTerms(std::string Hierarchy) const;

    /// used for tsurff computations, modify base disc operator strings for semibound regions (duplicate of the above?)
    static void dropTerms(std::vector<std::string> &terms, const Index *IIndex, const Index *JIndex);
    static void dropTerms(std::vector<std::string> &terms, std::string Hierarchy);
    static void dropTerms(std::vector<OperatorDefinition> &terms, std::string Hierarchy);

    /// modify for use in tSurff
    OperatorDefinition & tsurffDropTerms(const Index* IIndex, const Index* JIndex);

    //--------------- re-implemetations from OperatorData ------------------------

    /// keeps only those parts of Definition that start with <I,J> (or <1>), terminates if structure does not match
    OperatorDefinition extractBlock(size_t I, size_t J) const;

    ///expand such that first term is factor<def> or factor[[special]]
    std::vector<OperatorDefinition> singleTerms(const std::string & IHierarchy, const std::string & JHierarchy) const;

    /// return vector of terms extracted from definition
    std::vector<OperatorDefinition> terms() const;

    /// return vector of factors in a single term (abort if not single term)
    std::vector<OperatorDefinition> factors() const;

    /// return OperatorDefinition adjusted for given hierarchy
    OperatorDefinition expandStandard(std::string Hierarchy) const;

    OperatorDefinition first(bool Unsign=false) const; ///< return the first <..> factor of Def (if any)
    //-------------------------------------------------------------------------------

    /// operator string for Hierarchy given in the form X.Y.Z or Phi1.Eta1.Phi2.Eta2
    /// if duplicate coordinates or pairs specXXX - kXXX, first term will be ingnored,
    /// e.g specA.Z.kA -> Z.kA, or Y.B.Y -> B.Y
    static std::string get(std::string Name, std::string Hierarchy);

    std::string parameter(OperatorDefinition & Remainder) const;
    static void setParameters(const std::string Definition);
    /// return name without << >> and trailing parameters, set parameter value
    static std::string extractParameters(const std::string Name);
    /// resolve brackets (without expanding standard operators)
    static OperatorDefinition unBracket(std::string Term, std::string Front="",std::string Back="");

    static void truncationRadius(double Rsmooth, double Rtrunc); ///< too specialize for being here

    /// apply constraints to Mult, return remainder with adjusted contstraint, sanity checks on constraints
    void specialConstrain(UseMatrix & Mult,const Index * IIndex, const Index * JIndex) const;

    /// set those elements of Mult = 0 where this is specified in a definition
    /// (used to avoid computations where matrix elements are known to be =0)
    OperatorDefinition constrain(UseMatrix &Mult, const Index * IIndex, const Index * JIndex) const;

    /// true for non-special operators: containing only <singel> or <<predefined>> factors
    bool isStandard(const Index * IIndex, const Index * JIndex) const;
    /// does not contain any pre-defined factors
    bool isExpanded() const {return str().find("<<")==std::string::npos and str().find(">>")==std::string::npos;};
    /// local - no overlap different (FE) branches (implemented case by case)
    bool isLocal(const Index * IIndex, const Index * JIndex) const;
    /// indicates that operator has the format A (x) 1 + 1 (x) B
    bool isSeparable() const;
    /// basic syntax check (very rudimentary for now)
    void syntax() const;

    static std::string coordinates(const std::string  & Hierarchy);///< remove non-coordinate strings from hierarchy (i.e. remove FE and specXXX...)

    /// standardOperator[identifier][coordinates]...operator string for given coordinates, e.g. standardOperator["Laplacian"][Rn.Phi.Eta]
    static std::map<std::string,std::map<std::string,std::string> >standardOperators;
    static void setup(); ///< set standard operators for a range of coordinate systems
    static void Test();

    /// warn if time-dependent or otherwise non-constant parameters appears in operator
    void checkVariable(std::string Message, std::string Name) const;

};

/** @} */ // end group OperatorData

#endif // OPERATORDEFINITION_H
