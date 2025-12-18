// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATORNDIM_H
#define OPERATORNDIM_H

#include <vector>
#include <map>
#include <complex>

class OperatorDefinition;
class Algebra;
class Index;

/// \ingroup OperatorData
/// \brief convert an OperatorDefinition to an operator for use with BasisNdim
class OperatorNdim
{
    typedef std::complex<double> (*Kernel)(const std::vector<std::complex<double> > & Q);

    static std::map<std::string,Kernel> special;
    static std::complex<double> identity(const std::vector<std::complex<double> > & Coor){return 1.;} // complex scaling will produce complex potentials
    static void setDefaultKernels(){addPotNdim("Id",identity);}

    class Term {
        friend class OperatorNdim;
        int _ivd,_jvd; // index for value/partial derivative on either side
        std::complex<double> _mult;
        std::vector<Algebra*> alg; // one algebraic expression for each coordinate
        Kernel _potNdim;
    public:
        Term(std::string Single, std::complex<double> Multiplier);
        int ivd() const {return _ivd;} ///< rhs: 0...use value, 1,2,3... use partial derivatives
        int jvd() const {return _jvd;} ///< lhs: 0...use value, 1,2,3... use partial derivatives
        std::complex<double> kernel(std::vector<std::complex<double> > Point) const;
        std::string str() const;
    };
private:
    std::vector<Term>_terms;
public:
    OperatorNdim(const std::string &Def, const std::string & ICoor, const std::string & JCoor);
    unsigned int terms() const {return _terms.size();}
    const Term & term(unsigned int N) const {return _terms[N];}
    static void addPotNdim(std::string Name, Kernel Pot);
};

#endif // OPERATORNDIM_H
