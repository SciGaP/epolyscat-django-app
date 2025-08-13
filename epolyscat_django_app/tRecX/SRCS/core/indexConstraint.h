// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef INDEX_CONSTRAINT_H
#define INDEX_CONSTRAINT_H

#include <functional>
#include <vector>
#include <map>
#include <string>
#include <memory>

#include "str.h"
class ReadInput;
class Index;
class Coefficients;

/**
 * Impose constraints, like l1 >= 5 => l2 < 5, l2 >= 5 => l1 < 5
 */
class IndexConstraint{

    struct Constraint{
        virtual ~Constraint(){}
        const std::vector<std::string> required_physicals;
        Constraint(std::vector<std::string> RequiredPhysicals): required_physicals(RequiredPhysicals) {}
        virtual bool includes(std::map<std::string, double>& Physicals) =0;
        virtual std::string str() const=0;
    };

    struct RemoveOrbitalsConstraint: public Constraint{
        std::string _axes;
        std::string _orbitalDef;
        std::vector<int> _select;

        RemoveOrbitalsConstraint(std::string Axes,std::string Kind);
        bool includes(std::map<std::string, double>& Physicals) override {return true;}
        std::string str() const override {
            return "On "+_axes+" remove orbitals "+_orbitalDef+
                    (_select.size()?", select {" +tools::str(_select,3,"")+"}":"");
        }
    };

    struct LmMConstraint: public Constraint{
        const int threshold;
        const int excluded_l;
        const std::string l_axis;
        const std::string m_axis;

        LmMConstraint(std::string LAxis, std::string MAxis, int Threshold, int Excluded_l):
            Constraint({LAxis, MAxis}),
            threshold(Threshold),
            excluded_l(Excluded_l),
            l_axis(LAxis),
            m_axis(MAxis){}

        bool includes(std::map<std::string, double>& Physicals) override{
            return (Physicals[l_axis] + Physicals[m_axis]) < threshold
                    or Physicals[l_axis] < excluded_l;
        }
        std::string str() const override {return Str("L-M<","")+threshold+" for L>"+excluded_l+" at "+l_axis+"."+m_axis;}
    };

    struct MZeroConstraint: public Constraint{
        const std::string m1_axis;
        const std::string m2_axis;

        MZeroConstraint(std::string M1Axis, std::string M2Axis):
            Constraint({ M1Axis, M2Axis }),
            m1_axis(M1Axis),
            m2_axis(M2Axis){}

        bool includes(std::map<std::string, double>& Physicals) override{
            return Physicals[m1_axis] + Physicals[m2_axis] == 0.;
        }
        std::string str() const override {return Str("M1+M2=0 at ","")+m1_axis+"."+m2_axis;}
    };

    struct LShapeConstraint: public Constraint{
        int threshold;
        int sumL;
        const std::string l1_axis;
        const std::string l2_axis;

        LShapeConstraint(std::string L1Axis, std::string L2Axis, int Threshold, int SumL):
            Constraint({ L1Axis, L2Axis }),
            threshold(Threshold),
            sumL(SumL),
            l1_axis(L1Axis),
            l2_axis(L2Axis){}

        bool includes(std::map<std::string, double>& Physicals) override{
            return (Physicals[l1_axis] < threshold)
                    or (Physicals[l2_axis] < threshold)
                    or (Physicals[l1_axis] + Physicals[l2_axis] < sumL);
        }
        std::string str() const override {return Str("l1,l2<","")+threshold+" for l1+l2>"+sumL+" at "+l1_axis+"."+l2_axis;}
    };

    std::vector<std::unique_ptr<Constraint>> _constraints;

public:
    virtual ~IndexConstraint(){};
    static IndexConstraint main;
    static void readAll(ReadInput &Inp, std::string &Kind, std::string &Axes, int Line=1);
    void read(ReadInput& Inp);
    const std::vector<std::unique_ptr<Constraint>> & constraints() const {return _constraints;}
    bool includes(std::vector<const Index*> Path, std::vector<unsigned int> Pos) const;
    void apply(Index *Idx, std::vector<const Index*> Path, std::vector<unsigned int> Pos) const;
    std::string str() const;
    static std::vector<Coefficients *> vectors(const Index* Idx); ///< append constraint vectors to Vecs
};
#endif
