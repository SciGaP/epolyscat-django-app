// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "indexConstraint.h"
#include "readInput.h"
#include "index.h"
#include "tools.h"
#include "basisSub.h"
#include "basisOrbital.h"

#include <map>

IndexConstraint IndexConstraint::main;

std::string IndexConstraint::str() const {
    std::string s;
    for(size_t k=0;k<_constraints.size();k++)s+=_constraints[k]->str()+" ";
    return s;
}

void IndexConstraint::readAll(ReadInput &Inp, std::string &Kind, std::string &Axes, int Line){
    Inp.texdocuCategoryAdd("BasisConstraint","kind,axes",
                           R"tex(
                           Easy and flexible constraints on the basis functions are one of the essential
                           features of tRecX's discretization structure.
                           The exact functions used at any node of the hierarchy can be made dependent on the
                           higher levels in the hierarchy.

                           The list of constraints shown below can be easily extended (upon request or creating your own entry in
                           \lcode{class IndexConstraint})
                           )tex",
                           "22,23");
    Inp.read("BasisConstraint", "kind", Kind, "None",
             "Kinds: None...(default), L-M[c1;c2]...restrict to l-m<c1 except l<c2, M=0...total M,"
             "Lshape[c1;c2]...exclude (l1>=c1 && l2>=c1) except l1+l2<c2,"
             "Remove[orbitals:select]...remove selected orbitals", Line)
            .texdocu(R"tex(
                     Options to constrain admissible combinations of basis functions on two (or more) different coordinates.
                     \begin{itemize}
                     \item[L-M] Limit the difference between allowed total $L$ and $L_z$ quantum numbers.
                     This is useful in near-circular polarization, where there is a strong preponderance for
                     simultaneous changes $(L,m)\to (L\pm1,m\pm1)$. For $L<c2$ no constraint will be applied.
                     {\bf Note:} the rule depends on the chirality of the pulse,
                     the constraint will work only for one chirality.
                     Signs are not tracked carefully at the moment find corect chirality by experimentation!
                     \item[M=0] Select a given total M for the two-electron problem
                     (M=0 works for sure, non-zero values may work, try/inquire when needed)
                     \item[Lshape]
                     When a two-electron problem is dominated by single-electron processes,
                     only one $l_i$ out of angular momenta $(l_1,l_2)$ will become large. By symmtry, this designats
                     an L-shaped region in the $l_1,l_2$-plane. The region can be selected here as stated in the short info
                     \item[Remove]
                     Project onto an orthogonal subspace. \''orbitals\'' is a valid definition of an orbitital basis (see BasisNdim),
                     \''select\'' has the format \{1,2,0,7\} and is optional, default is all orbitals
                     \end{itemize}
                     )tex");
    Inp.read("BasisConstraint", "axes", Axes, "BLANK", "Axes used in constraint, ordering is important", Line)
            .texdocu(R"tex(
                     Specifies, which coordinate axes the constraint is to be applied to.
                     )tex");
}

void IndexConstraint::read(ReadInput& Inp){

    Inp.obsolete("IndexConstraint","name","use BasisContraint: kind instead");
    Inp.obsolete("IndexConstraint","axes","use BasisContraint: axes instead");

    _constraints.clear();
    for(int line=1;;line++){
        std::string kind, axes;
        readAll(Inp,kind,axes,line);
        if(kind=="None") break;

        if(kind.find("L-M") != std::string::npos){
            int threshold, excluded_l;
            std::string pars = tools::splitString(kind, '[')[1];
            pars = tools::splitString(pars,']')[0];

            if(pars.find(",") != std::string::npos){
                threshold = tools::string_to_int(tools::splitString(pars, ',')[0]);
                excluded_l = tools::string_to_int(tools::splitString(pars, ',')[1]);
            } else if(pars.find(";") != std::string::npos){
                threshold = tools::string_to_int(tools::splitString(pars, ';')[0]);
                excluded_l = tools::string_to_int(tools::splitString(pars, ';')[1]);
            }else{
                threshold = tools::string_to_int(pars);
                excluded_l = 1.66*threshold;
            }

            if(axes == "BLANK"){
                _constraints.push_back(std::unique_ptr<Constraint>(new LmMConstraint("Eta", "Phi", threshold, excluded_l)));
                _constraints.push_back(std::unique_ptr<Constraint>(new LmMConstraint("Eta1", "Phi1", threshold, excluded_l)));
                _constraints.push_back(std::unique_ptr<Constraint>(new LmMConstraint("Eta2", "Phi2", threshold, excluded_l)));
            }else{
                std::vector<std::string> ax = tools::splitString(axes, '.');
                if(ax.size() != 2) ABORT("Specify exactly two axes, like: Eta.Phi1");
                _constraints.push_back(std::unique_ptr<Constraint>(new LmMConstraint(ax[0], ax[1], threshold, excluded_l)));
            }

        }else if(kind.find("M=0") != std::string::npos){
            if(axes == "BLANK"){
                _constraints.push_back(std::unique_ptr<Constraint>(new MZeroConstraint("Phi1", "Phi2")));
            }else{
                std::vector<std::string> ax = tools::splitString(axes, '.');
                if(ax.size() != 2) ABORT("Specify exactly two axes, like: Phi1.Phi2");
                _constraints.push_back(std::unique_ptr<Constraint>(new MZeroConstraint(ax[0], ax[1])));
            }

        }else if(kind.find("Lshape") != std::string::npos){
            int threshold, sumL;
            std::string pars = tools::splitString(kind, '[')[1];
            pars = tools::splitString(pars,']')[0];
            if(pars.find(";") == std::string::npos){
                threshold = tools::string_to_int(pars);
                sumL = 3*threshold;
            }else{
                threshold = tools::string_to_int(tools::splitString(pars, ';')[0]);
                sumL = tools::string_to_int(tools::splitString(pars, ';')[1]);
            }

            if(axes == "BLANK"){
                _constraints.push_back(std::unique_ptr<Constraint>(new LShapeConstraint("Eta1", "Eta2", threshold, sumL)));
            }else{
                std::vector<std::string> ax = tools::splitString(axes, '.');
                if(ax.size() != 2) ABORT("Specify exactly two axes, like: Eta1.Eta2");
                _constraints.push_back(std::unique_ptr<Constraint>(new LShapeConstraint(ax[0], ax[1], threshold, sumL)));
            }

        }
        else if(kind.find("Remove[")!=std::string::npos){
            _constraints.push_back(std::unique_ptr<Constraint>(new RemoveOrbitalsConstraint(axes,kind)));
        }
        else
            ABORT("Unknown constraint: "+kind);

    }
}

bool IndexConstraint::includes(std::vector<const Index*> Path, std::vector<unsigned int> Pos) const{
    //not very beautiful...
    std::map<std::string, double> physicals;
    for(size_t i=0; i<Pos.size(); i++){
        const BasisAbstract* bas = Path[i]->basis();
        if(bas!=0){
            std::string ax=Path[i]->axisName();

            // set correct physical value where needed
            physicals[ax]=DBL_MAX;
            for(auto& c: _constraints)
                if(std::find(c->required_physicals.begin(),c->required_physicals.end(),ax)!=c->required_physicals.end()){
                    physicals[ax] = bas->physical(Pos[i]);
                    break;
                }
        }
    }

    for(auto& c: _constraints){
        for(auto v: c->required_physicals){
            if(physicals.find(v) == physicals.end()){
                goto Done;
            }
        }
        if(!c->includes(physicals)) return false;
Done:;
    }

    return true;
}



void IndexConstraint::apply(Index *Idx, std::vector<const Index*> Path, std::vector<unsigned int> Pos) const {
    if(_constraints.size()==0)return; // actually no constraints
    Path.push_back(0);
    Pos.push_back(-1);
    std::vector<int> subset;
    for(int k=Idx->childSize()-1;k>=0;k--){
        Path.back()=Idx->child(k);
        Pos.back()=k;
        apply(Idx->child(k),Path,Pos);
        if(Idx->child(k)->sizeCompute()==0 or not includes(Path,Pos))Idx->childErase(k);
        else                                                         subset.insert(subset.begin(),k);
    }
    Path.pop_back();
    Pos.pop_back();
    if(subset.size()!=Idx->childSize()){
        Idx->setBasis(BasisAbstract::factory(BasisSub::strDefinition(Idx->basis(),subset)));
    }
}

std::vector<Coefficients*> IndexConstraint::vectors(const Index* Idx){
    std::vector<Coefficients*> vecs;
    for(size_t k=0;k<main.constraints().size();k++){
        const RemoveOrbitalsConstraint* c;
        if((c=dynamic_cast<const RemoveOrbitalsConstraint*>(main._constraints[k].get()))){
            const BasisAbstract* bas=BasisAbstract::factory(c->_orbitalDef);
            const BasisOrbital* borb=dynamic_cast<const BasisOrbital*>(bas);
            if(not borb)ABORT("Not orbital basis "+c->_orbitalDef);
            // insert into matching place
            for(auto &b: borb->orbitals()){
                vecs.push_back(new Coefficients(Idx,0.));
                for(Coefficients* v=vecs.back();v;v=v->nodeNext())
                    if(v->idx()==b->idx())*v=*b;
                if(vecs.back()->isZero())
                    DEVABORT("could not locate orbital hierarchy "+b->idx()->hierarchy()+" in "+Idx->hierarchy());
            }
        }
    }
    return vecs;
}


IndexConstraint::RemoveOrbitalsConstraint::RemoveOrbitalsConstraint(std::string Axes, std::string Kind)
    :Constraint({}),_axes(Axes)
{
    // expects Remove[definition_of_orbitals:{1,2,4..}]
    if(Kind.find("Remove[")!=0 and Kind.back()!=']')
        ABORT("Incorrect constraint definition, need \"Remove[orbdef:select], got "+Kind);
    _orbitalDef=Kind.substr(std::string("Remove[").length());
    _orbitalDef.pop_back();
    size_t pos=tools::findFirstOutsideBrackets(_orbitalDef,":","([{",")]}");
    if(pos!=std::string::npos){
        std::string ssel=tools::cropString(_orbitalDef.substr(pos+1));
        _select=tools::string_to_intVec(ssel);
        _orbitalDef=_orbitalDef.substr(0,pos);
    }
}
