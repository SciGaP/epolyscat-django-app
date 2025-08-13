// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "basisMO.h"

#include "readInput.h"
#include "printOutput.h"
#ifdef _USE_HACC_
#include "quantumChemicalInput.h"
#include "mo.h"
#endif
#include "basisOrbitalNumerical.h"
#include "qtEigenDense.h"


void BasisMO::read(ReadInput & Inp){
    if(not Inp.found("Chemical","data"))return;

    std::string Name;
    Inp.obsolete("MolecularOrbitals","name","access from Axis definition, e.g.  as ChannelHF[Chemical{10,9,8,7}:single] ");
    VectorValuedFunction::add(Name,std::shared_ptr<VectorValuedFunction>(new BasisMO(Inp)));
}

void BasisMO::print(std::string RefIdx) const{
    //    BasisOrbitalNumerical orbs(_name+":"+RefIdx+":"+tools::str(_select,",",2));
    //    orbs.print("MolecularOrbitals");
    PrintOutput::title("Molecular orbitals: "+_name+":"+RefIdx);
}

BasisMO::BasisMO(ReadInput & Inp):BasisAbstract("MO")
{
#ifdef _USE_HACC_
    QuantumChemicalInput::read(Inp);
    Inp.obsolete("MolecularOrbitals","select","specify selection in basis for Axis, e.g. ChannelHF[Chemical{10,9,8,7}:single]");
    _vinayMO.reset(new mo(*QuantumChemicalInput::neutrals,0));
    _source=QuantumChemicalInput::neutrals->source();
    VectorValuedFunction::add(_name,std::shared_ptr<VectorValuedFunction>(this));
#else
    ABORT("for haCC, compile with -D_USE_HACC_");
#endif
}


BasisMO::BasisMO(const QuantumChemicalInput * System):BasisAbstract("Chemical"){

#ifdef _USE_HACC_
    // not trying to fix old mo(...)
    _vinayMO.reset(new mo(*const_cast<QuantumChemicalInput*>(System),0));
    _source=System->source();
#else
    ABORT("for haCC, compile with -D_USE_HACC_");
#endif
}

void BasisMO::add(const QuantumChemicalInput * System){
#ifdef _USE_HACC_
    std::string name=System->source();
    if(name.rfind("/")==name.length()-1)name=name.substr(0,name.length()-1);
    name=name.substr(name.rfind("/"));
    VectorValuedFunction::add(name,std::shared_ptr<VectorValuedFunction>(new BasisMO(System)));
#else
    ABORT("for haCC, compile with -D_USE_HACC_");
#endif
}

unsigned int BasisMO::size() const {
    if(_select.size()==0){
        Eigen::MatrixXd val;
        Eigen::Matrix<double, Eigen::Dynamic, 3> arg=Eigen::MatrixXd::Zero(1,3);
#ifdef _USE_HACC_
        const_cast<mo*>(vinayMO())->values(val,arg);
        for(int k=0;k<val.rows();k++)
            const_cast<BasisMO*>(this)->_select.push_back(k);
#endif

    }
    return _select.size();
}

std::string BasisMO::str(int Level) const {
    //    return Str("DVR(","")+_opol->name()+") ["+_lowBound+","+_upBound+"] ("+_shiftX+","+_scale+") "+_size+"["+_dvrX.size()+"]";

    Str strSel=_select.size()<10 ? Str(" {","")+_select+"} " : Str(" {")+_select.front()+"..."+_select.back()+"} ";
    return Str(VectorValuedFunction::name(),"")+strSel+"["+size()+"]"+info();
}
std::string BasisMO::info(int Number) const {
    Str strSel;
    if (Number>-1 or _select.size()==1)
        strSel=Str("{","")+_select[Number]+"}";
    else
        strSel=_select.size()<10 ? Str(" {","")+_select+"} " : Str(" {")+_select.front()+"..."+_select.back()+"} ["+size()+"]";
    return Str(VectorValuedFunction::name(),"")+strSel+" from "+_source;
}

std::vector<std::complex<double> > BasisMO::operator()(std::vector<double> X) const{
    Eigen::MatrixXd val;
    Eigen::Matrix<double, Eigen::Dynamic, 3> arg=Eigen::Map<Eigen::MatrixXd>(X.data(),1,3);
#ifdef _USE_HACC_
    _vinayMO->values(val,arg);
#else
    ABORT("for haCC, compile with -D_USE_HACC_");
#endif
    if(_select.size()==0)
        for(int k=0;k<val.rows();k++)
            const_cast<BasisMO*>(this)->_select.push_back(k);

    std::vector<std::complex<double> > res(_select.size());
    for(size_t k=0;k<_select.size();k++)res[k]=val(_select[k],0);
    return res;
}

