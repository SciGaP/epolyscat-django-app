// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "basisChannel.h"

#include "basisOrbitalNumerical.h"
#include "basisMatMatrix.h"

BasisChannel::BasisChannel(std::string Def,int Size,int First):BasisAbstract(Def){
    if(Def.find("ChannelHF")==0
            or Def.find("ChannelHole")==0
            or Def.find("ChannelDet")==0
            ){
        _orb=new BasisOrbitalNumerical(tools::stringInBetween(Def,"[","]")+":"+tools::str(First)+":"+tools::str(First+Size));
        _rho1.resize(Size);
        if(_orb->size()!=size())
            ABORT(Sstr+"in"+Def+"number of orbitals does not match channels: "+_orb->size()+"!="+size());
        for(size_t i=0;i<size();i++){
            for(size_t j=0;j<size();j++)
                _rho1[i].push_back(Eigen::MatrixXcd::Zero(size(),size()));
            if(Def.find("ChannelHF")==0 or Def.find("ChannelDet")==0){
                _rho1[i][i]=Eigen::MatrixXcd::Identity(size(),size());
                if(Def.find("ChannelHF")==0)_rho1[i][i](i,i)=0.;
            }
            else if(Def.find("ChannelHole")==0)_rho1[i][i](i,i)=-1.;
        }
    }
    else
        ABORT("cannot construct BasisChannel from "+Def);
    BasisMatMatrix::add("<1>:"+str()+"|"+str(),Eigen::MatrixXcd::Identity(size(),size()));
}

BasisChannel::BasisChannel(std::string Def, const BasisOrbital* Orbs):BasisAbstract(Def),_orb(Orbs){
    if(Def.find("ChannelHF")==0
            or Def.find("ChannelHole")==0
            or Def.find("ChannelDet")==0
            ){
        _rho1.resize(_orb->size());
        for(size_t i=0;i<size();i++){
            for(size_t j=0;j<size();j++)
                _rho1[i].push_back(Eigen::MatrixXcd::Zero(size(),size()));
            if(Def.find("ChannelHF")==0 or Def.find("ChannelDet")==0){
                _rho1[i][i]=Eigen::MatrixXcd::Identity(size(),size());
                if(i>0 and Def.find("ChannelHF")==0)_rho1[i][i](i,i)=0.;
            }
            else if(Def.find("ChannelHole")==0)_rho1[i][i](i,i)=-1.;
        }
    }
    else
        ABORT("cannot construct BasisChannel from "+Def);

    BasisMatMatrix::add("<1>:"+str()+"|"+str(),Eigen::MatrixXcd::Identity(size(),size()));
}

