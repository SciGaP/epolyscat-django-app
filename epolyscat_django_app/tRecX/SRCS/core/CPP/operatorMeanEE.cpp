// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorMeanEE.h"

/* Notes
 * 0) The matrix Mij is long along j-index - it will be multiplied with the vector of all coefficients C
 * 1) Since at each floor we integrate over the whole region,
 * there will be
 * that contain nodes: _nodeVal[0][j], _nodeVal[1][j] and
 * _nodeVal[2][j] for border and inner regions.
 * 2) _nBeg will be added to _dvrPoints and _dvrWeights
 * of each floor when the dvr property bi(a_m)=delta(m,i)n_i is used
 * 3) the last floor region is infinite and has its own dvrWeights
 * 5) E=-2.8616 ground state l=0
 * 6) for testing use Rn, 8, 0., 10., polynomial,4
 * 6.5) Use nodeRight and not nodeNext
 * 7) Works only for indexing Rn.Phi.Eta
 * 8) To compute polarizability get <Z> in eigenSolverNonLin at E=0.01
 * 9) To print norm of the Wf use timePropagatorOutput::print() and ::write() and just Sout
 * 10) 1.5 eV = 0.057 (frequency), 5ev = 247.97nm length Hochstuhl 27.55nm
 * 11) Field/Frequency = etas
 * 12) cut energy o speed up Operator: projection='(1/2<<Laplacian>>-1.5<<Coulomb>>):Rn.Phi.Eta'  and cutEnergy 30
 * 13) Hydrogen polarizability -4.5, Helium HF 1.322
*/


bool OperatorMeanEE::_iterations = false;
bool OperatorMeanEE::_noInteraction = false;
bool OperatorMeanEE::disableUpdate=false;
std::complex<double> OperatorMeanEE::ME=0;


void OperatorMeanEE::test(){

}
    
OperatorMeanEE::OperatorMeanEE(const std::vector<int> &Info, const std::vector<std::complex<double> > &Buf)
    :OperatorFloorNonLin("MeanEE")
{
    unpackBasic(Info,Buf);
}

OperatorMeanEE::OperatorMeanEE(std::string Pot, const Index* IIndex, const Index* JIndex, std::complex<double> Multiplier)
    : OperatorFloorNonLin("MeanEE")
{
    dat=0;
    oNorm=1.;
    if(IIndex->basis()->isAbsorptive() // zero norm in absorptive range or when non-local
            or IIndex->basis()->integrable()->lowBound()!=JIndex->basis()->integrable()->lowBound()
            or IIndex->basis()->integrable()->upBound()!=JIndex->basis()->integrable()->upBound()
            )
        oNorm=0.;
    _lsize=JIndex->root()->firstFloor()->parent()->childSize();
    _b=dynamic_cast<const BasisDVR*>(JIndex->basis());
    _R.assign(2*_lsize, std::vector<std::vector<double>>(JIndex->basis()->integrable()->size(), std::vector<double>()));
    _W.assign(_b->size(), 0.);

    _etaJ=JIndex->nSibling(); // current eta wrt J
    _etaI=IIndex->nSibling(); // current eta wrt I




    std::vector<double> dvrPoints, dvrWeights;
    _b->dvrRule(dvrPoints, dvrWeights);



    for(int i=0;i<_b->size();i++){    //needed for speed up for the mean field matrix element *2./(2.*_etaI+1.)
        _ovr.push_back(dvrWeights[i+_b->nBeg()]*std::norm(_b->valNodes()[i]));
    }


    _GauntMatrix.assign(_lsize,std::vector<std::vector<double>>(_lsize,std::vector<double>(2*_lsize,0.)));


    Gaunts GG;
    for(int i=0;i<_lsize;i++){
        for(int j=0;j<_lsize;j++){
            int l3min;
            std::vector<double> VecG=GG.vals(0,0,i,j,l3min);
            for(int k=l3min;k<i+j+1;k++){
                if((i+j+k)%2==1){
                    _GauntMatrix[i][j][k]=0;
                }else{
                    _GauntMatrix[i][j][k]=VecG[(k-l3min)/2];
                }
            }
        }
    }
    computeRadialMatrix(IIndex, JIndex);
}

void OperatorMeanEE::computeRadialMatrix(const Index *IIndex, const Index *JIndex){  //Works only for Rn.Phi.Eta
    for(const Index* jdxParent=JIndex->root()->child(0);jdxParent!=0;jdxParent=jdxParent->nodeRight()){
        const Index* jdx=jdxParent->firstFloor();
        const BasisDVR* bi=dynamic_cast<const BasisDVR*>(IIndex->basis());
        const BasisDVR* bj=dynamic_cast<const BasisDVR*>(jdx->basis());
        std::vector<double> dvrPointsI, dvrPointsJ, dvrWeightsI, dvrWeightsJ;
        bi->dvrRule(dvrPointsI, dvrWeightsI);
        bj->dvrRule(dvrPointsJ, dvrWeightsJ);
        MultipolePotential potEE(2*_lsize, "CoulombEE", IIndex->basis()->integrable(), jdx->basis()->integrable());
        for(int i=0; i<JIndex->basis()->size(); i++){
            for(int j=0; j<jdx->basis()->size(); j++){
                for(int l=0;l<2*_lsize;l++){
                    _R[l][i].push_back( (potEE.vals(l)(i,j)).real()*bj->valNodes()[j]*bi->valNodes()[i] );   // WITH Nj and Ni already!!
                }
            }
        }
    }
}

void OperatorMeanEE::updateNonLin(double Time, const Coefficients* C){

        _W.assign(_b->size(),0.);
        getMeanField(C);
}

void OperatorMeanEE::getMeanField(const Coefficients* C){

    if(true){
        for(int l=0;l<2*_lsize;l++){
            int offset=0;
            for(const Index* jdxParent=C->idx()->root()->child(0);jdxParent!=0;jdxParent=jdxParent->nodeRight()){ //wide horizontal loop for radial components
                const Index* jdx=jdxParent->child(0);
                const Coefficients* Crad=const_cast<Coefficients*>(C)->retrieve(jdx);
                for(int q=0;q<jdx->firstFloor()->childSize();q++){
                    std::complex<double> sumla=0;
                    for(int la=0;la<_lsize;la++){
                        std::complex<double> sumlb=0;
                        for(int lb=0;lb<_lsize;lb++){
                            if((la+lb+l)%2==0)sumlb+=Crad->child(lb)->data()[q]*_GauntMatrix[la][lb][l];
                        }
                        sumla+=std::conj(Crad->child(la)->data()[q])*sumlb;
                    }
                    std::complex<double> GRho=_GauntMatrix[_etaI][_etaJ][l]*sumla;
                    for(int a=0;a<_R[l].size();a++){
                        _W[a]+=_R[l][a][q+offset]*GRho;
                    }
                }
                offset+=jdx->firstFloor()->childSize();
            }
        }
    }else{
        for(int l=1;l<2;l++){
            int offset=0;
            for(const Index* jdxParent=C->idx()->root()->child(0);jdxParent!=0;jdxParent=jdxParent->nodeRight()){ //wide horizontal loop for radial components
                const Index* jdx=jdxParent->child(0);
                const Coefficients* Crad=const_cast<Coefficients*>(C)->retrieve(jdx);
                for(int q=0;q<jdx->firstFloor()->childSize();q++){
                    std::complex<double> sumla=0;
                    for(int la=0;la<_lsize;la++){
                        std::complex<double> sumlb=0;
                        for(int lb=0;lb<_lsize;lb++){
                            if((la+lb+l)%2==0)sumlb+=Crad->child(lb)->data()[q]*_GauntMatrix[la][lb][l];
                        }
                        sumla+=std::conj(Crad->child(la)->data()[q])*sumlb;
                    }
                    std::complex<double> GRho=_GauntMatrix[_etaI][_etaJ][l]*sumla;
                    for(int a=0;a<_R[l].size();a++){
                        _W[a]+=_R[l][a][q+offset]*GRho;
                    }
                }
                offset+=jdx->firstFloor()->childSize();
            }
        }
    }


    for(int i=0;i<_R[0].size();i++){
    }
}

std::complex<double> OperatorMeanEE::getC(const Coefficients* C, int l, int j){  //not needed anymore, can be deleted
    std::vector<std::complex<double>> Ceta;
    const Index* idx=C->idx();
    for(const Index* Parent=idx->firstFloor()->parent();Parent!=0;Parent=Parent->nodeRight()){
        const Index* Floor=Parent->child(l);
        for(int i=0;i<Floor->childSize();i++){
            Ceta.push_back(const_cast<Coefficients*>(C)->retrieve(Floor)->orderedData()[i]);
        }
    }
    return Ceta[j];
}

void OperatorMeanEE::axpy(const std::complex<double> & Alfa, const std::complex<double>* X, unsigned int SizX,
                          const std::complex<double> & Beta, std::complex<double>* Y, unsigned int SizY) const{
    if(_noInteraction){                 //not necessary, used for a hack to compute the non-modified matrix element <psi|H0|psi>
        for(int i=0;i<SizX;i++){
            Y[i]=Beta*Y[i];
        }
    }
    else{
        std::vector<double> dvrPoints, dvrWeights;
        _b->dvrRule(dvrPoints, dvrWeights);
        for(int i=0;i<SizX;i++){
            Y[i]=Beta*Y[i]+Alfa*(_W[i])*X[i];  //+_ovr[i]*ME
        }
    }
}

void OperatorMeanEE::pack(std::vector<int> &Info, std::vector<std::complex<double> > &Buf) const
{
    Buf.insert(Buf.end(),dat->begin(),dat->end());
    packBasic(Info,Buf);
}

