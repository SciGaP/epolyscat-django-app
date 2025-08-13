// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
// End of license
#include "parallel.h"
#include "operatorFloorXC.h"
#ifdef _USE_HACC_
#include "mo.h"
#endif
#include <memory>
#include "multipoleOrbital.h"
#include "multipolePotential.h"
#include "index.h"
#include "basisExpIm.h"
#include "basisOrbitalNumerical.h"
#include "operatorTree.h"
#include "gaunts.h"
#include "densityMatrix1.h"
#include "basisChannel.h"
#include "basisMO.h"
#include "patchRadial.h"
#include "algebra.h"
#include "basisDvr.h"
#include "string.h"
#include "eigenTools.h"
#include <valarray>
#include <bits/stdc++.h>

OperatorFloorXC::OperatorFloorXC(std::string Pot, const Index *AIndex, const Index *BIndex, std::complex<double> Multiplier)
{
    dat = 0;
    oNorm = 1.;
    if (AIndex->basis()->isAbsorptive() or BIndex->basis()->isAbsorptive())
        oNorm = 0.;
}

static int lMax(const Index *Idx)
{
    // HACK somewhat unsafe
    const Index *idx = Idx;
    while (idx != 0 and idx->axisName().substr(0, 3) != "Eta")
        idx = idx->nodeNext();
    if (idx == 0)
        DEVABORT("Index does not contain Eta-axis, hierarchy=" + Idx->hierarchy());
    int l = 0;
    for (size_t k = 0; k < idx->childSize(); k++)
        l = std::max(l, int(idx->basis()->physical(k)));
    return l;
}
static std::vector<const BasisIntegrable *> radialBases(const Index *Idx)
{

    std::vector<const BasisIntegrable *> res;
    // find first Rn floor
    const Index *idx = Idx->firstFloor();
    while (idx != 0 and not(idx->axisName() == "Rn" and idx->basis()->integrable()))
        idx = idx->nodeNext();
    if (idx == 0)
        DEVABORT("no Rn basis found in Idx hierarchy " + Idx->root()->hierarchy());

    for (; idx != 0; idx = idx->nodeRight(Idx))
    {
        if (idx->axisName() != "Rn")
            ABORT("need Rn radial axis, got: " + idx->strNode());
        size_t k;
        for (k = 0; k < res.size() and not(*res[k] == *idx->basis()->integrable()); k++)
            ;
        if (k == res.size())
            res.push_back(idx->basis()->integrable());
    }
    return res;
}
TIMER(XC, )
void OperatorFloorXC::postProcess(OperatorTree *Op, std::string Pot, const std::vector<std::vector<Eigen::MatrixXcd>> &RhoIJ,
                                  const BasisOrbitalNumerical *IOrb)
{

    DensityMatrix1 rho(RhoIJ); // initialize the density matrix
    // for(int i=0;i<RhoIJ.size();i++)
    //         for(int j=0;j<RhoIJ[0].size();j++)
    //         std::cout<<"RhoIJ("<<i<<","<<j<<") "<<std::endl<<RhoIJ[i][j]<<std::endl;

    std::vector<const BasisIntegrable *> radBasS(radialBases(IOrb->orbital(0)->idx()));
    std::vector<const BasisIntegrable *> radBasA(radialBases(Op->iIndex));

    int lpMax = lMax(IOrb->orbital(0)->idx()) + lMax(Op->iIndex);
    Eigen::VectorXcd mpPotAB = Eigen::VectorXcd::Zero(lpMax + 1);

    START(XC);
    //===========================================================================================================================================
    for (PatchRadial<OperatorFloorXC> patch(Op, rho, "full"); patch.leafs.size() > 0; patch = PatchRadial<OperatorFloorXC>(Op, rho, "full"))
    {
        // multipole potential
        std::cout << "patch.aBas()->lowBound(): " << patch.aBas()->lowBound() << " patch.bBas()->lowBound(): " << patch.bBas()->lowBound() << std::endl;

        MultipolePotential mpPot(lpMax, Pot, patch.aBas(), patch.bBas());
        rho.patch(*patch.aBas(), *patch.bBas(), *IOrb, *IOrb); // get spherical representation of density matrix for present patch
        Gaunts gaunts;
        const BasisDVR *aDvr = dynamic_cast<const BasisDVR *>(patch.aBas());
        const BasisDVR *bDvr = dynamic_cast<const BasisDVR *>(patch.bBas());
        std::vector<double> aVal(aDvr->valNodes());
        std::vector<double> bVal(bDvr->valNodes());

        // std::cout << "MPIwrapper::Rank(): " << MPIwrapper::Rank() << " MPIwrapper::Size(): " << MPIwrapper::Size() << std::endl;
        //============================================================
        int cntMll = 0;
        for (size_t a = 0; a < patch.aBas()->size(); a++)
        {
            for (size_t b = 0; b < patch.bBas()->size(); b++)
            {
                if (MPIwrapper::Rank() == (cntMll++) % MPIwrapper::Size())
                {
                    rho.setPoint(a, b);
                    for (unsigned int lp = 0; lp <= lpMax; lp++)
                        mpPotAB(lp) = mpPot.vals(lp)(a, b) * aVal[a] * bVal[b];

                    for (PatchRadial<OperatorFloorXC>::LeafM &leafM : patch.leafs)
                    {
                        // std::cout << " leafM.ma: " << leafM.ma << " leafM.mb: " << leafM.mb << std::endl;
                        int mpMax = std::max(leafM.ma - rho.jmMin(), leafM.mb - rho.imMin());
                        int mpMin = std::min(leafM.ma - rho.jmMax(), leafM.mb - rho.imMax());
                        // for (int mi = rho.imMin(); mi <= rho.imMax(); mi++)
                        //     for (int mj = rho.jmMin(); mj <= rho.jmMax(); mj++)
                        for (int mp = mpMin; mp < mpMax + 1; mp++)
                        {
                            // if(leafM.mb-mi!=leafM.ma-mj)continue;
                            int mi = leafM.mb - mp;
                            int mj = leafM.ma - mp;
                            DensityMatrix1::M *rhoM = rho(mi, mj);
                            // std::cout << " mi: " << mi << " rhoM->ilMin: " << rhoM->ilMin()<< " rhoM->ilMax: " << rhoM->ilMax() << std::endl;
                            // std::cout << " mj: " << mj << " rhoM->jlMin: " << rhoM->jlMin()<< " rhoM->jlMax: " << rhoM->jlMax() << std::endl;
                            for (PatchRadial<OperatorFloorXC>::LeafL &leafL : leafM.leafs)
                            {

                                for (int li : rhoM->listLi())
                                {
                                    int ilpMin;
                                    std::vector<double> &gI = gaunts.vals(-mi, leafM.mb, li, leafL.lb, ilpMin);
                                    
                                    // Eigen::MatrixXd GI=Eigen::Map<Eigen::MatrixXd>(aVal.data(), aVal.size(), 1);
                                    // std::cout<<GI<<std::endl;
                                    DensityMatrix1::L *rhoML = (*rhoM)(li);
                                    std::map<int, std::complex<double>> gauntPotI;
                                    // gauntPotI.clear();
                                    for (int gpI = 0; gpI < gI.size(); gpI++)
                                        gauntPotI[ilpMin + 2 * gpI] = gI[gpI] * mpPotAB[ilpMin + 2 * gpI];
                                    //  std::cout << " mi: " << mi << " leafM.mb: " << leafM.mb<< " li: " << li <<" leafL.lb: "<<leafL.lb<<" lpMin: "<<lpMin<<" gI.size(): "<< gI.size()<< std::endl;
                                    for (int lj : rhoM->listLj())
                                    {
                                        int jlpMin;
                                        std::vector<double> &gJ = gaunts.vals(mj, -leafM.ma, lj, leafL.la, jlpMin);
                                        
                                        // std::map<int, complex<double>> gauntPotJ;
                                        // for (int gpJ = 0; gpJ < gJ.size(); gpJ++)
                                        //     gauntPotJ[jlpMin + 2 * gpJ] = gJ[gpJ] ;
                                        int lpmin = std::max(ilpMin, jlpMin);
                                        int lpmax = std::min(ilpMin + 2 * (gI.size() - 1), jlpMin + 2 * (gJ.size() - 1));
                                        if (abs(ilpMin - jlpMin) % 2 != 0  or lpmax<lpmin)
                                            continue;
                                        
                                        
                                        // std::cout<<GJ.segment(int((lpmin-jlpMin) / 2),int((lpmax - lpmin) / 2))<<std::endl<<std::endl;
                                        
                                        
                                        // std::inner_product(gJ[slice((lpmin - jlpMin) / 2,(lpmax - lpmin) / 2)].beg)
                                        // complex<double> dotIJ=(Eigen::Map<Eigen::MatrixXd>(gJ.data(), 1, gJ.size()).block(0,(lpmin - jlpMin) / 2,0,(lpmax - lpmin) / 2)*Eigen::Map<Eigen::VectorXd>(gI.data(), gI.size(), 1).block(0,(lpmin - ilpMin) / 2,0,(lpmax - lpmin) / 2)).sum();
                                        
                                        
                                        
                                        
                                        DensityMatrix1::LL *rhoMLL = (*rhoML)(lj);
                                        
                                        std::complex<double> vLgauntIJ = 0;
                                        for (int lp = lpmin, Jp = (lpmin-jlpMin) / 2; lp <= lpmax; lp += 2, Jp++)
                                        {
                                            // if (gauntPotI.count(lp) == 0 or gauntPotJ.count(lp) == 0)
                                            //     continue;
                                            vLgauntIJ += gauntPotI[lp] * gJ[Jp]; // gJ[(lp - jlpMin) / 2];
                                            // std::cout << "gauntPot[" << lp << "] " << gauntPot[lp] <<" vLgauntI: "<<vLgauntIJ<< std::endl;
                                        
                                        }
                                        for (PatchRadial<OperatorFloorXC>::LeafC &leafC : leafL.leafs)
                                        {
                                            (*leafC.mat)(a, b) += vLgauntIJ * rhoMLL->densC[leafC.cc];
                                        }
                                        // std::cout << std::endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //================================================================================================================================================
    STOP(XC);
    OperatorFloorXC *floor;
    for (OperatorTree *op = Op->firstLeaf(); op != 0; op = op->nextLeaf())
    {
        if (0 != (floor = dynamic_cast<OperatorFloorXC *>(const_cast<OperatorFloor *>(op->floor()))) and floor->dat == 0)
        {

            floor->finalize(op); // this puts data into the floor
            if (floor->dat == 0)
                OperatorFloor::replace(op->floor()); // set OperatorDUM where no data
            // Parallel::moveFloor(op->floor(),Parallel::floorHost(op->idx(),op->jdx()));
        }
    }
    if (dynamic_cast<const BasisOrbitalNumerical *>(IOrb))
        test(Op, *IOrb);
}

void OperatorFloorXC::finalize(const OperatorTree *OpLeaf)
{

    if (_mat.size() == 0)
        DEVABORT("no matrix calculated")

    // post-multiply by function of r if desired
    std::string def = OpLeaf->def();
    if (def.find("*") != std::string::npos)
    {
        ABORT("Multiplication by a function of r not implemented in XC");
        // std::string algStr;
        // if(def.find("*XC")!=std::string::npos)
        //     algStr=tools::stringInBetween(OpLeaf->def(),"<","*XC",true);
        // else if(def.find("XC")>def.find("*"))
        //     ABORT("pre-multiply only as (Algebra of Q)*XC, found: "+def)
        //             else
        //             algStr=tools::stringInBetween(OpLeaf->def(),"*",">");
        // if(algStr.find("Q")!=std::string::npos)
        //     ABORT("must not use \"Q\" in prefactor, use axis name as the variable instead, here"+OpLeaf->iIndex->axisName())
        //             std::string algOrig(algStr);
        // algStr=tools::substringReplaceAll(algOrig,OpLeaf->iIndex->axisName(),"Q");

        // Algebra alg(algStr);
        // if(not alg.isAlgebra())
        //     ABORT("error in algebra string "+def+"\n"+Algebra::failures+"\nPrefactor after replacement: "+algStr);
        // std::vector<double> nodes=dynamic_cast<const BasisDVR*>(OpLeaf->iIndex->basis())->nodes();
        // for(size_t k=0;k<nodes.size();k++)_mat(k)*=alg.val(nodes[k]);
    }
    // this populates the operator floor
    MPIwrapper::AllreduceSUM(_mat.data(), _mat.size());
    if (MPIwrapper::Rank() == Parallel::floorHost(OpLeaf->idx(), OpLeaf->jdx()))
        construct(_mat, "XC");
    // data is stored in floor, release temporary memory
    _mat.resize(0, 0);
}
void OperatorFloorXC::test(const OperatorTree *Op, const BasisOrbitalNumerical &IOrb)
{
#ifndef _USE_HACC_
    return;
#else

    const std::vector<int> &s = IOrb.select();

    for (const OperatorTree *ppCC = Op->descend(); ppCC != 0; ppCC = ppCC->nodeNext())
    {
        if (ppCC->def().find("XC") != std::string::npos and ppCC->parent()->iIndex->axisName() == "Channel" and ppCC->iIndex->parent() == ppCC->jIndex->parent() and ppCC->iIndex->axisName() == "Phi")
        {

            // if(ppCC->def().find("HartreeRelative")!=string::npos)continue; // this cannot compare to full Hartree

            const BasisChannel *bChan = dynamic_cast<const BasisChannel *>(ppCC->parent()->iIndex->basis());
            // non-zero channel blocks
            int ci = ppCC->iIndex->nSibling();
            int cj = ppCC->jIndex->nSibling();
            // std::cout<<"ci="<<ci<<" cj="<<cj<<std::endl;
            // for(int i=0;i<bChan->rho1().size();i++)
            // for(int j=0;j<bChan->rho1()[0].size();j++)
            // std::cout<<"bChan->rho1("<<i<<","<<j<<") "<<std::endl<<bChan->rho(i,j)<<std::endl;
            // std::cout<<"bChan->rho1().size()"<<bChan->rho1().size()<<std::endl;
            Eigen::MatrixXcd rho(bChan->rho(ci, cj));
            // std::cout<<"rho"<<rho<<std::endl;
            Eigen::MatrixXcd potNum = Eigen::MatrixXcd::Zero(IOrb.size(), IOrb.size());
            Eigen::MatrixXcd potOrb = Eigen::MatrixXcd::Zero(IOrb.size(), IOrb.size());
            Coefficients tmp(ppCC->iIndex);
            for (size_t j = 0; j < IOrb.size(); j++)
            {
                *ppCC->tempRHS() = *IOrb.orbital(j);
                ppCC->apply(1., *ppCC->tempRHS(), 0., tmp);
                for (size_t i = 0; i < IOrb.size(); i++)
                {
                    potNum(i, j) = IOrb.orbital(i)->innerProduct(&tmp);
                    for (size_t k = 0; k < IOrb.size(); k++)
                    {
                        for (size_t l = 0; l < IOrb.size(); l++)
                        {
                            double vIJKL = dynamic_cast<const BasisMO *>(IOrb.orbitalFunctions())->vinayMO()->vEE(s[i], s[l], s[k], s[j]);
                            // std::cout << "The VIJKL are " << s[i] << " " << s[j] << " " << s[l] << " " << s[k] << std::endl;
                            potOrb(i, j) += rho(k, l) * vIJKL;
                        }
                    }
                }
            }

            MPIwrapper::AllreduceSUM(potNum.data(), potNum.size());
            EigenTools::purge(potNum, 1.e-10);
            EigenTools::purge(potOrb, 1.e-10);
            double err = 0.;
            double maxVal = potOrb.lpNorm<Eigen::Infinity>();
            for (int j = 0; j < potNum.cols(); j++)
                for (int i = 0; i < potNum.rows(); i++)
                {
                    if (potNum(i, j) + potOrb(i, j) != 0.)
                        err = max(err, abs(potNum(i, j) - potOrb(i, j)));
                }
            PrintOutput::message(Sstr + "verified XC for channel (" + ci + cj + ") against original potentials with relative error" + (tools::str(100 * err / maxVal, 2) + "% (abs=") + err + ")");
            if (err > maxVal * 1.e-2)
            {
                std::cout << "Calculated(XC): " << std::endl;
                std::cout << potNum << std::endl;
                std::cout << "From Columbus(XC): " << std::endl;
                std::cout << potOrb << std::endl;
                PrintOutput::warning("significant deviation -- may need larger basis");
            }
        }
    }
    // DEVABORT("test done")
#endif
}
