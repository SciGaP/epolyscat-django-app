// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "parallel.h"
#include "operatorHartree.h"
#include "basisOrbitalNumerical.h"
#include "multipolePotential.h"
#include "index.h"
#include "basisExpIm.h"
#include "basisAssocLeg.h"
#include "densityMatrix1.h"
#include "operatorFloorXC.h"
#include "gaunts.h"
#include "patchRadial.h"

#include "timer.h"
#include "printOutput.h"

//#define _USE_HACC_
#ifdef _USE_HACC_
#include "my4darray_shared.h"
#include "mo.h"
#endif

#include "eigenTools.h"
#include "basisChannel.h"
#include "gaunt.h"
#include "basisMO.h"
#include "basisDvr.h"
#include "vectorValuedFunction.h"
#include "basisMatMatrix.h"
#include "algebra.h"

#include "readInput.h"
#include "inverse.h"
#include "plot.h"

using namespace std;

std::string potName(std::string OpDef)
{
    std::string pot = OpDef.substr(OpDef.find("Hartree"));
    if (pot == OpDef)
        ABORT("no Hartree potential in " + OpDef);
    if (tools::findFirstOutsideBrackets(pot, "[", "{(", "})") == std::string::npos)
        return "CoulombEE";
    pot = tools::stringInBetween(OpDef, "[", "]");
    if (pot.find("Hartree") != string::npos)
        return "CoulombEE";
    return pot;
}

OperatorHartree::OperatorHartree(std::string Pot, const Index *IIndex, const Index *JIndex, std::complex<double> Multiplier)
{
    // matrix elements only appear within the same radial element
    // also, in the complex scaled region the mean field is ill-defined for now and will not be used
    // no other information is created here
    // these serve as indicatorst that operator blocks should be set up in postProcess(...)
    dat = 0;
    oNorm = 1.;
    if (IIndex->basis()->isAbsorptive() or IIndex->basis()->integrable()->lowBound() != JIndex->basis()->integrable()->lowBound() or IIndex->basis()->integrable()->upBound() != JIndex->basis()->integrable()->upBound())
        oNorm = 0.;
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

// in given patch, subtract from all (diagonal) channels the reference channel matrix
// this is auxiliary, not needed for standard Hartree setup
static void hartreeRelative(PatchRadial<OperatorHartree> &Patch, const DensityMatrix1 &Dens)
{
    if (Patch.firstC().op->def().find("HartreeRelative{") == string::npos)
        return; // do not form difference

    int RefChan = tools::string_to_int(tools::stringInBetween(Patch.firstC().op->name, "{", "}"));

    for (PatchRadial<OperatorHartree>::LeafM &leafM : Patch.leafs)
    {
        for (PatchRadial<OperatorHartree>::LeafL &leafL : leafM.leafs)
        {
            Eigen::MatrixXcd matRef = *leafL.leafs[Dens.cacb(RefChan, RefChan)].mat;
            for (int c = 0; c < Dens.nChan(); c++)
            {
                if (leafL.leafs[Dens.cacb(c, c)].mat->size() == 0)
                    DEVABORT("matrix no set up");
                *leafL.leafs[Dens.cacb(c, c)].mat -= matRef;
            }
        }
    }
}

/// maximum angular momentum l found in Idx
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

/// basA-wise multipole potentials
///
///     Wml[basA][M][L](a,cc)...int[ds] V[L](a,s) sum[mi,mj,li,lj] gaunt(li,mi,L,M,lj,mj) rho{cc}[mi,mj,li,lj](s,s)
///
/// (see notes)
static void getPotWml(std::string Pot, const std::vector<std::vector<Eigen::MatrixXcd>> &RhoIJ,
                      const BasisOrbitalNumerical &IOrb, const Index *IdxA,
                      std::map<double, std::map<int, std::vector<Eigen::MatrixXcd>>> &Wml)
{
    Gaunts gaunts;
    std::vector<const BasisIntegrable *> radBasS(radialBases(IOrb.orbital(0)->idx()));
    std::vector<const BasisIntegrable *> radBasA(radialBases(IdxA));

    DensityMatrix1 rho(RhoIJ);                                              // spherical repres. of density matrix rho for present radial basis
    int lpMax = std::min(2 * lMax(IdxA), 2 * lMax(IOrb.orbital(0)->idx())); // maximal possible lp

    // loop through radial basis
    for (const BasisIntegrable *basS : radBasS)
    {
        rho.patch(*basS, *basS, IOrb, IOrb); // set density to present (diagonal) x-patch

        // mpPots[basA][lp](a,s) = V[lp](a,s)*basA(a)*basS(s)
        // radial multipole pots multiplied by basis values at DVR points for given basS
        // below is very specific for DVR rep, see notes for explanations
        std::map<double, std::vector<Eigen::MatrixXcd>> mpPots;
        const BasisDVR *sDvr = dynamic_cast<const BasisDVR *>(basS);
        for (const BasisIntegrable *basA : radBasA)
        {
            // product of values ad DVR points
            const BasisDVR *aDvr = dynamic_cast<const BasisDVR *>(basA);
            if (aDvr == 0 or sDvr == 0)
                DEVABORT("for now, need dvr bases, got\n" + basA->str() + "\n" + basS->str());
            std::vector<double> aVal(aDvr->valNodes());
            std::vector<double> xVal(sDvr->valNodes());
            MultipolePotential mpPot(lpMax, Pot, basA, basS);
            Eigen::MatrixXd prodAS = Eigen::Map<Eigen::MatrixXd>(aVal.data(), aVal.size(), 1) * Eigen::Map<Eigen::MatrixXd>(xVal.data(), 1, xVal.size());
            for (int lp = 0; lp <= lpMax; lp++)
            {
                mpPots[basA->lowBound()].push_back(mpPot.vals(lp).cwiseProduct(prodAS));
            }
        }

        // the loop below tries to optimize memory access and avoid re-fetch and re-compute
        // as much as possible, while not bloating memory

        // loop through points in basS
        for (size_t s = 0; s < basS->size(); s++)
        {
            // set rho to the given radial values
            rho.setPoint(s, s);
            // sig[M][L](1,cc)=sum[li,mi,mj,lj] gaunt(li,mi|L,mi-mj,lj,mj) rho{cc}[mi,mj,li,lj]
            std::map<int, std::vector<Eigen::MatrixXcd>> sig;
            // rho is expanded into spherical basis (see notes), this runs through allowed im,jm, lhs/rhs m-values
            for (int mi = rho.imMin(); mi <= rho.imMax(); mi++)
                for (int mj = rho.jmMin(); mj <= rho.jmMax(); mj++)
                {
                    // the storage of the density matrix is structured hiearchically
                    // mi,mj is on top, get that block
                    DensityMatrix1::M *rhoM = rho(mi, mj);

                    // get storage where sig[M][L] is accumulated (depends only on M=mi-mj)
                    // sig is the sigma in notes
                    std::vector<Eigen::MatrixXcd> &sigM = sig[mi - mj];
                    if (sigM.size() == 0)
                        sigM.assign(rhoM->ilMax() + rhoM->jlMax() + 1, Eigen::MatrixXcd::Zero(1, rho.cacbMax() + 1));

                    // this runs through the li where gaunts may be non-zero
                    for (int li : rhoM->listLi())
                    {
                        DensityMatrix1::L *rhoML = (*rhoM)(li); // get the li-slice if the density matrix
                        //  run through lj
                        for (int lj : rhoM->listLj())
                        {
                            int lpMin;
                            const std::vector<double> &gauntIJ = gaunts.vals(-mi, mj, li, lj, lpMin);
                            DensityMatrix1::LL *rhoMLL = (*rhoML)(lj); // rho[mi,mj,li,lj](s,s;cc)
                            for (size_t gp = 0; gp < gauntIJ.size(); gp++)
                            {
                                // compressed storage for the contributions, only where non-zero
                                // sig is defined for all l, but gaunt is non-zero only between
                                // lpMin and lpMax, with step size 2: l=lpMin+2*gp
                                sig[mi - mj][lpMin + 2 * gp] += gauntIJ[gp] * Eigen::Map<Eigen::MatrixXcd>(rhoMLL->densC.data(), 1, rhoMLL->densC.size());
                            }
                        }
                    }
                }

            // at given s, loop through all radial bases
            for (const BasisIntegrable *basA : radBasA)
            {
                // at given s, loop  through all sig[M]
                for (auto &sigM : sig)
                {
                    // sig is a map, sigM.first = M
                    // Wml is a map of maps: first index the lower boundary of the interval that it refers to
                    std::vector<Eigen::MatrixXcd> &wmlA = Wml[basA->lowBound()][sigM.first];
                    // means wmlA starts at current element boundary and mi-mj=M
                    if (wmlA.size() == 0)
                        wmlA.assign(lpMax + 1, Eigen::MatrixXcd::Zero(basA->size(), sigM.second[0].cols()));
                    // sum all multipole terms that are non-zero at present patch
                    for (size_t lp = 0; lp < wmlA.size(); lp++)
                    {
                        Wml[basA->lowBound()][sigM.first][lp] += mpPots[basA->lowBound()][lp].col(s) * sigM.second[lp];
                    }
                }
            }
        }
    }
}

static void setOrbRho(std::string Def, const BasisOrbitalNumerical *&IOrb, std::vector<std::vector<Eigen::MatrixXcd>> &RhoIJ)
{
    std::vector<std::string> part = tools::splitString(Def, '@');
    if (IOrb != 0)
    {
        if (part.size() > 1 and part[1] != "")
            ABORT("orbitals defined internally as " + IOrb->str() + " do not specify in " + Def);
    }
    else
    {
        if (part.size() < 2)
            ABORT("Orbitals are not generated automatically, specify as " + Def + "@OrbitalDefinition");
        IOrb = dynamic_cast<const BasisOrbitalNumerical *>(BasisAbstract::factory("Orbital:" + part[1]));
        if (IOrb == 0)
            ABORT(part[1] + " in " + Def + " does not seem to define an orbital basis");
    }

    if (RhoIJ.size() != 0)
    {
        if (part.size() == 3)
            ABORT("rho defined internally - must not specify single particle density matrix in " + Def);
    }
    else
    {
        if (part.size() < 3)
            ABORT("rho not defined internally - must specify single particle density matrix in " + Def + " (after 2nd @)");
        const BasisMatMatrix *m = BasisMatMatrix::factory("<" + part[2] + ">", IOrb, IOrb);
        if (m->mat().size() == 0)
            ABORT("no single particle density matrix " + part[2] + " required in " + Def);
        RhoIJ = {{m->mat()}};
    }
}

TIMER(Hartree, )
TIMER(potwml, )
static bool firstPost = true;
void OperatorHartree::postProcess(OperatorTree *Op, const std::vector<std::vector<Eigen::MatrixXcd>> &RhoIJ, const BasisOrbitalNumerical *IOrb)
{

    std::vector<std::vector<Eigen::MatrixXcd>> rhoIJ(RhoIJ);
    std::string def = Op->def().substr(Op->def().find("Hartree"));
    setOrbRho(def.substr(0, def.find(">")), IOrb, rhoIJ);

    

    Gaunts gaunts;
    DensityMatrix1 rho(rhoIJ); // create the ml-structured density matrix
    // create first (non-finalized) radial patch
    // this logics is not great: patch is always the first non-finished
    PatchRadial<OperatorHartree> patch(Op, rho, "diagonal");
    if (patch.leafs.size() == 0)
        return; // already done

    IOrb->orbitals(); // make sure the orbitals are set up (a call to orbitals() finishes setup if needed)
    if (firstPost)
        IOrb->print("Setting up Hartree operator for " + IOrb->str());
    firstPost = false;

    // get mean field potentials in mp,lp-expansion
    std::map<double, std::map<int, std::vector<Eigen::MatrixXcd>>> wml;
    START(potwml);
    getPotWml(potName(Op->def()), rhoIJ, *IOrb, Op->iIndex, wml);
    STOP(potwml);
    // a patch, in principles, refers to a 2-dim area with aBas() x bBas
    // in Hartree, only aBas=bBas are non-zero

    //PARALL - parallelization of evaluation
    // assuming that all time and memory is in the loops below parallelize only that loop
    //

    // compute M{cc}[basA;ma,la,a;mb,lb] as sum over multipoles
    // PatchRadial logics is not great: patch is always the first non-finished in Op
    // this way we here have a loop through patches of Op until all are done
    START(Hartree);
    for (; patch.leafs.size() > 0; patch = PatchRadial<OperatorHartree>(Op, rho, "diagonal"))
    {
        // alternative
        // while(patch.next(Op,rho,diagonal){...}
        std::cout << "patch.aBas()->lowBound(): " << patch.aBas()->lowBound() << " patch.bBas()->lowBound(): " << patch.bBas()->lowBound() << std::endl;
        if (not(*patch.aBas() == *patch.bBas()))
            DEVABORT("non-diagonal radial term in Hartree");
        // loop through all pairs ma,mb
        for (PatchRadial<OperatorHartree>::LeafM &leafM : patch.leafs)
        {
            int mp = leafM.ma - leafM.mb;
            // loop through all la,lb
            for (PatchRadial<OperatorHartree>::LeafL &leafL : leafM.leafs)
            {
                int lpMin;
                std::vector<double> gAB = gaunts.vals(-leafM.ma, leafM.mb, leafL.la, leafL.lb, lpMin);
                // loop over channels cc
                //PARALL - in PatchRadial, we set only leafL.leafs on some threads, redistribute will be in a final step
                for(PatchRadial<OperatorHartree>::LeafC & leafC: leafL.leafs){
                    // sum over multipole L
                    for (size_t gp = 0; gp < gAB.size(); gp++)
                    {
                        *leafC.mat += wml[patch.aBas()->lowBound()][mp][lpMin + 2 * gp].col(leafC.cc) * gAB[gp];
                    }
                }
            }
        }
        hartreeRelative(patch, rho); // potential relative to one given channel (if desired)
    }
    
    // final setup of floors
    OperatorHartree* floor;
    for(OperatorTree* op=Op->firstLeaf();op!=0;op=op->nextLeaf()){

        if(0!=(floor=dynamic_cast<OperatorHartree*>(op->floor())) and floor->dat==0){
            floor->finalize(op); // put data into the floor
            // std::cout<<"Inside finalize"<<std::endl;
            //PARALL move to floor its host
            // if(floor->dat==0)OperatorFloor::replace(op->floor()); // set OperatorDUM where no data
            // Parallel::moveFloor(op->floor(),Parallel::floorHost(op->idx(),op->jdx()));
        }
    }
     STOP(Hartree);
    PrintOutput::timerWrite();
    
    if (dynamic_cast<const BasisOrbitalNumerical *>(IOrb))
        test(Op, *IOrb);
   
    if (Op->name.find("Exp") != 0)
        return;

    // compute matrix elements for checks (if operator is used for computing Expectation values
    if (dynamic_cast<const BasisOrbitalNumerical *>(IOrb))
        moMatrixElements(Op, dynamic_cast<const BasisOrbitalNumerical *>(IOrb), rhoIJ[0][0]);
}

void OperatorHartree::moMatrixElements(const OperatorTree *Op, const BasisOrbitalNumerical *IOrb, const Eigen::MatrixXcd &RhoIJ)
{
    if (Op->name.substr(0, 3) != "Exp")
        return;
    const std::vector<int> &s = IOrb->select();

    // compute matrix elements
    Eigen::MatrixXd exp(IOrb->size(), IOrb->size());
    const OperatorTree *op = Op;
    while (op and not(op->jdx() == IOrb->orbital(0)->idx() and op->idx() == op->jdx()))
        op = op->nodeNext();
    if (op == 0)
        DEVABORT("failed to find orbital block");
    Coefficients o(op->jdx());
    Coefficients oI(op->idx());
    for (size_t j = 0; j < IOrb->size(); j++)
    {
        o = *IOrb->orbital(j);
        op->apply(1., o, 0., oI);
        for (size_t i = 0; i < IOrb->size(); i++)
            exp(i, j) = IOrb->orbital(i)->innerProduct(&oI).real();
    }
    PrintOutput::matrix(exp, 3);

    Eigen::MatrixXcd orig(IOrb->size(), IOrb->size());
    orig.setZero();

#ifdef _USE_HACC_
    for (size_t j = 0; j < IOrb->size(); j++)
    {
        for (size_t i = 0; i < IOrb->size(); i++)
        {
            for (size_t k = 0; k < IOrb->size(); k++)
            {
                for (size_t l = 0; l < IOrb->size(); l++)
                {
                    double vIJKL = dynamic_cast<const BasisMO *>(IOrb->orbitalFunctions())->vinayMO()->vEE(s[i], s[j], s[k], s[l]);
                    orig(i, j) += RhoIJ(k, l) * vIJKL;
                }
            }
        }
    }
    PrintOutput::matrix(orig.real(), 3);
#endif
}

void OperatorHartree::finalize(const OperatorTree *OpLeaf)
{


    if(_mat.size()!=0){
        // post-multiply by function of r if desired
        std::string def=OpLeaf->def();
        if(def.find("*")!=std::string::npos){
            std::string algStr;
            if(def.find("*Hartree")!=std::string::npos)
                algStr=tools::stringInBetween(OpLeaf->def(),"<","*Hartree",true);
            else if(def.find("Hartree")>def.find("*"))
                ABORT("pre-multiply only as (Algebra of Q)*Hartree, found: "+def)
                        else
                        algStr=tools::stringInBetween(OpLeaf->def(),"*",">");
            if(algStr.find("Q")!=std::string::npos)
                ABORT("must not use \"Q\" in prefactor, use axis name as the variable instead, here"+OpLeaf->iIndex->axisName())
                        std::string algOrig(algStr);
            algStr=tools::substringReplaceAll(algOrig,OpLeaf->iIndex->axisName(),"Q");

            Algebra alg(algStr);
            if(not alg.isAlgebra())
                ABORT("error in algebra string "+def+"\n"+Algebra::failures+"\nPrefactor after replacement: "+algStr);
            std::vector<double> nodes=dynamic_cast<const BasisDVR*>(OpLeaf->iIndex->basis())->nodes();
            for(size_t k=0;k<nodes.size();k++)_mat(k)*=alg.val(nodes[k]);
        }
        // this populates the operator floor
        construct(_mat,"Hartree");
        // data is stored in floor, release temporary memory
        _mat.resize(0,0);
    }
}

// test against input MO integrals
void OperatorHartree::test(const OperatorTree *Op, const BasisOrbitalNumerical &IOrb)
{
#ifndef _USE_HACC_
    return;

#else
    const std::vector<int>& s = IOrb.select();

    for(const OperatorTree* ppCC=Op->descend();ppCC!=0;ppCC=ppCC->nodeNext()){
        if(ppCC->def().find("Hartree")!=string::npos
                and ppCC->parent()->iIndex->axisName()=="Channel"
                and ppCC->iIndex->parent()==ppCC->jIndex->parent()
                and ppCC->iIndex->axisName()=="Phi")
        {

            if (ppCC->def().find("HartreeRelative") != string::npos)
                continue; // this cannot compare to full Hartree

            const BasisChannel *bChan = dynamic_cast<const BasisChannel *>(ppCC->parent()->iIndex->basis());
            // non-zero channel blocks
            int ci = ppCC->iIndex->nSibling();
            int cj = ppCC->jIndex->nSibling();
            Eigen::MatrixXcd rho(bChan->rho(ci, cj));
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
                            double vIJKL = dynamic_cast<const BasisMO *>(IOrb.orbitalFunctions())->vinayMO()->vEE(s[i], s[j], s[k], s[l]);
                            potOrb(i, j) += rho(k, l) * vIJKL;
                        }
                    }
                }
            }

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
            PrintOutput::message(Sstr + "verified Hartree for channel (" + ci + cj + ") against original potentials with relative error" + (tools::str(100 * err / maxVal, 2) + "% (abs=") + err + ")");
            std::cout << "Calculated(Hartree): " << std::endl;
                std::cout << potNum << std::endl;
                std::cout << "From Columbus(Hartree): " << std::endl;
                std::cout << potOrb << std::endl;
            if (err > maxVal * 1.e-2)
            {
                
                PrintOutput::warning("significant deviation -- may need larger basis");
            }
        }
    }
    DEVABORT("test done")

        #endif
}
