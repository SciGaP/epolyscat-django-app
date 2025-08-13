// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "multipolePotential.h"
#include "basisIntegrable.h"

#include "timer.h"
#include "eigenTools.h"

#include "basisDvr.h"

static std::complex<double> radialCoulombEE(int Lambda, std::complex<double> R, std::complex<double> S)
{
    // Coulomb multipole term: min(r,s)^la / max(r,s)^(la+1)
    if (R.real() < S.real())
        std::swap(R, S);
    std::complex<double> Q = 1. / R;
    return Q * std::pow(S * Q, Lambda);
}
static std::complex<double> delta12(int Lambda, std::complex<double> R, std::complex<double> S)
{
    if (R != S)
        return 0.;
    if (Lambda != 0)
        return 0.;
    return 1.;
}
std::map<std::string, MultipolePotential::multipoleRadial> MultipolePotential::_listPots;

std::map<std::string, std::shared_ptr<std::vector<Eigen::MatrixXcd>>> MultipolePotential::_listVals;

MultipolePotential::MultipolePotential(int LMax, std::string PotRadial, const BasisIntegrable *IBas, const BasisIntegrable *JBas)
{
    _listPots["CoulombEE"] = radialCoulombEE;
    _listPots["Delta"] = delta12;

    if (_listPots.find(PotRadial) == _listPots.end())
        ABORT("Potential " + PotRadial + " not defined, available: " + tools::listMapKeys(_listPots));

    std::string hash = Str(PotRadial, "|") + IBas->label() + JBas->label();
    if (not _listPots.count(hash))
    {
        _listVals[hash] = integrate(LMax + 1, _listPots[PotRadial], IBas, JBas);
    }
    _vals = _listVals[hash];
}

TIMER(intAll, )
TIMER(intFunc, )
TIMER(intPot, )
TIMER(intMult, )
TIMER(intTens, )
TIMER(intVal, )
/// performes a two-dimensional integration over rectangle or lower triangle
static void integrateMulti(std::vector<Eigen::MatrixXcd> &Res,
                           MultipolePotential::multipoleRadial Pot,
                           const BasisIntegrable *IBas, const BasisIntegrable *JBas,
                           bool LowTriangle /** integrate over lower triangle only */)
{
    std::vector<double> qI, wI, qJ, wJ;
    IBas->quadRule(Res.size() + IBas->order(), qI, wI);
    JBas->quadRule(Res.size() + JBas->order(), qJ, wJ);

    for (size_t i = 0; i < qI.size(); i++)
    {
        std::vector<std::complex<double>> iVal(IBas->val(qI[i]));
        double scalI = 1.;
        if (LowTriangle)
            scalI = (qI[i] - IBas->lowBound()) / (IBas->upBound() - IBas->lowBound());

        for (size_t j = 0; j < qJ.size(); j++)
        {
            double pJ = JBas->lowBound() + (qJ[j] - JBas->lowBound()) * scalI;
            double qwiqwj = wI[i] * wJ[j] * scalI;

            // NOTE: at present most time is here - replace by more efficient basis,
            // e.g. Legendres, transform after integration
            std::vector<std::complex<double>> jVal(JBas->val(pJ));
            Eigen::MatrixXcd loTensorUp = qwiqwj *
                                          Eigen::Map<Eigen::MatrixXcd>(iVal.data(), IBas->size(), 1) *
                                          Eigen::Map<Eigen::MatrixXcd>(jVal.data(), JBas->size(), 1).adjoint();
            if (Pot == radialCoulombEE)
            {
                // special case: real-valued Coulomb
                double q, r;
                if (pJ > qI[i])
                {
                    q = 1. / pJ;
                    r = qI[i] * q;
                }
                else
                {
                    q = 1. / qI[i];
                    r = pJ * q;
                }
                for (size_t la = 0; la < Res.size(); la++, q *= r)
                    Res[la] += loTensorUp * (q * 4 * math::pi / double(2 * la + 1));
            }
            else if (Pot == delta12)
            {
                if (qI[i] == qJ[j])
                {
                    for (size_t la = 0; la < Res.size(); la++)
                        Res[la] += loTensorUp * (4 * math::pi);
                    //                    Res[0]+=loTensorUp*(4*math::pi);
                }
            }
            else
            {
                // the complex version of Pot is particularly inefficient, avoid whenever possible
                for (size_t la = 0; la < Res.size(); la++)
                    Res[la] += loTensorUp * Pot(la, qI[i], pJ);
                ;
            }
        }
    }
    for (size_t la = 0; la < Res.size(); la++)
        if (Res[la].lpNorm<Eigen::Infinity>() > 1.e10)
            Sstr + "multipole error at" + la + Res[la].lpNorm<Eigen::Infinity>() + Sendl;
}

static bool intervalOverlap(double A, double B, double C, double D)
{
    double eps = std::min(B - A, D - C) * 1.e-12;
    return (A - eps < C and B - eps > C) or (C - eps < A and D - eps > A);
}
static bool intervalSame(double A, double B, double C, double D)
{
    double eps = std::min(B - A, D - C) * 1.e-12;
    return (std::abs(A - C) < eps and std::abs(B - D) < eps);
}

std::shared_ptr<std::vector<Eigen::MatrixXcd>> MultipolePotential::integrate(int LMax, multipoleRadial Pot, const BasisIntegrable *IBas, const BasisIntegrable *JBas)
{

    std::shared_ptr<std::vector<Eigen::MatrixXcd>> res;
    if (intervalSame(IBas->lowBound(), IBas->upBound(), JBas->lowBound(), JBas->upBound()))
    {
        // create transposed, will be transposed below
        res.reset(new std::vector<Eigen::MatrixXcd>(LMax + 1, Eigen::MatrixXcd::Zero(JBas->size(), IBas->size())));
        // add upper triangle
        integrateMulti(*res, Pot, JBas, IBas, true);
        for (auto &m : *res)
            m.transposeInPlace();
        // add lower triangle
        integrateMulti(*res, Pot, IBas, JBas, true);
    }
    else
    {
        if (intervalOverlap(IBas->lowBound(), IBas->upBound(), JBas->lowBound(), JBas->upBound()))
            goto DomainError;
        res.reset(new std::vector<Eigen::MatrixXcd>(LMax, Eigen::MatrixXcd::Zero(IBas->size(), JBas->size())));
        integrateMulti(*res, Pot, IBas, JBas, false);
    }
    // here we may truncate in lambda
    return res;

DomainError:
    DEVABORT(Str("intervals must either be equal or disjoint, found: [", "") + IBas->lowBound() + "," + IBas->upBound() + "] and [" + JBas->lowBound() + "," + JBas->upBound() + "]");

    return res; // shut up the Clang compiler
}

#include "basisMonomial.h"
void MultipolePotential::Test()
{
    std::vector<double> radii = {0., 1.};
    std::vector<int> order(radii.size() - 1, 3);
    int lmax = 2;

    std::vector<BasisMonomial> bas;
    for (int k = 0; k < order.size(); k++)
    {
        bas.push_back(BasisMonomial(order[k], radii[k], radii[k + 1]));
    }
    for (int k = 0; k < bas.size(); k++)
        for (int l = 0; l < bas.size(); l++)
        {
            MultipolePotential pots(lmax, "CoulombEE", &bas[k], &bas[l]);
            Sstr + "val [k,l]" + k + l + Sendl;
            for (int l = 0; l < lmax + 1; l++)
                Sstr + EigenTools::str(pots.vals(l) / (4 * math::pi) * (2 * l + 1), 7) + Sendl;
        }
    Sstr + "end of tests" + Sendl;
    exit(0);
}
