// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorTreeVolkov.h"

#include "pulse.h"
#include "gaunt.h"
#include "tools.h"
#include "recursiveIntegrator.h"
#include "orthopol.h"
#include "index.h"
#include "basisIntegrable.h"
#include "basisGrid.h"
#include "timer.h"

#include "qtEigenDense.h"
#include "eigenTools.h"

OperatorTreeVolkov::IntegratedVectorPotential OperatorTreeVolkov::IntegratedVectorPotential::main;

void OperatorTreeVolkov::IntegratedVectorPotential::set(double time){
    static tools::RecursiveIntegrator<std::complex<double> > iAx(Pulse::iAx, 1.e-10);
    static tools::RecursiveIntegrator<std::complex<double> > iAy(Pulse::iAy, 1.e-10);
    static tools::RecursiveIntegrator<std::complex<double> > iAz(Pulse::iAz, 1.e-10);
    static tools::RecursiveIntegrator<std::complex<double> > Asq(Pulse::Asqu, 1.e-10);

    if(_t_last == -DBL_MAX) _t_last = Pulse::gettBegin();

    std::complex<double> int_iAx = iAx.integrate({ _t_last, time }, 10);
    std::complex<double> int_iAy = iAy.integrate({ _t_last, time }, 10);
    std::complex<double> int_iAz = iAz.integrate({ _t_last, time }, 10);
    std::complex<double> int_Asq = Asq.integrate({ _t_last, time }, 10);

    if(std::abs(int_iAx.real()) > 1.e-9) ABORT("Unexpected");
    if(std::abs(int_iAy.real()) > 1.e-9) ABORT("Unexpected");
    if(std::abs(int_iAz.real()) > 1.e-9) ABORT("Unexpected");
    if(std::abs(int_Asq.imag()) > 1.e-9) ABORT("Unexpected");

    _int_Ax += int_iAx.imag();
    _int_Ay += int_iAy.imag();
    _int_Az += int_iAz.imag();
    _int_Asq += int_Asq.real();

    _t_last = time;
}

double OperatorTreeVolkov::IntegratedVectorPotential::int_Asq() const{
    return _int_Asq;
}
double OperatorTreeVolkov::IntegratedVectorPotential::parallel_int_A() const{
    if(std::abs(_int_Ax) < 1.e-12 and std::abs(_int_Ay) < 1.e-12)
        return _int_Az;
    else if(std::abs(_int_Az) < 1.e-12)
        return std::sqrt(_int_Ax*_int_Ax + _int_Ay*_int_Ay + _int_Az*_int_Az);
    else
        ABORT("Not implemented");
}
double OperatorTreeVolkov::IntegratedVectorPotential::beta_int_A() const{
    if(std::abs(_int_Ax) < 1.e-12 and std::abs(_int_Ay) < 1.e-12)
        return 0.;
    else if(std::abs(_int_Az) < 1.e-12)
        return math::pi/2.;
    else
        ABORT("Not implemented");
}
double OperatorTreeVolkov::IntegratedVectorPotential::gamma_int_A() const{
    if(std::abs(_int_Ax) < 1.e-12 and std::abs(_int_Ay) < 1.e-12)
        return 0.;
    else if(std::abs(_int_Az) < 1.e-12)
        return -std::atan2(_int_Ay, _int_Ax);
    else
        ABORT("Not implemented");
}
double OperatorTreeVolkov::IntegratedVectorPotential::dt() const{
    return _t_last - Pulse::gettBegin();
}

OperatorTreeVolkov::UnrotatedMatrixElements OperatorTreeVolkov::UnrotatedMatrixElements::main;

void OperatorTreeVolkov::UnrotatedMatrixElements::initialize(const BasisGrid* Basis){
    if(Basis == 0 or Basis->size()==0) ABORT("Could not find ks");
    std::vector<double> ks;
    for(int k=0; k<Basis->size(); k++){
        ks.push_back(Basis->mesh()[k]);
    }

    if(momenta.size() == 0){
        momenta = ks;
    }else{
        for(int k=0;k<momenta.size();k++){
            if(momenta.size()!=ks.size() or abs(momenta[k]-ks[k])>1.e-12){
                Sstr+"momenta"+momenta+Sendl;
                Sstr+"     ks"+ks+Sendl;
                ABORT("Not compatible k basis");
            }
        }
    }

}

TIMER(volkovGetQuad,)
void OperatorTreeVolkov::UnrotatedMatrixElements::initialize(int lambda_max, int quad_order){
    OrthogonalLegendre pol;

    auto get_quad = [&pol](int lambda_max, int quad_order,
            Eigen::VectorXd& quad_points, std::vector<Eigen::VectorXd>& quad){
        std::vector<double> xs, ws;
        pol.quadratureGauss(quad_order, xs, ws);

        quad_points = Eigen::VectorXd(xs.size());
        for(int i=0; i<xs.size(); i++) quad_points(i) = xs[i];

        quad.clear();
        for(int i=0; i<=lambda_max; i++) quad.push_back(Eigen::VectorXd(xs.size()));

        for(int i=0; i<xs.size(); i++){
            std::vector<double> pol_vals = pol.val(lambda_max + 1, xs[i]);
            for(int l=0; l<=lambda_max; l++){
                quad[l](i) = pol_vals[l] * ws[i];
                quad[l](i) *= 2*math::pi; // Integral over phi
                quad[l](i) *= std::sqrt((2*l+1)/(4*math::pi)); //Normalization of Y_l,0
            }
        }
    };

    STARTDEBUG(volkovGetQuad);
    get_quad(lambda_max, quad_order, quad_points, quad);

    // Comment this line to disable checking
    get_quad(lambda_max, quad_order + 20, quad_points_check, quad_check);
    STOPDEBUG(volkovGetQuad);
}

TIMER(volkovMatrixElements,)
const Eigen::VectorXcd& OperatorTreeVolkov::UnrotatedMatrixElements::get(double a, int lambda){
    STARTDEBUG(volkovMatrixElements);
    // First try cache
    if(cache.size() > 0 and std::abs(cache_a-a) < 1.e-12){
        for(auto& v: cache){
            if(v.first == lambda){
                STOPDEBUG(volkovMatrixElements);
                return v.second;
            }
        }
    }

    // Start with 20 point quadrature
    // Alternatively pick big enough quadrature and disable checking, should be about twice as fast
    if(quad.size() <= lambda) initialize(lambda, quad_points.rows() > 0 ? quad_points.rows() : 20);

Restart:
    Eigen::VectorXcd result(momenta.size());

    for(int i=0; i<momenta.size(); i++){
        Eigen::VectorXcd func_vals(quad_points.rows());
        for(int j=0; j<quad_points.size(); j++){
            func_vals(j) = std::exp(std::complex<double>{0., quad_points(j) * a * momenta[i] });
        }

        result(i) = quad[lambda].dot(func_vals);

        // Check if bigger quadrature gives different result, if so go back
        if(quad_check.size() > lambda){
            Eigen::VectorXcd func_vals(quad_points_check.size());
            for(int j=0; j<quad_points_check.rows(); j++){
                func_vals(j) = std::exp(std::complex<double>{0., quad_points_check(j) * a * momenta[i] });
            }

            std::complex<double> result_check = quad_check[lambda].dot(func_vals);

            if(std::abs(result(i) - result_check) > 1.e-9){
                initialize(quad.size(), quad_points.rows() + 20);
                PrintOutput::DEVmessage("UnrotatedMatrixElements: quadrature order="+std::to_string(quad_points.rows()));
                goto Restart;
            }
        }
    }

    if(std::abs(cache_a - a) > 1.e-12) cache.clear();
    cache_a = a;
    cache.push_back({ lambda, result });

    STOPDEBUG(volkovMatrixElements);
    return cache.back().second;
}

TIMER(volkovExpkSq,)
const Eigen::VectorXcd& OperatorTreeVolkov::UnrotatedMatrixElements::getExp_iksq(double a){
    STARTDEBUG(volkovExpkSq);
    // First try cache
    if(cache_iksq.rows() > 0 and std::abs(cache_iksq_a - a) < 1.e-12){
        STOPDEBUG(volkovExpkSq);
        return cache_iksq;
    }

    Eigen::VectorXcd result(momenta.size());
    for(int i=0; i<momenta.size(); i++){
        result(i) = std::exp(std::complex<double>{0., 1.} * a * momenta[i] * momenta[i]);
    }

    cache_iksq = std::move(result);
    cache_iksq_a = a;
    STOPDEBUG(volkovExpkSq);
    return cache_iksq;
}

OperatorTreeVolkov::WignerMatrix OperatorTreeVolkov::WignerMatrix::main;

void OperatorTreeVolkov::WignerMatrix::initialize(int lambda_max){
    // https://arxiv.org/pdf/1403.7698.pdf

    std::map<int, const BasisAbstract*> assocLegendre;

    for(int lambda = d_lm_pihalf.size(); lambda<=lambda_max; lambda++){
        d_lm_pihalf.push_back(std::vector<double>());
        for(int mu = 0; mu<=lambda; mu++){
            if(assocLegendre.find(mu) == assocLegendre.end()){
                assocLegendre[mu] = BasisAbstract::factory(
                            BasisSetDef(
                                lambda_max + 1, -1., 2., "assocLegendre", true, true, true,
                                Coordinate::fromString("Eta"), {}, ComplexScaling(), false, { double(mu) }
                                )
                            );
            }

            UseMatrix x(1, 1);
            x(0) = 0.; // cos(pi/2)
            UseMatrix Y_lm = assocLegendre[mu]->integrable()->val(x);

            d_lm_pihalf.back().push_back(
                        std::sqrt(2./(2.*lambda+1.)) * Y_lm(0, lambda-mu).real()
                        );
        }
    }

    // Conflicts with timer
    // check();
}

TIMER(volkovWignerMatrix,)
std::complex<double> OperatorTreeVolkov::WignerMatrix::get(double beta, double gamma, int lambda, int mu){
    if(beta == 0.) return mu==0 ? 1. : 0.;
    else if(beta == math::pi/2.){
        STARTDEBUG(volkovWignerMatrix);
        if(d_lm_pihalf.size() <= lambda) initialize(lambda);

        double f = 1.;
        if(mu < 0 and (-mu)%2==1) f=-1.;

        auto result = std::exp(std::complex<double>{ 0., double(mu)*gamma }) * f * std::pow(-1., mu) *
                d_lm_pihalf[lambda][std::abs(mu)];
        STOPDEBUG(volkovWignerMatrix);
        return result;
    }else
        ABORT("Not implemented")
}

void OperatorTreeVolkov::WignerMatrix::check(){
    get(math::pi/2, 0, 3, 0);

    // Values from Mathematica
    if(std::abs(get(math::pi/2, 0, 0,  0) - (1.         )) > 1.e-7) ABORT("Error in D^0_{00}");
    if(std::abs(get(math::pi/2, 0, 1, -1) - (-1./sqrt(2))) > 1.e-7) ABORT("Error in D^1_{-10}");
    if(std::abs(get(math::pi/2, 0, 1,  0) - (0.         )) > 1.e-7) ABORT("Error in D^1_{00}");
    if(std::abs(get(math::pi/2, 0, 1,  1) - (1./sqrt(2) )) > 1.e-7) ABORT("Error in D^1_{10}");
    if(std::abs(get(math::pi/2, 0, 2, -2) - (sqrt(6)/4. )) > 1.e-7) ABORT("Error in D^2_{-20}");
    if(std::abs(get(math::pi/2, 0, 2, -1) - (0.         )) > 1.e-7) ABORT("Error in D^2_{-10}");
    if(std::abs(get(math::pi/2, 0, 2,  0) - (-.5        )) > 1.e-7) ABORT("Error in D^2_{00}");
    if(std::abs(get(math::pi/2, 0, 2,  1) - (0.         )) > 1.e-7) ABORT("Error in D^2_{10}");
    if(std::abs(get(math::pi/2, 0, 2,  2) - (sqrt(6)/4. )) > 1.e-7) ABORT("Error in D^2_{20}");

    for(int i=0; i<d_lm_pihalf.size(); i++){
        for(int j=0; j<d_lm_pihalf[i].size(); j++){
            if(d_lm_pihalf[i][j] != d_lm_pihalf[i][j])
                ABORT("d_lm(pi/2), l="+std::to_string(i)+", m="+std::to_string(j-i)+" is NaN");
        }
    }
}

OperatorTreeVolkov::OperatorFloorVolkov::OperatorFloorVolkov(const Index* IIndex, const Index* JIndex,
                                                             const AxisDescriptor& Axis):
    OperatorFloor(IIndex->sizeStored(), JIndex->sizeStored(), "Volkov"){

    /*
     * TODO!
     * For one-dimensional problems Axis.phi==Axis.eta=="NONE" and iM==iL==jM==jL==0
     */
    iM = 0;
    iL = 0;
    jM = 0;
    jL = 0;

    for(const Index* i=IIndex; i->parent()!=0; i=i->parent()){
        if(i->parent()->axisName()==Axis.eta) iL=i->physical();
        else if(i->parent()->axisName()==Axis.phi) iM=i->physical();
    }
    for(const Index* i=JIndex; i->parent()!=0; i=i->parent()){
        if(i->parent()->axisName()==Axis.eta) jL=i->physical();
        else if(i->parent()->axisName()==Axis.phi) jM=i->physical();
    }

    if(IIndex->axisName() == Axis.k){
        idSizeBefore = 1;
        idSizeAfter = IIndex->descend()->basis()->size();

        UnrotatedMatrixElements::main.initialize(IIndex->basis()->grid());
    }else if(IIndex->descend()->axisName() == Axis.k){
        idSizeBefore = IIndex->basis()->size();
        idSizeAfter = 1;
        UnrotatedMatrixElements::main.initialize(IIndex->descend()->basis()->grid());
    }


    if(IIndex->axisName() == Axis.k){
        idSizeBefore = 1;
        idSizeAfter = IIndex->descend()->basis()->size();
        UnrotatedMatrixElements::main.initialize(IIndex->basis()->grid());
    }else if(IIndex->descend()->axisName() == Axis.k){
        idSizeBefore = IIndex->basis()->size();
        idSizeAfter = 1;
        UnrotatedMatrixElements::main.initialize(IIndex->descend()->basis()->grid());
    }else{
        ABORT("Could not find axis name '"+Axis.k+"' in '"+IIndex->hierarchy()+"'");
    }

    // If we do not have two electrons, this happens
    if(idSizeAfter == 0) idSizeAfter=1;
    if(idSizeBefore == 0) idSizeBefore=1;
    
}

void OperatorTreeVolkov::OperatorFloorVolkov::axpy(const std::complex<double> & Alfa,
                                                   const std::complex<double>* X, unsigned int SizX, const std::complex<double> & Beta,
                                                   std::complex<double>* Y, unsigned int SizY) const{

    if(Beta!=1.){
        if(Beta==0.)for(unsigned int k=0;k<SizY;k++)Y[k]=0.;
        else        for(unsigned int k=0;k<SizY;k++)Y[k]*=Beta;
    }

    int mu = iM - jM;

    // Collect exp(.5i\int k^2)
    Eigen::VectorXcd volkov = UnrotatedMatrixElements::main.getExp_iksq(.5*IntegratedVectorPotential::main.dt());


    for(int lambda=std::abs(iL-jL); lambda<=iL+jL; lambda++){
        if(Gaunt::main.coeff_isZero(iL, lambda, jL, iM, mu, jM)) continue;

        double gaunt = Gaunt::main.coeff(iL, lambda, jL, iM, mu, jM);

        std::complex<double> rotation_factor = WignerMatrix::main.get(
                    IntegratedVectorPotential::main.beta_int_A(),
                    IntegratedVectorPotential::main.gamma_int_A(),
                    lambda, mu);

        // Corresponds to rotation by inverse angle
        rotation_factor *= std::pow(-1., mu);

        if(std::abs(gaunt * rotation_factor) < 1.e-12) continue;

        // Collect exp(-i\vec k * \int \vec A)
        Eigen::VectorXcd tmp = rotation_factor * volkov.cwiseProduct(
                    UnrotatedMatrixElements::main.get(-IntegratedVectorPotential::main.parallel_int_A(), lambda));

        if(idSizeAfter == 1){
            int quot = SizY / idSizeBefore;
            if(SizX != quot*idSizeBefore or SizY != quot*idSizeBefore) ABORT("Mismatch");

            for(int j=0; j<SizY; j+=quot){
                Eigen::Map<Eigen::VectorXcd>(Y+j, quot) +=
                        Alfa * gaunt * tmp.cwiseProduct(Eigen::Map<const Eigen::VectorXcd>(X+j, quot));

            }
        }else if(idSizeBefore == 1){
            int quot = SizY / idSizeAfter;
            if(SizX != quot*idSizeAfter or SizY != quot*idSizeAfter) ABORT("Mismatch");

            for(int j=0; j<idSizeAfter; j++){
                Eigen::Map<Eigen::VectorXcd, 0, Eigen::InnerStride<>>(Y+j, quot, Eigen::InnerStride<>(idSizeAfter)) +=
                        Alfa * gaunt * tmp.cwiseProduct(
                            Eigen::Map<const Eigen::VectorXcd, 0, Eigen::InnerStride<>>(
                                X+j, quot, Eigen::InnerStride<>(idSizeAfter)));
            }

        }else{
            ABORT("Unsupported: "+std::to_string(idSizeBefore)+", "+std::to_string(idSizeAfter));
        }
    }
}

void OperatorTreeVolkov::apply(std::complex<double> A, const Coefficients &Vec,
                               std::complex<double> B, Coefficients &Y) const{

    /**
     * Kind of ugly bypassing of ParallelOperator
     */
    _apply(A, Vec, B, Y);
}

OperatorTreeVolkov::OperatorTreeVolkov(const Index* IIndex, const Index* JIndex, const AxisDescriptor& Axis):
    OperatorTree("Volkov", IIndex, JIndex){

    if(IIndex->hasFloor()){
        floor() = new OperatorFloorVolkov(IIndex, JIndex, Axis);
        return;
    }

    bool diagonalAxis = true;
    if(IIndex->axisName() == Axis.phi) diagonalAxis=false;
    if(IIndex->axisName() == Axis.eta) diagonalAxis=false;

    for(int i=0; i<IIndex->childSize(); i++){
        for(int j=0; j<JIndex->childSize(); j++){
            if(diagonalAxis and std::abs(IIndex->child(i)->physical() - JIndex->child(j)->physical()) > 1.e-6)
                continue;

            childAdd(new OperatorTreeVolkov(IIndex->child(i), JIndex->child(j), Axis));
        }
    }
}

void OperatorTreeVolkov::update(double Time, const Coefficients* CurrentVec){
    IntegratedVectorPotential::main.set(Time);
}

double OperatorTreeVolkov::takeAsq(double Time){
    IntegratedVectorPotential::main.set(Time);
    return IntegratedVectorPotential::main.int_Asq();
}