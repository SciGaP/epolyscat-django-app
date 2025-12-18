// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef RECURSIVEINTEGRATOR_H
#define RECURSIVEINTEGRATOR_H

#include <cfloat> // for DBL_MIN
#include <vector>
#include <ostream>
#include <cmath>
#include <functional>

#include "qtAlglib.h"    // get the quadrature rules

#include "tools.h"
#include "abort.h"

namespace tools {

//==== Lots of vector mathematics =========================
// this could also be done with a more general template
template<typename vec_type>
double abs(const std::vector<vec_type>& vec) { // max norm
    double result=0.;
    using std::abs;
    for (size_t i=0; i!=vec.size(); ++i) {
        if (abs(vec[i])>result) {result=abs(vec[i]);}
    }
//    if (result==0) {return DBL_MIN;}
    return result;
}

template<typename vec_type>
std::vector<vec_type> operator+ (const std::vector<vec_type>& lhs, const std::vector<vec_type>& rhs) {
    // assume equal length of both vectors
    std::vector<vec_type> result(lhs);
    for (size_t i=0; i!=result.size(); ++i) {
        result.at(i)+=rhs.at(i);
    }
    return result;
}

template<typename vec_type>
std::vector<vec_type>& operator+= (std::vector<vec_type>& lhs, const std::vector<vec_type>& rhs) {
    // assume equal length of both vectors
    for (size_t i=0; i!=lhs.size(); ++i) {
        lhs.at(i)+=rhs.at(i);
    }
    return lhs;
}

template<typename vec_type>
std::vector<vec_type> operator- (const std::vector<vec_type>& lhs, const std::vector<vec_type>& rhs) {
    // assume equal length of both vectors
    std::vector<vec_type> result(lhs);
    for (size_t i=0; i!=result.size(); ++i) {
        result.at(i)-=rhs.at(i);
    }
    return result;
}

template<typename vec_type>
std::vector<vec_type>& operator-= (std::vector<vec_type>& lhs, const std::vector<vec_type>& rhs) {
    // assume equal length of both vectors
    for (int i=0; i!=lhs.size(); ++i) {
        lhs.at(i)-=rhs.at(i);
    }
    return lhs;
}

template<typename vec_type, typename scalar_type>
std::vector<vec_type> operator* (const std::vector<vec_type>& lhs, scalar_type rhs) {
    std::vector<vec_type> result(lhs.size());
    for (size_t i=0; i!=result.size(); ++i) {
        result[i]=lhs[i]*rhs;
    }
    return result;
}

template<typename vec_type, typename scalar_type>
std::vector<vec_type> operator* (scalar_type lhs, const std::vector<vec_type>& rhs) {
    return rhs*lhs;
}

template<typename vec_type, typename scalar_type>
std::vector<vec_type> operator/ (const std::vector<vec_type>& lhs, scalar_type rhs) {
    return lhs*pow(rhs, -1);
}



//==== the integrator template class ======================
template<typename ReturnType>
class RecursiveIntegrator
{
    int cnt=0;
public:
    RecursiveIntegrator(std::function<ReturnType(double)> func, double eps_rel=1.e-10, double eps_abs=1.e-10): func(func), eps_rel(eps_rel), eps_abs(eps_abs) {}
    ~RecursiveIntegrator() {}
    RecursiveIntegrator(const RecursiveIntegrator<ReturnType>& other): func(other.func), eps_rel(other.eps_rel), eps_abs(other.eps_abs) {}
    RecursiveIntegrator<ReturnType>& operator= (const RecursiveIntegrator<ReturnType>& rhs)
    {
        func = rhs.func;
        eps_rel = rhs.eps_rel;
        eps_abs = rhs.eps_abs;
    }

    /*!
     * \brief integrate
     * \param par integration boundaries
     * \param degree order of quadrature routine
     * \param kind sort of quadrature routine (default Gauss-Legendre quadrature)
     * \return
     */
    const ReturnType integrate(const std::vector<double>& par, int degree=5, const quadrature_kind& kind=gq_legendre)
    {
        ReturnType comparison;
        static_integration(par, comparison, degree, kind);
        return wrapped_integrate(par, comparison, degree, kind);
    }

    /*!
     * \brief wrapped_integrate recursively integrate function
     * \param par integration boundaries of current subinterval
     * \param comparison comparison value for quadrature on the whole interval
     * \param degree order of quadrature routine
     * \param kind sort of quadrature routine
     * \return
     */

    const ReturnType wrapped_integrate(std::vector<double> par, ReturnType comparison, int degree=5, const quadrature_kind& kind=gq_legendre)
    {
        // set subintervals
        static std::vector<double> sub_par0(par), sub_par1(par);
        sub_par0[0]=par[0]; sub_par1[1]=par[1];
        sub_par0[1]=(par[0]+par[1])/2.; sub_par1[0]=(par[0]+par[1])/2.;

        // integrate over subintervals
        static ReturnType sub_interval0; // instantiate to 0
        static ReturnType sub_interval1; // instantiate to 0
        static_integration(sub_par0, sub_interval0, degree, kind);
        static_integration(sub_par1, sub_interval1, degree, kind);

        // get the errors
        static ReturnType returnTypeError;
        static double doubleError;
        returnTypeError=sub_interval0;
        returnTypeError+=sub_interval1;
        returnTypeError-=comparison;
        using std::abs;
        using tools::abs;
        doubleError=abs(returnTypeError);

        //compare with quadrature over the entire interval (comparison); continue recursion if error is insufficient
        if (doubleError<eps_abs or doubleError/(abs(comparison)+DBL_MIN)<eps_rel) {
            return comparison;
        }
        else {
            // continue iteration
            if(++cnt>100)DEVABORT("");
            ReturnType intermediateResult(wrapped_integrate(sub_par0, sub_interval0, degree, kind));
            sub_par1[1]=par[1];
            sub_par1[0]=(par[0]+par[1])/2.;
            intermediateResult += wrapped_integrate(sub_par1, sub_interval1, degree, kind);
            return intermediateResult;
        }
    }

    /*!
     * \brief static_integration basic quadrature routine
     * \param par integration boundaries
     * \param result is overriden by the result of the integration
     * \param degree order of quadrature
     * \param kind sort of quadrature
     */
    void static_integration(const std::vector<double>& par, ReturnType& result, int degree=5, const quadrature_kind& kind=gq_legendre)
    {
        static alglib::real_1d_array point, weigh; // would be nicer with vector, but alglib insists on its real_1d_array
        get_quadrature_rule(degree, kind, point, weigh, par);

        result=func(point[0])*weigh[0];
        for (int i=1; i!=point.length(); ++i) {
            result += func(point[i])*weigh[i];
        }
    }

private:
    std::function<ReturnType(double)> func;
    double eps_rel, eps_abs;

}; // RecursiveIntegrator

void test_recursiveIntegrator_double();

} // tools

#endif // RECURSIVEINTEGRATOR_H
