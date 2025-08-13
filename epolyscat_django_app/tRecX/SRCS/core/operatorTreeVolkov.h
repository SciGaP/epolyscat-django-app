// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATOR_TREE_VOLKOV_H
#define OPERATOR_TREE_VOLKOV_H

#include <complex>

#include "operatorFloor.h"
#include "operatorTree.h"

class BasisGrid;


class OperatorTreeVolkov: public OperatorTree{
//public:
    class IntegratedVectorPotential{
        double _int_Ax;
        double _int_Ay;
        double _int_Az;
        double _int_Asq;

        double _t_last;

    public:
        IntegratedVectorPotential(): _t_last(-DBL_MAX), _int_Ax(0.), _int_Ay(0.), _int_Az(0.){}
        static IntegratedVectorPotential main;

        void set(double _time);

        double int_Asq() const; ///< integral of A*A
        double parallel_int_A() const; ///< component integral over vector A parallel to (rotated) z axis

        // Euler angles in z-y-z to rotate the z-axis parallel to \int \vec A
        // TODO! Fix active/passive interpretation
        // Euler angle \alpha can be chosen to be zero
        double beta_int_A() const; ///< Euler \beta angle of integral over vector A
        double gamma_int_A() const; ///< Euler \gamma angle of integral over vector A
        double dt() const; ///< Elapsed time
    };

    class UnrotatedMatrixElements{
        std::vector<Eigen::VectorXd> quad;
        Eigen::VectorXd quad_points;

        std::vector<Eigen::VectorXd> quad_check;
        Eigen::VectorXd quad_points_check;

        void initialize(int lambda_max, int quad_order);

        double cache_a;
        std::vector<std::pair<int, Eigen::VectorXcd>> cache;

        double cache_iksq_a;
        Eigen::VectorXcd cache_iksq;

        std::vector<double> momenta;
    public:
        static UnrotatedMatrixElements main;
        
        void initialize(const BasisGrid *Basis);
        const Eigen::VectorXcd& get(double a, int lambda);
        const Eigen::VectorXcd& getExp_iksq(double a);
    };

    class WignerMatrix{
        /// d_{m0}(\beta) = \sqrt{(l-|m|)!/(l+|m|)!} P_l^|m|(cos \beta)
        std::vector<std::vector<double>> d_lm_pihalf;

        void initialize(int lambda_max);
    public:
        static WignerMatrix main;

        /**
         * Get the matrix element D^\lambda_{\mu, 0}(0, \beta, \gamma),
         * where \beta and \gamma are z-y-z Euler angles
         */
        std::complex<double> get(double beta, double gamma, int lambda, int mu);

        void check();
    };

public:
    struct AxisDescriptor{
        std::string phi;
        std::string eta;
        std::string k;
    };

private:
    class OperatorFloorVolkov: public OperatorFloor{
        int iM;
        int iL;
        int jM;
        int jL;

        // Apply Id_{idSizeBefore} \otimes Volkov \otimes Id_{idSizeAfter}
        int idSizeBefore;
        int idSizeAfter;

    public:
        OperatorFloorVolkov(const Index* IIndex, const Index* JIndex, const AxisDescriptor& Axis);

        void axpy(const std::complex<double> & Alfa, const std::complex<double>* X, unsigned int SizX,
                  const std::complex<double> & Beta, std::complex<double>* Y, unsigned int SizY) const;
        void pack(std::vector<int> &Info, std::vector<std::complex<double> > &Buf) const{
            ABORT("Not implemented");
        }
    };


public:
    OperatorTreeVolkov(const Index* IIndex, const Index* JIndex, const AxisDescriptor& Axis);

    void update(double Time, const Coefficients* CurrentVec) override;
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const override;
    static  double  takeAsq(double Time) ; 
};



#endif
