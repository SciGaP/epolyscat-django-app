// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef SCHROEDINGERMACHINE_H
#define SCHROEDINGERMACHINE_H

#include <string>
#include <memory>
#include <vector>
#include "Eigen/Dense"


class BasisDVR;
class Index;
class OperatorTree;
class Coefficients;

// expose symbols by name
#pragma once
#ifdef __cplusplus
extern "C" {
#endif

/// structure information for the training data
void sm_StrucTrain(int *SizeElem /** points on each element, size=trainNElem() */,
                   double *BElem /** element boundaries, size=1+trainNElem() */);
void sm_SizesTrain(const char File[], int *NPsi, int *SizePsi);///< training wave functions from File, NPsi functions, each SizePsi
int  sm_NElemTrain();                                          ///< number of elements in training data
void sm_DataTrain(const char File[], double *TrainData);       ///< read into TrainData - buffer double[NPsi*SizePsi]
void sm_Quadrature(int Order, int Nelem, double *Node0, double *Weig0);  ///< Lobattto quadrature
void sm_FuncPsiC(int Order, int Nelem, const double* X,const double *PsiIn, double *PsiOut,const double *Grad);

void sm_LossL2(int Nelem, const double* Elem, const double* Psi, double *Res, const double *Grad=0);
double sm_ErrorL2(float * Psi, int NElem, float * Elem);                      ///< given FE-DVR values Psi, L2-error on new discretization
double sm_ErrorL2_grad(float * Psi, int NElem, float * Elem, float* GradIn=0, float *GradRes=0);    ///< given FE-DVR values Psi, grad wrt FE sizes
double sm_Loss(int BatchSize);                                 ///< loss function for a given batch (stub)


#ifdef __cplusplus
}
#endif

//------------------------------------------------------------------//

/// Machine-learning the (TD)SE
class SchroedingerMachine
{

    static std::string _task;
    static std::string _potential;
    static bool _exactIntegrals;

    std::shared_ptr<const OperatorTree> _ham;
    std::shared_ptr<const OperatorTree> _ovr;

public:
    static void read();
    SchroedingerMachine(){}
    SchroedingerMachine(const Index* Idx, std::string Potential);
    static bool exactIntegrals(bool ExactIntegrals);
    static void enter(const Index* Idx); ///< enter I/O cycle



    Eigen::MatrixXd LossL2(const Eigen::VectorXf & Elem, const Eigen::MatrixXcd & Psi, const Eigen::VectorXf Grad=Eigen::VectorXf()) const;
    Eigen::MatrixXcd funcB(const Eigen::VectorXf &Elem, const Eigen::MatrixXcd &Psi, const Eigen::VectorXf Grad=Eigen::VectorXf()) const;
    Eigen::MatrixXcd funcA(const Eigen::VectorXf &Elem, const Eigen::MatrixXcd &Psi, const Eigen::VectorXf Grad=Eigen::VectorXf()) const;
    Eigen::MatrixXd funcW(const Eigen::VectorXf &Elem, const Eigen::VectorXf Grad=Eigen::VectorXf()) const;
    Eigen::MatrixXd funcPhi(const Eigen::MatrixXd &X, const Eigen::MatrixXcd &Psi, const Eigen::MatrixXcd Grad) const;
    Eigen::VectorXd funcC(const Eigen::VectorXf &Elem, const Eigen::VectorXf Grad=Eigen::VectorXf()) const;
    Eigen::MatrixXd funcX(const Eigen::VectorXf &Elem, const Eigen::VectorXf Grad=Eigen::VectorXf()) const;

    void test(const Eigen::MatrixXd &X, const Eigen::MatrixXcd &Psi, const Eigen::MatrixXcd Grad) const;

};

#endif // SCHROEDINGERMACHINE_H
