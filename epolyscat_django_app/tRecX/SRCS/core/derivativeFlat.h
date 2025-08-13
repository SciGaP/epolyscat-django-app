// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef DERIVATIVEFLAT_H
#define DERIVATIVEFLAT_H

#include <vector>
#include <memory>

#include "parallel.h"
#include "coefficients.h"
#include "coefficientsLocal.h"
#include "coefficientsGlobal.h"
#include "derivativeBlock.h"
#include "arpack.h"
#include "operatorAbstract.h"
#include "operatorTree.h"


class OperatorTree;
class Wavefunction;
class DiscretizationSpectral;
class DiscretizationSpectralProduct;
class OperatorFloor;
class ParallelProcess;
class CoefficientsLocal;
class Inverse;
class ProjectSubspace;

#include "linSpaceMap.h"

/// derivative with operator tree flattened for performance
class DerivativeFlat : public OperatorAbstract
{
    friend class RungeKutta4;
    friend class Parallel;

    double _time;
    double start;
    bool precon; //temporary for debugging

    void setNonlinUpdate(); // setup non-linear update (not really used at present)
public:
    void updateNonLin(double Time, Coefficients *C);
    const OperatorAbstract *o;     // pointer to full operator
protected:
    Coefficients *setupXY;     // application is tied to a given input/output storage!
    CoefficientsGlobal *globXY;     // application is tied to a given input/output storage!

    CoefficientsLocal *localXY;     // application is tied to a given input/output storage!
    Coefficients *setupTemp;
    CoefficientsGlobal *globTemp;     // application is tied to a given input/output storage!
    CoefficientsLocal *localTemp;     // application is tied to a given input/output storage!

    const Inverse* inverseOverlap;
    std::shared_ptr<Coefficients> _lastUpdateVector;

    void test();
    void testProjection();

    ///////////////////////////////////
    class FlattenedOperatorTree{
        std::vector<DerivativeBlock> blocks;
        void add(std::complex<double> *Fac, const OperatorTree *Op, Coefficients* ICoeff, Coefficients* JCoeff,
                 bool AbsorbInverse, double ApplyEpsilon, const std::vector<unsigned int> &ISort,
                 const std::vector<unsigned int> &JSort);
    public:

        virtual ~FlattenedOperatorTree();//{delete par;}
        Parallel* par;
        FlattenedOperatorTree(const OperatorAbstract* Map, bool AbsorbInverse, double ApplyEpsilon, Coefficients* ICoeff,
                              Coefficients* JCoeff, std::string SendRecv);
    };
    FlattenedOperatorTree* op;
    ///////////////////////////////////
protected:
    class Projection{
    public:
        virtual ~Projection(){}
        virtual void apply() =0;
    };

    ///////////////////////////////////
    class ProjectionSingle: public Projection{
        Coefficients* cSpec;
        Coefficients* cXY;

        CoefficientsLocal* localXY;
        CoefficientsGlobal* globalXY;
        Coefficients* cTemp;
        CoefficientsLocal* localTemp;
        CoefficientsGlobal* globalTemp;
        CoefficientsLocal* localSpec;
        CoefficientsGlobal* globalSpec;

        FlattenedOperatorTree* mapFrom;
        FlattenedOperatorTree* mapTo;
//        std::shared_ptr<const Eigen::FullPivLU<Eigen::MatrixXcd>> _lu;
        std::shared_ptr<const Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>>> _lu;

    public:
        ~ProjectionSingle(){delete mapFrom; delete mapTo;}
        ProjectionSingle(const DiscretizationSpectral* ProjectionDisc, Coefficients* Coeff, double ApplyEpsilon);
        ProjectionSingle(const ProjectSubspace* Projection, Coefficients* Coeff, double ApplyEpsilon);
        void apply();
    };

    ///////////////////////////////////
    class ProjectionProduct: public Projection{
        Coefficients* cXY;
        CoefficientsLocal* localXY;
        CoefficientsGlobal* globalXY;
        Coefficients* cTemp;
        CoefficientsLocal* localTemp;
        CoefficientsGlobal* globalTemp;
        std::vector<Coefficients*> cSpec;
        std::vector<CoefficientsLocal*> localSpec;
        std::vector<CoefficientsGlobal*> globalSpec;
        std::vector<FlattenedOperatorTree*> mapFrom;
        std::vector<FlattenedOperatorTree*> mapTo;


    public:
        ProjectionProduct(const DiscretizationSpectralProduct* ProjectionDisc, Coefficients* Coeff, double ApplyEpsilon);
        void apply();
    };

    std::shared_ptr<ProjectSubspace> _projSub;
    std::shared_ptr<Projection> projection;
//    std::shared_ptr<Projection> subspace;
protected:
    double applyEpsilon;          ///< threshold for operator application
    void _construct(const OperatorTree* Op, const DiscretizationSpectral *ProjectionDisc,std::shared_ptr<ProjectSubspace> Project);
    void applyA(std::complex<double> A, CoefficientsLocal *localX) const;
    void applyB(std::complex<double> B, CoefficientsLocal &Y) const;
    void updateNonLin(const Coefficients* C, double Time);
public:
    static bool applyFlat;
    bool isEmpty() const {return op==0;}
    virtual const Index* idx() const{ return iIndex;}
    virtual ~DerivativeFlat();
    DerivativeFlat():o(0),setupXY(0),localXY(0),setupTemp(0),op(0),projection(0){}
    DerivativeFlat(const OperatorTree* Op, double ApplyThreshold, std::shared_ptr<ProjectSubspace> Project)
        :OperatorAbstract("der["+Op->name+"]",Op->iIndex,Op->jIndex),start(0.),o(Op),op(0),projection(0),applyEpsilon(ApplyThreshold)
    {_construct(Op,0,Project);}
    DerivativeFlat(const OperatorTree* Op, double ApplyThreshold, const DiscretizationSpectral *ProjectionDisc=0)
        :OperatorAbstract("der["+Op->name+"]",Op->iIndex,Op->jIndex),start(0.),o(Op),op(0),projection(0),applyEpsilon(ApplyThreshold)
    {_construct(Op,ProjectionDisc,0);}


    void eigen(std::vector<std::complex<double> > & Eval,std::vector<Coefficients*> &Evec, unsigned int NStat);
    void eigenValues(double Time, std::vector<std::complex<double> > Eval); ///< eigenvalues of the derivative operator

    const Coefficients & lhsVector() const;
    const Coefficients & rhsVector() const {return lhsVector();}
    void update(double Time, const Coefficients* CurrentVec=0);

    bool applyAlias() const {return true;} ///< input and output vector can be identical (will be copied internally)
    virtual void apply(std::complex<double> A, const Coefficients& X, std::complex<double> B, Coefficients &Y) const;
    virtual void apply(std::complex<double> A, CoefficientsLocal *localX, std::complex<double> B, CoefficientsLocal &Y) const;

    void project(Coefficients & C) const;

    class Arp:public Arpack{
        static std::complex<double> shift; // subtract shift*S from Hamiltonian to move desired eigenvalues to below 0
        static std::complex<double> shiftZero; // add shiftZero*P to Hamiltonian to push projected subspace to high eigenvalues
        const DerivativeFlat* der;
        std::vector<std::complex<double> *> pX,pY;
        CoefficientsGlobal x,y;
        CoefficientsLocal locX,locY;

        void apply(const std::complex<double> *X, std::complex<double> *Y);
    public:
        Arp(const DerivativeFlat* Der);

        /// return eigenvectors
        void eigen(std::vector<std::complex<double> > &Eval, std::vector<Coefficients* > &Rvec,
                   unsigned int Nvec, const std::string &Which, bool Restart);
    };

};


#endif // DERIVATIVEFLAT_H
