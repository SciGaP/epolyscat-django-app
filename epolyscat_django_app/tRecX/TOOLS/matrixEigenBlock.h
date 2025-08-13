#ifndef MATRIXEIGENBLOCK_H
#define MATRIXEIGENBLOCK_H

#include "matrixBlocking.h"
#include "printOutput.h"
#include "eigenTools.h"
#include "useMatrix.h"

#include <Sparse>

namespace tRecX{

//static int smallDimensions=100;

/// considering the matrix as a block of eigenvectors, (pseudo)-orthonormalizes Subset wrt Ovr
static void gramSchmidt(Eigen::MatrixXcd & Vec,const Eigen::MatrixXcd & Ovr,std::vector<unsigned int> Subset,bool Pseudo){
    if(Subset.size()==0)
        for(unsigned int k=0;k<Vec.cols();k++)Subset.push_back(k);
    Eigen::MatrixXcd OvrI(Vec.rows(),1);
    for (unsigned int i=0;i<Subset.size();i++){
        if(Ovr.size()!=0)OvrI=Ovr*Vec.col(Subset[i]);
        else             OvrI=Vec.col(Subset[i]);
        std::complex<double>ovrij;
        for(unsigned int j=0;j<i;j++){
            if(Pseudo) ovrij=(Vec.col(Subset[j]).transpose()*OvrI)(0,0);
            else       ovrij=(Vec.col(Subset[j]).adjoint()*OvrI)(0,0);
            Vec.col(Subset[i])-=Vec.col(Subset[j])*ovrij;
        }
        if(Pseudo) ovrij=(Vec.col(Subset[i]).transpose()*OvrI)(0,0);
        else       ovrij=(Vec.col(Subset[i]).adjoint()*OvrI)(0,0);
        if(abs(ovrij)<1.e-28)ABORT(Str("vectors are (pseudo-)linearly dependent")+ovrij+"pseudo="+Pseudo);
        Vec.col(Subset[i])*=1./sqrt(ovrij);
    }
}

/// return true if sucessfully orthonormalized or pseudo-orthonormalized
static bool eigenOrthonormalize(const Eigen::MatrixXcd &Val, Eigen::MatrixXcd &Vec, const Eigen::MatrixXcd &Ovr, bool Pseudo, double Eps=1e-12)
{
    // find blocks of near-degenerated values
    std::vector<bool> use(Val.size(),true);

    std::vector<unsigned int>block;
    do {
        block.clear();
        for(unsigned int k=0;k<Val.size();k++){
            if(use[k] and (block.size()==0 or
                           std::abs(Val(block[0])-Val(k))<Eps*std::max(1.,abs(Val(k))))){
                block.push_back(k);
                use[k]=false;
            }
        }
        // Schmidt-orthonormalize within block (may be single function!)
        if(block.size()>0)gramSchmidt(Vec,Ovr,block,Pseudo);

        // advance to next block
    } while(block.size()>0);
    return true;
}

static bool eigen(Eigen::MatrixXcd & Mat, Eigen::MatrixXcd & Eval, Eigen::MatrixXcd & Evec,
                  bool ComputeVectors, const Eigen::MatrixXcd& Ovr=Eigen::MatrixXcd(0,0))
{
    enum Eigen::DecompositionOptions comp=Eigen::EigenvaluesOnly;
    if(ComputeVectors)comp=Eigen::ComputeEigenvectors;

    if(EigenTools::isSelfadjoint(Mat,1e-12) and EigenTools::isSelfadjoint(Ovr,1e-12)){
        if(Ovr.size()==0){
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> slv(Mat,comp);
            Eval=slv.eigenvalues();
            if(ComputeVectors)Evec=slv.eigenvectors();

        }
        else if(Ovr.size()==Mat.size()){
            Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> slv(Mat,Ovr,comp);
            Eval=slv.eigenvalues();
            if(ComputeVectors)Evec=slv.eigenvectors();
        }
    }
    else {
        Eigen::SparseMatrix<std::complex<double>> sOvr;
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> slv;
        if(Ovr.size()!=0){
            sOvr=Ovr.sparseView();
            Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double> > > lu;
            lu.analyzePattern(sOvr);
            lu.factorize(sOvr);
            slv.compute(lu.solve(Mat),comp);
        }
        else
            slv.compute(Mat,comp);

        Eval=slv.eigenvalues();
        if(ComputeVectors)Evec=slv.eigenvectors();
    }
    return true;
}

static void eigenFull(Eigen::MatrixXcd &Mat, Eigen::MatrixXcd &Val,
                      bool RightVectors, Eigen::MatrixXcd &RVec,
                      bool DualVectors, Eigen::MatrixXcd &LVec,
                      const Eigen::MatrixXcd  &Ovr=Eigen::MatrixXcd(0,0))
{
    // solve sub-block eigenproblem
    eigen(Mat,Val,RVec,true,Ovr);
    double epsSymmetry=1e-7;
    if(DualVectors){
        std::string kind;
        if(EigenTools::isSelfadjoint(Mat,epsSymmetry) and EigenTools::isSelfadjoint(Ovr,epsSymmetry)){
            if(Ovr.size()>0)LVec=(Ovr*RVec).conjugate();
            else LVec=RVec.conjugate();
        }
        else if(EigenTools::isSymmetric(Mat,epsSymmetry) and EigenTools::isSymmetric(Ovr,epsSymmetry)){
            if(Ovr.size()>0)LVec=Ovr*RVec;
            else LVec=RVec;
        }
        else {
            PrintOutput::DEVwarning("non-s.a. and non-symm case - Duals are ill-defined");
            Eigen::MatrixXcd tMat=Mat.transpose();
            Eigen::MatrixXcd tOvr=Ovr.transpose();
            eigen(tMat,Val,LVec,true,tOvr);
        }
#ifdef _DEVELOP_
        if(not (EigenTools::isSelfadjoint(Mat,1.e-12) and EigenTools::isSelfadjoint(Ovr,1e-12)) and
                not  (EigenTools::isSymmetric(Mat,1.e-12) and EigenTools::isSymmetric(Ovr,1e-12)))
            PrintOutput::DEVwarning("Operators not s.a. or symmetric to 1e-12, duals may no be strictly orthogonal");
#endif
    }
    else if(RightVectors)
    {
        bool pseudo=not(EigenTools::isSelfadjoint(Mat,1e-12) and EigenTools::isSelfadjoint(Ovr,1e-12));
        if(pseudo and EigenTools::isSymmetric(Mat,1e-12) and EigenTools::isSymmetric(Ovr,1e-12)){
            if(not eigenOrthonormalize(Val,RVec,Ovr,pseudo)){
                std::string mess("encountered non-hermitian and non-symmetric matrix - Dual vectors may be wrong");
#ifdef _DEVELOP_
                PrintOutput::DEVwarning(mess);
#else
                DEVABORT(mess);
#endif
            }

        }
    }

    if(DualVectors and LVec.size()==0)DEVABORT("failed to compute dual vectors");
}

template<class MSparse>
void eigenBlock( const Eigen::SparseMatrix<std::complex<double>> &Mat, Eigen::MatrixXcd &Val,
                 bool RightVectors, MSparse &RVec,
                 bool DualVectors, MSparse &DVec,
                 const MSparse &Ovr=MSparse(0,0), int MinBlock=50){

    //HACK
    DualVectors=RightVectors;

    if(Mat.rows()<MinBlock){
        // do not block-decompose small matrices
        Eigen::MatrixXcd mat(Mat),ovr(Ovr),rvec,dvec;
        eigenFull(mat,Val,RightVectors,rvec,DualVectors,dvec,ovr);
        RVec=rvec.sparseView();
        DVec=dvec.sparseView();
    }

    else {

        // get blocking
        std::vector<const MSparse*> mats({&Mat,&Ovr});
        if(Ovr.size()==0)
            mats.pop_back();
        else if(Mat.rows()!=Ovr.rows() or Mat.cols()!=Ovr.cols())
            DEVABORT("dimensions of matrix and overlap do not match");

        std::vector<std::vector<unsigned int> > blocks=matrixBlocking(mats,1.e-12);

        unsigned int iVals=0; // block starting index
        // loop through blocks
        Val=Eigen::MatrixXcd(Mat.rows(),1);
        RVec=MSparse();
        DVec=MSparse();
        Eigen::VectorXi reserve(Mat.rows());

        size_t ib(0);
        std::vector<int> sizs;
        for(auto b: blocks){
            sizs.push_back(b.size());
            for(size_t i=0;i<b.size();i++,ib++)reserve(ib)=b.size();
        }


        PrintOutput::message(Sstr+"In sparse matrix, found "+blocks.size()+"independent subsets with size(s)"+sizs);

        if(RightVectors){RVec.resize(Mat.rows(),Mat.rows());RVec.reserve(reserve);}
        if( DualVectors){DVec.resize(Mat.rows(),Mat.rows());DVec.reserve(reserve);}

        for(auto block: blocks){
            if(block.size()>2000){
                if(block.size()>10000)ABORT("huge matrix in eigensolver");
                PrintOutput::warning(Sstr+"Large block in MatrixEigenBlock ",10,0,
                                     std::string("Possible fixes:")
                                     +"\n - try -eigenMethod=Arpack"
                                     +"\n - try reducing discretization size ");
            }

            // extract matrix and overlap
            Eigen::MatrixXcd bMat(block.size(),block.size());
            Eigen::MatrixXcd bVal(bMat.cols(),1),bRVec,bLVec;
            Eigen::MatrixXcd bOvr(0,0);
            if(Ovr.size())bOvr.resize(block.size(),block.size());

            for (unsigned int i=0;i<block.size();i++){
                for (unsigned int j=0;j<block.size();j++){
                    bMat(i,j)=coeff(Mat,block[i],block[j]);
                    if(bOvr.size())bOvr(i,j)=coeff(Ovr,block[i],block[j]);
                }
            }

            eigenFull(bMat,bVal,RightVectors,bRVec,DualVectors,bLVec,bOvr);


            // place block eigenvalues and std::vectors into global positions
            Val.block(iVals,0,bVal.rows(),1)=bVal;
            for(int j=0;j<bRVec.cols();j++)
                for(unsigned int i=0;i<bRVec.rows();i++){
                    if(RightVectors)RVec.insert(block[i],iVals+j)=bRVec(i,j);
                    if( DualVectors)DVec.insert(block[i],iVals+j)=bLVec(i,j);
                }
            iVals+=bVal.size();
        }
    }
    RVec.makeCompressed();
    DVec.makeCompressed();
}

}
#endif // MATRIXEIGENBLOCK_H
