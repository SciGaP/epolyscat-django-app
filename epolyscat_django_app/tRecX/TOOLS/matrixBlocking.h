#ifndef MATRIXBLOCKING_H
#define MATRIXBLOCKING_H

#include "abort.h"
#include "str.h"
#include <vector>


namespace tRecX{

template<class M> std::complex<double> coeff(const M& Mat, int I, int J){
    return std::complex<double>(Mat.coeff(I,J));
}

/// auxiliary routine for tRecX::matrixBlocking(...)
template<class M>
void connect(unsigned int From, std::vector<unsigned int> & Conns, std::vector<bool> & Use,
                        std::vector<const M*> Mat, double Eps) {
    if(Use.size()==0)Use.assign(Mat.front()->rows(),true);
    if(Use.size()!=Mat.front()->rows())ABORT("Use.size() must match rows()");

    // get all indices that connect From and that are not taken yet
    std::vector<unsigned int> newConns;
    for(unsigned int k=0;k<Mat.front()->rows();k++){
        if(Use[k]){

            // check whether k connects to present index
            bool connected=k==From;
            for(auto m: Mat){
                if(connected)break;
                connected=connected or
                        std::norm(coeff(*m,k,From))
                        >Eps*Eps*sqrt(std::norm(coeff(*m,k,k)*coeff(*m,From,From)));
            }

            if(connected){
                Use[k]=false;
                Conns.push_back(k);
                if(k!=From)newConns.push_back(k);
            }
        }
    }

    // add all connections from all new indices
    for (unsigned int k=0;k<newConns.size();k++)
        connect(newConns[k],Conns,Use,Mat,Eps);
}

/// find shared indepedent sub-blocks of a set of matrices
///
/// Matrix class must have functions rows(),cols(), element access by index through coeff<class M>(m,i,j) may need to be specialized
/// Element type must support std::norm(coeff(m,i,j))
template<class M>
std::vector<std::vector<unsigned int> > matrixBlocking(const std::vector<const M*> &Mat, double Eps){
    // check dimensions
    if(Mat.front()->rows()!=Mat.front()->cols())ABORT("implemented only for square matrices");
    for(auto m: Mat){
        if(Mat.front()->rows()!=m->rows())ABORT(Str("rows do not match"," ")+Mat.front()->rows()+m->rows());
        if(Mat.front()->cols()!=m->cols())ABORT(Str("columns do not match"," ")+Mat.front()->cols()+m->cols());
    }

    // list of block
    std::vector<std::vector<unsigned int> > blocks(0);

    // find connected subsets
    unsigned int next=0;
    std::vector<bool> Use(Mat.front()->rows(),true);
    while (next<Mat.front()->rows()) {
        blocks.push_back(std::vector<unsigned int>(1,next));
        Use[next]=false;
        connect(next,blocks.back(),Use,Mat,Eps);

        // find next unused
        for(;next<Use.size();next++)
            if(Use[next])
                break;
    }
    return blocks;
}
}





#endif // MATRIXBLOCKING_H
