#ifndef ROTATECOEFFICIENTS_H
#define ROTATECOEFFICIENTS_H

#include <memory>
#include <vector>
#include "qtEigenSparse.h"

#include "operatorTree.h"

#include "indexNew.h"

class OperatorRotate: public OperatorAbstract
{
    // recursive constructor for Index
    class IndexFull: public IndexNew {
    public:
        /// M=INT_MAX defines top of Idx, does NOT test Idx.isRoot()
        IndexFull(const Index* Idx, size_t Lmax, int M=INT_MAX);
        IndexNew* rotIndex() const; // extract rotational part of IndexFull
    };

    Eigen::SparseMatrix<std::complex<double>> _toSpecSparse,_fromSpecSparse;
    std::unique_ptr<OperatorTree> _toSpec,_fromSpec;
    std::unique_ptr<Coefficients> _cSpec;

    std::vector<size_t> _phaseK; // m-code vor _cSpec.data()
    mutable std::vector<std::complex<double>> _phases; // exp(+- i m Angle)

    // test is run for any new Idx appearing in apply
    static bool test(const Index* Idx); /// test sequence of rotations for C
public:
    static bool debug;
    /// rotates Coefficients of Jdx around oriented RotationAxis
    ///
    /// Jdx.hierarchy() must contain sequence Phi.Eta, i.e. Phi must be immediately followed by Eta
    ///<br> a new lhs Index (idx()) is created, which has complete spherical basis
    ///<br> orientation is given through the definition of <<AngularMomentumX>> in OperatorDefinition
    ///<br> rotation can be constrained to a Lmax (default: no all in basis)
    OperatorRotate(std::vector<double> RotationAxis, const Index* Jdx);

    /// set the rotation angle, use as Op.angle(0.123).apply(...)
    const OperatorRotate& angle(double Angle) const;
    void apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const override;

    static bool testCoriolis(const Index* Idx);

};


#endif // ROTATECOEFFICIENTS_H
