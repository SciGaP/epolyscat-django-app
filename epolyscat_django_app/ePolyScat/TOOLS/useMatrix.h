// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef USEMATRIX_H
#define USEMATRIX_H
#include "toolsHeader.h"
#include "qtEigenDense.h"
#include "qtEigenDense.h"
#ifdef _CBLAS_
#include "cblas.h"
#endif
//#include "arpack.h"

#define MapMatrix UseMatrix::UseMap
//#define MapVector UseMatrix::UseMap

#ifdef _WIN32
// for alternative tokens for windows ( e.g. and, or, not, ... ) 
#include <iso646.h>
#endif
class TriFactor;
class OperatorFloorInverse;
class BandedOverlap;

/** \brief
 \ingroup Linalg
 strongly limited matrix class for tRecX solver

 <br>treats double and complex<double> on equal footing
 <br>serves as an interface to linear algebra packages (mostly LAPACK)
 <br> matrices come in two categories
 <br> basic: full_normal, band_normal, band_symmetric (may be extended, e.g. sparse_normal)
 <br> derived: block,transpose
 <br> basic matrices own their data, derived point to the data of the basic type
 <br> derived matrices cannot be constructed, but are returned only by functions
 <br> upon assignment or copy, the derived types are converted to basic type,
 <br> i.e. operations like transposition, conjugation, blocking, and inversion performed explicitly
*/
class UseMatrix{
    friend class TriFactor;
    friend class TriFactorLU;
    friend class InverseIter;
    friend class OperatorFloorInverse;
    friend class BandedOverlap;
public:

    /// class for un-typed element-wise reference
    /// all element-wise operations are slow

    /// @cond DEV
    /// element arithmetic always goes through complex
    class Data{
        friend class UseMatrix;
        friend class TriFactor;
        friend class TriFactorLU;
    public:
        // only one of r,c can differ from zero, r==0 and c==0 is empty data (should always be handled)
        double * r,*rZero;
        std::complex<double> * c,*cZero;
        bool conjugated;
        Data(const UseMatrix * M, int I, int J);
        const UseMatrix * m;
    private:
        void changeData(); // check before changing a data item
    public:
        void operator=(double D);
        void operator=(Data Dat);
        bool operator==(const Data & Dat) const;
        inline bool operator!=(const Data & Dat) const {return(not operator==(Dat) );}

        /// key routine, used by all others: return value in complex variable
        inline std::complex<double> complex() const {if (r!=0) return *r;if(conjugated) return std::conj(*c); return *c; }
        inline double real() const {if(r!=0)return *r;return std::real(*c);}
        inline double imag() const {if(r!=0)return  0;if(conjugated) return -std::imag(*c); return std::imag(*c);}

        inline operator std::complex<double>(){return complex();}

        void operator =(std::complex<double> C);
        void operator+=(std::complex<double> C);
        void operator-=(std::complex<double> C);
        void operator*=(std::complex<double> C);
        void operator/=(std::complex<double> C);

        inline std::complex<double> operator-() const {return -complex();}

        inline std::complex<double> operator+(std::complex<double> B) const {return complex()+B;}
        inline std::complex<double> operator-(std::complex<double> B) const {return complex()-B;}
        inline std::complex<double> operator*(std::complex<double> B) const {return complex()*B;}
        inline std::complex<double> operator/(std::complex<double> B) const {return complex()/B;}

        inline std::complex<double> operator+(Data const & B) const {return complex()+B.complex();}
        inline std::complex<double> operator-(Data const & B) const {return complex()-B.complex();}
        inline std::complex<double> operator*(Data const & B) const {return complex()*B.complex();}
        inline std::complex<double> operator/(Data const & B) const {return complex()/B.complex();}

        inline bool operator==(std::complex<double> B){return complex()==B;}
        inline bool operator!=(std::complex<double> B){return complex()!=B;}
    };

    /// matrix kinds
    enum symmetry{unknown,unsymm,symmetric,hermitian};
    class Symmetry { public: symmetry symm; Symmetry(symmetry S):symm(S){}};

    /// @endcond

    enum EigenMethod {automatic,general,symmetric_full,symmetric_banded,hermitian_full,hermitian_banded,real_symmetric_full,real_symmetric_banded,arpack};
    static EigenMethod eigenMethod;

    /// frequently used matrices
    static UseMatrix Zero(unsigned int nrows, unsigned int ncols);      ///< zero matrix (for now, full)
    static UseMatrix Identity(unsigned int nrows, unsigned int ncols);  ///< identity matrix (for now, full)
    static UseMatrix Constant(unsigned int nrows, unsigned int ncols,const std::complex<double> C); ///< full constant matrix
    static UseMatrix Random(unsigned int nrows, unsigned int ncols); ///< random full matrix, (real, imag) in ([0,1],[0,1])
    static UseMatrix RandomReal(unsigned int nrows, unsigned int ncols); ///< random real full matrix, values in [0,1]
//    template<class T>
//    static UseMatrix FromVector(const T &vector);
    static UseMatrix FromVector(const std::vector<std::complex<double> >& vector);
    static UseMatrix FromVector(const std::vector<double>& vector);

    // NOTE: returns "derived" matrix
    static UseMatrix UseMap(std::complex<double> * Data, int Nrows, int Ncols); ///< construct with data (emulates Eigen::Map)
    static UseMatrix UseMap(std::complex<double> * Data, int Nrows);            ///< construct with data (emulates Eigen::Map)

    ~UseMatrix();
    UseMatrix();
    UseMatrix(const UseMatrix & other);

    UseMatrix(const unsigned int Rows, const unsigned int Cols); ///< generates a matrix
    UseMatrix(const Eigen::MatrixXcd & mat); ///< for legacy compatability
    /// generate basic banded matrix
    UseMatrix(const unsigned int Rows, const unsigned int Cols, const unsigned int SubDiag, const unsigned int SuperDiag, symmetry Symm=unknown);
    UseMatrix & assignMatrix(const UseMatrix & other, symmetry Sym=unknown, bool ForceFull=false); ///< A=B: overwrite A with a copy of B, converted to basic
    UseMatrix & operator=(const UseMatrix & other){return assignMatrix(other,unknown,false);} ///< A=B: overwrite A with a copy of B, converted to basic
    void swap(UseMatrix & B); ///< swap present with B

    // element-wise access is slow
    UseMatrix::Data operator()(unsigned int I,unsigned int J) const;  ///< access to element (SLOW!)
    UseMatrix::Data operator()(unsigned int I) const;                 ///< element from column/row vector (SLOW!)

    std::complex<double> *data() const {return cdata;}     ///< pointer to (complex) data
    std::complex<double> *data(unsigned int i, unsigned int j) const {return cdata+shape.locIJ(i,j);}   ///< pointer to (complex) data
    unsigned int rows() const {return shape.nrows;} ///< actual number of rows, taking into account transposition state
    unsigned int cols() const {return shape.ncols;} ///< actual number of cols, taking into account transposition state
    unsigned int size() const {return shape.nrows*shape.ncols;} ///< rows()*cols()
    unsigned int subD() const {return shape.lowBw-1;} ///< number of sub-diagonals (= rows()-1 for full matrix)
    unsigned int superD() const {return shape.upBw-1;} ///< number of super-diagonals (=cols()-1 for full matrix)
    unsigned int leadDim() const {return shape.leadDim;} ///< leading dimension of column-wise matrix storag

    // views of the matrix
    UseMatrix transpose() const; ///< return transpose view on the data (without changing the data)
    UseMatrix adjoint() const;   ///< return adjoint view on the data (without changing the data)
    UseMatrix block(unsigned int I, unsigned int J, unsigned int M, unsigned int N) const;///< M x N sub-block with upper left corner at I,J
    UseMatrix band(unsigned int SubD, unsigned int SuperD) const;///< subD/superD band of matrix
    UseMatrix col(unsigned int I) const; ///< I'th column
    UseMatrix row(unsigned int I) const; ///< I'th row
    UseMatrix leftCols(unsigned int N); ///< return view on the N leftmost columns
    UseMatrix rightCols(unsigned int N);///< return view on the N rightmost columns
    UseMatrix asDiagonal() const; ///< interprete single row or column matrix as diatonal matrix

    std::vector<std::complex<double> > extractSubmatrix(std::vector<unsigned int> Rows, std::vector<unsigned int> Cols) const;
    std::vector<std::complex<double> > extractDiagonal(std::vector<unsigned int> Rows, std::vector<unsigned int> Cols) const;

    // linear system solving
    UseMatrix inverse() const;                  ///< return inverse matrix
    UseMatrix & solve(UseMatrix & Rhs) const;   ///< replace Rhs with M^-1 Rhs

    std::vector<std::complex<double> > extractDiagonal(); ///< vector of diagonal values
    std::vector<std::complex<double> > extractSubmatrix(const std::vector<unsigned int> &Rows,const std::vector<unsigned int> &Cols); ///< column major storage of submatrix

    void eigenSetNext(EigenMethod Method);

    void eigenValues(UseMatrix & Val, const UseMatrix & Ovr=UseMatrix()) const; ///< eigenvalues and -vectors of generalized eigenproblem
    void      eigen(UseMatrix &Val,UseMatrix &Vec, const UseMatrix &Ovr=UseMatrix(), bool Vectors=true) const; ///< eigenvalues and -vectors of generalized eigenproblem
    void eigenBlock(UseMatrix &Val,UseMatrix &Vec, UseMatrix &Ovr, bool Vectors); ///<eigenvalues and -vectors with automatic blocking detection

    /// \brief eigenvalues and -vectors with automatic blocking detection, include dual eigen vectors
    ///
    /// dual: dual_i (Mat - e_i Ovr) = 0, i.e. dual_i.innerProduct(eigen_j)=delta_ij * norm
    void      eigen(UseMatrix & Val, const UseMatrix &Ovr, bool RightVectors, UseMatrix &RVec, bool DualVectors, UseMatrix & DVec) const;
    ///\brief eigenvalues and -vectors with automatic blocking detection, include dual eigen vectors
    void eigenBlock(UseMatrix &Val, UseMatrix &Ovr, bool RightVectors, UseMatrix &RVec,bool DualVectors, UseMatrix & DVec);
    /// @brief (pseudo-)orthonormalize degenerate sub-spaces
    bool eigenOrthonormalize(UseMatrix & Val, UseMatrix & Vec, const UseMatrix &Ovr, double eps=1.e-8) const;
    /// @brief Gram-Schmidt (pseudo-)orthonormalize subset of vectors
    void gramSchmidt(const UseMatrix & Ovr,std::vector<unsigned int> Subset, bool Pseudo);
    /// find unconnected sets of indices
    std::vector<std::vector<unsigned int> > blocking(const std::vector<const UseMatrix *> &Mat=std::vector<const UseMatrix*>(0), /**< additional connection matrices */
            double Eps=0. /**< row-wise zero threshold */) const;

    /// list of all indices that a given index connects to
    void connect(unsigned int From, /**< starting index */
                 std::vector<unsigned int> & Conn, /**< append all connected indices to this vector */
                 std::vector<bool> &Use, /**< subset of all indices to be used (Use.size()==0 initializes to all true */
                 std::vector<const UseMatrix*> Other, /**< check also these matrices for connection */
                 double Eps=0. /**< threshold for considering a matrix element =0 */) const;

    /// compress data to banded and real where possible
    UseMatrix & compress(double Eps=0.   /**< [in] tolerance for conjugation, zero etc. */,
                  bool Symm=true  /**< [in] true: exploit symmetry if possible*/);
    UseMatrix & expand();   ///< expand to full storage
    UseMatrix expandConst() const;   ///< expand to full storage
    UseMatrix & purge(double Eps=1.e-12, double EpsAbs=1.e-14); ///< set near-zeros to exactly =0

    // matrix multiplications
    UseMatrix & operator*=(const UseMatrix & B);      ///< in-place matrix multiplication
    UseMatrix   operator* (const UseMatrix & B) const;///< out-of-place matrix multiplication
    Eigen::VectorXcd operator*(const Eigen::VectorXcd & V); ///< interfacing to Eigen

    /// generic matrix multiplication C=c * C + a * this * B
    void multiply(const UseMatrix & B, UseMatrix & C,
                  const std::complex<double> & a, const std::complex<double> & c) const;

    /// element-wise operations, e.g. return this + B;
    UseMatrix     operator+(const UseMatrix & B) const {UseMatrix C;ternaryOp(Add,     B,C,false);return C;}
    UseMatrix     operator-(const UseMatrix & B) const {UseMatrix C;ternaryOp(Subtract,B,C,false);return C;}
    UseMatrix  cwiseProduct(const UseMatrix & B) const {UseMatrix C;ternaryOp(Product, B,C,false);return C;}
    UseMatrix cwiseQuotient(const UseMatrix & B) const {UseMatrix C;ternaryOp(Quotient,B,C,false);return C;}

    /// in-place element-wise matrix-matrix operations, e.g. A+=B;
    UseMatrix &      operator+=(const UseMatrix & B){binaryOp(BinaryInAdd,     B); return *this;}
    UseMatrix &      operator-=(const UseMatrix & B){binaryOp(BinaryInSubtract,B); return *this;}
    UseMatrix &  cwiseProductIn(const UseMatrix & B){binaryOp(BinaryInProduct, B); return *this;}
    UseMatrix & cwiseQuotientIn(const UseMatrix & B){binaryOp(BinaryInQuotient,B); return *this;}

    /// in-place matrix-scalar operations, e.g. A+=a;
    UseMatrix &  operator=(std::complex<double> a){scalarOp(ScalarInAssign, a); return *this;}
    UseMatrix & operator+=(std::complex<double> a){scalarOp(ScalarInAdd,    a); return *this;}
    UseMatrix & operator-=(std::complex<double> a){scalarOp(ScalarInAdd,   -a); return *this;}
    UseMatrix & operator*=(std::complex<double> a){scalarOp(ScalarInProduct,std::complex<double>(a)); return *this;}
    UseMatrix & operator/=(std::complex<double> a){scalarOp(ScalarInProduct,std::complex<double>(1./a)); return *this;}

    /// out-of-place matrix-scalar operations, e.g. return this+a or this*a
    UseMatrix operator+(std::complex<double> a){UseMatrix C;ternaryOp(ScalarOutAdd     ,*this,C,false,1.,1.,std::complex<double>( a)); return C;}
    UseMatrix operator+(             double  a){UseMatrix C;ternaryOp(ScalarOutAdd     ,*this,C,false,1.,1.,std::complex<double>( a)); return C;}
    UseMatrix operator-(std::complex<double> a){UseMatrix C;ternaryOp(ScalarOutAdd     ,*this,C,false,1.,1.,std::complex<double>(-a)); return C;}
    template<class T> UseMatrix operator/(T a){UseMatrix C;ternaryOp(ScalarOutProduct,*this,C,false,1.,1.,std::complex<double>(1./a)); return C;}
    UseMatrix operator*(std::complex<double> a){UseMatrix C;ternaryOp(ScalarOutProduct,*this,C,false,1.,1.,std::complex<double>(a)); return C;}
    UseMatrix operator*(             double  a){UseMatrix C;ternaryOp(ScalarOutProduct,*this,C,false,1.,1.,std::complex<double>(a)); return C;}


    /// dot product of two column vectors (=single-row matrices)
    std::complex<double> dot(const UseMatrix & B) const; // matrix.transpose * matrix

    /// queries
    bool isIdentity(double eps=0.) const;
    bool isDiagonal(double eps=0.) const;
    bool isZero(double eps=0.) const;
    bool operator==(const UseMatrix &rhs) const;
    bool compare(const UseMatrix &rhs, double Eps) const;
    inline bool operator!=(const UseMatrix& rhs) const {return !this->operator==(rhs);}

    unsigned int nonZeros(double eps=0.) const; ///< count of entries > eps
    std::vector<unsigned int> locNonZero(double Eps=0.) const; ///< first index pair where matrix-element=Val

    /// reduce operations
    double maxAbsVal() const; //!< max_{(i,j)} |A(i,j)|

    /// printing, testing, etc.
    void show(std::string text="row\\col  ") const; //!< show the structure of a matrix
    std::string str(std::string text="row\\col  ", unsigned int Digits=3) const; ///< printable string, Digits==0 prints structure
    void print(std::string text="row\\col  ", unsigned int Digits=3) const; ///< print, Digits==0 prints structure
    void printShape(std::string text) const {shape.print(text);} ///< print matrix shape and storage info
    std::string strShape() const {return shape.str();} ///< print matrix shape and storage info
    static void Usage(); ///< a few usage examples
    static void Test();  ///< various tests (mostly for development)


protected:

//    class Arp:public Arpack{
//        friend class UseMatrix;
//        virtual void apply(const UseMatrix X, UseMatrix Y);
//        void eigen(UseMatrix &Eval, UseMatrix &Rvec, unsigned int Nvec, const std::string &Which, bool Restart);
//        Arp(const UseMatrix & O);
//        // internal data
//        const UseMatrix* o;
//    };

/// @cond DEV
    class Shape;
public:
    typedef unsigned int(*Location)(const unsigned int I,const unsigned int J, const Shape & shape);
protected:
    class Shape{
        /** complete info about matrix shape
         */
        friend class TriFactor;
        friend class TriFactorLU;
        friend class UseMatrix;
    public:

    private:
        unsigned int nrows,ncols;      // rows and columns
        Location locationIJ;           // compute location in data storage
        unsigned int leadDim,total;    // leading dimensions and total size of storage (total==0: emtpy matrix)
        unsigned int lowBw,upBw;       // lower and upper band width, i.e. 1 + # of sub/super diagonals
        unsigned int row0,col0;        // starting indices of sub-block
        const Shape* stored;           // all non-basic matrices point to their parent's shape

        Shape() : nrows(0),ncols(0),locationIJ(full_normal),
            leadDim(1),total(0),
            lowBw(0),upBw(0),
            row0(0),col0(0),
            stored(0){}

        /// standard full storage shape
        Shape(unsigned int Nrows, unsigned int Ncols) :
            nrows(Nrows),ncols(Ncols),locationIJ(full_normal),
            leadDim(nrows),total(leadDim*Ncols+1),// NOTE: use last data point to store 0.
            lowBw(Nrows),upBw(Ncols),
            row0(0),col0(0),
            stored(0){}

        /// standard banded storage shape
        Shape(unsigned int Nrows, unsigned int Ncols,unsigned int SubDiag,unsigned int SuperDiag,symmetry Symm=unknown) :
            nrows(Nrows),ncols(Ncols),locationIJ(band_normal),
            leadDim(SubDiag+SuperDiag+1),total(leadDim*Ncols+1), // NOTE: use last data point to store 0.
            lowBw(SubDiag+1),upBw(SuperDiag+1),row0(0),col0(0),stored(0)
        {
            if(lowBw==nrows and upBw==ncols){locationIJ=full_normal;leadDim=nrows;total=leadDim*ncols+1;};
            if(Symm==symmetric or Symm==hermitian){locationIJ=band_symmetric;leadDim=SubDiag+1;total=leadDim*Ncols+1;}
        }

        /// complete contructor (assigning all members)
        Shape(unsigned int nrows, unsigned int ncols,
              Location locationIJ,
              unsigned int leadDim,unsigned int total,
              unsigned int lowBw,unsigned int upBw,
              unsigned int row0,unsigned int col0,
              const Shape* stored) :
            nrows(nrows),ncols(ncols),
            locationIJ(locationIJ),
            leadDim(leadDim),total(total),
            lowBw(lowBw),upBw(upBw),
            row0(row0),col0(col0),
            stored(stored){}

    public:
        bool storeBasic() const; ///< criterion for basic shape
        bool storeContiguous() const; ///< true for matrices with contiguous data storage
        bool storeFull() const;  ///< true when full matrix storage (not banded or sparse)
        bool storeBand() const; ///< true when banded matrix storage (not full or sparse)
        bool storeTranspose() const;  ///< true when transposed of matrix is stored
        bool storeSymmetric() const;  ///< true when symmetric storage (i.e. i <= j)
        bool operator==(const Shape & other) const;
        void print(std::string text="") const;
        std::string str() const;
        unsigned int locIJ(unsigned int I,unsigned int J) const {return locationIJ(I,J,*this);}
        std::string strLocIJ() const;

        inline unsigned int begCol(unsigned int J) const {
            if(stored==0)return locIJ(std::max(0,int(J)-int(upBw)+1),J);
            return stored->locIJ(row0,col0+J);}
        inline unsigned int endCol(unsigned int J) const {
            if(stored==0)return locIJ(std::min(nrows,J+lowBw)-1   ,J)+1;
        return stored->locIJ(row0+std::min(nrows,J+lowBw)-1,col0+J)+1;} // AFTER last element in column
    };
    /// @endcond

public:
    // various forms of accessing the data
    static unsigned int full_normal(unsigned int I,unsigned int J, const Shape & shape);
    static unsigned int band_normal(unsigned int I,unsigned int J, const Shape & shape);
    static unsigned int band_symmetric(unsigned int I,unsigned int J, const Shape & shape);
    // views on the previous
    static unsigned int full_normal_map(unsigned int I,unsigned int J, const Shape & shape);
    static unsigned int asDiagonal_map (unsigned int I,unsigned int J, const Shape & shape);
    static unsigned int any_transpose(unsigned int I,unsigned int J, const Shape & shape);
    static unsigned int any_block(unsigned int I,unsigned int J, const Shape & shape);
    static unsigned int any_band (unsigned int I,unsigned int J, const Shape & shape);
protected:
    unsigned int storeIncI() const; ///< increment in storage corresponding to I++

    bool isShape(Location loc) const {return shape.locationIJ==loc;} ///< true if given shape
public:
    bool isFull() const;  ///< true if all elements of each column are stored
    bool isBand() const;  ///< true if banded storage
    bool isBandNormal() const;  ///< true if band_normal storage
protected:
    bool isBasic() const; ///< criterion for basic matrix
    bool sharesData(const UseMatrix & other) ///< true if two matrices share their data storage
    const {return (cdata!=0 and cdata==other.cdata) or  (rdata!=0 and rdata==other.rdata);}

    // copy matrix data into standard shape
    void copyData(const UseMatrix & other); ///< copy data only (double or complex)
    void freeData();                        ///< free data if owner
    void allocateData(bool Complex);        ///< allocate data storage for standard shape
    template<class D>void cornerData(const D & data) const; ///< clean the corners of banded data storage

    template<class D> inline D begCol(D data,unsigned int J) const {return data+shape.begCol(J);} ///< pointer to first data of column
    template<class D> inline D endCol(D data,unsigned int J) const {return data+shape.endCol(J);} ///< pointer to 1 beyond last of column

public:

    // parameters for all matrix operations
    enum operation {assign,add,subtract,quotient,product};
    enum     terms {unary,scalar,binary,ternary}; ///< connect one, two, or three matrices or matrix and a scalar
    /// @cond DEV
    struct Operation{
        const operation kind; // add, multiply, etc.
        const terms term;     // which and how many terms are involved
        const bool inPlace;   // output storage in the same place as input?
        Operation(operation Kind,terms Term, bool InPlace):kind(Kind),term(Term),inPlace(InPlace){}
    };
    /// @endcond

    bool isHermitian(double Eps=1e-12, double Threshold=1e-12) const;
    bool isSymmetric( double Eps=1.e-12) const;
    bool isReal( double Eps=1.e-12) const;
    std::string type(double BandRatio=1.) const; ///< Lapack style strings to characterize matrices

protected:
    ///< in place scalar operations, e.g. A*=a;
    static const Operation ScalarInAssign;
    static const Operation ScalarInAdd;
    static const Operation ScalarInProduct;

    ///< out-of-place scalar operations, e.g. A=B*a;
    static const Operation ScalarOutAdd;
    static const Operation ScalarOutProduct;

    ///< in place element-wise binary operations, e.g. A+=B;
    static const Operation BinaryInAssign;
    static const Operation BinaryInAdd;
    static const Operation BinaryInSubtract;
    static const Operation BinaryInProduct;
    static const Operation BinaryInQuotient;

    ///< out-of-place ternary operations, e.g. C=A+B
    static const Operation Add;
    static const Operation Subtract;
    static const Operation Product;
    static const Operation Quotient;

    void  scalarOp(Operation Op, std::complex<double> a=1.);
    void  binaryOp(Operation Op, const UseMatrix &B, std::complex<double> a=1.,std::complex<double> b=1.);
    void ternaryOp(Operation Op,const UseMatrix & B, UseMatrix & C, bool keepC, std::complex<double> a=1., std::complex<double> b=1., std::complex<double> c=1.) const;

    ///< perform operations involving 3 matrices
    template<class Ad,class Bd,class Cd,class As,class Bs,class Cs>
    void ternaryData(Operation Op,        const As a, const Ad & Adata,
                     const UseMatrix & B, const Bs b, const Bd & Bdata,
                     const UseMatrix & C, const Cs c, const Cd & Cdata) const;

    ///< apply operation column-ranges
    template<class Ad,class Bd,class Cd,class As,class Bs,class Cs>
    void rangeOp(operation Kind, terms Term, bool InPlace, bool aconjg, bool bconjg,
                 Ad Abeg,Ad Aend,unsigned int Ainc,
                 Bd Bbeg,Bd Bend,unsigned int Binc,
                 Cd Cbeg,Cd Cend,unsigned int Cinc,
                 As as, Bs bs, Cs cs) const;

    // auxiliary routines

    /// conjugate a stretch of data
    inline void conjg(std::complex<double>*&a,std::complex<double>*&end, int & inc, std::complex<double>*& temp) const {
        temp=new std::complex<double>[end-a];
        std::complex<double>* c=temp;
        for(;a<end;a+=inc,c++)*c=std::conj(*a);
        a=temp;
        end=c;
        inc=1;
    }
    void conjg(double* &a,double *&end, int & inc, double*& temp) const {}// dummy for overloading only

    /// change the band width of a matrix without changing its entry values
    UseMatrix  reband(unsigned int SubDiagonal, unsigned int SuperDiagonal, bool Crop=false) const;
    void triFactor() const; // (re)-calculate the triangular factorization

    // structure analysise
    bool zeroData(const unsigned int Loc, double Eps=0.) const; ///< abs(data+Loc)<Eps
    bool realData(const unsigned int Loc, double Eps=0.) const; ///< abs(imag(data+Loc))<Eps
    bool equalData(const unsigned int LocA, const unsigned int LocB, double Eps=0.) const; ///< abs(dataA-dataB)<=Eps
    bool conjgData(const unsigned int LocA, const unsigned int LocB, double Eps=0.) const; ///< abs(dataA-conj(dataB)com)<Eps

public:
    double maxNorm(std::vector<double> &RowMax, std::vector<double> &ColMax) const; ///< maximum norm by row, column
    double maxRealImag(std::vector<double> &RowMax, std::vector<double> &ColMax) const; ///< max(|real|,|imag|) by row, column
    bool polar(std::vector<std::complex<double> > &Phase,std::vector<double>& RealMat) const; ///< true if Matrix = Lphas*RealMat
    void diagnose(double Eps, unsigned int &NonZero, unsigned int &TrueSub, unsigned int &TrueSuper, char &DataType) const;

protected:
    unsigned int trueSubD(double Eps=0.) const;     ///< number of non-zero sub-diagonals.
    unsigned int trueSuperD(double Eps=0.) const;   ///< number of non-zero super-diagonals.

    void storeSymmetricBand(double Eps);
    void storeBand(double Eps);

public:
#ifdef _CBLAS_
    CBLAS_TRANSPOSE cblasTrans() const; ///< translate normal,transpose,hermitian into corresponding CBLAS values
#endif
    char lapackTrans() const; ///< return Lapack transposition 'n','t','h'
    //===============================================================================================
protected:
    Shape shape;                ///< all shape information
    // only one of rdata and cdata can differ from zero, both==0 indicates empty matrix (must be handled correctly)
    double * rdata;               ///< data for real matrix
    std::complex<double> * cdata; ///< data for complex matrix
    bool conjugated;              ///< interprete (complex) data as conjugated
    Symmetry *symm;               ///< store information about symmetry of the matrix (pointer construction for using in const functions)
public:
    TriFactor * trf;              ///< factorization into triangular matrices
};

inline std::complex<double> operator+(std::complex<double> A, const UseMatrix::Data B){return A+B.complex();}
inline std::complex<double> operator-(std::complex<double> A, const UseMatrix::Data B){return A-B.complex();}
inline std::complex<double> operator*(std::complex<double> A, const UseMatrix::Data B){return A*B.complex();}
inline std::complex<double> operator/(std::complex<double> A, const UseMatrix::Data B){return A/B.complex();}

inline void operator+=(std::complex<double> A, const UseMatrix::Data B){A+=B.complex();}
inline void operator-=(std::complex<double> A, const UseMatrix::Data B){A-=B.complex();}
inline void operator*=(std::complex<double> A, const UseMatrix::Data B){A*=B.complex();}
inline void operator/=(std::complex<double> A, const UseMatrix::Data B){A/=B.complex();}

inline double abs(const UseMatrix::Data &B){return std::abs(B.complex());}

#endif // USEMATRIX_H
