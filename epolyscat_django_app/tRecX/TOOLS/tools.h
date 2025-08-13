// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef __TOOLS__
#define __TOOLS__

#include "abort.h"
#include "stringTools.h"
#include "units.h"

class UseMatrix;


#ifdef _WIN32
// for alternative tokens for windows ( e.g. and, or, not, ... ) 
#include <iso646.h>
#endif


class Discretization;

namespace alglib {class real_1d_array;}



/** @defgroup Tools Tools
 * \brief code general use: string manipulation, I/O, physical constants and units, many more...
 */

namespace tools {

/** @defgroup Misc Miscellaneous
 *  @ingroup Tools
 *  \brief Mixed bag of utilities
 *  @{ */

//void lapack_zggev(UseMatrix &A, UseMatrix &B, std::vector<std::complex<double> > & Eval, UseMatrix & vr, bool Vectors=true);
std::vector<double> insertInterval(double A, double B, vector<double> Split);


/// create a new file name by appending numbers, if file exists
std::string newFile(const std::string File);

/// convenient interface to sleep
void sleepSeconds(double Secs);

enum compare {equal,smaller,larger,differs};
/// true if any elements of the vector fulfills comparison (equal,smaller,larger)
template<typename T>
bool anyElement(std::vector<T> v,compare comp,T s){
    switch (comp){
    case equal:   for (unsigned int i=0;i<v.size();i++) if(v[i]==s)return true; return false;
    case differs: for (unsigned int i=0;i<v.size();i++) if(v[i]!=s)return true; return false;
    case smaller: for (unsigned int i=0;i<v.size();i++) if(v[i] <s)return true; return false;
    case larger:  for (unsigned int i=0;i<v.size();i++) if(v[i] >s)return true; return false;
    default: abort();
    }
}
template<typename T>
bool anyElementAbs(std::vector<T> v,compare comp,const double s){
    switch (comp){
    case equal:   for (unsigned int i=0;i<v.size();i++) if(abs(v[i])==s)return true; return false;
    case differs: for (unsigned int i=0;i<v.size();i++) if(abs(v[i])!=s)return true; return false;
    case smaller: for (unsigned int i=0;i<v.size();i++) if(abs(v[i]) <s)return true; return false;
    case larger:  for (unsigned int i=0;i<v.size();i++) if(abs(v[i]) >s)return true; return false;
    default: abort();
    }
}
/// first occurance of element in vector (use stride Stride and start from From)
/// if not found, returns vector size
template<typename T>
unsigned int locateElement(std::vector<T> v,T s, unsigned int Stride=1, unsigned int From=0){
    for (unsigned int i=From;i<v.size();i+=Stride) if(v[i]==s)return i;
    return v.size();
}

/// use the following routinely for float comparisons
inline bool doubleLess(double x,double y,double eps){return x<y-abs(eps);}
inline bool doubleMore(double x,double y,double eps){return x>y+abs(eps);}
inline bool doubleSame(double x,double y,double eps){return (not doubleMore(x,y,eps)) and (not doubleLess(x,y,eps));}
inline bool doubleBelow(double x,double a,double b){if(a==-DBL_MAX)return false; return x-a< -std::min((b-a)*1.e-14,1.);}
inline bool doubleAbove(double x,double a,double b){if(b== DBL_MAX)return false; return x-b>  std::min((b-a)*1.e-14,1.);}
inline bool doubleInside(double x, double a, double b){return (not doubleBelow(x,a,b)) and (not doubleAbove(x,a,b));}

/// write/read contiguous data
template<typename T> void write(std::ostream & stream, const T c){stream.write(reinterpret_cast<const char *>(&c),sizeof(c));}
template<typename T> void write(std::ostream & stream, const T*c, size_t CountT){for(unsigned k=0;k<CountT;k++)stream.write(reinterpret_cast<const char *>(c+k),sizeof(*c));}
template<typename T> void write(std::ostream & stream, std::vector<T> & v){if(v.size()>0)stream.write(reinterpret_cast<char * >(v.data()),v.size()*sizeof(v[0]));}

template<typename T> void read(std::istream & stream, T & c){stream.read(reinterpret_cast<char *>(&c),sizeof(c));}
template<typename T> void read(std::istream & stream, T*c, size_t CountT){for(unsigned k=0;k<CountT;k++)stream.read (reinterpret_cast<char *>(c+k),sizeof(*c));}
template<typename T> void read(std::istream & stream, std::vector<T> & v){if(v.size()>0)stream.read( reinterpret_cast<char * >(v.data()),v.size()*sizeof(v[0]));}


typedef bool(*CompareComplex)(const complex<double> & A,const complex<double> & B);
bool lessReal(const complex<double> & A, const complex<double> & B);
bool lessImag(const complex<double> & A, const complex<double> & B);
bool lessAbs( const complex<double> & A, const complex<double> & B);

/// auxiliary class for sorting of one array by criteria on a different array
template <class T>
struct ComPair{
    complex<double> key;
    T val;
    CompareComplex comp;
    ComPair(const complex<double> Key,const T Val, const CompareComplex Comp):key(Key),val(Val),comp(Comp){}
    bool operator<(const ComPair & other) const{return comp(key,other.key);}
};

/// @brief sort complex Key's and, optionally, an arbitrary type vector of equal length by various criteria
template<typename T>
void sortKey(const std::string Comp,/**< SmallReal,SmallImag,SmallAbs */
             vector<complex<double> > & Key, vector<T> & Val=vector<T>()){
    if(Comp=="NoSort")return;
    CompareComplex comp;
    if     (Comp=="SmallReal")comp=tools::lessReal;
    else if(Comp=="SmallAbs")comp=tools::lessAbs;
    else if(Comp=="SmallImag")comp=tools::lessImag;
    else ABORT("no such sorting criterion implemented: "+Comp);

    // Keys only
    if(Val.size()==0)
        sort(Key.begin(),Key.end(),comp);

    // Keys and values
    else
    {   if(Val.size()!=Key.size())ABORT("Val.size() and Key.size() do not match");
        vector<ComPair<T> > v;
        for(unsigned int k=0;k<Key.size();k++)v.push_back(ComPair<T>(Key[k],Val[k],comp));
        std::sort(v.begin(),v.end());
        for(unsigned int k=0;k<Key.size();k++){
            Key[k]=v[k].key;
            Val[k]=v[k].val;
        }
    }
}

/// @brief auxiliary type for sortByKey
template<typename K, typename V>
struct KeyVal{
    K key;
    V val;
    bool operator<(const KeyVal & B) const {return key<B.key;}
    KeyVal(const K & Key, const V & Val):key(Key),val(Val){}
};


/// @brief sort key/value data for keys in ascending order
template<typename K, typename V>
void sortByKey(vector<K> & Key, vector<V> & Val){

    if(Val.size()!=Key.size())ABORT("Val.size() and Key.size() do not match");
    vector<KeyVal<K,V> > v;
    for(unsigned int k=0;k<Key.size();k++)
        v.push_back(KeyVal<K,V>(Key[k],Val[k]));

    std::stable_sort(v.begin(),v.end());

    for(unsigned int k=0;k<Key.size();k++){
        Key[k]=v[k].key;
        Val[k]=v[k].val;
    }
}

/// true if all Vec[k-1]-Vec[k] are equal (within tolerance)
bool isEquidistant(const std::vector<double> & Vec);


/// Cholesky style decomposition for a tri-diagnal complex symmetric (not hermitian!) matrix
void pseudoCholeskyTri(vector<complex<double> > & Diag, vector<complex<double> > & Super);

/// recursively generates a Gram-Schmidt transformation for a Metric matrix
void gramSchmidtTrans(const vector<vector<double> > & Metric, vector<vector<double > > & Trans);


/// character indicates magnitude of Val (., o, x, or X for increasing size)
char zero(const std::complex<double> & val);

/// returns {Num,Denom} for |Value-Num/Denom| < Eps*|Value|, Denom <= if no rational representation
///
/// only small prime-factors are used
vector<int> fractionSmallPrimes(double Value,
                                double Eps=1.e-12 /** accuracy of rational representation */,
                                double MaxDen=100 /** largest guaranteed denominator */);

//==== get quadrature rules ===============================
enum quadrature_kind { // to be extended
    gq_legendre,
    gq_hermite,
    gq_laguerre,
     q_equidistant
};

void get_quadrature_rule(int order, quadrature_kind kind, alglib::real_1d_array& points, alglib::real_1d_array& weights,
                         const std::vector<double>& endpoints=std::vector<double>(0) );

/// return P such that Perm[k]=Model[P[k]], empty if not permutation
template<typename T>
std::vector<unsigned int> permutation(const vector<T> & Model, const std::vector<T> & Perm)
{
    std::vector<unsigned int> p;
    std::vector<bool> avail(Model.size(),true);
    if(Model.size()!=Perm.size())return p;
    for(unsigned int k=0;k<Perm.size();k++)
        for(unsigned int l=0;l<Model.size();l++){
            if(avail[l] and Model[l]==Perm[k]){
                avail[l]=false;
                p.push_back(l);
                break;
            }
        }
    if(p.size()!=Model.size())p.clear();
    return p;
}
/** @} */

}



#endif
