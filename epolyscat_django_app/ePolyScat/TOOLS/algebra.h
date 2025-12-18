// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef ALGEBRA_H
#define ALGEBRA_H
#include "tools.h"
#include "stdio.h"
#include <iostream>
#include <fstream>
#include <map>
#include "printOutput.h"
#include "integrate.h"
class ReadInput;


/** @defgroup Functions
 * @ingroup Tools
 * \brief string-controled algebra, orthogonal polynomials, etc.
 * @{
 */

/// \brief Algebra - a tree of operations -, +, /, or *, as well as single-argument functions
///
/// a valid Algebra expression is composed of terms
/// <br> related by the standard operatins +,* and their inverses -,/
/// <br>  each term is a constant (1,3.1415,137, etc.), the linear function Q or
/// <br>  a set of functions sqrt(Q), exp(Q) etc. (incomplete list)
/// <br>  the arguments of functions can be general Algebra strings, e.g. sqrt(1-Q*Q)
/// <br>  round brackets (...) can be used to group terms as usual
/// <br>  see Test() for examples
/// <br>  an effort is made to detect malformed strings
class Algebra
{

    // suppress copy and asignement constructors
    Algebra(const Algebra&){}
    Algebra & operator=(const Algebra&){return *this;}

    static std::map<std::string,std::complex<double> > specialConstants;
    typedef const Algebra* (*algebraFactory)(std::string Term);
    static std::map<std::string, algebraFactory > externalFactory;
 public:
    static const std::string Infty;

    /// add external factory to list
    ///
    /// factories in list will be tried to for creating an Algebra, if standard definitions fail
    /// <br> Factory takes a string as an argument and returns an Algebra
    static void addExternalFactory(algebraFactory Factory, std::string Name="from FactoryName");
    static void addSpecialConstant(std::string Name, std::complex<double> Value);
    static void addUpdatableConstant(std::string Name, std::complex<double> Value);
protected:
    static std::map<std::string,std::complex<double> > updatableConstants;
    std::vector<const Algebra*> A,I; ///< direct and inverse contributions
    std::string argument(std::string Term) const; ///< get argument of function
    bool add; ///< if true, consider terms as to be added (or negatively added), else to be multiplied (or divided)
    void addFailure(std::string Message){A.push_back(0);failures+=" ... "+Message;} ///< push 0 to A and append message to Algebra::failures
    const Algebra* factory(const std::string Term,bool Add); ///< return pointers to single term Algebra's
    std::string _definition;
public:
    virtual ~Algebra();
    Algebra():add(true){}
    std::string definition() const {return _definition;};
    Algebra(const std::string Definition,const bool Add=true);

    /// evaluate for argument Q
    virtual std::complex<double> val(const std::complex<double> Q) const;

    /// evaluate algebra of constants
    std::complex<double> constVal() const{
        if(not isAlgebraOfConsts())ABORT("not an algebra of constants: "+definition());
        return val(0);
    }

    /// full precision numerical integral
    virtual std::complex<double> integral(const std::complex<double> Q0,const std::complex<double> Q1) const;

    void plot(std::string File,double From, double To, int Points) const;

    std::string str() const;

    bool isAlgebra() const; ///< true for correctly formed Algebra instance
    bool isAlgebraOfConsts() const; ///< algebraic expression of all constants
    bool isAlgebraOfUpdatables() const; ///< algebraic expression of all constants and updatable constants

    static bool isAlgebra(std::string Definition); ///< check string for well-formed algebra
    static std::string failures; ///< list of places and reasons of malformation (if any)

    virtual std::vector<std::complex<double>> nonAnalyticQ() const; ///< function is non-analytic at these Q's
    /// locate sign change of algebra in [Low,Up] to within Epsilon (abort if there is no sign change)
    double signChangeQ(double Low, double Up, double ValueLow, double ValueUp, double Epsilon) const;
    /// zeros in [Low,Up], resolved to Delta, located with Epsilon
    virtual std::vector<double> zeros(double Low, double Up, double Delta, double Epsilon=DBL_EPSILON) const;
    /// poles in [Low,Up], resolved to Delta, located with Epsilon
    virtual std::vector<double> poles(double Low, double Up, double Delta, double Epsilon=DBL_EPSILON) const;

    /// check whether Term is a numerical constant
    static bool isConstant(std::string Term);

    /// check whether Term is a specially defined constant
    static bool isSpecialConstant(std::string Term);

    /// true if not dependent on algebra argument, but updatable
    static bool isUpdatableConstant(std::string Term);

    /// optionally input named constants
    static void readConstants(ReadInput & Inp);

    static std::vector<const Algebra*> integrand; // point to desired algebra
    static double getParameter(unsigned int k,std::string Definition); ///< evaluate k'th parameter in [alg0,alg1,...]
    static double smooth3rdDegree(double x, double a, double b, int der=0); ///< smooth 3rd degree polynomial with p(a)=0, p(b)=1
    static double valUp1(double x, double a, double b); ///< polynomial: value=1 at b, derivative=0, val(a)=0
    static double derUp1(double x, double a, double b); ///< polynomial: derivative=1 at b, value=0, val(a)=0
    static void Test(); ///< usage examples
    static std::complex<double> constantValue(std::string Term); ///< evaluate Algebra of constants
    static double realConstant(std::string Term); ///< real valued Algebra of constants
    static int integerConstant(std::string Term); ///< integer valued Algebra of constants
    static std::string listConstants(){return tools::listMapKeys(specialConstants," ");}
    static std::string listStandard();

    class Int:public Integrate::Tools
    {
    public:
        Int( double AccRel=1.e-12, double AccAbs=1.e-12,
             std::vector<std::vector<unsigned int> > NQuad=std::vector<std::vector<unsigned int> >(0),
             std::string Kind="GaussLegendre",
             std::string KindInf="GaussLaguerre")
            : Integrate::Tools(AccRel,AccAbs,NQuad,Kind,KindInf){}

        std::complex<double> nDim(const std::vector<std::vector<double > > Vol,
                                  const std::function<std::complex<double>(const std::vector<double> &)> Func,
                                  const std::vector<double> Params=std::vector<double>(0)
                ) {return Integrate::NDim<std::complex<double>,double,Int>(*this,Vol,Func,Params);}

        std::complex<double> recursive(const std::vector<std::vector<double> > Vol,
                                       const std::function<std::complex<double>(const std::vector<double> &)> Func,
                                       const std::vector<double> Params=std::vector<double>(0)
                ) {return Integrate::Recursive<std::complex<double>,double,Int>(*this,Vol,Func,Params);}
    };
};

class AlgebraConstant:public Algebra{
    std::complex<double> c;
public:
    AlgebraConstant(const std::string Term, const std::complex<double> Val):c(Val){_definition=Term;}
    AlgebraConstant(const std::string Term)
        :c(tools::string_to_complex(tools::stringInBetween(Term,"const[","]",true)))
    {_definition=Term;}
    inline std::complex<double> val(std::complex<double> Q) const {return c;}
};

class  AlgebraExternalFunction:public Algebra {
    std::function<std::complex<double>(std::complex<double>)> _funcQ;
public:
    AlgebraExternalFunction(const std::string Term,std::function<std::complex<double>(std::complex<double>)> FuncQ)
        :_funcQ(FuncQ){_definition=Term;}
    inline std::complex<double> val(std::complex<double> Q) const {return _funcQ(Q);}
};

class AlgebraConstantUpdatable:public Algebra {
public:
    AlgebraConstantUpdatable(const std::string Term){_definition=Term;}
    inline std::complex<double> val(std::complex<double> Q) const {return updatableConstants[_definition];}
};

class AlgebraQ:public Algebra{
public:
    AlgebraQ(){_definition="Q";}
    inline std::complex<double> val(std::complex<double> Q) const {return Q;}
};

class AlgebraSqrt:public Algebra{
public:
    AlgebraSqrt(const std::string Term){_definition="sqrt";A.push_back(new Algebra(argument(Term)));}
    inline std::complex<double> val(const std::complex<double> Q) const {return sqrt(A[0]->val(Q));}
};

class AlgebraSin:public Algebra{
public:
    AlgebraSin(const std::string Term){_definition="sin";A.push_back(new Algebra(argument(Term)));}
    inline std::complex<double> val(const std::complex<double> Q) const {return sin(A[0]->val(Q));}
};

class AlgebraCos:public Algebra{
public:
    AlgebraCos(const std::string Term){_definition="cos";A.push_back(new Algebra(argument(Term)));}
    inline std::complex<double> val(const std::complex<double> Q) const {
//        std::ofstream fil;
//        fil.open((ReadInput::main.output()+"cos").c_str(), std::ios_base::app);
//        fil<<Q.real()<<", "<<cos(A[0]->val(Q)).real()<<", "<<cos(A[0]->val(Q)).imag()<<std::endl;
        return cos(A[0]->val(Q));}
};

class AlgebraExp:public Algebra{
public:
    AlgebraExp(const std::string Term){_definition="exp";A.push_back(new Algebra(argument(Term)));}
    inline std::complex<double> val(const std::complex<double> Q) const {return exp(A[0]->val(Q));}
};

/// power
class AlgebraPow:public Algebra{
public:
    AlgebraPow(const std::string Term){
        _definition=Term.substr(0,Term.find("]")+1);
        p=getParameter(0,_definition);
        A.push_back(new Algebra(argument(Term)));}
    inline std::complex<double> val(const std::complex<double> Q) const {
        return pow(A[0]->val(Q),p);
    }
private:
    std::complex<double> p;
};

/// characteristic function of the real part
class AlgebraChi:public Algebra{
    double x0,x1;
public:
    AlgebraChi(const std::string Term){
        _definition=Term.substr(0,Term.find("]")+1);
        x0=getParameter(0,_definition);
        x1=getParameter(1,_definition);
        A.push_back(new Algebra(argument(Term)));}
    inline std::complex<double> val(const std::complex<double> Q) const {
        std::complex<double> q=A[0]->val(Q);
        if(x0>q.real() or q.real()>x1)return 0.;
        return 1.;
    }
    std::vector<std::complex<double>> nonAnalyticQ() const{return {x0,x1};}
};

/// morse potential, parameters: strength,minimum position; asymptotically zero, near minimum ~ x^2/2
class AlgebraMorse:public Algebra{
public:
    AlgebraMorse(const std::string Term){
        _definition=Term.substr(0,Term.find("]")+1);
        stiffness=getParameter(0,_definition);
        minpos=getParameter(1,_definition);
        A.push_back(new Algebra(argument(Term)));}
    inline std::complex<double> val(const std::complex<double> Q) const {
        std::complex<double> arg=exp(stiffness*(minpos-A[0]->val(Q)))-1.;

        return (arg*arg-1.)/(2.*stiffness*stiffness);
    }
private:
    double stiffness,minpos;
};

/// truncation function for smooth cutoff of potentials
class AlgebraTrunc:public Algebra{
public:
    AlgebraTrunc(const std::string Term){
        _definition=Term.substr(0,Term.find("]")+1);
        from=getParameter(0,_definition);
        to=getParameter(1,_definition);
        der=0;
        if(tools::subStringCount(Term,",")==2)der=getParameter(2,_definition);
        if(to<0. or to<=from)ABORT("defined smooth interval above 0, is: "+Term);
        A.push_back(new Algebra(argument(Term)));}
    inline std::complex<double> val(const std::complex<double> Q) const {
        std::complex<double> c=A[0]->val(Q);
        double q=c.real();
        if (q<=-to or to<=q)return 0.;
        if(abs(c.imag())>1.e-14)ABORT("cannot call "+_definition+" with complex argument in non-zero range, is "+tools::str(c));
        if(-from<=q and  q<=from){if(der==0)return 1.;return 0.;}
        if( from<q  and  q<to)return smooth3rdDegree( q,from,to,der);
        if( from<-q and -q<to)return smooth3rdDegree(-q,from,to,der)*pow(-1,der);
        return 0.;
    }
    std::vector<std::complex<double>> nonAnalyticQ() const;//{return {-from,-to,from,to};}
private:
    double from,to;
    int der;
};

/// value=1 at upper boundary, other boundary values an derivatives=0
class AlgebraValUp1:public Algebra{
public:
    AlgebraValUp1(const std::string Term){
        _definition=Term.substr(0,Term.find("]")+1);
        from=getParameter(0,_definition);
        to=getParameter(1,_definition);
        if(to<0. or to<=from)ABORT("defined smooth interval above 0, is: "+Term);
        A.push_back(new Algebra(argument(Term)));}
    inline std::complex<double> val(const std::complex<double> Q) const {
        std::complex<double> c=A[0]->val(Q);
        double q=c.real();
        if( from<q  and  q<to)return smooth3rdDegree(q,from,to);
        if(abs(c.imag())>1.e-14)ABORT("cannot call "+_definition+" with complex argument in non-zero range, is "+tools::str(c));
        return 0.;
    }
private:
    double from,to;
};
/// derviative=1 at upper boundary, other boundary values an derivatives=0
class AlgebraDerUp1:public Algebra{
public:
    AlgebraDerUp1(const std::string Term){
        _definition=Term.substr(0,Term.find("]")+1);
        from=getParameter(0,_definition);
        to=getParameter(1,_definition);
        if(to<0. or to<=from)ABORT("defined smooth interval above 0, is: "+Term);
        A.push_back(new Algebra(argument(Term)));}
    inline std::complex<double> val(const std::complex<double> Q) const {
        std::complex<double> c=A[0]->val(Q);
        double q=c.real();
        if( from<q  and  q<to)return derUp1(q,from,to);
        if(abs(c.imag())>1.e-14)ABORT("cannot call "+_definition+" with complex argument in non-zero range, is "+tools::str(c));
        return 0.;
    }
private:
    double from,to;
};

class AlgebraSpherBessel:public Algebra{
public:
    AlgebraSpherBessel(const std::string Term){
        _definition=Term.substr(0,Term.find("]")+1);
        if(std::count(_definition.begin()+_definition.find('['),_definition.begin()+_definition.find(']'),',')!=2){
            DEVABORT("SpherBessel needs three arguments [der,l,R], found "+Term);
        }
        der=getParameter(0,_definition);
        order=int(getParameter(1,_definition));
        surface=getParameter(2,_definition);
        if(der>1)ABORT("first argument 'der' must be 0 or 1, found: "+Term);
        if(order>200)PrintOutput::warning("found absurdly high order for spherical basis, something wrong? "+Term);
        if(surface<=0)ABORT("last argument 'surface' must be >0, found: "+Term);

        A.push_back(new Algebra(argument(Term)));
    }

    std::complex<double> val(const std::complex<double> Q) const;

private:
    unsigned int der;     /// 0 or 1
    unsigned int order;   /// l
    double surface;       /// evaluate at surface*Q
};

class AlgebraExpI:public Algebra{
public:
    AlgebraExpI(const std::string& Term) {
        _definition=Term.substr(0,Term.find("]")+1);
        der=getParameter(0,_definition);
        surface=getParameter(1,_definition);

        A.push_back(new Algebra(argument(Term)));
    }

    std::complex<double> val(const std::complex<double> Q) const;
private:
    int der; //!< can be positive (derivative), zero (value) or negative (antiderivative)
    double surface; //!< evaluate at surface*Q
};

/** @} */
#endif // ALGEBRA_H
