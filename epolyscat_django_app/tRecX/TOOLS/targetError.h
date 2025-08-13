#ifndef TARGETERROR_H
#define TARGETERROR_H

#include <vector>
#include <string>
#include <math.h>
#include <memory>

#include "tools.h"
#include "str.h"

/// estimates error of convergence measure under a given hypothesis, e.g. exponential or power law
class TargetError{
protected:
    double _fac; // err(X)=_fac*func(X)
    double _par; // current value of parameter
    double _parMin,_parMax; // admissible range
    TargetError & setPar(double Par){_par=Par;return *this;}
    void fitDifference(double RefPar, std::vector<double> P, std::vector<double> D ,double FitMin,double FitMax,const double FitEps);
public:
    /// try various error behaviors, return best fit near kup
    static std::shared_ptr<TargetError>  factory(int ref, int klo,int kup, double ParMin, double ParMax,
                                                 const std::vector<double> &CPars,const std::vector<double> &Deltas,double FitEps);
    virtual void fitDifference(double RefPar, std::vector<double> P, std::vector<double> D,const double FitEps){
        fitDifference(RefPar,P,D,_parMin,_parMax,FitEps);}
    TargetError(double FitMin, double FitMax):_fac(1),_parMin(FitMin),_parMax(FitMax){}
    virtual double err(double X) const =0;
    virtual double inv(double Err) const =0;
    virtual std::string algebra() const=0; ///< returns algebra string for function at present parameters
    std::vector<double> parameters(){return {_par,_fac};}
    virtual std::vector<double> estimates(std::vector<double> &Pars, std::vector<double> &Deltas) const;

};

/// assume error fac*exp(-kappa*Par), fit kappa and fac
class ErrorExp: public TargetError{
public:
    ErrorExp(double FitMin,double FitMax):TargetError(FitMin,FitMax){}
    double err(double X) const {return _fac*std::exp(-_par*X);}
    double inv(double Err) const {return (std::log(_fac)-std::log(Err))/_par;}
    std::string algebra() const {return tools::str(_fac,4)+"*exp[-"+tools::str(_par,4)+"](Q)";}
};

/// assume error fac/(Par^kappa), fit kappa and fac
class ErrorPolynom: public TargetError{
public:
    ErrorPolynom(double FitMin,double FitMax):TargetError(FitMin,FitMax){}
    double err(double X) const {return _fac*std::pow(X,-_par);}
    double inv(double Err) const {return std::pow(Err/_fac,-1./_par);}
    std::string algebra() const {return tools::str(_fac,4)+"*pow["+tools::str(-_par,4)+"](Q)";}
};

/// assume error fac/(Par^kappa), fit kappa and fac
class ErrorNoise: public TargetError{
    virtual void fitDifference(double RefPar, std::vector<double> P, std::vector<double> D , const double FitEps){
        size_t firstInc=std::is_sorted_until(D.begin(),D.end(),[](double a,double b){return a>b;})-D.begin();
        _fac=firstInc==D.size()?D.back():*std::max_element(D.begin()+firstInc,D.end());
        _par=P.back();
    }
public:
    ErrorNoise(double CurPar,double ParMax):TargetError(CurPar,ParMax){}
    double err(double X) const {X=_par;return _fac;}
    double inv(double Err) const {Err=_fac;return std::min(_parMax,std::max(_par*1.5,_par+1.));}
    std::string algebra() const {return "noise("+tools::str(_fac,2)+")";}
    virtual std::vector<double> estimates(std::vector<double> &Pars, std::vector<double> &Deltas) const;
};

/// no meaningful convergence extracted
class ErrorNone: public TargetError{
public:
    ErrorNone(double ParCur,double FitMax):TargetError(ParCur,FitMax){}
    double err(double X) const {X=_par;return DBL_MAX;}
    double inv(double Err) const {Err=DBL_MAX;return std::min(_parMax,std::max(_par*1.5,_par+1.));}
    std::string algebra() const {return " no estimate ";}
};

#endif // TARGETERROR_H
