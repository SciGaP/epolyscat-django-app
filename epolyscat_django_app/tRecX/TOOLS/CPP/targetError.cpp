#include "targetError.h"

#include "abort.h"
#include "str.h"

// simple binary search
void TargetError::fitDifference(double RefPar, std::vector<double> P, std::vector<double> D, double FitMin, double FitMax, const double FitEps){
//    if(P.size()!=2 or D.size()!=2)DEVABORT(Sstr+"expect sizes 2, got "+P.size()+D.size());
    if(RefPar<P[1] or P[1]<P[0])DEVABORT(Sstr+"need strictly increasing parameters, got: "+P+RefPar);

    if(D[0]<D[1] or P.size()>2){
        // no convergent behavior
        _par=0.;
        _fac=DBL_MAX/10.;
        return;
    }

    double res(FitMin);
    double emin=DBL_MAX,rat=D[0]/D[1];
    double fitEps=std::min(0.5,0.10001*(FitMax-FitMin));

    for(_par=FitMin;_par+fitEps<FitMax;_par+=fitEps){
        double est=(err(RefPar)-err(P[0]))/(err(RefPar)-err(P[1]));
        double ecur=rat-est;
        if(std::abs(ecur)<emin){
            emin=ecur;res=_par;
        }
    }
    if(FitEps<fitEps)fitDifference(RefPar,P,D,res,res+fitEps,FitEps);
    _par=std::min(FitMax,res+fitEps*0.5);

    // recalculate factor
    _fac=1.;
    _fac=D[1]/std::abs(err(RefPar)-err(P[1]));
};

// select function that best fits behavior of CPars in the vicinity of CPars[ref]
std::shared_ptr<TargetError> TargetError::factory(int ref,int klo, int kup, double ParMin, double ParMax, const std::vector<double> &CPars,
                                                  const std::vector<double> &Deltas, double FitEps){

    // no error estimate - changes do not decrease well enough above present
    if(*std::min_element(Deltas.begin()+klo,Deltas.end())>0.3)
        return std::shared_ptr<ErrorNone>(new ErrorNone(ParMin,ParMax));

    // determine suitable vicinty
    // (currently done outside)

    // select parameters to use for fitting
    std::vector<double> fPars({CPars[klo],CPars[kup]}),fDels({Deltas[klo],Deltas[kup]});

    // detect onset of noise
    int noise=std::is_sorted_until(Deltas.begin()+klo,Deltas.end(),[](double a,double b){return a>b;})-Deltas.begin();
    if(noise<ref){
        // not strictly decreasing errors - assume noise and estimate width
        fPars={CPars.begin()+noise,CPars.end()};
        fDels={Deltas.begin()+noise,Deltas.end()};
    }

    std::vector<double> dev;
    std::vector<std::shared_ptr<TargetError>> funcs({std::shared_ptr<ErrorPolynom>(new ErrorPolynom(ParMin,ParMax))
                                                     ,std::shared_ptr<ErrorExp>(new ErrorExp(ParMin,ParMax))
                                                     ,std::shared_ptr<ErrorNoise>(new ErrorNoise(CPars.back(),ParMax))
                                                    });
    for(auto f: funcs){
        f->fitDifference(CPars.back(),fPars,fDels,FitEps);
        dev.push_back(0.);
        for(size_t k=0;k<Deltas.size();k++){
            double d=Deltas[k],e=f->err(CPars[k]);
            if(e*100>Deltas[kup]){
                // consider only errors that are on the scale of relevant Deltas
                dev.back()+=std::pow((e-d)/(e+d),2);
            }
        }
        if(CPars.size()<4)break; // assume first hypothesis, if only 3 datapoints
    }
    auto f=funcs[std::min_element(dev.begin(),dev.end())-dev.begin()];

    // return best hypothesis
    return funcs[std::min_element(dev.begin(),dev.end())-dev.begin()];
}

std::vector<double> TargetError::estimates(std::vector<double> &Pars, std::vector<double> &Deltas) const {
    std::vector<double> res(Deltas);
    res.push_back(err(Pars.back()));
    for(size_t k=0;k<Deltas.size();k++)
        if(Deltas[k]>100*Deltas.back())res[k]=Deltas[k];
    return res;
}

std::vector<double> ErrorNoise::estimates(std::vector<double> &Pars, std::vector<double> &Deltas) const {
    std::vector<double> res(Deltas);
    res.push_back(err(Pars.back()));
    for(size_t k=0;k<Pars.size();k++){
        res[k]=sqrt(std::pow(Deltas[k],2)+std::pow(err(Pars[k])-Deltas[k],2));
    }
    return res;
}
