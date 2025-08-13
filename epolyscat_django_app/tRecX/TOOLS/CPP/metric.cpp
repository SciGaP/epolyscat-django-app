#include "metric.h"
#include "qtEigenDense.h"
#include "abort.h"
#include "printOutput.h"

RelativeDistanceWeighted::RelativeDistanceWeighted(double Weight):Metric("RelativeDistanceWeight"),_weight(Weight){
    if(_weight >= 1.)DEVABORT("must have weight < 1 when computing relative distance");
    if(_weight < 1e-2)PrintOutput::DEVwarning(Sstr+"small weight of"+_weight+"for relative distance, recommended value>0.01");
}
double RelativeDistanceWeighted::_distance(const Eigen::VectorXd &A, const Eigen::VectorXd &B) const{
    if(_weight==0.)DEVABORT("cannot have weight=0");
    Eigen::VectorXd s=Eigen::Map<const Eigen::VectorXd>(A.data(),A.size()).cwiseAbs()+Eigen::Map<const Eigen::VectorXd>(B.data(),B.size()).cwiseAbs();
    // simple metric L2-size of relative error
    double eps=s.lpNorm<Eigen::Infinity>()*_weight;
    for(int i=0;i<s.size();i++)s(i)=2*(A(i)-B(i))/std::max(s(i),eps);
    return s.lpNorm<2>();
}

std::string RelativeDistanceWeighted::definition() const{
    return name()+","+tools::str(_weight);
}
