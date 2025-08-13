#ifndef METRIC_H
#define METRIC_H

#include <vector>
#include "qtEigenDense.h"

class Metric
{
    std::string _name;
protected:
    virtual double _distance(const Eigen::VectorXd &A, const Eigen::VectorXd &B) const =0;
public:
    Metric(std::string Name):_name(Name){}
    std::string name() const {return _name;}
    double distance(const std::vector<double> &A, const std::vector<double> &B){
        return _distance(Eigen::Map<const Eigen::VectorXd>(A.data(),A.size()),Eigen::Map<const Eigen::VectorXd>(B.data(),B.size()));
    }
    double distance(const Eigen::VectorXd &A, const Eigen::VectorXd &B){return _distance(A,B);}
    virtual std::string definition() const = 0;
};

/// reltive distance of two vectors
class RelativeDistanceWeighted: public Metric{
    const double _weight;
    double _distance(const Eigen::VectorXd &A, const Eigen::VectorXd &B) const;
public:
    /// L2 norm of d: d(i)=2*(A(i)-B(i))/max(eps,|A(i)|+|B(i)|), eps=max[i](|A|(i)+|B|(i)) / Weight
    RelativeDistanceWeighted(double Weight);
    std::string definition() const;
};

#endif // METRIC_H
