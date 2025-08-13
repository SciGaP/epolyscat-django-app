#ifndef CONVERGENCETARGET_H
#define CONVERGENCETARGET_H

#include <string>
#include <vector>

/// supplies columns of data for convergence
class ConvergenceTarget{
    std::string _name;
protected:
    std::string _def;
    virtual std::string file() const =0; ///< file from where data was drawn
public:
    ConvergenceTarget(std::string Name):_name(Name){}
    std::string name() const {return _name;}
    /// columns of target - empty if file does not exist (yet)
    virtual std::vector<std::vector<double>> cols(std::string Dir)=0;
    virtual std::string str() const {return " --- no target info ---";}
    virtual std::string definition() const=0;
};


#endif // CONVERGENCETARGET_H
