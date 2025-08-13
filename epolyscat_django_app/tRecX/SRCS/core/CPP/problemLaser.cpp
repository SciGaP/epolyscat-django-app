#include "problemLaser.h"

#include "pulse.h"
#include "constants.h"

ProblemLaser::ProblemLaser(ReadInputList &Inp):ProblemLaser()
{
    Pulse laser(Inp,"Laser");
    std::vector<double> lambdas=laser.wavelengths();
    std::vector<double> intens=laser.intensities();
    tools::sortByKey(lambdas,intens);
    for(size_t k=lambdas.size();k>0;k--){
        if(lambdas[k]==lambdas[k-1]){
            intens[k-1]+=intens[k];
            lambdas.erase(lambdas.begin()+k);
            intens.erase(intens.begin()+k);
        }
    }
    _intensWcmSq=0.;
    if(lambdas.size()==0)return;

    _lambda_nm=lambdas.front();
    _intensWcmSq=intens.front();
    _polarization=laser.polarization();
    double photonEnergy=Units::convert(physics::planck_constant*physics::speed_of_light/(_lambda_nm*1.e-9),"SI|energy","au");

    _durationOptCyc=laser.duration()*photonEnergy/(2*math::pi);

    // get actual laser input lines
    allLinesInCategory(Inp.outputTopDir()+"inpc","Laser",_inputLines);
}

void ProblemLaser::addToInput(std::ofstream &Stream) const{

    Stream<<std::endl;
    if(_intensWcmSq==0)
        Stream<<"# --- no laser in this problem ------ "<<std::endl;
    else{
        Stream<<"# --- laser characteristics, see input below for full definition ------ "<<std::endl;
        Stream<<"# axWavelength "+tools::str(_lambda_nm,4)+" nm @ "+tools::str(_intensWcmSq,3)+" W/cm2, duration "
                +tools::str(_durationOptCyc,3)+" OptCyc, "+_polarization+" polarization"<<std::endl;
    }
    Stream<<std::endl;

    // write the problem's actual laser input, not just key parameters
    writeInputLines(Stream);
}

std::string ProblemLaser::str() const{
    std::string s="-- no Laser --";
    if(_intensWcmSq>0.){
        s=Sstr+_intensWcmSq+" W/cm2, "+_lambda_nm+" nm, "+_durationOptCyc+" OptCyc, "+_polarization+" polarization";
    }
    return "Laser: "+s;
}
