#include "randomPotential.h"
#include "readInput.h"

#include "axisTree.h"

static double random01(){return (std::rand())/double(RAND_MAX);}

std::map<std::string,RandomPotential> RandomPotential::_list;

void RandomPotential::read(ReadInput &Inp){

    std::string name;
    Inp.read("RandomPotential","name",name,"randomPotential","name by which potential can be used");

    if(name.find("random")==std::string::npos)
        ABORT("random potential name must contain string \"random\", suggested name: random"+name);

    std::string ax;
    std::string kind;
    Inp.read("RandomPotential","kind",kind,"","type of random potential: step[n]... n steps on elements of axis");

    double vMin,vMax;
    Inp.read("RandomPotential","vMin",vMin,"0","lowest value of potential, 0=-infty");
    Inp.read("RandomPotential","vMax",vMax,"0","highest value of potential, 0=+infty");
    Inp.obsolete("RandomPotential","integral","removed, as it produces unphysical artefacts");
    //    Inp.read("RandomPotential","integral",integ,"0","integral over potential - scale such that int pot = this value, 0..do not scale");
    Inp.read("RandomPotential","axis",ax,"","name of axis for 1d random potential - must match one of the axis names");

    int seed;
    Inp.read("RandomPotential","seed",seed,"0","if !=0 reset seed, create reproducible sequence of potentials on a given machine");
    if(not Inp.found("RandomPotential"))return;

    if(seed!=0)std::srand(seed);

    std::vector<double> qs;
    int line(1);
    AxisTree axes(Inp,line);
    std::vector<std::string> avail;
    double lRange(0),uRange(0); // plot range
    for(const AxisTree* a=&axes;a!=0;a=a->nodeNext()){
        avail.push_back(a->name);
        if(a->name==ax){
            qs=a->elementBoundaries();
            lRange=a->lowerRange();
            uRange=a->upperRange();
        }
    }
    if(qs.size()==0)ABORT(Sstr+ax+" does not match any axis, available: "+avail);

    // remove end points (which may be infinite)
    if(qs.back()>= DBL_MAX/2)qs.pop_back();
    if(qs.front()<=-DBL_MAX/2)qs.erase(qs.begin());

    std::string s=tools::stringInBetween(kind,"[","]");
    if(s==kind)ABORT("need kind specification by [n], got: "+kind);
    int n=tools::string_to_int(s);

    if(n+1>int(qs.size()))
        ABORT(Sstr+"must have "+(n+1)+" finite interval boundaries for "+n+"pieces on"+ax+", got:"+qs.size());

    // randomly remove interval boundaries until desired number of intervals
    while(int(qs.size())>n+1)qs.erase(qs.begin()+(std::rand()%qs.size()));

    std::vector<std::string> defs(1,"0");
    if(kind.find("step[")==0){
        for(size_t k=0;k<qs.size()-1;k++)
            defs.push_back(tools::str(vMin)+"+Random01*"+tools::str((vMax-vMin)));
    }
    else
        ABORT("unknown kind "+kind);
    defs.push_back("0");

    RandomPotential rpot(qs,defs);
    if (_list.count(name) and not (rpot==_list[name]))
        DEVABORT("RandomPotential "+name+" was generated before with differend parameters");
    _list[name]=rpot;

    // make factory visible to Algebra
    Algebra::addExternalFactory(RandomPotential::factory,"AlgebraPiecewise");

    if(lRange<uRange){
        auto alg=factory(name);
        alg->plot(Inp.output()+"randomPotential",lRange,uRange,200);
    }
}

RandomPotential::RandomPotential(std::vector<double>Qs, std::vector<std::string> Defs)
    :_qs(Qs),_defs(Defs){
    for(auto d: _defs)_random01.push_back(random01());
}

const Algebra *RandomPotential::factory(std::string Name){
    if(not _list.count(Name))return nullptr;

    auto defs=_list[Name]._defs;
    for(size_t k=0;k<defs.size();k++){
        defs[k]=tools::substringReplace(defs[k],"Random01",tools::str(_list[Name]._random01[k]));
    }

    AlgebraPiecewise* res=new AlgebraPiecewise(_list[Name]._qs,defs);
    return res;
}
