#include "basisIntegrableRpow.h"
#include <complex>

#include "printOutput.h"
#include "basisSetDef.h"
#include "basisIntegrable.h"
#include "basisDvr.h"
#include "basisAssocLeg.h"
#include "basisMonomial.h"
#include "basisMat1D.h"
#include "eigenTools.h"
#include "index.h"

// this should got to basis integrable
int defOrder(std::string Def){
    std::string sOrder=tools::cropString(tools::splitString(Def,',')[2]);
    if(sOrder.find_first_of("123456789")!=0)DEVABORT("3rd entry (order) must be integer, found: "+Def);
    return tools::string_to_int(sOrder);
}

// Schmidt-orthogonalize all wrt basis function Ref
static Eigen::MatrixXd transSchmidt(const BasisIntegrable& Bas, int Ref){
    BasisMat1D ovr("<1>",&Bas,&Bas);
    if(not ovr.mat().imag().isZero(1.e-12))DEVABORT("non-real overlap for "+Bas.name());

    Eigen::MatrixXd o(ovr.mat().real());

    // move Ref to first index
    o.col(Ref).swap(o.col(o.cols()-1));
    o.row(Ref).swap(o.row(o.rows()-1));

    std::vector<std::vector<double>> vo,vt;
    for(int j=0,j0=0;j<o.cols();j++,j0+=o.rows())vo.push_back(std::vector<double>(o.data()+j0,o.data()+j0+o.rows()));

    tools::gramSchmidtTrans(vo,vt);

    //NOTE: likely vt[n] is the n-th COLUMN of the transformation
    o=Eigen::MatrixXd::Zero(o.rows(),o.cols());
    for(int j=0;j<o.cols();j++)
        for(size_t i=0;i<vt[j].size();i++)o(i,j)=vt[j][i];

    // move first to Ref
    o.col(Ref).swap(o.col(o.cols()-1));
    o.row(Ref).swap(o.row(o.rows()-1));
    return o;
}

static std::string joinString(const std::vector<std::string> & VecStrings, std::string Deliminator){
    std::string res;
    for(auto s: VecStrings)res+=s+Deliminator;
    return res.substr(0,res.length()-Deliminator.length());
}

static int power(std::string Def){
    std::string sPow=tools::stringInBetween(Def.substr(0,Def.find("*")),"[","]");
    sPow=tools::cropString(tools::stringInBetween(sPow,"l{","}",true));
    if(sPow.find_first_of("0123456789")!=0)DEVABORT("need integer power, got: "+Def);
    int p=tools::string_to_int(sPow);
    if(p<0 or p>INT_MAX/2)ABORT(Sstr+"illegal power"+p+"from"+Def);
    return p;
}

BasisIntegrableRpow::BasisIntegrableRpow(const std::string Def)
    :BasisIntegrable(Def.substr(Def.find("*")+1)){

    if(lowBound()*upBound()<0.)DEVABORT("cannot use Rpow across 0");
    _name=tools::splitString(Def,',',"[(","])")[0];
    _powQ=lowBound()==0.?power(Def):0;

    std::string defIntegrable=Def.substr(Def.find("*")+1);
    std::vector<std::string> def=tools::splitString(defIntegrable,',',"[(","])");
    // replace size for integrable basis
    if(def[1].find_first_of("123456789")!=0)DEVABORT("2nd entry (size) must be non-zero integer, found: "+Def);
    int siz=tools::string_to_int(def[1]);
    if(siz==1){
        _basInt=new BasisMonomial(_powQ,lowBound(),upBound(),_powQ);
    }
    else{
        def[0]=def[0].substr(def[0].find("*")+1);
        def.push_back(", noDVR");
        _basInt=dynamic_cast<const BasisIntegrable*>(BasisAbstract::factory(joinString(def,", ")));
        if(_basInt==nullptr)DEVABORT("basis after * not BasisIntegrable: "+Def);
    }

    if(lowBound()==0.){

        // orthonormalize againts last
        _trans=transSchmidt(*this,size()-1);
        // verify
        BasisMat1D ovr("<1>",this,this);
        if(not ovr.mat().isIdentity(1.e-12))
            ABORT(Sstr+"transform to orthonormal failed for "+_basInt->isDVR()+str()+"\n"+EigenTools::str(ovr.mat(),2,'R'));

        // ensure value of last == 1:
        std::vector<std::complex<double>>val,der;

        valDer({upBound()},val,der);
        if(std::abs(val[size()-1])!=0.)_trans.col(size()-1)/=val[size()-1].real();
    }
}



void BasisIntegrableRpow::valDer(const std::vector<std::complex<double> > &X,
                                 std::vector<std::complex<double> > &Val, std::vector<std::complex<double> > &Der, bool ZeroOutside) const {
    _basInt->valDer(X,Val,Der,ZeroOutside);

    if(_trans.size()>0){;
        Eigen::Map<Eigen::MatrixXcd>(Der.data(),X.size(),size())
                =Eigen::Map<Eigen::MatrixXcd>(Der.data(),X.size(),size())*_trans;
        Eigen::Map<Eigen::MatrixXcd>(Val.data(),X.size(),size())
                =Eigen::Map<Eigen::MatrixXcd>(Val.data(),X.size(),size())*_trans;
    }
}

void BasisIntegrableRpow::resolvePower(BasisSetDef &Def, std::vector<unsigned int> Branch, const std::vector<const Index*> Path){

    std::string pathHier;
    for(auto p: Path)pathHier+=p->axisName()+".";
    pathHier.pop_back();

    // not Rpow
    if(Def.funcs.substr(0,4)!="Rpow")return;

    // already resolved
    if(Def.funcs.substr(0,5)=="Rpow[" and Def.funcs.find_first_of("123456789")==5)return;

    // not starting from R=0
    if(Def.lowBound()!=0.){
        Def.funcs=Def.funcs.substr(Def.funcs.find("*")+1);
        return;
    }
    if(Def.funcs.substr(Def.funcs.find("*")+1)!="polynomial")
        ABORT("for now, Rpow only for polynomial basis, got: "+Def.funcs);


    // determine axis
    std::string ax;
    if(Def.funcs.substr(0,6)=="RpowL*"){
        // seek uniqe Eta-axis
        std::vector<const Index*> iEta;
        for(auto ix: Path)
            if(ix->axisName().substr(0,3)=="Eta")iEta.push_back(ix);
        if(iEta.size()==0)ABORT("need  Eta-axis above "+Path.back()->axisName()+", found: "+pathHier);
        if(iEta.size()>1)ABORT("multiple Eta-axis above "+Path.back()->axisName()+": "+pathHier+", specify as RpowL{Eta..}");
        ax=iEta[0]->axisName();
    }
    else if(Def.funcs.substr(0,5)=="Rpow["){
        std::string powDef=tools::stringInBetween(Def.funcs,"Rpow[","]",true);
        if(powDef.find("l{")!=std::string::npos)ax=tools::stringInBetween(powDef,"l{","}",true);
    }
    if(ax=="")ABORT("admissible Rpow formats Rpow[l{Eta..}], RpowL, got: "+Def.funcs);

    std::string res=Def.funcs;

    int pow=INT_MAX;
    for(size_t k=0;k<Path.size();k++){
        if(Path[k]->axisName()==ax and dynamic_cast<const BasisAssocLeg*>(Path[k]->basis())){
            pow=dynamic_cast<const BasisAssocLeg*>(Path[k]->basis())->lValue(Branch[k]);
            break;
        }
    }
    if(pow==INT_MAX)ABORT("axis "+ax+" does not have assocLegendre basis for defining L");
    if(pow>80){
        PrintOutput::warning(Sstr+"power>80, may cause numerical instability, set power=80",1);
        pow=80;
    }
    Def.funcs=pow>0?(Str("Rpow[","")+pow+"]*DVR: jacobi[0,"+2*pow+"]"):std::string("polynomial");
    Def.order=std::max(3,int(Def.order-pow));
}

std::string BasisIntegrableRpow::inputToStrDefinition(BasisSetDef Def){
    if(Def.funcs.find("Rpow")>2)return "";

    std::string powDef=Def.funcs.substr(0,Def.funcs.find("*"));
    Def.funcs=Def.funcs.substr(Def.funcs.find("*")+1);
    const BasisIntegrable* bas=dynamic_cast<const BasisIntegrable*>(BasisAbstract::factory(Def));
    if(not bas)ABORT("cannot interprete "+Def.str());

    return Def.lowBound()==0.?powDef+"*"+bas->strDefinition():bas->strDefinition();
}

std::string BasisIntegrableRpow::strDefinition() const{
    return lowBound()==0.?"Rpow["+tools::str(_powQ)+"]*"+_basInt->strDefinition():_basInt->strDefinition();
}























