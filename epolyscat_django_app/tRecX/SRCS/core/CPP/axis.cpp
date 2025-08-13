// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "axis.h"
#include "readInput.h"
#include "basisAbstract.h"
#include "basisIntegrable.h"
#include "basisNdim.h"
#include "printOutput.h"
#include "str.h"
#include "axisTree.h"

using namespace std;
using namespace tools;

Axis::Axis(string Coor, unsigned int Ncoefs, double Qlow, double Qup, string Functions, unsigned int Order, ComplexScaling Comsca)
    : comsca(Comsca) {
    construct(Coor,Ncoefs,Qlow,Qup,Functions,Order);
}

Axis::Axis(const string Name, const ComplexScaling Comsca, const BasisSetDef &Def):
    name(Name),comsca(Comsca)
{
    coor = Def.coor;
    basDef.assign(1,Def);
}

Axis::Axis(const std::string Name,const std::vector<const BasisAbstract*> & Bases)
    :name(Name),coor(Coordinate::fromString(Name)),bases(Bases),comsca(ComplexScaling())
{}


// separate this for use by other routines;
std::string Axis::readName(ReadInput & Inp, int Line, string Default, std::string InputName){
    string coor; //
    Inp.read(InputName,"name",coor,Default,
             "coordinate string: X,Y,Z,R,Rn,Rho,Eta,Phi,Ndim,Vec,Channel,Ion,Neutral,Orbital,Xperiodic,Yperiodic,PXi,PEta,BLANK",Line)
            .texdocu(R"tex(
                     Name of the discrete or continnous ``coordinate''.
                     Continuous axes can be specified section-wise by repeating the same axis name on multiple lines
                     (or leaving the axis name blank, which defaults to the preceding name. In that case \lcode{upper end} of one section
                     must match the \lcode{lower end} of the subsequent one, see \nameref{docu:Axis:lower end})
                     \begin{itemize}
                     \item[X,Y,Z] cartesian coordinates
                     \item[Rn] radial spherical coordinate $\l \phi|\chi\r=\int_0^\infty dr \phi^*(r)\chi(r)$
                     \item[Rho] polar coordinate $\l \phi|\chi\r=\int_0^\infty d\rho \rho \phi^*(\rho)\chi(\rho)$
                     \item[Eta] $\eta=\cos\th$ $\l \phi|\chi\r=\int_{-1}^1 d\eta \phi^*(\eta)\chi(\eta)$
                     \item[Eta] azimutal angle, $\l \phi|\chi\r=\int_{0}^{2\pi} d\varphi \phi^*(\varphi)\chi(\varphi)$
                     \item[Ndim] index count in multi-dimensional basis, see \nameref{???}
                     \item[Vec] general vector index
                     \item[Channel] channel indes for a multi-channel calculation (haCC, see \nameref{???})
                     \item[BLANK] (internal, not to be input by user) indicates end of coordinate list
                     \item[Ion] CI molecular ion basis (multi-electron functions)
                     \item[Neutral] CI neutral molecule basis (multi-electron functions)
                     \item[Xperiodic] (also Yperiodic) cartesian with periodic boundary condistions
                     \item[PXi] and PEta $\xi,\eta$ for parabolic coordinates:
                     $\l \phi|\chi\r=\int_{0}^\infty \int_{0}^\infty d\xi d\eta \phi^*(\xi,\eta)(\xi+\eta)\chi(\xi,\eta)$
                     \end{itemize}
                     )tex");
    return coor;
}
std::string Axis::readAlternate(ReadInput & Inp, int Line, std::string InputName,std::string Modify, std::string AlternateModify)
{
    std::string amod,res(Modify);
    std::string altDefault="alternate";
    std::string nam="alternateAxis"+(AlternateModify=="alternate"?"":"["+AlternateModify+"]");
    Inp.read(InputName,nam,amod,"","as Axis:modify, but for idxName in alternateAxis[idxName] (defaults to idxName=\""+altDefault+"\") {allow all}",Line);
    if(AlternateModify!=""){
        if(amod=="" and Inp.found(InputName,nam))
        ABORT(Sstr+"for Axis:"+nam+"at line"+Line+"specify keyword: ommit, keep, only, order= etc. (Axis:modify for full list)" );
        if(amod!="keep")res=amod;
    }
    return res;
}

Axis::Axis(ReadInput & in,unsigned int Line,string InputName,std::string AlternateModify){
    if(Line==0)DEVABORT("cannot read axis from Line=0");
    string Coor,functions;
    double qlow,qup;
    int ncoef,order;
    bool exact;

    // read axis name backward from current line until non-blank: we need not repeat axis name on each line
    for (int l=Line;l>0;l--){
        if("BLANK"!=(Coor=readName(in,l,"BLANK",InputName)))break;
    }

    if(Line!=0 and functions.find("grid")!=string::npos)
        ABORT("finite difference grid are only supported in 1d with a single Axis input line");

    Coordinate coor(Coor);
    in.read(InputName,"functions",functions,coor.defaultFunction(),"which functions to use",Line)
            .texdocu(R"tex(
                     Which function to use the on the given axis(-section) Axis:\nameref{docu:Axis:name}s, defaults exist for many,
                     examples see in Tutorials \nameref{docu:tutorials}
                     )tex");
    int nElem(-1);
    in.exclude(InputName,InputName,"nElem","nCoefficients");
    if(in.found(InputName,"nElem")){
        ncoef=-1;
        // new input style, using number of elements
        in.read(InputName,"nElem",nElem,"1","number of elements on axis, defaults to 1",Line);

        in.read(InputName,"order",order,tools::str(Coordinate::automaticOrder(Coor,ncoef)),"number of functions on element (default - adjust to 4 elements)",Line)
                .texdocu(R"tex(
                         Order of the discretization on the \lcode{Axis}.
                         Defaults to Axis:\nameref{docu:Axis:nCoefficients} for non-finite-element (FE) discetretizations.
                         For FE-discretization, the number of FE's is determinde as \lcode{int(nCoefficients/order)}.
                         )tex");
        if(order==-1)ABORT("must specify "+InputName+":order when using "+InputName+":nElem");
        ncoef=nElem*order;
    }
    else {
        string scoef=BasisNdim::allNdim(functions).size()?"1":"-1";
        in.read(InputName,"nCoefficients",ncoef,scoef,"(approximate) number of coefficients",Line)
                .texdocu(R"tex(
                         Number of linear coefficients for discretization of a given Axis.
                         Depending on the discretization, functions, and boundary conditions
                         the actual number of linearly independent coefficients may be smaller than this, see also Axis:\nameref{docu:Axis:order}
                         )tex");
        if(not in.noInputFile() and ncoef==-1)return;
        in.read(InputName,"order",order,tools::str(Coordinate::automaticOrder(Coor,ncoef)),"number of functions on element (default - adjust to 4 elements)",Line)
                .texdocu(R"tex(
                         Order of the discretization on the \lcode{Axis}.
                         Defaults to Axis:\nameref{docu:Axis:nCoefficients} for non-finite-element (FE) discetretizations.
                         For FE-discretization, the number of FE's is determinde as \lcode{int(nCoefficients/order)}.
                         )tex");
    }
    in.read(InputName,"lower end",qlow,tools::str(coor.qmin),"lower boundary of the axis",Line)
            .texdocu(R"tex(
                     Lower end of the axis(-section). For some Axis:\nameref{docu:Axis:name} there are defaults (discrete starts at 1).
                     If preceede  by an \lcode{Axis} with identical \lcode{name} must be  = \lcode{upper end} of preceding (else abort).
                     )tex");
    in.read(InputName,"upper end",qup ,tools::str(coor.qmax), "upper boundary of the axis",Line)
            .texdocu(R"tex(
                     Upper end of the axis(-section). For some Axis:\nameref{docu:Axis:name} there are defaults (discrete starts at 1).
                     If followed by an \lcode{Axis} with identical \lcode{name} must be  = \lcode{lower end} of next (else abort).
                     )tex");
    std::string modify;
    in.read(InputName,"modify",modify,"keep","path-dependence of basis, e.g. \"order=10 where Eta L<4\" or \"only Eta L<10 \" {allow all}",Line)
            .texdocu(R"tex(
                     examples:
                     \lcode{only Eta L<4} use axis section only where Eta<4
                     \lcode{order=20} sets order to 20 for L < 10, else as main
                     \lcode{order=20 & only Eta L<4} combines modifiers
                     )tex");

    // read alternate modification
    modify=readAlternate(in,Line,InputName,modify,AlternateModify);

    string strExact="false";
    if(coor.exactIntegration)strExact="true";
    in.read(InputName,"exactIntegral",exact,strExact,"use exaxt integrations for axis, even if overall femDVR",Line)
            .texdocu(R"tex([EXPERT USE])tex");

    if(Coor=="Orbital"
            and functions.find("Orbital[")>Coor.length()
            and functions.find("Eigenbasis[")>Coor.length() )
        functions="Orbital["+functions+"]";

    if(functions=="laguerreExp" or functions=="legendre" or functions=="lobatto"  or functions=="radauExp"  )
        ABORT("functions \""+functions+"\" obsolete, use \"polynomial\" or \"polExp[...]\" instead");

    if(functions.find("polExp")==0){
        if(functions.find("polExp[")!=0)ABORT("specify exponent in the form \"polExp[0.123]\"");
        if(order!=-1 and order!=ncoef)ABORT("cannot specify order!=nCoefficients for basis polExp");
        order=ncoef;
    }
    else if(functions.find("besselCoulomb")==0){
        //        if(not functions.find("besselCoulomb[")==0)ABORT("specify exponent in the form \"besselCoulomb[0.123]\"");
        if(order!=-1)ABORT("do not specify order for basis besselCoulomb (is automatic)");
        order=ncoef;
    }

    if(Line!=1 and Coor=="Vec")ABORT("component axis Comp must be first in hierarchy");
    comsca=ComplexScaling(in,Coor);

    if(not in.noInputFile())construct(Coor,ncoef,qlow,qup,functions,order,exact);

    for(auto &d: basDef)d._modify=modify;
}

void Axis::extendBox(std::vector<double> Box){
    if(name!="Rn")DEVABORT("for now only for Rn-axis, is: "+name);
    if(basDef.back().funcs.find("polExp")!=0)DEVABORT("assumes last is polExp, found: "+basDef.back().str());

    BasisSetDef ext=basDef[basDef.size()-2];
    double box=Box[1],len=ext.upBound()-ext.lowBound();
    if(abs(box-boxsize())<1.e-7*len)return; // do not extend, same size

    if(box-boxsize()<0.11*len)DEVABORT("extend size differst too little from present");

    int addElem=int((box-boxsize())/len+0.9);
    ext.scale=ext.scale*(box-boxsize())/(addElem*len);
    for(int k=0;k<addElem;k++){
        ext.shift+=len;
        basDef.insert(basDef.end()-1,1,ext);
    }
    basDef.back().shift=ext.upBound();
    // comsca info is duplicated - NOT NICE!
    comsca._r0upper=ext.upBound();

    for(auto &b: basDef){
        b.comSca=comsca;
    }
}

void Axis::truncateBox(std::vector<double> Box){

    double eps=(Box[1]-Box[0])*1.e-10;
    auto ele=basDef.begin();
    while(ele!=basDef.end()){
        if(ele->upBound()<Box[0]+eps or ele->lowBound()>Box[1]-eps)
            ele++; // keep element
        else if(ele->lowBound()>Box[0]-eps and ele->upBound()<Box[1]+eps)
            basDef.erase(ele,ele+1);
        else
            DEVABORT(Str("truncation boundary does not coincide with element boundary: [")+Box+"] vs ["+ele->lowBound()+ele->upBound()+"]");
    }
}

std::vector<double> Axis::elementBoundaries() const{
    std::vector<double> res;
    if(bases.size()){
        for(auto b: bases)if(b->integrable())res.push_back(b->integrable()->lowBound());
        res.push_back(bases.back()->integrable()->upBound());
    }
    else if(basDef.size()){
        for(auto d: basDef)res.push_back(d.lowBound());
        res.push_back(basDef.back().upBound());
    }
    return res;
}

bool Axis::isElementBoundary(double Val) const{
    //    double eps=1.e-12;
    //    for(unsigned int k=0;k<basDef.size();k++){
    //        if(((basDef[k].lowBound()+eps<Val and Val<basDef[k].upBound()-eps) or
    //            (basDef[k].lowBound()+eps<Val and Val<basDef[k].upBound()-eps)))return false;
    //    }
    auto bs(elementBoundaries());
    double eps=std::abs(bs.back()-bs.front())*1e-12;
    for(auto b: bs)
        if(std::abs(Val-b)<eps)return true;
    return false;
}

void Axis::fromFile(ReadInput &In, std::vector<Axis> & Ax, string Subset, string ReadCategory){
    Ax.clear();
    int l=0;
    string sub="";
    vector<Axis> allAx;
    while(not In.endCategory(ReadCategory,++l) ){
        Axis ax=Axis(In,l,ReadCategory);
        if(ax.basDef.size()==0)break; // empty axis

        // keep track of all axes on input
        allAx.push_back(ax);

        // check for subset
        if(Subset!=""){
            //Note: using previous sub as default - need to specify only for first axis of block
            //            In.read("Axis","subset",sub,sub,"axis will be read into named subset",l);
            sub=AxisTree::readSubset(In,ReadCategory,l,sub);
            if(sub!=Subset)continue; // ignore axes from non-matching subset
        }
        // compose subsequent sections of the same axis
        if(Ax.size()>0 and ax.name==Ax.back().name){
            Ax.back()=Ax.back().append(ax);
        } else {
            Ax.push_back(ax);
        }
    }
    for(unsigned int i=0;i<Ax.size();i++)Ax[i].check();
    for(unsigned int i=0;i<ComplexScaling::names.size();i++)
        if(not In.noInputFile())isNameOfAxis(ComplexScaling::names[i],allAx,"Absorption: axis");

    // quick error check for special case of spherical harmonics input
    for(unsigned int k=0;k<Ax.size();k++){
        if(Ax[k].name.find("Phi")==0){
            for(unsigned int l=0;l<Ax.size();l++){
                if(Ax[l].name.find("Eta")==0){
                    if(Ax[l].basString().find(Ax[k].name)!=string::npos){
                        // related axes
                        if(l<k or Ax[l].maxOrder()*2<Ax[k].maxOrder()+1){
                            Axis::print(Ax);
                            ABORT("polar coordinates: Phi must be above Eta and must have (Phi order) < 2*(Eta order)");
                        }
                    }

                }
            }
        }
    }

}

bool Axis::isNameOfAxis(const string Name, const std::vector<Axis> Ax, const string Mess){
    string s;
    for(unsigned int l=0;l<Ax.size();l++){
        s+=" "+Ax[l].name;
        if(Name==Ax[l].name)return true;
    }
    if(Mess!="")ABORT(Mess+" \""+Name+"\"  does not name an axis:"+s);
    return false;
}

void Axis::construct(string Coor, int Ncoefs, double Qlow, double Qup, string Functions, int Order, bool ExactIntegral){
    if(Coor=="BLANK")return; // do not try to construct blank axis
    // many more consistency checks and default supplementation should go here
    name=Coor;
    coor=Coordinate(Coor);
    if (Ncoefs<1)error("need at least one coefficient");
    if (Order<1)error("order must be > 0");
    if (Order>Ncoefs)error("fewer coefficients than orders");
    if (Qlow >= Qup and Ncoefs>1)ABORT(Str("lower axis boundary must be below upper: ")+name+Qlow+Qup+Ncoefs);
    if (Functions=="automatic")Functions = coor.defaultFunction();

    // correct for tiny floating differences
    if(abs(Qlow-coor.qmin)<=1.e-10*abs(coor.qmin))Qlow=coor.qmin;
    if(abs(Qup -coor.qmax)<=1.e-10*abs(coor.qmax))Qup =coor.qmax;

    vector<int> marg(2);
    marg[0] = 0;       // left margin is first function of element
    marg[1] = Order - 1; // right margin is last function of element
    basDef.clear();
    double x0=Qlow,scal;

    string func=Functions.substr(0,Functions.find('['));
    if(func=="grid"){
        // for grids, make one single element, communicate fd order separately
        if(Ncoefs<1)ABORT("grid axis must have at least two points, found "+tools::str(Ncoefs));
        if(Order<3 or Order%2==0)ABORT("grid axis must have odd order >=3, found "+tools::str(Order));

        vector<double> parms;
        parms.push_back(double(Order));

        // exponential grid scaling damping factor
        parms.push_back(0.);
        if(Functions.find('[')!=string::npos)
            parms.back()=tools::string_to_double(tools::stringInBetween(Functions,"[","]"));
        basDef.push_back(BasisSetDef(Ncoefs-1,x0,Qup-Qlow,func,ExactIntegral,true,true,coor,marg,comsca,false,parms));
    }
    else if(BasisNdim::factory(func)){
        for(auto b: BasisNdim::allNdim(func))
            basDef.push_back(BasisSetDef(b->size(),0.,1.,b->BasisAbstract::name(),false,true,true,Coordinate("Ndim")));
    }
    else {
        int nElem = Ncoefs / Order;
        for (int n = 0; n < nElem; n++){
            scal=(Qup-Qlow)/nElem;
            if((Functions.find("polExp[")==0)){
                //scaling is by inverse of exponent!
                scal=0.5/abs(tools::string_to_double(tools::stringInBetween(Functions,"[","]")));
                if(Qlow<-DBL_MAX/2){
                    scal=-scal; // polyExp on first element is reverted
                    x0=Qup;
                }
                else if(Qup<DBL_MAX/2)ABORT(Functions+" must have lower or upper boundary \"Infty\", is: ["+tools::str(Qlow,4)+","+tools::str(Qup,4)+"]");
            }
            string functions=Functions;
            if(Functions.find("Rl*")==0 and x0!=0.)functions=Functions.substr(3,Functions.find("{")-3);
            //if(Functions.find("sqrt*")==0 and x0!=0.)functions=Functions.substr(5,Functions.find("{")-3);
            basDef.push_back(BasisSetDef(Order,x0,scal,functions,ExactIntegral,n==0,n==nElem-1,coor,marg,comsca));
            x0 +=(Qup-Qlow)/nElem;
        }
    }
}

void Axis::check(){
    if(basDef.size()==0)ABORT("no basis defined");
    if(Coordinate::isDiscrete(coor.cString))return; // no checks on discrete coordinates
    //    if(coor.jaco!=Coordinate::J_one)ABORT("coordinates with non-trivial Jacobian temporarily disabled");

    // check for suspicioius boundaries
    if (lowerEnd()<coor.qmin)error("axis boundary below minimum by: "+tools::str(lowerEnd()-coor.qmin)+" min: "+tools::str(coor.qmin));
    if (upperEnd()>coor.qmax)error("axis boundary above maximum by: "+tools::str(upperEnd()-coor.qmax)+" max: "+tools::str(coor.qmax));
    if (abs(coor.qmin) != DBL_MAX and lowerEnd()!= coor.qmin)
        PrintOutput::DEVwarning(Sstr+"lower boundary "+coor.cString+"="+lowerEnd()+"above minimum="+coor.qmin);
    if (abs(coor.qmax) != DBL_MAX and upperEnd()!= coor.qmax)
        PrintOutput::DEVwarning(Sstr+"upper boundary "+coor.cString+"="+upperEnd()+"below maximum="+coor.qmax);
    string forbid=".R.Rho.";
    if(forbid.find(coor.name())!=string::npos and comsca.eta!=1.)ABORT("complex scaling on R, Rho temporarily not admitted");
}

/// re-do the axis setup, using present data
void Axis::remake(bool Deriv){
    /// recalculate basis, possibly with different Deriv
    vector<int> marg(2);
    for (unsigned int i = 0; i < basDef.size(); i++){
        basDef[i].margin[0] = 0;					// left margin is first function of element
        basDef[i].margin[1] = basDef[i].order - 1;	// right margin is last function of element
        basDef[i].first=(i==0);
        basDef[i].last=(i==basDef.size()-1);
        basDef[i].deriv=Deriv;
    }

}

/// restrict range an order of an axis
void Axis::constrain(double Lower, double Upper, unsigned int Order){
    double lower=max(Lower,coor.qmin),upper=min(Upper,coor.qmax);
    vector<BasisSetDef> con;
    for (unsigned int i = 0; i < basDef.size(); i++){
        if(basDef[i].upBound()-lower<1e-12*abs(lower) or basDef[i].lowBound()-upper>-1.e-12*abs(upper))continue; // outside range
        con.push_back(basDef[i]);
        con.back().order=min(basDef[i].order,Order);
        con.back().margin[1]=con.back().order-1;
    }
    if(con.size()==0)ABORT("constrained axis has zero size: "+tools::str(Lower)+", "+tools::str(Upper)+" ("+tools::str(Order)+") ");
    con[0].first=true;
    con.back().last=true;
    basDef=con;
    if(lowerEnd()!=lower)PrintOutput::warning("constraint of axis boundary moved from "+tools::str(lower)+" to "+tools::str(lowerEnd()));
    if(upperEnd()!=upper)PrintOutput::warning("constraint of axis boundary moved from "+tools::str(upper)+" to "+tools::str(upperEnd()));
}

/// re-do the axis setup, using present data for overlap matrices for prolate spheroidal coordinates
void Axis::setupXiBasis(std::complex<double> s_k, std::complex<double> q_k){
    /// recalculate basis, possibly with some given overlap matrices
    vector<int> marg(2);
    for (unsigned int i = 0; i < basDef.size(); i++){
        basDef[i].margin[0] = 0;					// left margin is first function of element
        basDef[i].margin[1] = basDef[i].order - 1;	// right margin is last function of element
        basDef[i].first=i==0;
        basDef[i].last=(i==basDef.size()-1);
        basDef[i].deriv=false;
    }
}

/// append an axis to present
Axis Axis::append(const Axis & app) const {
    if(basDef.size()==0)return app;
    if(app.basDef.size()==0)return Axis(*this);
    if(coor.name()!=app.coor.name())
        ABORT("cannot append axis "+app.coor.name()+" to "+coor.name());
    if(basDef.size()>1 and basDef.back().upBound()==DBL_MAX){
        show("axis");
        ABORT("cannot append to infinitely long axis");
    }

    // re-create axis with present comsca
    Axis sum((*this));
    sum.basDef.clear();
    for (unsigned int k=0;k<basDef.size();k++){
        sum.basDef.push_back(basDef[k]);
        sum.basDef.back().comSca=sum.comsca;
        sum.basDef.back().last=false;
    }

    // check whether we can append
    if(!(sum.comsca==app.comsca))ABORT("cannot append: complex scaling parameters do not match");
    if(abs(sum.upperEnd()-app.lowerEnd())>1.e-12*abs(basDef.back().scale)){
        sum.show("sum");
        app.show("app");
        ABORT("ends of axes must match: "+tools::str(sum.upperEnd())+" != "+tools::str(app.lowerEnd()));
    }

    // append further basis functions
    for (unsigned int k=0;k<app.basDef.size();k++){
        sum.basDef.push_back(app.basDef[k]);
        sum.basDef.back().comSca=sum.comsca;
        sum.basDef.back().first=false;
    }
    return sum;
}

void Axis::appendInPlace(const Axis & app){
    if(app.basDef.size()==0)return;
    if(coor.name()!=app.coor.name())
        error("cannot append axis "+app.coor.name()+" to "+coor.name());
    if(basDef.size()>1 and basDef.back().upBound()==DBL_MAX){
        show("axis");
        ABORT("cannot append to infinitely long axis");
    }

    // re-create axis with present comsca

    // check whether we can append
    if(!(comsca==app.comsca))ABORT("cannot append: complex scaling parameters do not match");
    if(abs(upperEnd()-app.lowerEnd())>1.e-12*abs(basDef.back().scale)){
        show("sum");
        app.show("app");
        ABORT("ends of axes must match: "+tools::str(upperEnd())+" != "+tools::str(app.lowerEnd()));
    }

    // append further basis functions
    basDef.back().last=false;
    for (unsigned int k=0;k<app.basDef.size();k++){
        basDef.push_back(app.basDef[k]);
        basDef.back().comSca=comsca;
        basDef.back().first=false;
    }
}

unsigned int Axis::maxSize() const {
    if(     basDef[0].funcs.find("polExp")    ==std::string::npos and
            basDef[0].funcs.find("expIm")    ==std::string::npos and
            basDef[0].funcs.find("polynomial")==std::string::npos
            )
        return basDef[0].order;

    unsigned int sum=0;
    for(unsigned int k=0;k<basDef.size();k++)sum+=basDef[k].order-1;
    if(not BasisFunction::asympZero(basDef[0].funcs) and basDef[0].coor.zeroLow)sum--;      // subtract if first coefficient is discarded
    if(BasisFunction::asympZero(basDef.back().funcs) or not basDef.back().coor.zeroUp)sum++; // add if last coefficient is not discarded
    return sum;
}
unsigned int Axis::maxOrder() const {
    unsigned int maxorder=0;
    for(unsigned int k=0;k<basDef.size();k++)maxorder=max(maxorder,basDef[k].order);
    return maxorder;
}
unsigned int Axis::minOrder() const {
    unsigned int minorder=INT_MAX;
    for(unsigned int k=0;k<basDef.size();k++)minorder=min(minorder,basDef[k].order);
    return minorder;
}

string Axis::basString(string Kind) const {
    string s=basDef[0].funcs;
    for(unsigned int n=0;n<basDef.size();n++)
        if(basDef[n].funcs!=s and s.find(basDef[n].funcs)==string::npos)
            s+="/"+basDef[n].funcs;
    if(s=="grid")s+="["+tools::str(basDef[0].par[1])+"]";
    if(Kind=="brief" and s.length()>20)s=s.substr(0,17)+"...";
    return s;
}

void Axis::show(string Text) const {
    if(Text!="")cout<<Text<<endl;
    cout<<"Name: "+str();
}

string Axis::str(int Brief) const {
    string s,b;
    s=name;
    int ncoef=0,maxOrder=0,minOrder=INT_MAX;
    if(bases.size()){
        for(unsigned int n=0;n<bases.size();n++){
            ncoef+=bases[n]->size();
            if(n==0)ncoef++;
            maxOrder=max(maxOrder,(int) bases[n]->size());
            minOrder=min(minOrder,(int) bases[n]->size());
        }
    } else {
        for(unsigned int n=0;n<basDef.size();n++){
            ncoef+=basDef[n].order-1;
            if(n==0)ncoef++;
            maxOrder=max(maxOrder,(int) basDef[n].order);
            minOrder=min(minOrder,(int) basDef[n].order);
        }
    }
    b+=tools::str(ncoef);
    if(maxOrder==minOrder){
        s+=" Approx. size (elements, order) "+tools::str(ncoef);
        s+=" ("+tools::str(bases.size()?bases.size():basDef.size())+", "+tools::str(minOrder)+")";
        if(ncoef>maxOrder)b+="("+tools::str(maxOrder)+")";
    } else {
        s+=" Approx. size (elements, min/max order) "+tools::str(ncoef);
        s+=" ("+tools::str(bases.size()?bases.size():basDef.size())+", "+tools::str(minOrder)+"/"+tools::str(maxOrder)+")";
        b+="("+tools::str(minOrder)+"/"+tools::str(maxOrder)+")";
    }
    // add interval, unless full coordinate axis
    if(lowerEnd()!=coor.qmin or upperEnd()!=coor.qmax or arg(comsca.eta)!=0.){
        b+="["+tools::str(lowerEnd(),3,DBL_MAX/2)+",";
        if(arg(comsca.eta)!=0.)b+=tools::str(int(comsca.r0up()))+"/"+tools::str(upperEnd(),3,DBL_MAX/2)+"@"+tools::str(int(arg(comsca.eta)/(math::pi)*180.))+"]";
        else b+=tools::str(upperEnd(),3,DBL_MAX/2)+"]";
    }
    if(bases.size()){
        for(unsigned int n=0;n<bases.size();n++)s+="  "+bases[n]->str();
    }else {
        for(unsigned int n=0;n<basDef.size();n++)s+="  "+basDef[n].str();
    }
    if(Brief>0)return b;
    return s;
}

void Axis::print(const vector<Axis> & Ax,string File) {
    PrintOutput::end();
    PrintOutput::paragraph();
    PrintOutput::newRow();
    PrintOutput::rowItem("Axis");
    PrintOutput::rowItem("maxCoeff");
    PrintOutput::rowItem("from");
    PrintOutput::rowItem("  to");
    PrintOutput::rowItem("order");
    PrintOutput::rowItem("basis");
    PrintOutput::rowItem("absorption: eta,inner");
    vector<string>basisDetails;
    for(unsigned int n=0;n<Ax.size();n++){
        PrintOutput::newRow();
        PrintOutput::rowItem(Ax[n].name);
        PrintOutput::rowItem(Ax[n].maxSize());
        if(BasisNdim::factory(Ax[n].basDef[0].funcs)!=0)continue;
        if(Ax[n].lowerEnd()>-DBL_MAX/2)PrintOutput::rowItem(Ax[n].lowerEnd());
        else  PrintOutput::rowItem("-Infty|"+tools::str(Ax[n].basDef[0].upBound(),3));
        if(Ax[n].upperEnd()< DBL_MAX/2)PrintOutput::rowItem(Ax[n].upperEnd());
        else  PrintOutput::rowItem(tools::str(Ax[n].basDef[Ax[n].basDef.size()-1].lowBound(),3)+"|Infty");
        if(Ax[n].minOrder()+1>=Ax[n].maxOrder())PrintOutput::rowItem(Ax[n].maxOrder());
        else     PrintOutput::rowItem(tools::str(Ax[n].minOrder())+"-"+tools::str(Ax[n].maxOrder()));
        PrintOutput::rowItem(Ax[n].basString());
        PrintOutput::rowItem(Ax[n].comsca.str());

        if(Ax[n].basString()!=Ax[n].basDef[0].funcs)
            basisDetails.push_back(Ax[n].basString("full"));

    }
    PrintOutput::end();

    if(basisDetails.size()>0){
        PrintOutput::paragraph();
        PrintOutput::lineItem("Full basis name(s)",basisDetails[0]);
        for(size_t k=1;k<basisDetails.size();k++)
            PrintOutput::lineItem("",basisDetails[k]);
        PrintOutput::end();
    }
    PrintOutput::paragraph();
    if(File=="")return;
}

double Axis::lowerEnd() const {return basDef.front().lowBound();}
double Axis::upperEnd() const {return basDef.back().upBound();}
//double Axis::lowerRange() const {return basDef.front().lowBound()>-DBL_MAX/2?basDef.front().lowBound():basDef.front().upBound()*1.1-basDef.back().lowBound()*0.1;}
//double Axis::upperRange() const {return basDef.back().upBound()<DBL_MAX/2?basDef.back().upBound():basDef.back().lowBound()*1.1-basDef.front().upBound()*0.1;}
double Axis::lowerRange() const {return basDef.front().lowBound()>-DBL_MAX/2?basDef.front().lowBound():basDef.front().upBound()-3*std::abs(basDef.front().scale);}
double Axis::upperRange() const {return basDef.back().upBound()<DBL_MAX/2?basDef.back().upBound():basDef.back().lowBound()+3*std::abs(basDef.back().scale);}

double Axis::boxsize() const{
    DEVABORT("re-implement");
    return 0;
    //    double boxSize=lowerEnd();
    //    for(unsigned int k=0;k<basDef.size();k++){
    //        BasisSet* B = BasisSet::get(basDef[k]);
    //        if(basDef[k].funcs.find("polExp[")==0 or B->isAbsorptive())break;
    //        boxSize=B->upBound();
    //    }
    //    return boxSize;
}

// plot all basis functions on axis
void Axis::plot(string File, int Points, double QLow, double QUp) const {
    double plow=max(QLow,lowerEnd()),pup=min(QUp,upperEnd());

    // maximum number of basis functions on axis
    unsigned int maxcols=0;
    for(unsigned int n=0;n<basDef.size();n++)maxcols=max(maxcols,basDef[n].order);
    ofstream out;
    out.open(File.c_str());

    // loop through elements on axis
    for(double xi=plow;xi<=pup;xi+=(pup-plow)*(1.-1.e-15)/(max(1,Points-1)))
    {
        out<<setprecision(5);
        out<<setw(10)<<xi;
        for (unsigned int n=0;n<basDef.size();n++){
            if(basDef[n].lowBound()>xi or basDef[n].upBound()<xi)continue;
            // evaluate basis functions on this grid
            ABORT("implement this");
            UseMatrix vals;//=basDef[n].val(UseMatrix::UseMap(&cxi,1,1),true);
            // write to file
            out<<setprecision(5);
            for(unsigned int k=0;k<vals.cols();k++)out<<" "<<setw(10)<<vals(0,k).real();
            for(unsigned int k=vals.cols();k<maxcols;k++)out<<" "<<setw(10)<<0.;
            out<<endl;
        }
    }

    out.close();
}

void Axis::error(string Message){ABORT("Axis "+str()+"\n"+Message);}

