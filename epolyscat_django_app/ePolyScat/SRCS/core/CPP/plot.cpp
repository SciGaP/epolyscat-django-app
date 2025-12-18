// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "plot.h"
#include "coefficients.h"
#include "indexGrid.h"
#include "inverse.h"
#include "multiIndex.h"
#include "coefficientsFloor.h"
#include "index.h"
#include "readInput.h"
#include "printOutput.h"
#include "asciiFile.h"
#include "discretizationGrid.h"
#include "algebra.h"
#include "sphericalHarmonicReal.h"
#include "operatorMap.h"
#include "basisGrid.h"

#include "log.h"
#include "treeEmpty.h"
#include "operatorDefinition.h"
#include "basisMat1D.h"
#include "algebra.h"
#include "operatorMap.h"
#include "indexGridHybrid.h"
#include "mapGridHybrid.h"

using namespace std;
using namespace tools;

Plot::~Plot(){
    //    delete plotStore;
}

Plot::Plot():_realValues(false),discIndex(0),gIndex(0),gMap(0),nDigits(8),_append(false),plotStore(0){}

Plot::Plot(const Index *Idx, vector<std::string> Axis, vector<string> Use,vector<unsigned int> GridPoints,vector<vector<double> > Bounds)
    :_realValues(false),discIndex(Idx),nDigits(8),_append(false),plotStore(0)
{
    if(GridPoints.size()==0)GridPoints.resize(Axis.size(),0);
    if(Bounds.size()==0)Bounds.resize(Axis.size(),{0.,0.});
    std::vector<std::vector<double>> g,w;
//    DiscretizationGrid::gridWeight(g,w,Axis,GridPoints,Bounds);
    IndexGrid::gridWeight(g,w,Axis,GridPoints,Bounds);
    construct(Axis,Use,g,w);
}
Plot::Plot(const Index *Idx, const PlotKind* Kind)
    :_realValues(false),discIndex(Idx),nDigits(8),_append(false),plotStore(0)
{
    construct(Kind->axes(),Kind->use(),Kind->grid(),Kind->weig());
}

double Plot::integral() const{
    std::string ovrDef;
    for(const Index* ix=plotStore->idx();ix;ix=ix->descend()){
        if(ix->basis()->grid())ovrDef+="<GridWeight>";
        else if(ix->axisName()!="NONE")ovrDef+="<1>";
    }
    Coefficients one(plotStore->idx(),1.);
    OperatorTree ovr("overlap",ovrDef,plotStore->idx(),plotStore->idx());
    std::complex<double> integ=ovr.matrixElementUnscaled(one,*plotStore);
    if(std::abs(integ.imag())>std::abs(integ.real())*1.e-7)
        PrintOutput::warning(Sstr+"found non-zero imaginary part in integral of spectrum - set=0: "+integ);
    return integ.real();
}

Plot::Plot(const Index *Idx, ReadInput &Inp, bool WarnIfEmpty)
    :_realValues(false),discIndex(Idx),_append(false),plotStore(0){
    /// Plot: axis, points, lowerBounds, upperBounds, usage
    /// - axis...       name of axis, e.g., R,X, R1,Eta, Phi,Phi1, Phi2 (see class Axis)
    /// - points..      number of (equidistant) points for plot (=0: do not transform to grid)
    /// - lowerBound... start plot
    /// - upperBound... end plot
    /// - usage...      grid/separate/sum the axis points

    vector<string>gUsage,gAxis;
    vector<unsigned int> gPoints;
    vector<vector<double> > gBounds;

    Inp.obsolete("Plot","currents","dipole densities will automatically be included with density plots (if defined)",1);

    Inp.read("Plot","onlyMatching",_onlyMatching,"false","true = abort rather than skip if plot axes do not match index");
    Inp.read("Plot","append",_append,"false","append to file");
    Inp.read("Plot","tagMin",_tagMin,"-"+Algebra::Infty,"plot only tags >= tagMin");
    Inp.read("Plot","tagMax",_tagMax,Algebra::Infty,"plot only tags <= tagMax");
    Inp.read("Plot","interval",_tagInterval,tools::str(_tagIntervalDefault),
             "minimal interval for plots (<=0 - all w/o duplicates, <0 - can be overwritten by setDefaultInterval())");
    unsigned int nDig;
    Inp.read("Plot","digits",nDig,"8","number of digits for plot values");
    digits(nDig);

    if(_tagInterval>0.){
        if(Units::isDefined("OptCyc"))
            _tagInterval=Units::convert(_tagInterval,"DEFAULT_SYSTEM","OptCyc");
    }

    PlotKind kind("",Inp);
    if(kind.axes().size()==0){
        if(WarnIfEmpty)PrintOutput::warning("empty plot definition on file "+Inp.file()+" - see docu in "+Inp.docFile());
        return;
    }

    construct(kind.axes(),kind.use(),kind.grid(),kind.weig());
}

void Plot::construct(vector<std::string> Axis, vector<string> Use,
                     std::vector<std::vector<double>> Grid, std::vector<std::vector<double>> Weig)
{

    std::vector<std::string> hier=tools::splitString(discIndex->hierarchy(),'.');
    for(auto a: Axis){
        if(std::find(hier.begin(),hier.end(),a)==hier.end()){
            std::string mess="No plot generated: Axis "+a+" not in hierarchy "+discIndex->hierarchy();
            if(_onlyMatching)PrintOutput::message(mess);
            else PrintOutput::warning(mess,50,0,"specify Plot: onlyMatching=true to convert warning to message" );
            return;
        }

    }
    // axes to be converted to grid
    vector<string> gAxis,jAxis;
    vector<vector<double> > gGrid,gWeig;
    if(Use.size()!=Axis.size() or Grid.size()!=Axis.size() or Weig.size()!=Axis.size())
        ABORT("must specify Use, Grid, Weig for each Axis");
    for(unsigned int k=0;k<Axis.size();k++){
        // Axis to convert to grid (check IndexGrid for usage rules for gGrid and gWeig)
        if(Use[k]=="g" or Grid[k].size()!=0){
            gAxis.push_back(Axis[k]);
            gGrid.push_back(Grid[k]);
            gWeig.push_back(Weig[k]);
        }
    }

    Index::build=true;
    if(discIndex->isHybrid() and discIndex->hierarchy().find("Ndim")==std::string::npos){
        IndexGridHybrid* gHyb=new IndexGridHybrid(discIndex,gAxis,gGrid);
        gIndex.reset(gHyb);
        gMap.reset(new MapGridHybrid(gHyb,discIndex));
    } else {
        gIndex.reset(new IndexGrid(discIndex,gAxis,gGrid,gWeig));
        gMap.reset(new OperatorMap(gIndex.get(),discIndex));
    }
    gvec.reset(new Coefficients(gIndex.get()));
    lvec.reset(new Coefficients(gIndex.get()));
    iFull.reset(new Index(*gIndex));

    // need to expand for permutations
    iFull->bottomExpandAll();
    gvec->idx()->bottomExpandAll();
    iFull->resetFloor(iFull->firstLeaf()->depth());

    // determine permutation such that top is "s" (summation), "p" (separate data colums), lowest "g"
    vector<unsigned int > plotUse(iFull->heightAboveFloor(),0); // default is 0=sum over index
    for(unsigned int k=0;k<Axis.size();k++){
        const Index* lev;
        if(0==(lev=iFull->axisIndex(Axis[k]))){
            ABORT("could not find axis in Index: "+Axis[k]+" not in hierarchy "+iFull->hierarchy());
        }
        //   0- 99: summation
        // 100-199: separate columns
        // 200-299: actual grid
        unsigned int kAx=lev->depth();
        // with 'p' keep axes in the order of the discretization (there may be a necessary hierarchy)
        if(Use[k]=="p")     plotUse[kAx]=100+2*kAx;
        // within 'g', use sorting of Axes[k] (i.e. 2*k, not 2*kAx)
        else if(Use[k]=="g")plotUse[kAx]=200+2*k;

        // place level corresponding to continuity axis just below continuity level
        if(0!=(lev=lev->descend()->axisIndex(Axis[k])))plotUse[lev->depth()]=plotUse[kAx]+1;
    }

    unsigned int sumDepth=0;
    plotDepth=0;
    for(unsigned int k=0;k<plotUse.size();k++){
        if     (plotUse[k]<100)sumDepth++;
        else if(plotUse[k]<200)plotDepth++;
    }

    vector<unsigned int> perm;
    for(unsigned int k=0;k<plotUse.size();k++)perm.push_back(k);
    tools::sortByKey(plotUse,perm);

    // what is this supposed to do?
    while(plotUse.front()<100)plotUse.erase(plotUse.begin());

    // get re-sorted index for plotting: grid axes at deepest level
    // note: this should be cast into a constructor
    iPlot.reset(new Index());
    iFull->permute(perm,*iPlot);
    iPlot->sizeCompute();
    // create a view on gvec with floor lowered
    LOG_PUSH("rightView");

    // left and right view down to single grid values, and desired permutations
    rightView.reset(new Coefficients());
    rightFullView=new Coefficients(iFull.get(),gvec.get());
    rightFullView->permute(perm,*rightView,true);
    leftView.reset(new Coefficients());
    leftFullView=new Coefficients(iFull.get(),lvec.get());
    leftFullView->permute(perm,*leftView,true);
    LOG_POP();

    // set up plot storage
    // WARNING: we assume that all levels for plotting are equivalent
    //          as of now, this is not checked
    iStor.reset(new Index(*iPlot->descend(sumDepth)));
    plotStore=new Coefficients(iStor.get());
    plotStore->idx()->descend(plotDepth)->sizeCompute();
    plotStore->treeOrderStorage();

    generateAxes(plotStore->idx()->descend(plotDepth),axCols);

    // axis coordinate columns
    LOG_PUSH("plotStore");
    Coefficients* p=plotStore->descend(plotDepth);
    while(p->height()>0){
        if(p->idx()->depthOfDuplicate()==Index::npos)columnHeader+=p->idx()->axisName()+", ";
        p=p->descend();
    }
    columnHeader=columnHeader.substr(0,columnHeader.rfind(","))+":";
    LOG_POP();

    // append: repeated tags will be written into an axis column
    if(_append)append(_append);

    if(axCols.size()>2)
        PrintOutput::warning(Str("there are")+axCols.size()+"axes ("+columnHeader+") - plotting may be ambigous");

    // axes that are summed over
    Index * l=iPlot.get();
    string sum;
    while(l->depth()<sumDepth){
        if(l->depthOfDuplicate()==Index::npos)sum+=l->axisName()+",";
        l=l->descend();
    }
    if(sum.length()>0)columnHeader+=" sum["+sum.substr(0,sum.length()-1)+"]";

    // data columns
    string col;
    p=plotStore;
    while(p->depth()<plotDepth){
        // avoid double-counting of continuity levels
        if(p->idx()->depthOfDuplicate()==Index::npos)col+=p->idx()->axisName()+",";
        p=p->descend();
    }
    if(col.length()>0){
        col.resize(col.length()-1);
        col+=") =";
        while(p!=0){
            vector<unsigned int> pi(p->index());
            col+=" (";
            string pars=") ";
            const Index * up=p->idx()->parent();
            for(int k=pi.size();k>0;k--)
            {
                pars=","+tools::str(up->basis()->physical(pi[k-1]),3)+pars;
                up=up->parent();
            }
            col+=pars.substr(1);
            p=p->nodeRight(plotStore);
        }

        columnHeader+=" ("+col;
    }
    else {
        _header.push_back("");
        columnHeader+=" density";
    }

    // overlap/quadrature for axes that will be summed over
    OperatorDefinition opDef;
    std::vector<std::string> idxAx=tools::splitString(discIndex->hierarchy(),'.');
    if(idxAx.back()=="NONE")idxAx.pop_back();
    for(size_t k=0;k<idxAx.size();k++){
        const Index* ix=gvec->idx()->descend(k);
        if(ix->axisName().find("spec")==0)
            continue;
        else if(std::find(idxAx.begin()+k+1,idxAx.end(),ix->axisName())!=idxAx.end())
            continue; // FE=axis
        else if(std::find(Axis.begin(),Axis.end(),ix->axisName())!=Axis.end() and
                Use[std::find(Axis.begin(),Axis.end(),ix->axisName())-Axis.begin()]!="s")
            opDef+="<Id>"; // not a sum-axis
        else if(ix->basis()->integrable() )
            opDef+="<1>";
        else if(ix->basis()->grid())
            opDef+="<GridWeight>";
        else
            opDef+="<Id>"; // direct sum
    }
    _densityOp.reset(new OperatorTree("Quadratures",opDef,gvec->idx(),gvec->idx()));
    if(_densityOp->isIdentity())_densityOp=0;

    const Index*ix;
    for(ix=iFull.get();not ix->isBottom();ix=ix->descend());
    for(;ix!=0 and ix->isBottom();ix=ix->nodeRight())ix->bottomUnexpand();
    for(ix=gvec->idx();not ix->isBottom();ix=ix->descend());
    for(;ix!=0 and ix->isBottom();ix=ix->nodeRight())ix->bottomUnexpand();
}

void Plot::addOperator(const OperatorAbstract * Op){
    _ops.push_back(Op);
    _header.push_back(Op->name+"...."+Op->def());
    columnHeader+=",   "+Op->name;
}

void Plot::generateAxes(const Index *Idx,vector<vector<double> > &Cols,unsigned int IAx) const{

    if(Cols.size()==0)Idx->sizeCompute();
    if(not Idx->isLeaf() and Idx->basis()->isGrid()){
        if(Cols.size()<IAx+1)Cols.resize(IAx+1,vector<double>());
        for(unsigned int k=0;k<Idx->childSize();k++)
            Cols[IAx].resize(Cols[IAx].size()+Idx->child(k)->sizeStored(),Idx->basis()->grid()->mesh()[k]);
        IAx++;
    }
    for(unsigned int k=0;k<Idx->childSize();k++)generateAxes(Idx->child(k),Cols,IAx);
    if(Idx->parent()==0 and Cols.size()==0)PrintOutput::DEVwarning(Sstr+"no axes - empty plot");
}

double Plot::sumLeftConjRight(Coefficients *LVec,Coefficients *RVec, Coefficients &Store) const{

    double norm(0.);
    for(unsigned int k=0;k<LVec->childSize();k++)
        if(LVec->height()>Store.height())norm+=sumLeftConjRight(LVec->child(k),RVec->child(k),Store);
        else                             norm+=sumLeftConjRight(LVec->child(k),RVec->child(k),*Store.child(k));

    if(Store.isLeaf()){
        for(unsigned int k=0;k<Store.idx()->sizeCompute();k++){
            double dens=(conj(LVec->floorData()[k])*RVec->floorData()[k]).real();
            Store.floorData()[k]+=dens;
            norm+=dens;
        }
    }
    return norm;
}

void Plot::intoColumns(const Coefficients &Store, std::vector<std::vector<double> > &Cols) const{
    // new data column for this index
    if(Store.depth()==plotDepth)Cols.push_back(vector<double>());
    if(Store.isLeaf()){
        // append data to current column
        for(complex<double>* a=Store.floorData();a!=Store.floorData()+Store.idx()->sizeCompute();a++)
            Cols.back().push_back(a->real());
    } else {
        // descend further
        for(unsigned int k=0;k<Store.childSize();k++)intoColumns(*Store.child(k),Cols);
    }
}

string Plot::str() const {

    string s=discIndex->hierarchy();

    // which axes will be plotted
    s+="->axes: "+iPlot->hierarchy();
    return s;
}

void Plot::print() const {
    if(plotStore==0)return;

    // write header if plot-definition is given plot definition
    PrintOutput::title("GRID AND PLOT");
    PrintOutput::paragraph();
    if(_tagInterval>=0.)
    {
        string prtInt="all";
        if(_tagInterval>0.){
            prtInt=tools::str(_tagInterval,3);
            if(Units::isDefined("OptCyc"))prtInt+="(OptCyc)";
        }
        PrintOutput::lineItem("interval",prtInt);
    }
    PrintOutput::lineItem("digits",tools::str(nDigits));
    PrintOutput::paragraph();
    PrintOutput::paragraph();
    PrintOutput::newRow();
    PrintOutput::rowItem("Axis");
    PrintOutput::rowItem("  usage");
    PrintOutput::rowItem(" points");
    PrintOutput::rowItem("   from");
    PrintOutput::rowItem("     to");
    PrintOutput::newRow();

    Coefficients* c=plotStore;
    do {
        int siz=c->idx()->basis()->size();
        if(not c->isLeaf() and c->idx()->axisName()==c->child(0)->idx()->axisName()){
            siz=0;
            for(size_t k=0;k<c->childSize();k++)siz+=c->child(k)->idx()->basis()->size();
            c=c->child(0);
        }
        if(siz>1){
            PrintOutput::rowItem(c->idx()->axisName());
            string usage="separate";
            if(c->depth()>=plotDepth)usage="grid";
            PrintOutput::rowItem(usage);
            PrintOutput::rowItem(siz);
            if(c->idx()->basis()->grid()){
                PrintOutput::rowItem(c->idx()->basis()->grid()->lowBound());
                PrintOutput::rowItem(c->idx()->basis()->grid()->upBound());
            }
            PrintOutput::newRow();
        }
        if(c->isLeaf())break;
        c=c->child(0);
    } while(true);

}

void Plot::plotJAD(const Coefficients *C, double E1, double E2, string fileName, unsigned int thetaNum){
    unsigned int k1Idx = 0, k2Idx = 0;
    for(unsigned int i = 0; i < C->firstFloor()->idx()->basis()->grid()->mesh().size(); i++){
        if(C->firstFloor()->idx()->basis()->grid()->mesh()[i] < sqrt(2 * E1))k1Idx = i;
        if(C->firstFloor()->idx()->basis()->grid()->mesh()[i] < sqrt(2 * E2))k2Idx = i;
    }
    SphericalHarmonicReal spherHarm;
    // get the depth of the l1, l2, m1, m2 level, usually, m1 = -m2 and the maximum size
    unsigned int m1Depth=0, m2Depth=0, l1Depth=0, l2Depth=0;
    size_t m1Max(0),l1Max(0),m2Max(0),l2Max(0);
    unsigned int i = 0;
    while(!C->descend(i)->isLeaf()){
        if(C->descend(i)->idx()->axisName() == "Phi1"){
            m1Depth = i;
            m1Max = C->descend(i)->childSize();
            const Coefficients * test = C->descend(i);
            while(test != 0){
                if(m1Max < test->childSize())
                    m1Max = test->childSize();
                test = test->nodeRight();
            }
        }else if(C->descend(i)->idx()->axisName() == "Eta1"){
            l1Depth = i;
            l1Max = C->descend(i)->childSize();
            const Coefficients * test = C->descend(i);
            while(test != 0){
                if(l1Max < test->childSize())
                    l1Max = test->childSize();
                test = test->nodeRight();
            }
        }else if(C->descend(i)->idx()->axisName() == "Phi2"){
            m2Depth = i;
            m2Max = C->descend(i)->childSize();
            const Coefficients * test = C->descend(i);
            while(test != 0){
                if(m2Max < test->childSize())
                    m2Max = test->childSize();
                test = test->nodeRight();
            }
        }else if(C->descend(i)->idx()->axisName() == "Eta2"){
            l2Depth = i;
            l2Max = C->descend(i)->childSize();
            const Coefficients * test = C->descend(i);
            while(test != 0){
                if(l2Max < test->childSize())
                    l2Max = test->childSize();
                test = test->nodeRight();
            }
        }
        i++;
    }
    if(m1Depth > l1Depth or m2Depth > l2Depth)ABORT("The phi axsis is lower than eta axsis, please change it");
    if(m1Depth > m2Depth)ABORT("The phi1 is lower than phi2, please permute it");
    if(m1Max != m2Max and l1Max != l2Max)ABORT("The code currently only deals with m1 = -m2, l1max=l2max case");
    //Collect all the scattering amplitute
    vector<UseMatrix> ampl;

    for(size_t m1 = 0; m1 < m1Max; m1++){
        UseMatrix mat(l1Max - abs(int(m1-m1Max)/ 2), l1Max - abs(int(m1-m1Max)/ 2));
        Coefficients * test = C->descend(m1Depth)->child(m1)->firstLeaf();
        unsigned int Nk = test->idx()->childSize();
        while(test){
            mat(test->index()[l1Depth], test->index()[l2Depth]) = *(test->data() + k1Idx * Nk + k2Idx);
            test = test->nodeRight();
        }
        ampl.push_back(mat);
    }
    // The hierachy of Ylm vector m1(-m2), l1, l2
    double dTheta = math::pi / (thetaNum - 1);
    vector<vector<UseMatrix> > Ylm;
    Ylm.resize(m1Max);
    for(size_t m1 = 0; m1 < m1Max; m1++){
        for(unsigned int i = 0; i < 4; i++){
            UseMatrix mat(l1Max - std::abs(int(m1 - m1Max / 2)), thetaNum - 1);
            for(unsigned int row = 0; row < mat.rows(); row++)
                for(unsigned int col = 0; col < mat.cols(); col++)
                    //i =:  0 ==> 1 puls, 1 ==> 2 puls, 2 ==> 1 minus, 3 ==> 2 minus
                    if(i == 0)mat(row, col) = spherHarm(row + std::abs(int(m1 - m1Max / 2)), m1 - m1Max / 2, dTheta * col, 0.);
                    else if(i == 1)mat(row, col) = spherHarm(row + std::abs(int(m1 - m1Max / 2)), m1 - m1Max / 2, dTheta * col, math::pi);
                    else if(i == 2)mat(row, col) = spherHarm(row + std::abs(int(m1 - m1Max / 2)), -(m1 - m1Max / 2), dTheta * col, 0.);
                    else mat(row, col) = spherHarm(row + std::abs(int(m1 - m1Max / 2)), -(m1 - m1Max / 2), dTheta * col, math::pi);
            Ylm[m1].push_back(mat);
            if(i == 0)
                mat.print();
        }
    }
    UseMatrix totalJAD(2 * thetaNum - 3, 2 * thetaNum - 3);
    for(unsigned int m1 = 0; m1 < m1Max; m1++){
        UseMatrix temJAD(2 * thetaNum - 3, 2 * thetaNum - 3);
        temJAD.block(0, 0, thetaNum - 1, thetaNum - 1) = Ylm[m1][0].transpose() * ampl[m1] *Ylm[m1][0];
        UseMatrix tem2 = Ylm[m1][1].transpose() * ampl[m1] * Ylm[m1][0];
        for(unsigned int i = thetaNum - 1; i < 2 * thetaNum - 3; i++)
            temJAD.block(i, 0, 1, thetaNum - 1) = tem2.row(2 * thetaNum - 3 - i);
        tem2 = Ylm[m1][0].transpose() * ampl[m1] * Ylm[m1][3];
        for(unsigned int i = thetaNum - 1; i < 2 * thetaNum - 3; i++)
            temJAD.block(0, i, thetaNum - 1, 1) = tem2.col(2 * thetaNum - 3 - i);
        tem2 = Ylm[m1][1].transpose() * ampl[m1] * Ylm[m1][3];
        for(unsigned int i = thetaNum - 1; i < 2 * thetaNum - 3; i++)
            for(unsigned int j = thetaNum - 1; j < 2 * thetaNum - 3; j++)
                temJAD(i,j) = tem2(2 * thetaNum - 3 - i, 2 * thetaNum - 3 - j);
        totalJAD += temJAD;
    }
    ofstream JAD(fileName);
    for(unsigned int i = 0; i < totalJAD.rows(); i++){
        for(unsigned int j = 0; j < totalJAD.cols(); j++)
            JAD << (totalJAD(i,j).real() * totalJAD(i,j).real() + totalJAD(i,j).imag() * totalJAD(i,j).imag())<< " ";
        JAD<< endl;
    }
    JAD.close();
}

const Plot &Plot::generate(const Coefficients &C, const OperatorAbstract * Op) const {

    Coefficients PlotC(plotStore->idx());
    PlotC.treeOrderStorage();
    gMap->apply(1,C,0.,*gvec);
    *lvec=*gvec;

    if(_realValues){
        double imMax=0,reMax=0;
        for(size_t k=0;k<lvec->size();k++){
            imMax=max(imMax,std::abs(lvec->data()[k].imag()));
            reMax=max(reMax,std::abs(lvec->data()[k].real()));
            lvec->data()[k]=1.;
        }
        if(imMax/(max(reMax,1.))>1.e-12)
            PrintOutput::warning(Sstr+"Asking to plot real, but input has imMax/reMax="+(imMax/(max(reMax,1.))));
    }

    const OperatorAbstract* op=Op;
    if(op){
        // apply Op and convert to grid
        op->apply(1.,C,0.,*op->tempLHS());
        gMap->apply(1.,*op->tempLHS(),0.,*gvec);
    }
    else{
        // apply density operators
        _densityOp->apply(1,*lvec,0.,*gvec);
    }

    // summation levels
    PlotC.setToZero();
    sumLeftConjRight(leftView.get(),rightView.get(),PlotC);
    *const_cast<Plot*>(this)->plotStore=PlotC; // there is a design inconsistency
    return *this;
}

void Plot::plot(const Coefficients &Cin, const std::string &File, const vector<string> &Head, const string Tag, bool OverWrite) const {

    //temporary: for transition to new hybrid grid
    bool oldHyb=Cin.idx()->isHybrid() and gMap and 0==dynamic_cast<const IndexGridHybrid*>(gMap->idx());
    const Coefficients *C=oldHyb?Cin.child(1):&Cin;

    if(dimension()==0)return;
    if(not acceptPlot(Tag))return;

    double inf=DBL_MAX; // here should be the general "infinity"
    if(tagMin()>-inf/2 or tagMax()<inf/2){
        double dTag=tools::string_to_double(Tag);
        if(dTag<tagMin()-1.e-10*max(1.,abs(tagMin())) or
                tagMax()+1.e-10*max(1.,abs(tagMax()))< dTag)return;
    }

    if(C->idx()!=discIndex)DEVABORT("wrong plot: indices do not match");

    // separation of files by tag
    string file=File;
    std::vector<std::vector<double>> dum=sum(*C,file,1.,Tag);

    // write to (new) file
    write(file,Head,OverWrite);

    const_cast<Plot*>(this)->clear(); // plot() erases after write
}

void Plot::append(bool Append){
    if(Append and not _append){
        axCols.insert(axCols.begin(),vector<double>(axCols[0].size(),0.));
        columnHeader="tag, "+columnHeader;
        _append=Append;
    }
    else if(_append and not Append)
        PrintOutput::warning("append mode for plot has been set earlier, cannot switch off");
}

vector<vector<double> > Plot::sum(const Coefficients &C, std::string &File, double Scale, const std::string & Tag) const
{
    vector<vector<double> > cols(axCols);
    for(size_t k=0;k<_ops.size()+1;k++){
        if(k>0)generate(C,_ops[k-1]);
        else   generate(C);
        // write into new columns
        intoColumns(*plotStore,cols);
    }

    double tagVal=0.;
    if(_append){
        if(Tag=="")ABORT("for appending to plot file, call sum with Tag!=\"\"");
        tagVal=tools::string_to_double(Tag);
        cols[0]=vector<double>(cols[0].size(),tagVal);
    }
    else if(Tag!="")
        PrintOutput::DEVwarning("Tag specified, but not append mode - Tag will be ignored");

    if(_cols.size()==0){
        // if new Tag, set to 0
        _cols=vector<vector<double> >(cols.begin(),cols.begin()+dimension());
        _cols.resize(cols.size(),vector<double>(cols[0].size(),0.));
    }

    if(_append and _cols[0].back()!=tagVal){
        // new tag - need to extend columns
        _cols[0].resize(_cols[0].size()+cols[0].size(),tagVal);
        _cols[1].insert(_cols[1].end(),_cols[1].begin(),_cols[1].begin()+cols[1].size());
        for(size_t k=2;k<_cols.size();k++)_cols[k].resize(_cols[k].size()+cols[0].size(),0.);
    }

    int row0=_cols[0].size()-cols[0].size();
    for(size_t k=dimension();k<cols.size();k++)
        for(size_t l=0;l<cols[k].size();l++)
            _cols[k][row0+l]+=cols[k][l]*Scale;
    return _cols;
}

void Plot::write(const string &File, const vector<string> &Head, bool OverWrite) const{
    unsigned int blank=_cols.size();
    if(axCols.size()>1)blank=axCols.size()-2;

    AsciiFile plotFile(File,blank,", ");
    if(plotFile.empty() or OverWrite or _overWrite){
        vector<string>comm(Head);
        comm.insert(comm.begin()+Head.size(),_header.begin(),_header.end());
        comm.push_back("");
        comm.push_back(columnHeader);
        plotFile.writeComments(comm);
    }
    else
        plotFile.writeBlankRow();
    plotFile.writeCols(_cols,nDigits);
}

bool Plot::plotable(const Index *Idx, std::vector<string> Axes, std::vector<unsigned int> Points){
    for(size_t k=0;k<Axes.size();k++){
        const Index* iAx=Idx->axisIndex(Axes[k]);
        if(iAx==0)ABORT("plot axis "+Axes[k]+" not in index "+Idx->hierarchy());

        if(not iAx->basis()->isGrid())return false;
        if(Points.size()>0 and Points[k]>0){
            if(iAx->basis()->size()!=Points[k])
                ABORT(Str("axis")+Axes[k]+"with size"+iAx->basis()->size()+"grid does not match desired number of"
                      +Points[k]+"points:"+iAx->basis()->str());
        }
    }
    return true;
}
