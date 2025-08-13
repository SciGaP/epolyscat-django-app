// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorGradient.h"

#include "readInput.h"
#include "printOutput.h"
#include "discretizationGrid.h"
#include "index.h"
#include "coefficients.h"
#include "basisAbstract.h"

using namespace std;

OperatorGradient::OperatorGradient(const Discretization *Parent,std::vector<string> &Ax,
                                   std::vector<unsigned int> &Points, std::vector<std::vector<double> > &Bounds)
{ 
    /// Gradient:
    /// - Ax...     names of axes as specified for basic discretization
    /// - Points..  number of (equidistant) points for each axis (=0: do not transform)
    /// - Bounds... first and last points for each axis

    // transformation to values and derivatives on grid
    grid=new DiscretizationGrid(Parent->idx(),Ax,Points,Bounds,true);

    name="grad_";
    for(unsigned int k=0;k<Ax.size();k++)name+=Ax[k]+".";
    name.resize(name.length()-1);

    jIndex=Parent->idx();

    // lump together sub-indices for values and derivatives
    Index* permIndex=new Index();
    permIndex->childAdd(new Index(*grid->mapFromParent()->iIndex));
    for(unsigned int k=0;k<grid->mapDerivative.size();k++){
        permIndex->childAdd(new Index(*grid->mapDerivative[k]->iIndex));
    }
//    permIndex->setBasis(BasisSet::getDummy(grid->mapDerivative.size()+1));
    permIndex->setBasis(BasisAbstract::factory("Vector:"+tools::str(grid->mapDerivative.size()+1)));

    // permutation: top level to bottom
    std::vector<unsigned int> perm;
    for(unsigned int k=1;k<Parent->idx()->height();k++)perm.push_back(k);
    perm.push_back(0);

    // operator left hand index with val/der index at floor
    Index* tIndex=new Index();
    permIndex->permute(perm,*tIndex,false);
    tIndex->purge(permIndex->height());
    iIndex=const_cast<const Index*>(tIndex);
    delete permIndex;

    // Coefficients to accept grid maps output
    permModel=new Coefficients();
    permModel->setIdx(new Index()); // provide idx() (for completeness)
    Coefficients * tmp=new Coefficients(grid->mapFromParent()->iIndex);
    tmp->treeOrderStorage();
    permModel->childAdd(tmp);
    for(unsigned int k=0;k<grid->mapDerivative.size();k++){
        tmp=new Coefficients(grid->mapDerivative[k]->iIndex);
        tmp->treeOrderStorage();
        permModel->childAdd(tmp);
    }
    const_cast<Index*>(permModel->idx())->setBasis(BasisAbstract::factory("Vector:"+tools::str(permModel->childSize())));

    // view on permModel in ordering of OperatorGradient::iIndex
    model=new Coefficients();
    permModel->permute(perm,*model,true);
}

void OperatorGradient::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    // map into permModel
    grid->mapFromParent()->apply(A,Vec,B,*permModel->child(0));
    for(unsigned int k=1;k<permModel->childSize();k++)
        grid->mapDerivative[k-1]->apply(1.,Vec,0.,*permModel->child(k));
    // add properly permuted view into output
    Y.axpy(A,*model,B);
}

OperatorGradient* OperatorGradient::read(const Discretization *Parent, ReadInput &Inp){

    std::vector<std::string> axis;
    std::vector<unsigned int> points;
    std::vector<std::vector<double> > bounds;

    int p;
    double lb,ub;
    string axStr;
    vector<Axis> ax=Parent->getAxis();


    Inp.read("Gradient","axis",axStr,"BLANK","which axis to convert",1);
    if(axStr!="BLANK"){
        // write header if plot-definition is given plot definition
        PrintOutput::title("GRID FOR GRADIENT");
        PrintOutput::paragraph();
        PrintOutput::newRow();
        PrintOutput::rowItem("Axis");
        PrintOutput::rowItem("nPoints");
        PrintOutput::rowItem("   from");
        PrintOutput::rowItem("     to");
    }

    int line=0;
    axStr="";
    while (axStr!="BLANK") {
        line++;

        Inp.read("Gradient","axis",axStr,"BLANK","which axis to convert",line);
        Inp.read("Gradient","points",p,"0","number of points",line);
        Inp.read("Gradient","lowerBound",lb,"0","lower grid boundary",line);
        Inp.read("Gradient","upperBound",ub,tools::str(lb),"upper grid boundary; defaults to =lower (only for single point)",line);

        if(axStr=="BLANK")break;

        // find axis number in present discretization
        int n=0;
        for(;n<ax.size();n++)if(Parent->hierarchy[n]==axStr)break;
        if(n==Parent->idx()->height()){
            string message;
            message="no coordinate = '"+axStr+"': ";
            for(unsigned int i=0;i<Parent->hierarchy.size();i++)message+=" "+tools::str(i)+"="+Parent->hierarchy[i];
            ABORT(message);
        }

        axis.push_back(axStr);
        points.push_back(p);
        bounds.push_back(vector<double>(2,lb));
        bounds.back()[1]=ub;

        PrintOutput::newRow();
        PrintOutput::rowItem(axStr);
        PrintOutput::rowItem(p);
        PrintOutput::rowItem(lb);
        PrintOutput::rowItem(ub);

        // a few sanity checks
        if(p==0){PrintOutput::paragraph();ABORT(axStr+": need at least one point on axis");}
        if(p==1 and ub!=lb){PrintOutput::paragraph();ABORT(axStr+": for single point need upperBound==lowerBound");}
        if(p>1 and ub<=lb){PrintOutput::paragraph();ABORT(axStr+": for multiple points need lowerBound < upperBound");}

    }
    PrintOutput::paragraph();

    if(axis.size()==0)return 0;
    return new OperatorGradient(Parent,axis,points,bounds);

}
