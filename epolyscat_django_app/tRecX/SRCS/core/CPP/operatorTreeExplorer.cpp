// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "operatorTreeExplorer.h"

#ifdef _OPERATOR_TREE_EXPLORER_

#include <string>
#include <iostream>
#include <sstream>

#include "tools.h"
#include "str.h"
#include "readInput.h"

#include "hMatrix.h"
#include "operatorSVD.h"
#include "operatorTucker.h"


/*
 * Static helpers
 */

static void printIndex(const Index* idx, ViNCurses::Buffer& target){
    std::vector<const Index*> path=idx->path();
    path.push_back(idx);

    if(path.size()==1){
        target()<<"ROOT";
    }

    for(unsigned int i=0; i<path.size()-1; i++){
        if(i!=0) target()<<"-";
        target()<<path[i]->axisName()<<"("<<path[i+1]->physical()<<")";
    }

    if(idx->hasFloor()){
        target()<<" F";
    }
}


static void average(UseMatrix& src, UseMatrix*& target, unsigned int iWidth, unsigned int jWidth){
    target = new UseMatrix(src.rows()/iWidth, src.cols()/jWidth);
    for(unsigned int i=0; i<target->rows(); i++){
        for(unsigned int j=0; j<target->cols(); j++){
            *(target->data(i,j))=0;
            for(unsigned int n=0; n<iWidth; n++){
                for(unsigned int m=0; m<jWidth; m++){
                    *(target->data(i,j))+=src(i*iWidth+n,j*jWidth+m)/(iWidth*jWidth);
                }
            }
        }
    }
}

static const OperatorTree* findChildRec(const OperatorTree* optree, const Index* iIndex, const Index* jIndex){
    for(unsigned int i=0; i<optree->childSize(); i++){
        if(optree->child(i)->iIndex==iIndex and optree->child(i)->jIndex==jIndex) return optree->child(i);
        if(         (optree->child(i)->iIndex==optree->iIndex or optree->child(i)->iIndex==iIndex) 
                and (optree->child(i)->jIndex==optree->jIndex or optree->child(i)->jIndex==jIndex)){
            return findChildRec(optree->child(i),iIndex,jIndex);
        }
    }

    return 0;
}

/*
 * Basic stuff and helpers
 */

UseMatrix* OperatorTreeExplorer::matrix(const OperatorAbstract* op){ 
    if(mat.find(op)==mat.end()){

        //TODO: Possible memeory leak
        std::vector<std::complex<double> >* m=new std::vector<std::complex<double> >();
        status("Creating matrix... ");
        op->matrix(*m);
        status("Creating matrix... done");
        UseMatrix* m_ = new UseMatrix(UseMatrix::UseMap(m->data(), optree->iIndex->sizeCompute(), optree->jIndex->sizeCompute()));
        mat[op]=m_;
        return m_;
    }else{
        return mat[op];
    }
}


void OperatorTreeExplorer::descend(int nChild){

    if(nChild>=optree->childSize()){
        status("Index too large!");
        return;
    }

    traverse(optree->child(nChild));
}

void OperatorTreeExplorer::descend(int nI, int nJ){

    if(optree->iIndex->childSize()<=nI or optree->jIndex->childSize()<=nJ){
        status("Indices too large!");
        return;
    }

    const OperatorTree* _optree = findChildRec(optree, optree->iIndex->child(nI), optree->jIndex->child(nJ));
    if(_optree!=0){
        traverse(_optree);
    }else{
        status("Child not found!");
    }
}

void OperatorTreeExplorer::ascend(){

    if(optree->parent()==0){
        status("Can't ascend further!");
        return;
    }

    traverse(optree->parent());
}

void OperatorTreeExplorer::traverse(const OperatorTree* _optree){
    optree=_optree;
    w_main->onTraverse();
    w_nav->onTraverse();
    w_appc->onTraverse();
    w_sing->onTraverse();
    w_play->onTraverse();
}

/*
 * ViNCurses stuff
 */
void OperatorTreeExplorer::init_windows(){
    w_main=new MainWindow();
    w_nav=new NavigatorWindow();
    w_appc=new ApplicationCostWindow();
    w_sing=new SVDWindow();
    w_play=new PlaygroundWindow();

    add_window(w_main);
    add_window(w_nav, w_main->box()->split('L', 0.33));
}

bool OperatorTreeExplorer::command(std::string command, bool before_windows){
    if(command=="a"){
        ascend();
        return true;
    }else if(command==":appc"){
        if(w_appc->assigned()) remove_window(w_appc);
        else add_window(w_appc, w_nav->box()->split('J', 0.2));
        return true;
    }else if(command==":sing"){
        if(w_sing->assigned()) remove_window(w_sing);
        else add_window(w_sing, w_main->box()->split('J', 0.2));
        return true;
    }else if(command==":playground"){
        if(w_play->assigned()) remove_window(w_play);
        else add_window(w_play, w_main->box()->split('J', 0.5));
        return true;
    }
    
    return false;
}

/*
 * Windows
 */
void OperatorTreeExplorer::MainWindow::render(){
    buffer.clear();
    const OperatorTree* optree = dynamic_cast<OperatorTreeExplorer*>(parent())->getOptree();

    if(state=='i'){
        buffer(0,0)<<"OPTREE: "<<optree->name;
        buffer(1,0)<<"definition: "<<optree->getDefinition();
        buffer(2,0)<<"dim: "<<optree->iIndex->sizeCompute()<<"x"<<optree->jIndex->sizeCompute();

        buffer(4,0)<<"children: "<<optree->childSize();
        buffer(5,0)<<"floor: ";
        if(optree->floor()==0) buffer()<<"0";
        else{ 
            buffer()<<optree->floor()->strInfo();
            // buffer()<<", tFac: "<<optree->floor()->timeDepFac;
        }

        buffer(6,0)<<"blockdiagonal: "<<optree->isBlockDiagonal();

    }else if(state=='s' or state=='m'){
        UseMatrix* mat = dynamic_cast<OperatorTreeExplorer*>(parent())->getMatrix();
        int iWidth = std::max((int)std::ceil(double(optree->iIndex->sizeStored())/100.), 1);
        int jWidth = std::max((int)std::ceil(double(optree->jIndex->sizeStored())/100.), 1);
        UseMatrix* mat_=mat;
        if(iWidth>1 || jWidth>1) average(*mat_,mat_,iWidth,jWidth);

        buffer.stream(&std::cout, 0,0);

        std::cout<<"Averaging over "<<iWidth<<"x"<<jWidth<<" blocks"<<std::endl;
        if(state=='s') mat_->show();
        else if(state=='m') mat_->print();

        if(mat_!=mat) delete mat_;

    }
}

bool OperatorTreeExplorer::MainWindow::command(std::string command){
    command_move_buffer(command);
    
    if(command=="s"){
        state='s';
        parent()->status("Structure");
    }else if(command=="m"){
        state='m';
        parent()->status("Matrix");
    }else if(command=="i"){
        state='i';
        parent()->status("Info");
    }else{
        return false;
    }

    stale();
    return true;
}

void OperatorTreeExplorer::MainWindow::onTraverse(){
    set_offset(0,0);
    stale();
}

void OperatorTreeExplorer::NavigatorWindow::render(){ 
    buffer.clear();

    const OperatorTree* optree = dynamic_cast<OperatorTreeExplorer*>(parent())->getOptree();
    buffer(0,0)<<"I: ";
    printIndex(optree->iIndex, buffer);
    buffer()<<", J: ";
    printIndex(optree->jIndex, buffer);
    for(unsigned int i=0; i<optree->childSize(); i++){
        if(i==selectedChild) buffer(5+i,4,VIN_SELECTION);
        else buffer(5+i,4);
        
        buffer()<<"+I: ";
        printIndex(optree->child(i)->iIndex, buffer);
        buffer()<<", J: ";
        printIndex(optree->child(i)->jIndex, buffer);
    }

}

bool OperatorTreeExplorer::NavigatorWindow::command(std::string command){ 
    command_move_buffer(command);
    
    if(command=="j"){
        if(dynamic_cast<OperatorTreeExplorer*>(parent())->getOptree()->childSize()>selectedChild+1) selectedChild++;
    }else if(command=="k"){
        if(selectedChild>0) selectedChild--;
    }else if(command=="\n"){
        dynamic_cast<OperatorTreeExplorer*>(parent())->descend(selectedChild);
    }else{
        return false;
    }

    stale();
    return true;
}

void OperatorTreeExplorer::NavigatorWindow::onTraverse(){
    set_offset(0,0);
    selectedChild=0;
    stale();
}

void OperatorTreeExplorer::ApplicationCostWindow::render(){
    buffer.clear();

    const OperatorTree* optree=dynamic_cast<OperatorTreeExplorer*>(parent())->getOptree();

    buffer()<<"COST: "<<100*optree->applicationCost(); //<<", THEO: "<<optree->applyCount();

    buffer(1,0)<<"i\\j";

    double appc[optree->iIndex->childSize()*optree->jIndex->childSize()];
    double total=0.;
    int count=0;

    for(unsigned int i=0; i<optree->iIndex->childSize(); i++){
        for(unsigned int j=0; j<optree->jIndex->childSize(); j++){
            const OperatorTree* child = findChildRec(optree, optree->iIndex->child(i), optree->jIndex->child(j));
            if(child!=0) appc[i*optree->jIndex->childSize()+j]=child->applicationCost();
            else         appc[i*optree->jIndex->childSize()+j]=0;

            if(child!=0) count++;
            total+=      appc[i*optree->jIndex->childSize()+j];
        }
    }

    if(total==0.) total=1.;
    if(count==0)  count=1;

    for(unsigned int i=0; i<optree->iIndex->childSize(); i++) buffer(2+i, 0     )<<i;
    for(unsigned int j=0; j<optree->jIndex->childSize(); j++) buffer(1,   5+10*j)<<j;
    
    for(unsigned int i=0; i<optree->iIndex->childSize(); i++){
        for(unsigned int j=0; j<optree->jIndex->childSize(); j++){

            double val = 100.*appc[i*optree->jIndex->childSize()+j]/total;
            if(val>100./count) buffer(2+i, 5+10*j, VIN_HIGHLIGHT)<<int(val);
            else               buffer(2+i, 5+10*j               )<<int(val);
        }
    }
    
}

bool OperatorTreeExplorer::ApplicationCostWindow::command(std::string command){
    return command_move_buffer(command);
}

void OperatorTreeExplorer::ApplicationCostWindow::onTraverse(){
    set_offset(0,0);
    stale();
}

void OperatorTreeExplorer::SVDWindow::render(){
    buffer.clear();

    // OperatorSVD svd(dynamic_cast<OperatorTreeExplorer*>(parent())->getOptree());
    //
    // std::vector<double> sv = svd.getSingularValues();
    //
    // buffer(0,0)<<"SVD";
    // for(unsigned int i=0; i<sv.size(); i++){
    //     buffer(1+i/8, (i%8)*15)<<sv[i];
    // } 
}

bool OperatorTreeExplorer::SVDWindow::command(std::string command){
    return command_move_buffer(command);
}

void OperatorTreeExplorer::SVDWindow::onTraverse(){
    stale();
}


void OperatorTreeExplorer::PlaygroundWindow::render(){}

bool OperatorTreeExplorer::PlaygroundWindow::command(std::string command){
    command_move_buffer(command);

    const OperatorTree* op = dynamic_cast<OperatorTreeExplorer*>(parent())->getOptree();
    if(command==":hmat"){
        buffer.clear();
        buffer.stream(&std::cout);
        
        std::cout<<"Starting H Matrix creation, nr: "<<hmat_counter<<std::endl;
        HMatrix* hmat = HMatrix::truncate(op);
        std::cout<<"Done!"<<std::endl;
        std::cout<<"H Matrix"<<std::endl;
         
        for(int k=0; k<10; k++){
            hmat->rank(k);
            std::cout<<"k="<<k<<std::endl;
            
            double rel, abs;
            hmat->accuracy(rel, abs);
            std::cout<<"Cost: "<<100*hmat->asOperatorAbstract()->applicationCost()<<" Theo: "<<hmat->applyCount()<<" Accuracy: "<<rel<<"/"<<abs<<std::endl;
        }

        std::cout<<"Accuracy 1.e-6"<<std::endl;
        hmat->rankForAccuracy(1.e-6, 0.);
        std::cout<<"Cost: "<<100*hmat->asOperatorAbstract()->applicationCost()<<" Theo: "<<hmat->applyCount()<<std::endl;
       
        std::ofstream output;
        output.open((ReadInput::main.output()+"hmat_"+std::to_string(hmat_counter++)+".tex").c_str());
        hmat->writeStructureTikZ(&output);
        output.close();

        buffer.free_stream(&std::cout);
        buffer.flush();

        return true;
    }else if(command==":hmatad"){
        buffer.clear();
        buffer.stream(&std::cout);
        
        std::cout<<"Starting H Matrix creation, nr: "<<hmat_counter<<std::endl;
        HMatrix* hmat = HMatrix::truncate(op, new HMatrix::AdaptiveTruncationStrategy());
        std::cout<<"Done!"<<std::endl;
        std::cout<<"H Matrix"<<std::endl;
         
        std::cout<<"Accuracy 1.e-6"<<std::endl;
        hmat->rankForAccuracy(1.e-6, 0.);
        std::cout<<"Cost: "<<100*hmat->asOperatorAbstract()->applicationCost()<<" Theo: "<<hmat->applyCount()<<std::endl;

        std::ofstream output;
        output.open((ReadInput::main.output()+"hmat_"+std::to_string(hmat_counter++)+".tex").c_str());
        hmat->writeStructureTikZ(&output);
        output.close();

        buffer.free_stream(&std::cout);
        buffer.flush();

        return true;
    }else if(command==":hosvd"){
        buffer.clear();
        buffer.stream(&std::cout);

        std::cout<<"Starting HOSVD";
        OperatorTucker* opTucker = OperatorTucker::truncate(op);

    }

    return false;
}

void OperatorTreeExplorer::PlaygroundWindow::onTraverse(){
    buffer.clear();
}

#endif //_OPERATOR_TREE_EXPLORER_
