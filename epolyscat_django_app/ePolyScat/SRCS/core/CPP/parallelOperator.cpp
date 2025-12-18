// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "parallelOperator.h"

#include "operatorFloor.h"
#include "parallel.h"
#include "parallelProcess.h"
#include "str.h"
#include "parameters.h"
#include "printOutput.h"

#include "log.h"

TIMER(sync,ParOp)
TIMER(fuse1,ParOp)
TIMER(fuse2,ParOp)
TIMER(fuse3,ParOp)
TIMER(fuse4,ParOp)

using namespace std;

void ParallelOperator::unsetHost(const OperatorTree *Op){
    // apply only in parallel mode
    if(MPIwrapper::isMPI()){
        auto pHost=_host.find(Op->hash());
        if(pHost!=_host.end())_host.erase(pHost);
    }
}

void ParallelOperator::sync(OperatorTree* Op){
    if(MPIwrapper::Size()==1)return;
    ParallelOperator par(Op);
    par.syncNormCost();
}
void ParallelOperator::bcast(OperatorTree* Op){
    if(MPIwrapper::Size()==1)return;
    ParallelOperator par(Op);
    par.bcast();
}
void ParallelOperator::setDistribution(const OperatorAbstract* Op, const Host &Hosts){
    if(MPIwrapper::Size()==1)return;
    const OperatorTree * treeOp=dynamic_cast<const OperatorTree*>(Op);
    if(treeOp==0)
        _host[Op->hash()]=all; // for now, only trees can be distributed
    else{
        ParallelOperator par(treeOp);
        if(dynamic_cast<const PresentHost*>(&Hosts))
            par.setDistribution();
        else
            par.reDistribute(Hosts);
    }
}

std::map<std::string,int> ParallelOperator::_host;

/// return host as string
string ParallelOperator::strHost(const OperatorAbstract* Op){
    auto pos=_host.find(Op->hash());
    if(pos==_host.end())return "?";
    if(pos->second==ParallelOperator::none)return "n";
    if(pos->second==ParallelOperator::all)return "a";
    if(pos->second==ParallelOperator::undefined)return "u";
    if(pos->second==ParallelOperator::distributed)return "d";
    return tools::str(pos->second);
}

int ParallelOperator::getHost(const OperatorAbstract* Op,bool Warn)
{
    if(MPIwrapper::Size()==1)return all;

    auto pos=_host.find(Op->hash());
    if(pos==_host.end()){

        PresentHost pres;
        if(dynamic_cast<const OperatorTree*>(Op))pos->second=pres.redetermine(dynamic_cast<const OperatorTree*>(Op));
    }
    return pos->second;
}
//TIMER(allreducePres,)
int ParallelOperator::PresentHost::operator()(const OperatorTree* Op) const {

    auto ihost=_host.find(Op->hash());
    if(ihost==_host.end()){
        vector<int>posCnt(2,0);
        // if not dummy, count 1 and communicate rank
        if(dynamic_cast<const OperatorDUM*>(Op->floor())==0){
            posCnt={MPIwrapper::Rank(),1};
        }
        if(Op->name == "Commutator" and /*ReadInput::main.flag("DEBUGCommutator","suppress parallel of commutator")*/ MPIwrapper::Size() > 1)
            _host[Op->hash()]=all;

        else if(MPIwrapper::Size() > 1){
            //            STARTDEBUG(allreducePres);
            MPIwrapper::AllreduceSUM(posCnt.data(),posCnt.size()); // can be very time-consuming when using multiple nodes
            //            STOPDEBUG(allreducePres);
        }
        else
            return Op->iIndex->index()[0];

        if(posCnt[1]==0)_host[Op->hash()]=none;
        else if(posCnt[1]==MPIwrapper::Size())_host[Op->hash()]=all;
        else if(posCnt[1]==1)_host[Op->hash()]=posCnt[0];
        else ABORT(Sstr+"floor defined on several, but not on all processes:"+posCnt+" oper: "+Op->name+Op->index());
        ihost=_host.find(Op->hash());
    }
    return ihost->second;
}

int ParallelOperator::PresentHost::redetermine(const OperatorTree* Op) const {

    MPIwrapper::Barrier();
    auto ihost=_host.find(Op->hash());
    int shost=ihost==_host.end()?undefined:_host[Op->hash()];

    // if not dummy, count 1 and communicate rank
    vector<int>posCnt(2,0);
    if(!dynamic_cast<const OperatorDUM*>(Op->floor()))posCnt={MPIwrapper::Rank(),1};

    if(Op->name=="Commutator" and MPIwrapper::Size() > 1)
        _host[Op->hash()]=all;

    else if(MPIwrapper::Size() > 1){
        MPIwrapper::Barrier();
        MPIwrapper::AllreduceSUM(posCnt.data(),posCnt.size()); // can be very time-consuming when using multiple nodes
    }
    else
        return Op->iIndex->index()[0];

    if(posCnt[1]==0)_host[Op->hash()]=none;
    else if(posCnt[1]==MPIwrapper::Size())_host[Op->hash()]=all;
    else if(posCnt[1]==1)_host[Op->hash()]=posCnt[0];
    else ABORT(Sstr+"floor defined on several, but not on all processes:"+posCnt+" oper: "+Op->name+Op->index());
    if(shost!=undefined and shost!=_host[Op->hash()])
        PrintOutput::DEVwarning(Sstr+"wrong host was stored"+shost+"actual"+_host[Op->hash()]+Op->strNode());
    MPIwrapper::Barrier();
    return _host[Op->hash()];
}

ParallelOperator::ProcessHost::ProcessHost(const OperatorTree* Op, int Host){
    for(const OperatorTree * b=Op;b!=0;b=b->nodeNext(Op))_proc[b]=Host;
}

ParallelOperator::ProcessHost::ProcessHost(const vector<ParallelProcess*> & Proc){
    for(size_t k=0;k<Proc.size();k++){
        for(size_t l=0;l<Proc[k]->localB.size();l++)
            _proc[Proc[k]->localB[l]->oLeaf]=k;

        for(size_t l=0;l<Proc[k]->sendB.size();l++)
            for(size_t m=0;m<Proc[k]->sendB[l].size();m++)
                _proc[Proc[k]->sendB[l][m]->oLeaf]=k;

        for(size_t l=0;l<Proc[k]->recvB.size();l++)
            for(size_t m=0;m<Proc[k]->recvB[l].size();m++)
                _proc[Proc[k]->recvB[l][m]->oLeaf]=k;
    }
}

int ParallelOperator::DiagonalBlockHost::operator()(const OperatorTree* Op) const {
    const OperatorTree* block=Op;
    while(block!=_root and block!=0)block=block->parent();
    if(block==0)ABORT("block does not belong to root");

    // descend from root towards present until not block-diagonal
    while(block->isBlockDiagonal() and block->iIndex->overlap()->isBlockDiagonal()){
        block=block->child(Op->index()[block->depth()]);
    }

    // Store in _host - required by EigenSolver
    if(_host.count(block->hash())
            and block->parent()
            and _host[block->hash()]!=int(block->levelRank(_root)%MPIwrapper::Size()))
        PrintOutput::DEVwarning(Sstr+block->name+"illegal host reset\n"+
                                _host[block->hash()]+"->"+block->levelRank(_root)
                +block->iIndex->strNode()+block->jIndex->strNode(),10);

    _host[block->hash()]=block->levelRank(_root)%MPIwrapper::Size();
    return _host[block->hash()];
}


int ParallelOperator::ProcessHost::operator ()(const OperatorTree* Op) const {
    if(_proc.count(Op)==0)return none;
    return const_cast<ProcessHost*>(this)->_proc[Op];
}
ParallelOperator::SavedHost::SavedHost(const ParallelOperator *Par){
    for(int k=0;k<Par->_view.triplets();k++)
        for (int l=0;l<Par->_view.blocks(k);l++)
            savedHost[Par->_view.block(k,l)]=_host[Par->_view.block(k,l)->hash()];
}
int ParallelOperator::SavedHost::operator()(const OperatorTree *Op) const{
    map<const OperatorAbstract*,int>::const_iterator sh;
    sh=savedHost.find(Op);
    if(sh==savedHost.end())ABORT("no host saved for operator");
    return sh->second;
}

ParallelOperator::ParallelOperator(const OperatorAbstract *Op)
{
    _op=const_cast<OperatorAbstract*>(Op);

    OperatorTree* treeOp=dynamic_cast<OperatorTree*>(_op);
    if(treeOp!=0)_view=BlockView(treeOp,BlockView::depthFloor,BlockView::depthFloor);
    else         _view=BlockView(Op);
}


void ParallelOperator::purge(){
    Parameters::updateToOne(); // make sure all parameters are non-zero
    syncNormCost();
    OperatorTree* treeOp=dynamic_cast<OperatorTree*>(_op);
    treeOp->purge(emptyLeafOrDummy);
    treeOp->purge();
    ParallelOperator fused(treeOp);
    swap(_view,fused._view);
    syncNormCost();
    Parameters::restoreToTime();
}


void ParallelOperator::fuse() {
    syncNormCost();
    OperatorTree* treeOp=dynamic_cast<OperatorTree*>(_op);
    if(treeOp==0)ABORT("cannot fuse non-tree operator "+_op->name);


    treeOp->fuseBottomUp();
    ParallelOperator fused(treeOp);
    swap(_view,fused._view);
    syncNormCost();
}

void ParallelOperator::syncNormCost(){
    if(MPIwrapper::Size()==1)return; // single process

    for(int tripl=0;tripl<_view.triplets();tripl++){
        for(int k=0;k<_view.blocks(tripl);k++){
            OperatorTree* op=const_cast<OperatorTree*>(dynamic_cast<const OperatorTree*>(_view.block(tripl,k)));
            if(op!=0 and dynamic_cast<OperatorDUM*>(op->floor()) == 0){
                op->floor()->norm();
                op->floor()->applicationCost(false);
            }
        }
    }

    STARTDEBUG(sync);
    PresentHost presHost;
    for(int tripl=0;tripl<_view.triplets();tripl++){
        for(int k=0;k<_view.blocks(tripl);k++){
            OperatorTree* op=const_cast<OperatorTree*>(dynamic_cast<const OperatorTree*>(_view.block(tripl,k)));
            if(op==0){
                Str s("block is not OperatorTree, orginal pointer=");
                s=s+_view.block(tripl,k);
                if(_view.block(tripl,k)!=0)s=s+" operator="+_view.block(tripl,k)->name;
                ABORT(s);
            }
            int host=presHost(op);
            if(host==none)
                if(op->floor())DEVABORT("floor does not have any owner: "+op->str());
            if(host==all)host=MPIwrapper::master();
            vector<double> nrmCst(2,0.);
            if(host==MPIwrapper::Rank())
                nrmCst={op->floor()->norm(),op->floor()->applicationCost(false)};
            MPIwrapper::Bcast(nrmCst.data(),nrmCst.size(),host); // can be time-consuming when setting up, however, seems unavoidable
            op->floor()->forceNorm(nrmCst[0]);
            op->floor()->forceCost(nrmCst[1]);
        }
    }
    STOPDEBUG(sync);

}

void ParallelOperator::bcast(){
    if(MPIwrapper::Size()==1)return; // single process


    PresentHost presHost;
    for(int tripl=0;tripl<_view.triplets();tripl++){
        for(int k=0;k<_view.blocks(tripl);k++){
            OperatorTree* op=const_cast<OperatorTree*>(dynamic_cast<const OperatorTree*>(_view.block(tripl,k)));
            if(op==0)ABORT("block is not OperatorTree");
            int host=presHost(op);
            if(host==none)continue;//ABORT("no host data found");
            if(host!=all)
            {
                vector<int> info(5,0);
                vector<complex<double> > buf;
                if(host==MPIwrapper::Rank())op->floor()->pack(info,buf);

                MPIwrapper::Bcast(info.data(),info.size(),host);

                if(info[4]>0){
                    buf.resize(info[4]);
                    MPIwrapper::Bcast(buf.data(),buf.size(),host);
                }

                if(host!=MPIwrapper::Rank()){
                    //even dummy floors carry the factors
                    complex<double>*tFac=op->oFloor->factor();
                    delete op->oFloor;
                    op->oFloor=OperatorFloor::unpackFactory(info,buf);
                    op->oFloor->setFactor(tFac);
                }
            }
        }
    }
    _view=BlockView(dynamic_cast<const OperatorTree*>(_op),BlockView::depthFloor,BlockView::depthFloor);

    // register all pieces as all
    _host[_op->hash()]=all;
    OperatorTree * opTree=dynamic_cast<OperatorTree*>(_op);
    while (opTree!=0){
        _host[opTree->hash()]=all;
        opTree=opTree->nodeNext(dynamic_cast<OperatorTree*>(_op));
    }

}


TIMER(reDistribute,)
void ParallelOperator::reDistribute(const Host& NewHosts){
    const SingleHost* singHo=dynamic_cast<const SingleHost*>(&NewHosts);
    const OperatorTree* opTree=dynamic_cast<const OperatorTree*>(_op);
    if(not opTree)DEVABORT("only OperatorTree can be distributed");
    if(not singHo){
        if(not opTree->floor())_host[_op->hash()]=distributed;
    }
    else {
        // assign all, except leaf:
        const OperatorTree* branch(opTree);
        while(branch){
            if(not branch->isLeaf())_host[branch->hash()]=NewHosts(branch);
            branch=branch->nodeNext(opTree);
        }
    }
    if(MPIwrapper::Size()>1){
        PresentHost presHost;
        int totalBlocks = 0;
        for(int tripl=0;tripl<_view.triplets();tripl++){
            for(int k=0;k<_view.blocks(tripl);k++){
                totalBlocks ++;
                OperatorTree* op=const_cast<OperatorTree*>(dynamic_cast<const OperatorTree*>(_view.block(tripl,k)));
                if(op==0)ABORT("block is not OperatorTree");
                // Get current host - not stored in _host!
                presHost.redetermine(op);
                int host=presHost(op);
                if(host==all){
                    // Store in _host
                    _host[_op->hash()]=all;
                    _host[op->hash()]=all;

                }else if(host!=none and _host[_op->hash()]==all){
                    // If the root of the operator tree has "all", all children need to be "all" or "none"
                    ABORT(op->root()->str()+"\ninconsistent operator distribution for "+op->name);

                }else{
                    // Get new host
                    int newHost=NewHosts(op); //NOTE: must be after the re-definition of _host[_op->hash()]
                    if(host!=Parallel::moveFloor(op->floor(),newHost))DEVABORT(Sstr+"location does not match expected "+host);
                    // Store new host
                    _host[op->hash()]=newHost;
                }
            }
        }
    }
}

std::string ParallelOperator::str() const {
    Str s("Parallel: ","",4); // set default width for number conversion
    PresentHost presHost;
    const OperatorTree * op;
    for(int k=0;k<_view.triplets();k++)
    {
        if(k==0)s=s+_view.block(0,0)->name+"\n";

        op=dynamic_cast<const OperatorTree*>(_view.block(k,0));
        int host=presHost(op);
        if(host==all)s=s+" all";
        else if(host==none)s=s+"none";
        else s=s+"<"+host+">";
        s=s+_view.block(k,0)->str()+" #"+_view.blocks(k);
        s=s+"\n";
    }
    return string(s);
}

bool ParallelOperator::emptyLeafOrDummy(const OperatorTree* Op){
    if(not Op->isLeaf())return false;
    if(Op->floor()==0)return true; // empty leaf
    int dummy=(dynamic_cast<const OperatorDUM*>(Op->floor())!=0);
    MPIwrapper::AllreduceSUM(&dummy,1); // not very time-consuming
    return dummy==MPIwrapper::Size(); // dummy on all
}
