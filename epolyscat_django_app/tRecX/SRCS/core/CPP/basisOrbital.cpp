// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "basisOrbital.h"

#include "mapGridHybrid.h"
#include "indexGridHybrid.h"
#include "discretizationGrid.h"
#include "operatorTree.h"
#include "operatorFloor.h"
#include "index.h"
#include "operatorDefinition.h"
#include "readInput.h"
#include "printOutput.h"
#include "operatorHartree.h"
#include "parallelOperator.h"
#include "mpiWrapper.h"
#include "plot.h"

#include "basisMatMatrix.h"
#include "eigenTools.h"
#include "basisIntegrable.h"
#include "operatorMap.h"
#include "index.h"

using namespace std;

// this should actually go into Index
std::map<std::string,const Index*> BasisOrbital::referenceIndex;
void BasisOrbital::addIndex(std::string Name, const Index* Idx){
    // add if named
    if(tools::cropString(Name)!="")
        referenceIndex[Name]=Idx;
}

void BasisOrbital::generate() const {

    if(_orb.size()==0 or _orb[0].idx()==0)
        const_cast<BasisOrbital*>(this)->generateOrbitals(0);
    if(_orb.size()==0 or _orb[0].idx()==0)
        DEVABORT("default generation of orbitals failed, run generateOrbitals(Index) before accessing orbitals");
}

Eigen::MatrixXd BasisOrbital::overlap() const{
    OperatorTree op("<<Overlap>>",string("<<Overlap>>"),_orb[0].idx(),_orb[0].idx());
    Coefficients d(_orb[0].idx());
    Eigen::MatrixXd res(size(),size());
    for(size_t j=0;j<size();j++){
        op.apply(1.,_orb[j],0.,d);
        for(size_t i=0;i<size();i++){
            res(i,j)=_orb[i].innerProduct(&d).real();
        }
    }
    return res;
}

std::vector<double> BasisOrbital::expectationValues(string OpDef) const{
    if(_orb.size()==0)return vector<double>();
    if(_orb.size()>0 and _orb[0].idx()==0)ABORT("run generateOrbitals() before computing expectationValues");
    OperatorTree op(OpDef,OpDef,_orb[0].idx(),_orb[0].idx());
    vector<double> expec;
    for(Coefficients c: _orb)expec.push_back(op.matrixElement(c,c).real());
    return expec;
}


BasisOrbital::BasisOrbital(string Name, ReadInput *Inp):BasisAbstract(Name),_isOrthonormal(false){
    if(Inp==0)return;
    // read plot definition
    int line(0);
    std::string plt;
    do{
        Inp->read(Name,"plot",plt,"","specify AxName:{lowbound}:upboud{:npts}, lowbound defaults to -upbound, npts defaults to 100",++line);
        _plotDef.push_back(plt);
    } while (plt!="");
    _plotDef.pop_back();
}

void BasisOrbital::plot() const {
    if(size()==0)return;
    Plot plt(_orb[0].idx(),ReadInput::main);
    if(plt.isEmpty())return;
    std::string pltFile=ReadInput::main.output()+tools::stringInBetween(name(),"[","]");
    tools::substringReplaceAll(pltFile,":","_");

    int wid= size()>100 ? 3 : 2;
    // (overwrite) orbitals
    for(int k=0;k<int(size());k++)
        plt.plot(_orb[k],pltFile+tools::str(k,wid,'0'),{},"",true);
    PrintOutput::message("Orbital plots on "+pltFile+tools::str(0,wid,'0')+" through ..."+tools::str(size()-1,wid,'0'));
    PrintOutput::paragraph();

}

std::vector<size_t> BasisOrbital::inScaledRegion() const{
    std::vector<size_t> res;
    for(size_t k=0;k<orbitals().size();k++){
        for(const Coefficients* c=orbital(k)->firstLeaf();c!=0;c=c->nextLeaf()){
            if(c->idx()->basis()->isAbsorptive() and c->norm()>1.e-12){
                res.push_back(k);
                break;
            }
        }
    }
    return res;
}


const Coefficients *BasisOrbital::orbital(int N) const{
    generate();
    if(int(size())<N+1)
        ABORT(Str("there are only")+size()+"orbitals, cannot return N="+N);
    return &_orb[N];
}

std::vector<const Coefficients*> BasisOrbital::orbitals() const {
    generate();
    vector<const Coefficients*> res;
    for(const Coefficients &c: _orb)res.push_back(&c);
    for(auto c: _orb)if(std::abs(c.cMaxNorm())>1.e5)ABORT("too large");
    return res;
}

void BasisOrbital::orthonormalize(bool Warn){
    if(_orb.size()==0)return;

    const OperatorAbstract * overlap=_orb[0].idx()->overlap();
    Coefficients *sCk=overlap->tempLHS();
    int k=0;
    double errOrtho=0.;
    std::vector<double> norms;
    while(k<int(_orb.size())){
        overlap->apply(1.,_orb[k],0.,*sCk);
        norms.push_back(std::abs(_orb[k].innerProduct(sCk)));
        complex<double> a;
        for(int l=0;l<k;l++){
            a=_orb[l].innerProduct(sCk);
            _orb[k].axpy(-a,&_orb[l]);
            errOrtho=max(errOrtho,abs(a));
        }
        a=_orb[k].innerProduct(sCk);
        if(abs(a)<1.e-20)ABORT("vectors linearly dependent");
        std::complex<double>phas=_orb[k].cMaxNorm();
        _orb[k]*=std::abs(phas)/(phas*sqrt(a));
        k++;
    }
    double errNorm=0.;
    for(double &n: norms)n-=1.;
    for(double n: norms)errNorm=std::max(errNorm,std::abs(n));
    if(errNorm>1.e-12){
        PrintOutput::warning(Str("orbitals "+name()+" not normalized by maximally ","")+errNorm+" - may need more accurate basis",10,0,
                             "               norms-1 are "+tools::str(norms));
    }
    if(Warn and errOrtho>1.e-12)
        PrintOutput::warning(Str("orthonormalizing: expanded "+name()+" not orthogonal by maximally"," ")+errOrtho+"- may need more accurate basis");
    _isOrthonormal=true;
}

static Eigen::MatrixXcd allCoefs( vector<const Coefficients*> ICoefs){
    if(ICoefs.size()==0)DEVABORT("emtpy vector of Coefficients");
    Eigen::MatrixXcd allMat(ICoefs[0]->size(),ICoefs.size());
    for(size_t k=0;k<ICoefs.size();k++)allMat.col(k)=Eigen::Map<Eigen::MatrixXcd>(ICoefs[k]->orderedData(),ICoefs[k]->size(),1);
    return allMat;
}

// true if time-depedent factors agree (or none on both)
bool equalTimeDep(const OperatorTree* A, const OperatorTree* B){
    if(A->floor()==0)return B->floor()==0;
    if(B->floor()==0)return false;
    return A->floor()->factor()==B->floor()->factor();
}


static void contract(OperatorTree* &Node, OperatorTree * OpFull, vector<const Coefficients*> ICoefs, vector<const Coefficients*> JCoefs
                     ,std::complex<double> Multiplier, std::vector<std::complex<double> *> TFac){

    ParallelOperator::bcast(OpFull);
    if(OpFull->floor()){
        if(not MPIwrapper::isMaster()){
            Node->floor()=new OperatorDUM();
        } else {
            Eigen::MatrixXcd ijmat;
            if(OpFull->floor()->rows()==0){
                // some parts of the full operator may not have been set up - defer to post-processing
                if(dynamic_cast<const OperatorHartree*>(OpFull->floor())){
                    Node->floor()=new OperatorDUM(1.);
                    Node->floor()->setFactor(OpFull->floor()->factor());
                    return;
                }
                else
                    DEVABORT("empty floor: "+OpFull->strNode(-1)+" - cannot contract");
            }

            Eigen::MatrixXcd fmat=OpFull->floor()->matrix();
            //NOTE: we use standard scalar product, in case of complex scaled region, this is not correct
            //      orbitals should not be non-zero in that region
            Eigen::MatrixXcd imat =ICoefs.size()>0 ? allCoefs(ICoefs).adjoint()*fmat : fmat;
            ijmat=JCoefs.size()>0 ? imat*allCoefs(JCoefs) : imat;


            if(Node->floor()){
                if(Node->floor()->factor()!=OpFull->floor()->factor())
                    DEVABORT("cannot contract, components have different time-depenendent factors:\n"
                             +OpFull->def());
                ijmat+=Node->floor()->matrix();
                delete Node->floor();
            }
            string hash(OpFull->def()+Node->iIndex->hash()+Node->jIndex->hash());
            ijmat*=Multiplier;
            Node->floor()=OperatorFloor::factory(vector<const Eigen::MatrixXcd*>(1,&ijmat),hash);
        }
        if(OpFull->floor()->factor() and TFac.size()>1)DEVABORT("can only have single time-dependente factor in a single operator term");
        Node->floor()->setFactor(OpFull->floor()->factor());
    }

    for (size_t n=0;n<OpFull->childSize();n++){
        vector<const Coefficients*>iCoefs;
        if(OpFull->child(n)->iIndex==OpFull->iIndex)
            iCoefs=ICoefs; // level does not change - keep coefs
        else for(const Coefficients*c: ICoefs)
            iCoefs.push_back(c->child(OpFull->child(n)->iIndex->nSibling()));

        vector<const Coefficients*>jCoefs;
        if(OpFull->child(n)->jIndex==OpFull->jIndex)
            jCoefs=JCoefs; // level does not change - keep coefs
        else for(const Coefficients*c: JCoefs)
            jCoefs.push_back(c->child(OpFull->child(n)->jIndex->nSibling()));

        const Index* iChild= ICoefs.size()>0 ? Node->iIndex : OpFull->child(n)->iIndex;
        const Index* jChild= JCoefs.size()>0 ? Node->jIndex : OpFull->child(n)->jIndex;

        // seek match: node or child (indices and possible time-dependent parameter)
        OperatorTree* nodeChild=0;
        if(Node->iIndex==iChild and Node->jIndex==jChild and equalTimeDep(Node,OpFull->child(n)))
            nodeChild=Node;
        else for(size_t k=0;k<Node->childSize();k++)
            if(Node->child(k)->iIndex==iChild and Node->child(k)->jIndex==jChild and equalTimeDep(Node->child(k),OpFull->child(n)))
                nodeChild=Node->child(k);

        // no match - create new
        if(nodeChild==0){
            Node->childAdd(new OperatorTree(OpFull->name,iChild,jChild));
            nodeChild=Node->childBack();
        }
        contract(nodeChild,OpFull->child(n),iCoefs,jCoefs,Multiplier,TFac);
    }

}

void orbitalMatrix(const Eigen::MatrixXcd & Mat, OperatorTree* Node){
    if(not MPIwrapper::isMaster()){
        Node->floor()=new OperatorDUM();
    } else {
        string hash(Node->def()+Node->iIndex->hash()+Node->jIndex->hash());
        Node->floor()=OperatorFloor::factory(vector<const Eigen::MatrixXcd*>(1,&Mat),hash);
    }
}

static std::map<std::string, std::shared_ptr<OperatorTree>> _orbitalOp;

void BasisOrbital::attachOperator(OperatorTree* Node, std::string Name, const OperatorDefinition &Def, const Index* IIndex, const Index* JIndex,
                                  std::complex<double> Multiplier, std::vector<std::complex<double>*> TFac){

    const BasisMatMatrix* m(0);
    if(Def.str().find("[[")<Def.str().find("<")){
        std::string mDef="<"+tools::stringInBetween(Def.str(),"[[","]]",true)+">"; // matrices are defined as <NameOfMatrix>
        m=BasisMatMatrix::factory(mDef,IIndex->basis(),JIndex->basis());
        if(not m)DEVABORT("did not find special"+mDef+" | "+Def.str());
    }
    if(m!=0){
        orbitalMatrix(m->mat(),Node);
    }
    else {
        // get full indices on both sides
        const Index *iFull=IIndex,*jFull=JIndex;
        const BasisOrbital *ib,*jb;
        if(0!=(ib=dynamic_cast<const BasisOrbital*>(IIndex->basis())))iFull=ib->orbital(0)->idx();
        if(0!=(jb=dynamic_cast<const BasisOrbital*>(JIndex->basis())))jFull=jb->orbital(0)->idx();

        if(ib and ib->inScaledRegion().size())
            PrintOutput::warning(Sstr+"Orbitals"+ib->inScaledRegion()+"extend into scaled region - operator may be inaccurate");
        if(jb and jb!=ib and jb->inScaledRegion().size())
            PrintOutput::warning(Sstr+"Orbitals"+jb->inScaledRegion()+"extend into scaled region - operator may be inaccurate");

        // get full operator
        std::string hash=iFull->hash()+jFull->hash()+Def.str();
        //        OperatorDefinition def=Def.expandStandard(iFull->hierarchy());
        if(!_orbitalOp.count(hash))_orbitalOp[hash].reset(new OperatorTree(Name,Def,iFull,jFull));

        // recursively build
        vector<const Coefficients*> iCoefs,jCoefs;
        if(ib)iCoefs=ib->orbitals();
        if(jb)jCoefs=jb->orbitals();

        // here we need to have Multiplier and Tfac
        contract(Node,_orbitalOp[hash].get(),iCoefs,jCoefs,Multiplier,TFac);
    }
}
void BasisOrbital::clear(){_orbitalOp.clear();}

void BasisOrbital::setMap(const Index *Grid, const Index* OrbJdx) const{
    if(not _mapOrbs or _mapOrbs->idx()!=Grid){
        Grid->sizeCompute();
        _mapOrbs.reset(new Map(Grid,OrbJdx,this));
    }
}

BasisOrbital::Map::Map(const Index* Grid, const Index* OrbJdx, const BasisOrbital* Orb)
    :OperatorAbstract(Orb->name(),Grid,OrbJdx),_orb(Orb)
{
    if(_mapV.size()==0 or _mapV[0].idx()!=Grid){
        std::unique_ptr<OperatorAbstract> map;
        if(dynamic_cast<const IndexGridHybrid*>(Grid)){
            map.reset(new MapGridHybrid(dynamic_cast<const IndexGridHybrid*>(Grid),Orb->orbital(0)->idx()));
        }
        else {
            map.reset(new OperatorMap(Grid,Orb->orbital(0)->idx()));
        }

        for(size_t k=0;k<_orb->orbitals().size();k++){
            _mapV.push_back(Coefficients(Grid,0.));
            map->apply(1.,*_orb->orbital(k),0,_mapV.back());
        }
    }
}

void BasisOrbital::Map::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{
    Y.scale(B);
    for(size_t j=0;j<jdx()->size();j++){
        Y.axpy(A*Vec.data()[j],_mapV[j]);
    }
}


