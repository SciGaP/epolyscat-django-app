#include "mapGridHybrid.h"
#include "indexGrid.h"
#include "indexGridHybrid.h"
#include "basisGrid.h"
#include "basisIntegrable.h"
#include "basisOrbital.h"
#include "discretizationGrid.h"

MapGridHybrid::MapGridHybrid(const IndexGridHybrid* I, const Index* J)
    :OperatorAbstract("MapFromHybrid",I,J){

    if(J->axisName().find('&')!=std::string::npos){
        for(size_t k=0;k<J->childSize();k++){
            childAdd(new MapGridHybrid(I,J->child(k)));
        }
    }
    else if(I->basis()==J->basis()){
        // no transformation on this axis
        for(size_t k=0;k<I->childSize();k++)
            childAdd(new MapGridHybrid(I->childHybrid(k),J->child(k)));
    }
    else if(J->depthOfDuplicate()!=Index::npos){
        // FE axis
        for(size_t k=0;k<J->childSize();k++){
            childAdd(new MapGridHybrid(I,J->child(k)));
        }
    }
    else if (not J->isHybrid()){
        std::vector<std::vector<double>>grid;
        std::vector<std::string> ax;
        I->axesAndGrids(ax,grid);

        if(J->basis()->orbital()){
            J->basis()->orbital()->setMap(I,J);
            _map=J->basis()->orbital()->map();
        }
        else {
            _map.reset(new OperatorMap(new IndexGrid(J,ax,grid,{}),J));
        }
        _auxC.reset(new Coefficients(_map->idx()));
        I->imap(_auxC->idx(),_imap);
    }
    else if (I->basis()->grid() and J->basis()->integrable()){
        _toGrid=OperatorMap::basisMap(I->basis(),J->basis());
        std::vector<double>mesh(I->basis()->grid()->mesh());
        const BasisIntegrable* bi=J->basis()->integrable();
        for(size_t i=0;i<I->childSize();i++){
            if((mesh[i]<bi->upBound() or (mesh[i]==bi->upBound() and J->upperNeighbor()==0))
                    and bi->lowBound()<=mesh[i]){
                _imap.push_back(childSize());
                if(not _auxC)_auxC.reset(new Coefficients(I->childHybrid(i)));
            }
        }
        if(_auxC){
            // create map to auxiliary coefficient
            for(size_t j=0;j<J->childSize();j++)
                childAdd(new MapGridHybrid(dynamic_cast<const IndexGridHybrid*>(_auxC->idx()),J->child(j)));
        }
        // remove missing grid points from _toGrid
        if(int(_imap.size())!=_toGrid.rows()){
            for(size_t i=0;i<_imap.size();i++){
                _toGrid.row(i)=_toGrid.row(_imap[i]);
            }
            _toGrid.resize(_imap.size(),_toGrid.cols());
        }
    }
    else
        DEVABORT("map not covered"+I->strNode()+" <-- "+J->strNode());
}

std::string MapGridHybrid::strNode(int Level) const{
    if(dynamic_cast<const OperatorMap*>(_map.get()))
        return dynamic_cast<const OperatorMap*>(_map.get())->str(Level);
    else
        return Level==Tree_noNodeInfo?"":idx()->strNode()+" <-- "+jdx()->strNode();
}

void MapGridHybrid::apply(std::complex<double> A, const Coefficients& X, std::complex<double> B, Coefficients& Y) const{
    const Index* I=Y.idx();
    const Index* J=X.idx();

    Y.scale(B);
    if(J->axisName().find('&')!=std::string::npos){
        // hybrid axis - all hybrid components map into same Y
        for(size_t k=0;k<childSize();k++)
            child(k)->apply(A,*X.child(k),1.,Y);
    }
    else if(I->basis()==J->basis()){
        // no transformation: each section maps separately
        for(size_t k=0;k<childSize();k++)
            child(k)->apply(A,*X.child(k),1.,*Y.child(k));
    }
    else if(J->depthOfDuplicate()!=Index::npos){
        // FE axis: ignored here, all FE's map into overall grid
        for(size_t k=0;k<childSize();k++){
            child(k)->apply(A,*X.child(k),1.,Y);
        }
    }
    else if (not J->isHybrid()){
        // standard hierarchy: use standard map to auxiliary grid and distribute into hybrid grid
        _map->apply(A,X,0.,*_auxC);
        for(size_t k=0;k<_imap.size();k++)Y.data()[_imap[k]]+=_auxC->data()[k];
        if(X.norm()!=0. and _auxC->norm()==0.)ABORT("zero map");
    }
    else if (I->basis()->grid() and J->basis()->integrable()){
        // convert axis, but admit hybrid below present axis (cannot use
        for(size_t j=0;J->childSize();j++){
            child(j)->apply(A,*X.child(j),0,*_auxC);
            for(int i=0;i<_toGrid.rows();i++){
                Y.child(_imap[i])->axpy(_toGrid(i,j),*_auxC);
            }
        }
    }
}
