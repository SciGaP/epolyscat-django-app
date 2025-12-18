#include "indexGridHybrid.h"

#include "discretizationGrid.h"
#include "basisGrid.h"
#include "basisSub.h"
#include "basisIntegrable.h"
#include "basisOrbital.h"

IndexGridHybrid::IndexGridHybrid(const Index* Idx, const std::vector<std::string> &Ax,const std::vector<std::vector<double>> &Grid)
    :_gdx(0)
{
    std::vector<std::string> ax(Ax);
    std::vector<std::vector<double>> grid(Grid);
    const Index* idx=Idx;
    //    // keep only 12 digits for safe float comparisons
    //    for(auto &g: grid){
    //        for(auto &p: g)
    //            p=tools::string_to_double(tools::str(p,12));
    //    }

    if(idx->depthOfDuplicate()!=Index::npos){
        for(size_t k=1;k<idx->childSize();k++)
            if(idx->child(0)->hierarchy()!=idx->child(k)->hierarchy())
                DEVABORT("must have equivalent hierarchies on all FEs");
        idx=idx->childBack();
    }
    while(idx->axisName().find("&")!=std::string::npos)
        idx=idx->childBack();

    nodeCopy(idx,false);
    size_t kAx=std::find(ax.begin(),ax.end(),idx->axisName())-ax.begin();
    if(kAx<Ax.size()){
        // axis is to be converted
        if(idx->depthOfDuplicate()==Index::npos){
            setBasis(BasisAbstract::factory(BasisGrid::factory(grid[kAx])->strDefinition()));
            // present axis is grid - set grid basis
            ax.erase(ax.begin()+kAx);
            grid.erase(grid.begin()+kAx);
        }
    }
    if(not idx->isLeaf()){
        if(nodeEquivalent(idx)){
            for(size_t k=0;k<childSize();k++){
                childAdd(new IndexGridHybrid(idx->child(k),ax,grid));
            }
        }
        else {
            for(size_t k=0;k<basis()->size();k++){
                childAdd(new IndexGridHybrid(idx->childBack(),ax,grid));
            }
        }
    }
    sizeCompute();
}

void IndexGridHybrid::axesAndGrids(std::vector<std::string> &Ax, std::vector<std::vector<double> > &Grid) const{
    if(basis()->grid()){
        Ax.push_back(axisName());
        Grid.push_back(basis()->grid()->mesh());
    }
    if(not isLeaf())childHybrid(0)->axesAndGrids(Ax,Grid);
}

void IndexGridHybrid::imap(const Index* Gdx, std::vector<int> &I, std::vector<double> Ipt, std::vector<double> Gpt) const{
    if(Gdx->axisName().find("&")!=std::string::npos)DEVABORT("must not be called with hybrid Gdx");

    if(nodeEquivalent(Gdx) and not Gdx->isLeaf()){
        // equivalent nodes - consecutive indices
        for(size_t k=0;k<childSize();k++){
            childHybrid(k)->imap(Gdx->child(k),I,Ipt,Gpt);
        }
    }
    else if(Gdx->depthOfDuplicate()!=Index::npos){
        for(size_t k=0;k<Gdx->childSize();k++)imap(Gdx->child(k),I,Ipt,Gpt);
    }
    else if(basis()->grid()){
        std::vector<double> imesh=(basis()->grid()->mesh());
        if(not (Gdx->basis()->grid()))DEVABORT("additional grid "+Gdx->strNode()
                                               +" does not match present "+strNode());
        std::vector<double> gmesh=Gdx->basis()->grid()->mesh();

        // loop to matching points
        size_t l=0;
        Ipt.push_back(0.);
        Gpt.push_back(0.);

        for(size_t k=0;k<gmesh.size();k++){
            while(l<imesh.size() and imesh[l]<gmesh[k])l++;
            if(l==imesh.size())DEVABORT(Sstr+"additional grid point"
                                        +gmesh[k]+" not contained in grid"+imesh);
            Ipt.back()=imesh[l];
            Gpt.back()=gmesh[k];
            if(isLeaf() or Gdx->isLeaf()){
                if(Ipt!=Gpt)DEVABORT(Sstr+"points do not match: "+Ipt+" != "+Gpt);
                I.push_back(posIndex()+l);
            }
            else {
                childHybrid(l)->imap(Gdx->child(k),I,Ipt,Gpt);
            }
        }
    }
}
