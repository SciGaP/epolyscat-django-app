#include "operatorSubindex.h"

#include "index.h"
#include "basisSub.h"
#include "basisExpIm.h"


/// list of functions of Sub in Bas
static std::vector<int> subset(const BasisAbstract* Bas, const BasisAbstract* Sub){
    std::vector<int> res;
    if(*Bas==*Sub){
        for(size_t k=0;k<Bas->size();k++)res.push_back(k);
    }
    else if(*BasisSub::superBas(Sub)==*Bas){
        res=BasisSub::subset(Sub);
    }
    else if(Bas->name()==Sub->name()){
        for(size_t k=0;k<Sub->size();k++){
            for(size_t l=0;l<Bas->size();l++){
                if(Bas->physical(l)==Sub->physical(k)
                        and std::find(res.begin(),res.end(),l)==res.end()){
                    res.push_back(l);
                    break;
                }
            }
        }
    }
    return res;
}

static const Index* mdx(const OperatorTree* Op,bool MatchI){
    return MatchI? Op->idx():Op->jdx();
}

OperatorSubindex::OperatorSubindex(const OperatorTree*Op,const Index* SubIdx, bool MatchI)
{
    name=Op->name+"_Subindex";
    nodeCopy(Op,false);
    if(MatchI)iIndex=SubIdx;
    else      jIndex=SubIdx;

    if(not Op->isLeaf()){
        if(mdx(Op,MatchI)==mdx(Op->child(0),MatchI)){
            for(size_t k=0;k<Op->childSize();k++)
                childAdd(new OperatorSubindex(Op->child(k),SubIdx,MatchI));
        }
        else {
            std::vector<int> sub(subset(mdx(Op,MatchI)->basis(),SubIdx->basis()));
            if(sub.size()!=SubIdx->childSize())
                DEVABORT("SubIdx basis "+SubIdx->basis()->str()
                         +"\n not contained in Op mdx basis "+mdx(Op,MatchI)->basis()->str());
            for(size_t k=0;k<Op->childSize();k++){
                // SubIdx.child(subp) matches Op.child(k) index
                size_t matchK=std::find(sub.begin(),sub.end(),mdx(Op->child(k),MatchI)->nSibling())-sub.begin();
                if(matchK!=sub.size())
                    childAdd(new OperatorSubindex(Op->child(k),SubIdx->child(matchK),MatchI));
            }
        }
    }
}
