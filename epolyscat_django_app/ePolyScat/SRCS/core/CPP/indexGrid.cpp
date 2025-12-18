#include "indexGrid.h"

#include "qtEigenDense.h"
#include "printOutput.h"
#include "log.h"

#include "basisDvr.h"
#include "basisSub.h"
#include "basisGridQuad.h"

// convert selected levels to grid
IndexGrid::IndexGrid(const Index *I, std::vector<unsigned int> Level,
                                   std::vector<unsigned int> Point, std::vector<double> Limit)
{
    setAxisName(I->axisName());
    setBasis(I->basis());

    setKind(I->indexKind());

    unsigned int kTrans=tools::locateElement(Level,I->depth()); // locate level
    if(kTrans>=Level.size() or I->basis()->isGrid()){
        // not in conversion list or is already grid
        for(unsigned int k=0;k<I->childSize();k++)
            Tree<Index>::childAdd(new IndexGrid(I->child(k),Level,Point,Limit));
    }
    else {
        // number of points can be specified (default = order underlying the basis)
        int nPoint=I->basis()->integrable()->order();
        const BasisIntegrable* superB=dynamic_cast<const BasisIntegrable*>(BasisSub::superBas(I->basis()));
        //HACK
        if(I->axisName().find("Eta")==0){

            //Eta axes usually vary for different m's, need grid suitable for all m's
            std::string depAx=tools::stringInBetween(superB->name(),"{","}");
            depAx=depAx.substr(0,depAx.find(".")); //remove possible constraint specification
            // now we should check whether that axis is to be converted

            nPoint=0;
            for(const Index* idx=I;idx!=0;idx=idx->nodeRight())
                nPoint=std::max(nPoint,int(idx->basis()->integrable()->order()));
        }
        if(kTrans<Point.size())nPoint=Point[kTrans];



        UseMatrix points,weights;
        if(2*kTrans+1>Limit.size()){
            // no Limit specified, use quadrature grid
            if(superB->isDVR() and nPoint==int(superB->order()))
                dynamic_cast<const BasisDVR*>(superB)->dvrRule(points,weights);
            else
                superB->quadRule(nPoint,points,weights);
        }
        else {
            // equidistant points between limits
            double dx=0;
            if(nPoint>1){
                if(Limit[2*kTrans]>=Limit[2*kTrans+1])
                    ABORT(Str("need lowerLimit<upperLimit for multiple grid points, axis and interval:")
                          +I->axisName()+Limit[2*kTrans]+Limit[2*kTrans+1]);
                dx=(Limit[2*kTrans+1]-Limit[2*kTrans])/(nPoint-1);
            }

            // determine points falling into [low,upp), [low,upp] for last element;
            double upp=superB->integrable()->upBound();
            double low=superB->lowBound();
            double eps=(upp-low)*1.e-10;
            const Index* uppN=I->upperNeighbor();
            if(uppN==0)upp+=eps;   // raise by epsilon for last element

            const Index* lowN=I->lowerNeighbor();
            if(lowN==0)low-=eps;   // lower by epsilon for first element
            else       low=lowN->basis()->integrable()->upBound(); // use upper boundary of lower neighbor

            std::vector<double> pts;
            for(int k=0;k<nPoint;k++){
                if(upp<=Limit[2*kTrans]+k*dx)break;
                if(low<=Limit[2*kTrans]+k*dx)pts.push_back(Limit[2*kTrans]+k*dx);
            }

            if(pts.size()!=0){
                points=UseMatrix(pts.size(),1),weights=UseMatrix(pts.size(),1);
                for (unsigned int n=0;n<pts.size();n++){
                    points(n)=pts[n];
                    weights(n)=dx;
                }
            }
        }

        if(points.size()>0){
            std::vector<double> pts,wgs;
            for(size_t k=0;k<points.size();k++){
                pts.push_back(points(k).real());
                wgs.push_back(weights(k).real());
            }
            setBasis(BasisGridQuad::factory(pts,wgs));
        }
        if(not I->isLeaf())
            for(unsigned int k=0;k<points.size();k++)
                Tree<Index>::childAdd(new IndexGrid(I->child(0),Level,Point,Limit));


        // converted lower indices must be equivalent
        for(size_t k=1;k<I->childSize() and childSize()>0;k++)
            if(not IndexGrid(I->child(k),Level,Point,Limit).treeEquivalent(child(0))){
                Sstr+I->child(0)->str()+Sendl;
                Sstr+I->child(k)->str()+Sendl;
                DEVABORT("branches are not equivalent after conversion - cannot convert to grid on present level "+I->axisName()+"\n"
                         +IndexGrid(I->child(k),Level,Point,Limit).str()
                         +"\n"+child(0)->str());
            }
    }
    sizeCompute();

}
/// grid index from Index
///
/// Rules:
///<br> axisName==Axis[k] will be converted to grid given in Grid[k] with quadrature weights Weig[k]
///<br> if Grid[k].size()==0, use Weigh[k].size() quadrature points, if Weigh[k].size()==0, use Basis order() points
///<br> BasisIntegrable is evaluated, BasisGrid interpolated for Grid points, other bases cannot be converted
///<br> on finite element axes, only Grid[k] points falling into given element will be converted (-> piece-wise grid)
IndexGrid::IndexGrid(const Index *I, const std::vector<std::string> & Ax,
                                   std::vector<std::vector<double> > Grid,
                                   std::vector<std::vector<double> > Weig)
{
    if(Weig.size()==0)Weig.resize(Ax.size());
    if(Grid.size()!=Ax.size() or Weig.size()!=Ax.size())DEVABORT("must specify grid points and weights for each axis");


    setAxisName(I->axisName());
    setBasis(I->basis());
    setKind(I->indexKind());

    size_t nAx=(std::find(Ax.begin(),Ax.end(),axisName())-Ax.begin());
    const BasisAbstract* supBas=BasisSub::superBas(I->basis());

    // not in list index or identical grid
    if(nAx==Ax.size() or (supBas->isGrid() and Grid[nAx].size()==0)){
        if(nAx<Ax.size() and not supBas->isIndex() and (supBas->isGrid()!=(Grid[nAx].size()==0)))
            DEVABORT("must give Grid points if on FE axis: "+axisName());
        // not in conversion list or not integrable
        for(unsigned int k=0;k<I->childSize();k++)
            Tree<Index>::childAdd(new IndexGrid(I->child(k),Ax,Grid,Weig));
    }
    // FE-level
    else if(supBas->isIndex()){
        // check for subbasis
        std::vector<int> mSub=BasisSub::subset(I->basis());
        if(true or supBas!=I->basis()){
            std::map<int,const Index*> mdx;
            for(size_t k=0;k<mSub.size();k++)mdx[mSub[k]]=I->child(k);

            // find largest subset shared among all, add subindices where missing
            for(const Index* ix=I->root()->descend(I->depth());ix!=0;ix=ix->nodeRight()){
                if(*BasisSub::superBas(I->basis())==*BasisSub::superBas(ix->basis())){
                    std::vector<int> iSub=BasisSub::subset(ix->basis());
                    for(size_t k=0;k<iSub.size();k++){
                        // value at it will be less or equal to i, if not equal, insert to mSub and add to mdx
                        auto it=std::lower_bound(mSub.begin(),mSub.end(),iSub[k]);
                        if(it==mSub.end() or *it!=iSub[k]){
                            mSub.insert(it,iSub[k]);
                            mdx[iSub[k]]=ix->child(k);
                        }
                    }
                }
            }
            // create grid for maximally extended index
            setBasis(BasisAbstract::factory(BasisSub::strDefinition(supBas,mSub)));

            //CHECK: tmp may not be needed after all...
            Index tmp;
            tmp.setBasis(basis());
            tmp.setAxisName(axisName());
            for(auto k: mSub)tmp.childAdd(new IndexNew(*mdx[k]));
            tmp.sizeCompute();
            for(auto k: mSub){
                Tree<Index>::childAdd(new IndexGrid(tmp.child(k),Ax,Grid,Weig));
            }
        }
        else {
            // not in conversion list or not integrable
            for(unsigned int k=0;k<I->childSize();k++)
                Tree<Index>::childAdd(new IndexGrid(I->child(k),Ax,Grid,Weig));
        }
    }

    // interpolation from grid
    else if (supBas->isGrid()) {
        const BasisGrid* grdBas=dynamic_cast<const BasisGrid*>(supBas);
        for(double g: Grid[nAx])
            if(g<grdBas->mesh()[0] or g>grdBas->mesh().back())
                ABORT(Sstr+"cannot interpolate"+g+"is outside grid ["+grdBas->mesh()[0]+","+grdBas->mesh().back()+"], Index"+I->strNode());
        setBasis(BasisGrid::factory(Grid[nAx]));
        if(not I->isLeaf())
            for(unsigned int k=0;k<basis()->size();k++)
                Tree<Index>::childAdd(new IndexGrid(I->child(0),Ax,Grid,Weig));
    }
    else if (dynamic_cast<const BasisIntegrable*>(supBas)){
        const BasisIntegrable* intBas=dynamic_cast<const BasisIntegrable*>(supBas);
        if(Grid[nAx].size()==0){
            if(Weig[nAx].size()>0)intBas->quadRule(Weig[nAx].size(),Grid[nAx],Weig[nAx]);
            else                  intBas->quadRule(  intBas->order(),Grid[nAx],Weig[nAx]);
        }
        if(Weig[nAx].size() and Weig[nAx].size()!=Grid[nAx].size())DEVABORT("grid points and weights do not match");


        // determine points falling into [low,upp), [low,upp] for last element;
        double upp=intBas->upBound();
        double low=intBas->lowBound();
        double eps=(upp-low)*1.e-10;
        const Index* uppN=I->upperNeighbor();
        if(uppN==0)upp+=eps;   // raise by epsilon for last element

        const Index* lowN=I->lowerNeighbor();
        if(lowN==0)low-=eps;   // lower by epsilon for first element
        else       low=lowN->basis()->integrable()->upBound(); // use upper boundary of lower neighbor

        // select grid points and weights falling into present basis
        std::vector<double> pts,wgs;
        for(double g: Grid[nAx]){
            if(upp<=g)break;
            if(low<=g){
                pts.push_back(g);
                if(Weig[nAx].size())wgs.push_back(Weig[nAx][std::find(Grid[nAx].begin(),Grid[nAx].end(),g)-Grid[nAx].begin()]);
            }
        }
        if(pts.size()>0){
            if(wgs.size())setBasis(BasisGridQuad::factory(pts,wgs));
            else          setBasis(BasisGrid::factory(pts));
            if(not I->isLeaf())
                for(unsigned int k=0;k<pts.size();k++)
                    Tree<Index>::childAdd(new IndexGrid(I->child(0),Ax,Grid,Weig));

            // converted lower indices must be equivalent
            for(size_t k=1;k<I->childSize() and childSize()>0;k++)
                if(not IndexGrid(I->child(k),Ax,Grid,Weig).treeEquivalent(child(0))){
                    Sstr+I->child(0)->strNode()+Sendl;
                    Sstr+I->child(k)->strNode()+Sendl;
                    DEVABORT("branches are not equivalent after conversion - cannot convert to grid on present level "+I->axisName()+"\nchild(k)\n"
                             +IndexGrid(I->child(k),Ax,Grid,Weig).str()
                             +"\nchild(0)\n"+child(0)->str());
                }
        }
        else
            _indexBas=0; // no basis
    }
    else
        PrintOutput::warning("basis on axis is not integrable - will not convert: "+basis()->str(),1);

    // remove empty branches
    std::vector<int> subs;
    for(size_t k=0;k<childSize();k++)subs.push_back(k);
    for(int k=childSize()-1;k>=0;k--){
        if(not child(k)->basis()){
            childErase(k);
            subs.erase(subs.begin()+k);
        }
    }

    // sub-basis if branches were removed
    if(not I->isLeaf() and basis()->size()>subs.size()){
        setBasis(BasisAbstract::factory(BasisSub::strDefinition(BasisSub::superBas(basis()),subs)));
    }
    sizeCompute();
}

void IndexGrid::gridWeight(std::vector<std::vector<double>> &Grid, std::vector<std::vector<double>> &Weig,
                                    std::vector<Integrate::string> Axes, std::vector<unsigned int> Points, std::vector<std::vector<double> > Bounds) {

    if(Points.size()==0)Points.assign(Axes.size(),0);
    if(Bounds.size()==0)Bounds.assign(Axes.size(),{0.,-1.});
    if(Points.size()!=Axes.size())DEVABORT(Sstr+"must supply Points for all Axes"+Axes+"is:"+Points);
    if(Bounds.size()!=Axes.size())DEVABORT(Sstr+"must supply Bounds for all Axes"+Axes+"got only"+Bounds.size());

    Grid.assign(Axes.size(),{});
    Weig.assign(Axes.size(),{});
    for(size_t k=0;k<Axes.size();k++){
        if(Points[k]==1){
            if(Bounds[k][0]!=Bounds[k][1])ABORT(Sstr+"for single point, lower and upper boundary must be equal, is:"+Bounds[k]);
            Grid[k]={Bounds[k][0]};
        }
        else if(Points[k]>1) {
            if(Bounds[k][0]<Bounds[k][1]){
                for(size_t l=0;l<Points[k];l++)
                    Grid[k].push_back(Bounds[k][0]+l*(Bounds[k][1]-Bounds[k][0])/(Points[k]-1));
            } else {
                Weig[k].assign(Points[k],0.);
            }
        }
    }
}
