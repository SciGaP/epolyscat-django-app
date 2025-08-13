#include "operatorRotate.h"

#include "qtEigenSparse.h"
#include "index.h"
#include "basisAbstract.h"
#include "basisExpIm.h"
#include "operatorDefinition.h"
#include "discretizationSpectral.h"
#include "operatorMap.h"
#include "operatorExpandIndex.h"
#include "operatorSubindex.h"
#include "operatorDiagonal.h"

#include "eigenTools.h"
#include "printOutput.h"
#include "timeCritical.h"
#include "parameters.h"
#include "pulse.h"
#include "operatorRotate.h"
#include "parallelOperator.h"
#include "operatorFloor.h"

#include "matrixEigenBlock.h"
#include "basisExpIm.h"
#include "basisAssocLeg.h"

bool OperatorRotate::debug=false;

static double normSq(const Coefficients & C){
    double nrm=0;
    for(size_t k=0;k<C.size();k++)nrm+=std::norm(C.data()[k]);
    return nrm;
}

static void buildMap(OperatorTree* Map, const Eigen::SparseMatrix<std::complex<double>> & Vecs){

    if(Map->idx()->hasFloor() and Map->jdx()->hasFloor()){
        int iF=Map->idx()->posIndex();
        int jF=Map->jdx()->posIndex();
        Eigen::MatrixXcd f(Vecs.block(iF,jF,Map->idx()->size(),Map->jdx()->size()));
        if(f.isZero())Map->floor()=0;
        else  Map->floor()=OperatorFloor::factory({&f},Sstr+"rot"+iF+"."+jF);
    }
    else if(Map->jdx()->hasFloor()){
        for(size_t i=0;i<Map->idx()->childSize();i++){
            Map->childAdd(new OperatorTree("map",Map->idx()->child(i),Map->jdx()));
            buildMap(Map->childBack(),Vecs);
            if(Map->childBack()->isLeaf() and Map->childBack()->floor()==0)Map->childPop();
        }
    }
    else {
        for(size_t j=0;j<Map->jdx()->childSize();j++){
            Map->childAdd(new OperatorTree("map",Map->idx(),Map->jdx()->child(j)));
            buildMap(Map->childBack(),Vecs);
            if(Map->childBack()->isLeaf() and Map->childBack()->floor()==0)Map->childPop();
        }
    }
}


static void rotationMaps(const Index* RotI,
                         const Eigen::SparseMatrix<std::complex<double>> & Duals,
                         const Eigen::SparseMatrix<std::complex<double>> & Vecs,
                         OperatorTree* &fromSpec, OperatorTree* &toSpec)
{
    int lMax=int(sqrt(double(Vecs.cols()))+1.e-5);
    if(lMax*lMax!=Vecs.cols())DEVABORT(Sstr+"number of eigenvalues does not seem to square of integer"+lMax);

    // build index for spectral values Subspace0.Subspace.SubSubspace
    Index* iSpec=new Index({BasisAbstract::factory(Sstr+"Vector: "+lMax)},{"SubSpace"});
    for(size_t k=0;k<iSpec->childSize();k++){
        iSpec->childReplace(k,new Index({BasisAbstract::factory(Sstr+"Vector: "+int(2*k+1))},{"SubSubspace"}));
    }
    iSpec->resetFloor(1);
    iSpec->sizeCompute();

    if(RotI->size()!=Duals.rows()  or  RotI->size()!=Vecs.rows())DEVABORT("length of eigenvector does not match matrices");
    if(iSpec->size()!=Duals.cols() or iSpec->size()!=Vecs.cols())DEVABORT("incorrect spectral index or eigenvalues");

    fromSpec=new OperatorTree("rotFromSpec",RotI,iSpec);
    toSpec=new OperatorTree("rotToSpec",iSpec,RotI);
    buildMap(fromSpec,Vecs);
    buildMap(toSpec,Duals.transpose());

    std::complex<double> one(1.);
    auto m=fromSpec->matrixSparse();
    m-=Vecs;
    m.prune(one);
    if(m.nonZeros()>0)DEVABORT("construction of sparse operator Vecs failed");
    m=Duals.transpose();
    m-=toSpec->matrixSparse();
    m.prune(one);
    if(m.nonZeros()>0)DEVABORT("construction of sparse operator Duals failed");
}

//static Eigen::SparseMatrix<std::complex<double>> extendSparseToSpec(const Eigen::SparseMatrix<std::complex<double>> & To,
//                                                                    const Index* Rot, const Index* In){
//    Eigen::SparseMatrix<std::complex<double>> res(To.rows(),In->size());
//    const BasisTrigon* rPhi=dynamic_cast<const BasisTrigon*>(Rot->basis());
//    const BasisTrigon* iPhi=dynamic_cast<const BasisTrigon*>(In->basis());
//    if(rPhi==0 or iPhi==0)DEVABORT("does not start with Trigon type basis");

//    const Index* iEta0=In->axisIndex("Eta");
//    for(const Index* iEta=iEta0;iEta;iEta=iEta->nodeRight()){
//        if(iEta->size()!=iEta0->size())DEVABORT("only for strict tensor products");


////        const BasisAssocLeg* rEta=dynamic_cast<const BasisAssocLeg*>(Rot->descend()->basis());
//        // m-range of In will be smaller or equal as Rot
//        for(size_t m=0;m<In->childSize();m++){
//            size_t cLI=In->child(m)->posIndex(In);
//            size_t cLR=Rot->child(m)->posIndex(Rot);
//            // l-range of Rot will be smaller or equal as In, positions increment by respective branch sizes
//            for(size_t l=0;l<Rot->child(m)->childSize();cLR+=Rot->child(m)->child(l)->size(),l++)
//            {
//                // expand for columns of In
//                for(size_t k=0;k<In->child(m)->child(l)->size();k++){
//                    // descend non-zero parts of row
//                    for(Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(To,cLR+l);it;++it){
//                        res.insert(it.row()*iEta0->size()+k,cLI+k)=it.value();
//                    }
//                }

//            }
//        }
//    }

//    return res;
//}


//static Eigen::SparseMatrix<std::complex<double>> sparseToSpec(const Eigen::SparseMatrix<std::complex<double>> & DVecT, std::vector<int> Subset){

//    Eigen::SparseMatrix<std::complex<double>> res(DVecT.rows()*2,Subset.size());

//    // storage for transformation to Spec (with constraint on maximal eigenvalue)
//    std::vector<int> inner;
//    for(int j=0;j<DVecT.cols();j++){
//        inner.push_back(DVecT.col(j).nonZeros());
//        inner.push_back(inner.back());
//    }
//    res.reserve(inner);

//    int pSub=0;
//    for (int k=0; k<DVecT.outerSize(); ++k){
//        Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(DVecT,k);
//        while(pSub<int(Subset.size()) and it.col()<Subset[pSub])pSub++;
//        if(pSub!=int(Subset.size())){
//            for(Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(DVecT,k);it;++it)
//            {
//                res.insert(it.row(),Subset[pSub]*2  )=it.value();
//                res.insert(it.row(),Subset[pSub]*2+1)=it.value();
//            }
//        }
//    }
//    return res;
//}


static std::vector<const Index*> _tested;
bool OperatorRotate::test(const Index *Idx){
    return true;
    if(std::find(_tested.begin(),_tested.end(),Idx)!=_tested.end())return true;
    _tested.push_back(Idx);

    timeCritical::suspend();
    // x,y,z axes
    std::vector<double> axX({1.,0.,0.});
    std::vector<double> axY({0.,1.,0.});
    std::vector<double> axZ({0.,0.,1.});

    // alternative axes with same handedness as x,y,z
    axX={1,1,1};
    axY={1,-1,0};
    axZ={1,1,-2};

    // rotate by pi/2 around y
    Coefficients c0(Idx);
    OperatorRotate rotY(axY,c0.idx());
    OperatorRotate rotZ(axZ,rotY.idx());
    OperatorRotate rotX(axX,rotY.idx());

    Coefficients cR(rotY.idx()),cS(rotY.idx());
    c0.setToRandom();
    std::vector<int> s;
    for(int p=3;p<5;p++){
        s={2*(p%2)-1,2*((p/2)%2)-1,2*((p/4)%2)-1};
        rotY.angle(math::pi/2).apply(1,c0,0,cR);
        rotZ.angle(s[0]*math::pi/2).apply(1.,cR,0.,cS);
        rotX.angle(s[1]*math::pi/2).apply(1.,cS,0.,cR);
        rotZ.angle(s[2]*math::pi/2).apply(1.,cR,0.,cS);
        cR=cS;
        rotY.angle(0).apply(-1,c0,1.,cR);
        if(not cR.isZero(1.e-12) or std::abs(normSq(cS)-normSq(c0))>normSq(c0)*1e-12){
#ifdef _DEVELOP_
            Sstr+s+tools::str(cR.values(),0,"      ")+Sendl;
            Sstr+"norms"+normSq(cS)+normSq(c0)+Sendl;
#endif
            DEVABORT("OperatorRotation failed test for "+Idx->hierarchy());
        }
    }
    timeCritical::resume();
    return cR.isZero(1e-10);
}

bool OperatorRotate::testCoriolis(const Index* Idx) {
    return true;
    if(std::find(_tested.begin(),_tested.end(),Idx)!=_tested.end())return true;
    _tested.push_back(Idx);

    // verify for time-constant Psi with U(theta)=exp(-i theta Ly), U^*=exp(i theta Ly)
    // U=U(-theta) is defined such that
    // i|A| d/dz U(-theta) Psi = U [iA(t).Nabla] Psi
    // d/dt U(-theta) Psi = d/dt exp(i*theta*Ly) Psi = U(-theta) (idTheta/dt Ly) Psi

    // set up rotation around y
    OperatorRotate rotY({0.,1.,0.},Idx);
    OperatorTree angLy("angLy",OperatorDefinition("<<AngularMomentumY>>",Idx->hierarchy()),Idx,Idx);

    // set up fixed derivative
    std::string defFix("iLaserAx[t]<<D/DX>>+iLaserAz[t]<<D/DZ>>");
    std::string defRot("iLaserAabs[t]<<D/DZ>>");

    Parameters::updateToOne();
    OperatorTree dipFix("dipFix",defFix,Idx,Idx);
    OperatorTree dipRot("dipRot",defRot,Idx,Idx);

    Coefficients psi(Idx),uPsi(Idx),dipPsi(Idx),dChi(Idx),udChi(Idx);
    psi.setToRandom();
    psi.makeContinuous();
    for(double t=0.;t<100.;t+=10.){
        // i|A| d/dz U Psi = U [iA(t).Nabla] Psi
        rotY.angle(-Pulse::thetaA(t)).apply(1.,psi,0.,uPsi);
        dipRot.time(t).apply(1.,uPsi,0.,dChi);
        dipFix.time(t).apply(1.,psi,0.,dipPsi);
        rotY.angle(-Pulse::thetaA(t)).apply(1.,dipPsi,-1,dChi);
        if(not dChi.isZero(1e-12))ABORT("not zero");

        // i d/dt U^* Psi = U^* dTheta/dt Ly) Psi
        dipPsi.setToConstant(1.);
        double dTheta;
        for(double h=0.01;h>1e-7;h*=0.1){
            // numerical differentiation
            std::vector<double> derT({-h,h}),derW({-0.5/h,0.5/h});

            // check consistency of thetaA
            dTheta=0.;
            for(size_t d=0;d<derT.size();d++)
                dTheta+=Pulse::thetaA(t+derT[d])*derW[d];

            uPsi.setToZero();
            for(size_t d=0;d<derT.size();d++)
                rotY.angle(-Pulse::thetaA(t+derT[d])).apply(derW[d],psi,1.,uPsi);
            dipPsi-=uPsi;
            if(dipPsi.normL2sq()<1.e-19)break;
            dipPsi=uPsi;

        }
        rotY.angle(-Pulse::thetaA(t)).apply(1,psi,0,dipPsi);
        angLy.apply(Pulse::current.idThetaA(t),dipPsi,0,udChi);

        if(((dChi=udChi)-=uPsi).normL2sq()>1.e-12)
            DEVABORT(Sstr+"failed test i d/dt U^* Psi = U^* dTheta/dt Ly) Psi  by ||Deltat(Psi)||="
                     +dChi.normL2sq());
        if(std::abs(dTheta-Pulse::current.dThetaA(t))>1e-7*std::max(1.,std::abs(dTheta)))
            DEVABORT(Sstr+"theta derivative "+Pulse::current.dThetaA(t)
                     +"inconsistent with numerical value"+dTheta);
    }
    PrintOutput::DEVmessage("pass testCoriolis");
    return true;
}

static Eigen::SparseMatrix<std::complex<double>> putSparse(const OperatorTree & Op){
    std::vector<Eigen::Triplet<std::complex<double>>> list;
    for(const OperatorTree* op=Op.firstLeaf();op!=0;op=op->nextLeaf()){
        Eigen::MatrixXcd mat=op->matrix();
        int i0=op->idx()->posIndex();
        int j0=op->jdx()->posIndex();
        for(int j=0;j<mat.cols();j++)
            for(int i=0;i<mat.rows();i++){
                if(std::abs(mat(i,j))>1.e-12)
                    list.push_back(Eigen::Triplet< std::complex<double> >(i0+i,j0+j,mat(i,j)));
            }
    }
    Eigen::SparseMatrix<std::complex<double>> res(Op.idx()->size(),Op.jdx()->size());
    res.setFromTriplets(list.begin(),list.end());
    PrintOutput::DEVmessage(Sstr+"non-zeros"+Op.name+":"+list.size()+"of"+(Op.idx()->size()*Op.jdx()->size()));
    return res;
}

TIMER(rotate,);
TIMER(rotate1,);
TIMER(rotate2,);
TIMER(rotate3,);
TIMER(rotate4,);
TIMER(rotate5,);
OperatorRotate::OperatorRotate(std::vector<double> RotationAxis, const Index* Jdx)
    :OperatorAbstract("Rot("+tools::str(RotationAxis,3,",")+")",0,Jdx)
{
    // seek Phi and Eta axes
    size_t lMax=0;
    for(const Index* idx=Jdx;idx!=0;idx=idx->descend()){
        if(idx->axisName()=="Eta"){
            if(idx->basis()->name().find("assocLegendre[")!=0)
                DEVABORT("need AssocLeg, got: "+idx->basis()->name());
            for(const Index* edx=idx;edx!=0;edx=edx->nodeRight()){
                for(size_t k=0;k<edx->basis()->size();k++)
                    lMax=std::max(lMax,size_t(edx->basis()->physical(k)));
            }
            break;
        }
    }

    // full rotation index
    IndexFull* idxF=new IndexFull(Jdx,lMax);
    iIndex=idxF->treeEquivalent(jIndex)?jIndex:idxF;
    Index* rotI(idxF->rotIndex());

    // get generator of rotation around RotationAxis wrt full basis
    double nrm=0.;
    for(auto a: RotationAxis)nrm+=a*a;
    std::string rotDef,rotDef2;
    std::vector<std::string> coor({"X","Y","Z"});
    std::vector<std::string> defs;
    defs.push_back(OperatorDefinition::standardOperators["AngularMomentumX"]["Phi.Eta"]);
    defs.push_back(OperatorDefinition::standardOperators["AngularMomentumY"]["Phi.Eta"]);
    defs.push_back(OperatorDefinition::standardOperators["AngularMomentumZ"]["Phi.Eta"]);
    for(size_t k=0;k<RotationAxis.size();k++){
        double a=RotationAxis[k]/sqrt(nrm);
        if(a!=0.)rotDef +=(a>0?"+":"")+tools::str(a)+defs[k]+"<1>";
    }
    LOG_PUSH("Rot1");
    OperatorTree rotOp("Rot",rotDef,rotI,rotI);
    LOG_POP();
    LOG_PUSH("Rot2");
    ParallelOperator pRot(&rotOp);
    LOG_POP();
    LOG_PUSH("Rot3");
    pRot.bcast();
    DiscretizationSpectral rotD(&rotOp);
    LOG_POP();
    LOG_PUSH("Rot4");


    Eigen::MatrixXcd vals;
    Eigen::SparseMatrix<std::complex<double>> vecs,duals;
    tRecX::eigenBlock(rotOp.matrixSparse(true),vals,true,vecs,true,duals,
                      Eigen::SparseMatrix<std::complex<double>>(0,0),1);
    OperatorTree* _from,*_to;
    rotationMaps(rotI,duals,vecs,_from,_to);

    // map from spectral to full index
   _fromSpec.reset(new OperatorExpandIndex(dynamic_cast<const OperatorTree*>(rotD.mapToParent()),iIndex,true));
//    _fromSpec.reset(new OperatorExpandIndex(dynamic_cast<const OperatorTree*>(_from),iIndex,true));

    LOG_POP();
    LOG_PUSH("Rot5");
    // extend spectral map to full index and restrict map to original index
    OperatorExpandIndex fullToSpec(dynamic_cast<const OperatorTree*>(rotD.mapFromParent()),idxF,false);
//    OperatorExpandIndex fullToSpec(dynamic_cast<const OperatorTree*>(_to),idxF,false);
    _toSpec.reset(new OperatorSubindex(&fullToSpec,Jdx,false));
    _toSpec->replaceIndex(_fromSpec->jdx(),jIndex); // indices were constructed twice
    delete fullToSpec.idx(); // this index was replaced
    if(iIndex==jIndex)delete idxF; // no expanded Index needed to be constructed

    _cSpec.reset(new Coefficients(_toSpec->idx()));

    // expand the vector of eigenvalues onto product space
    Coefficients specVal(rotD.idx(),rotD.eigenvalues());
    Coefficients specX(fullToSpec.coefficients(specVal,_toSpec->idx()));
    LOG_POP();
    LOG_PUSH("Rot6");

    // set pointers to phases storage
    _phases.resize(rotI->basis()->size(),1.);
    for(size_t k=0;k<specX.size();k++){
        int m=std::round(specX.data()[k].real());
        if(std::abs(double(m)-specX.data()[k].real())>1.e-10){
            DEVABORT(Sstr+"non-integer eigenvalues at k="+k+m+specX.data()[k].real()+specX.values());
        }
        else
            specX.data()[k]=double(m);
        _phaseK.push_back((m<=0?-2*m:2*m-1));
    }
    LOG_POP();
    LOG_PUSH("Rot7");

    _fromSpecSparse=putSparse(*_fromSpec);
    _toSpecSparse  =putSparse(*_toSpec);
    _fromSpec=0;
    _toSpec=0;
    LOG_POP();
}

const OperatorRotate& OperatorRotate::angle(double Angle) const {
    //  update phases
    _phases[0]=1.;
    _phases[1]=std::exp(std::complex<double>(0.,-Angle));
    _phases[2]=1./_phases[1];
    for(size_t k=3;k<_phases.size();k=k+2){
        _phases[k  ]=_phases[k-2]*_phases[1];
        _phases[k+1]=_phases[k-1]*_phases[2];
    }
    return *this;
}
void OperatorRotate::apply(std::complex<double> A, const Coefficients &Vec, std::complex<double> B, Coefficients &Y) const{

#ifdef _DEVELOP_
    if(not test(Vec.idx()))ABORT("OperatorRotate: test failed for hierarchy "+Vec.idx()->hierarchy());
#endif

    // map to spectral rep of rotation
    if(_toSpec!=0)
        _toSpec->apply(A,Vec,0,*_cSpec);
    else
        Eigen::Map<Eigen::VectorXcd>(_cSpec->data(),_cSpec->size())=
                A*_toSpecSparse*Eigen::Map<const Eigen::VectorXcd>(Vec.data(),Vec.size());

    // multiply by rotation phases
    for(size_t k=0;k<_cSpec->size();k++)_cSpec->data()[k]*=_phases[_phaseK[k]];

    // map back to spherical basis
    if(_fromSpec!=0){
        _fromSpec->apply(1,*_cSpec,B,Y);
    }
    else {
        Y.scale(B);
        Eigen::Map<Eigen::VectorXcd>(Y.data(),Y.size())+=
                _fromSpecSparse*Eigen::Map<const Eigen::VectorXcd>(_cSpec->data(),_cSpec->size());
    }
}

OperatorRotate::IndexFull::IndexFull(const Index* Idx, size_t Lmax, int M){
    // expand to full rotational basis, RotIdx contains only rotational parts
    nodeCopy(Idx,false);

    if(Idx->axisName()=="Phi"){
        std::string def=Idx->basis()->name()+": { 0";
        for(int m=1;m<=int(Lmax);m++)def+=","+tools::str(m)+","+tools::str(-m);
        def+=" } ["+tools::str(2*Lmax+1)+"]";
        setBasis(BasisAbstract::factory(def));
        for(size_t k=0;k<basis()->size();k++)
            childAdd(new IndexFull(Idx->child(0),Lmax,basis()->physical(k)));
    }

    else if (Idx->axisName()=="Eta"){
        //        if(Idx->parent()->axisName()!="Phi")ABORT("RotateCoefficients: need Phi.Eta sequence in hierarchy, have: "+Idx->root()->hierarchy());
        if(not Idx->subEquivalent())DEVABORT("non-equivalent sub-indices at Eta - cannot expand to full rotational space");

        std::string def="AssociatedLegendre: "+tools::str(std::abs(M))+","+tools::str(1+Lmax-std::abs(M));
        setBasis(BasisAbstract::factory(def));
        for(size_t k=0;k<basis()->size();k++)
            childAdd(new IndexFull(Idx->child(0),Lmax,M));
    }
    else
        for(size_t k=0;k<Idx->childSize();k++)
            childAdd(new IndexFull(Idx->child(k),Lmax,M));

    if(M==INT_MAX){
        //        if(hierarchy().find("Phi.Eta")==std::string::npos)
        //            ABORT("cannot rotate basis, need axes Phi.Eta in hierarchy, got: "+Idx->hierarchy());
        sizeCompute();
        resetFloor(Idx->heightAboveFloor());
        buildOverlap();
    }
}

// attach single-coefficient floor level
static void addVec1(Index* Idx){
    for(size_t k=0;k<Idx->basis()->size();k++){
        Idx->childAdd(new Index());
        Idx->childBack()->setAxisName("Vec");
        Idx->childBack()->setBasis(BasisAbstract::factory("Vector: 1"));
        Idx->childBack()->setFloor();
    }
}
// extract the rotational part from the full Index
IndexNew* OperatorRotate::IndexFull::rotIndex() const{
    IndexNew* res=new IndexNew();
    const Index* idx=this;
    while(idx->axisName()!="Phi")idx=idx->descend();
    res->nodeCopy(this,false);
    size_t siz=idx->childSize();
    while(idx->axisName()!="Eta")idx=idx->descend();
    while(idx and res->childSize()!=siz){
        res->childAdd(new IndexNew());
        res->childBack()->nodeCopy(idx,false);
        addVec1(res->childBack());
        idx=idx->nodeRight();
    }
    res->sizeCompute();
    res->buildOverlap();
    return res;
}


