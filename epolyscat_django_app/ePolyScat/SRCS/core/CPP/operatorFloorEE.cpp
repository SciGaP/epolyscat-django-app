// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "operatorFloorEE.h"
#include "printOutput.h"
#include "algebra.h"
#include "basisAbstract.h"
#include "gaunt.h"
#include "coefficientsFloor.h"
#include "index.h"
#include "qtAlglib.h"
#include "basisDvr.h"
#include "inverseDvr.h"
#include "operatorFloorEE.h"
#include "basisMat1D.h"
using namespace std;
#include "eigenNames.h"
#include "indexNew.h"
#include "overlapDVR.h"


std::vector<std::vector<std::vector<UseMatrix > > > OperatorFloorEE::Helpers_e_e::TwoElecOverlap(0);
std::vector<std::vector<std::vector<double > > > OperatorFloorEE::Helpers_e_e::TwoElecOverlapMaxValue(0);
double OperatorFloorEE::epsGaunt=1.e-14;
double OperatorFloorEE::epsTwoElec=1.e-14;
int OperatorFloorEE::lambda_upper_limit;
static int lambda_upper_limit_default=4;

OperatorFloorEE::BasicDiscOnlyRadial::BasicDiscOnlyRadial(ReadInput &In)
{
    bool oldDVR=BasisDVR::femDVR;
    BasisDVR::femDVR=true;
    PrintOutput::DEVwarning("BasisDVR::femDVR forced = true in BasicDiscOnlyRadial -- really bad hack");
    setAxis(In);

    lmax = 0;

    vector<Axis> ax;
    std::vector<std::string> hier;

    // Remove particle defined by remove
    for(unsigned int i=0;i<axis.size();i++){
        if(axis[i].name.find("Eta")!=std::string::npos){
            if(axis[i].basDef.size()!=1) ABORT("Size not 1");
            lmax = max(lmax,(int)(axis[i].basDef[0].order-1));
        }
        if(axis[i].name.find("Eta")!=std::string::npos or axis[i].name.find("Phi")!=std::string::npos){}
        else{
            ax.push_back(axis[i]);
            hier.push_back(hierarchy[i]);
        }
    }
    axis = ax;

    hierarchy = hier;
    if(axis.size()!=2)ABORT("Axis size is not 2: "+tools::str(axis.size()));
    if(axis[0].basDef.size() != axis[1].basDef.size()) ABORT("Only equal radial bases supported");
    for(size_t i=0; i<axis[0].basDef.size(); i++){
        if(axis[0].basDef[i].order != axis[1].basDef[i].order) ABORT("Only equal radial bases supported");
    }
    continuityLevel.clear();
    continuityLevel.push_back(0);
    continuityLevel.push_back(1);

    bool previousUse=IndexNew::doNotUse;
    IndexNew::doNotUse=true;

    if(axis.size()>0){
        _axisTree.reset(new AxisTree(axis));
        construct();
    }
    IndexNew::doNotUse=previousUse;
    BasisDVR::femDVR=oldDVR;
}

void OperatorFloorEE::Helpers_e_e::absorbInverse(const Index *Root){

    // locate element levels
    unsigned int lev0=Index::npos,lev1=Index::npos;
    const Index* lev=Root;
    while(not lev->isLeaf()){
        if(lev->axisName()=="Rn1")lev0=lev->depth()-1;
        if(lev->axisName()=="Rn2")lev1=lev->depth()-1;
        lev=lev->descend();
    }
    if(lev0==Index::npos)ABORT(Root->str()+"/nRn1 or Rn2 axis missing");

    // multiply weight on helpers
    vector<unsigned int> idx(6,0);
    for(unsigned int n0=0;n0<Helpers_e_e::TwoElecOverlap.size();n0++)
        for(unsigned int n1=0;n1<Helpers_e_e::TwoElecOverlap[n0].size();n1++)
        {
            idx[lev0]=n0;
            idx[lev1]=n1;
            const InverseDVR* invDVR=dynamic_cast<const InverseDVR*>(Root->nodeAt(idx)->inverseOverlap());
            if(invDVR->iIndex->sizeStored()!=TwoElecOverlap[n0][n1][0].size()){
                ABORT("sizes do not match: "+tools::str(invDVR->iIndex->sizeStored())
                      +" "+tools::str(int(TwoElecOverlap[n0][n1][0].size())));
            }
            for(unsigned int la=0;la<Helpers_e_e::TwoElecOverlap[n0][n1].size();la++)
                for(unsigned int k=0;k<invDVR->iIndex->sizeStored();k++)
                    Helpers_e_e::TwoElecOverlap[n0][n1][la]*=invDVR->invWeig(k);
        }
}

int OperatorFloorEE::read(ReadInput& Inp){
    Inp.obsolete("OperatorFloorEE", "lambdaMax","use Operator: multipoleOrder=VALUE instead");
    Inp.read("Operator", "multipoleOrder", lambda_upper_limit, "-1", "Maximum number of l-terms in multipole expansion of ee interaction");
    return lambda_upper_limit==-1?lambda_upper_limit_default:lambda_upper_limit;
}

void OperatorFloorEE::Helpers_e_e::initialize(const Index *Root, int& lambda_upper_limit)
{
    if(TwoElecOverlap.size()!=0)  return;      // Already initialized, should only be done once

    BasicDiscOnlyRadial *D = new BasicDiscOnlyRadial(ReadInput::main);

    if(lambda_upper_limit == -1){
        lambda_upper_limit =std::min(lambda_upper_limit_default,2*D->lmax);
        PrintOutput::message(Sstr+"Truncate ee-Interaction at default multipoleOrder="+lambda_upper_limit);
        PrintOutput::message(" default usually suffices, change by Operator: multipoleOrder=VALUE");
    }
    else if(lambda_upper_limit < 2*D->lmax){
        PrintOutput::message(Sstr+"Truncate ee-Interaction at orderMultipole="+lambda_upper_limit);
    }
    else {
        PrintOutput::message(Sstr+"Exact ee-Interaction: save compute time by setting Operator:orderMultipole");
    }
    lambda_upper_limit=std::min(lambda_upper_limit,2*D->lmax);

    if(lambda_upper_limit>2*lambda_upper_limit_default){
        PrintOutput::warning(Sstr+"Large multipoleOrder="+lambda_upper_limit,1,0,
                             Sstr+"for most problems default multipoleOrder="+lambda_upper_limit_default+"suffices"
                             +"\nif in doubt, de- or increase value and check convergence");
    }
    if(lambda_upper_limit<std::min(4,2*D->lmax))
        PrintOutput::warning(Sstr+"Small Operator:multipoleOrder="+lambda_upper_limit+" ee-interaction may be inaccurate",1,0,
                             "too small multipole order may give seemingly good ground state energies"
                             "\nor energies below the variational limit");

    if(TwoElecOverlap.size()!=0)  return;      // Already initialized, should only be done once

    vector<UseMatrix > Kin_mat,Q2_mat;
    vector<vector<bool> > WgtIncludedInFunction;
    vector<VectorXd> points,weig;

    int ne = -1;  // Number of unscaled elements set below
    double r_max=-1;
    int nInterval;
    int NR1(0), NR2(0),k = 0;
    while(not Root->descend(k)->hasFloor()){
        if(Root->descend(k)->axisName() == "Rn1")NR1 = Root->descend(k)->childSize();
        else if(Root->descend(k)->axisName() == "Rn2")NR2 = Root->descend(k)->childSize();
        k++;
    }

    if(NR1 != NR2)ABORT("The EE interaction only supports the same interval numbers of the two electron")
            nInterval = NR1;
    Gaunt::main.initialize(2*D->lmax+1);          // Initialize gaunt

    // determine possible truncation
    double rTrunc=DBL_MAX;
    if(Algebra::isSpecialConstant("Rtrunc"))rTrunc=Algebra("Rtrunc").val(0.).real();

    for(ne=0;ne<nInterval;ne++){
        vector<bool> wgs_inc;
        UseMatrix pts1,wgs1;
        const Index * lev = Root;
        //n1,n2,m1,m2,l1,l2
        while(not lev->hasFloor()){
            if(lev->axisName()=="Rn1")lev = lev->child(ne);
            else{lev = lev->descend();}
        }
        const BasisAbstract * B = lev->basis()->integrable();
        if(B->integrable()->lowBound()+1.e-7>rTrunc)break; // outside truncation
        if(B->isAbsorptive() or B->integrable()->upBound()>DBL_MAX/2.) break;  // at absorption or infinite interval

        analyseDVR(B,pts1,wgs1,wgs_inc);
        VectorXd pts = VectorXd::Zero(wgs_inc.size());
        VectorXd wgs = VectorXd::Zero(wgs_inc.size());
        for(unsigned int k=0;k<wgs_inc.size();k++) {pts(k) = pts1(k).real(); wgs(k) = wgs1(k).real();}
        points.push_back(pts);
        weig.push_back(wgs);
        UseMatrix kin,Q2;
        BasisMat1D * b1 = new BasisMat1D("<d_1_d>", B, B);
        kin = b1->useMat();
        BasisMat1D * b2 = new BasisMat1D("<1/(Q*Q)>", B, B);
        Q2 = b2->useMat();
        Kin_mat.push_back(kin);
        Q2_mat.push_back(Q2);
        WgtIncludedInFunction.push_back(wgs_inc);

        r_max = B->integrable()->upBound();
    }

    // Set size
    TwoElecOverlap.resize(ne);
    TwoElecOverlapMaxValue.resize(ne);
    for(int n=0;n<ne;n++){
        TwoElecOverlap[n].resize(ne);
        TwoElecOverlapMaxValue[n].resize(ne);
        for(int n2=0;n2<ne;n2++){
            TwoElecOverlap[n][n2].resize(lambda_upper_limit+1);
            TwoElecOverlapMaxValue[n][n2].resize(lambda_upper_limit+1);
        }
    }

    // Smooth Coulomb
    AlgebraTrunc* smooth=0;
    if(Algebra::isSpecialConstant("Rtrunc")){
        if(not Algebra::isSpecialConstant("Rsmooth"))ABORT("for using potential cutoff radius Rcutoff, must also define Rsmooth")
                string smoothFunction  ="trunc[Rsmooth,Rtrunc]";
        smooth = new AlgebraTrunc(smoothFunction);
    }

    //Setting up the overall T_Vee matrix
    for(int lambda=0;lambda<=lambda_upper_limit;lambda++){

        int size=0;
        for(int n=0;n<ne;n++)
            size+=Q2_mat[n].rows();
        size=size-ne+1;
        MatrixXcd T_Vee = MatrixXcd::Zero(size,size);

        int start=0;
        for(int n=0;n<ne;n++){
            MatrixXcd tp = Map<MatrixXcd>(Kin_mat[n].data(),Kin_mat[n].rows(),Kin_mat[n].cols())
                    +(double)lambda*((double)lambda+1)*Map<MatrixXcd>(Q2_mat[n].data(),Q2_mat[n].rows(),Q2_mat[n].cols());

            T_Vee.block(start,start,tp.rows(),tp.cols())+=tp;
            start=start+tp.rows()-1;
        }

        MatrixXcd temp = T_Vee.block(0,0,T_Vee.rows()-1,T_Vee.cols()-1);   //Removing the last function; boundary conditions
        T_Vee = MatrixXcd::Zero(T_Vee.rows(), T_Vee.cols());
        T_Vee.block(0,0,T_Vee.rows()-1,T_Vee.cols()-1) = temp.inverse();

        //distrubute it back to each of the voxels
        int start1 =0;
        for(int n1=0;n1<ne;n1++){
            int start2=0;
            for(int n2=0;n2<ne;n2++){
                MatrixXcd T_Vee_all_parts = T_Vee.block(start1,start2,Q2_mat[n1].rows(),Q2_mat[n2].cols());
                TwoElecOverlap[n1][n2][lambda] = UseMatrix::Zero(T_Vee_all_parts.rows()*T_Vee_all_parts.cols(),1);
                MatrixXcd temp2e = MatrixXcd::Zero(T_Vee_all_parts.rows(),T_Vee_all_parts.cols());
                unsigned int kk12=0;
                for(int k1=0;k1<T_Vee_all_parts.rows();k1++){
                    for(int k2=0;k2<T_Vee_all_parts.cols();k2++){
                        double r_j = points[n2](k2);
                        double r_i = points[n1](k1);
                        double w_j = weig[n2](k2);
                        double w_i = weig[n1](k1);

                        kk12++;

                        if(n1==ne-1 and k1==T_Vee_all_parts.rows()-1) continue;
                        if(n2==ne-1 and k2==T_Vee_all_parts.cols()-1) continue;

                        complex<double> term1 = pow(r_j*r_i,lambda)/pow(r_max,2*lambda+1);
                        complex<double> term2 = (2*(double)lambda+1)*T_Vee_all_parts(k1,k2)/(r_j*r_i);

                        if(WgtIncludedInFunction[n1][k1] and WgtIncludedInFunction[n2][k2])
                            temp2e(k1,k2) = term1 + term2/sqrt(w_i*w_j);
                        else if(WgtIncludedInFunction[n1][k1] and not WgtIncludedInFunction[n2][k2])
                            temp2e(k1,k2) = term1 * w_j + term2 *w_j/sqrt(w_i);
                        else if(not WgtIncludedInFunction[n1][k1] and WgtIncludedInFunction[n2][k2])
                            temp2e(k1,k2) = term1 * w_i + term2*w_i/sqrt(w_j);
                        else
                            temp2e(k1,k2) = (term1 + term2) * w_i*w_j;

                        if(smooth)
                            TwoElecOverlap[n1][n2][lambda](k1*T_Vee_all_parts.cols()+k2) = temp2e(k1,k2)*smooth->val(r_i)*smooth->val(r_j);
                        else
                            TwoElecOverlap[n1][n2][lambda](k1*T_Vee_all_parts.cols()+k2) = temp2e(k1,k2);

                    }

                    TwoElecOverlapMaxValue[n1][n2][lambda] = TwoElecOverlap[n1][n2][lambda].maxAbsVal();
                }
                start2 = start2+Q2_mat[n2].cols()-1;
            }
            start1 = start1+Q2_mat[n1].rows()-1;
        }
    }

    if(smooth!=0) delete smooth;
}

const std::vector<std::complex<double>> overlapWeights(const Index* Idx){
    if(!Idx->hasFloor())DEVABORT("only on floor level")
            const Index* idx=Idx;
    // ascend index until overlap is found, record path
    const OperatorTree*oDvr;
    std::vector<unsigned int> pos;
    while(0!=idx and 0==(oDvr=dynamic_cast<const OperatorTree*>(idx->overlap()))){
        pos.insert(pos.begin(),idx->nSibling());
        idx=idx->parent();
    }
    if(idx==0)DEVABORT(Sstr+"no OverlapDVR found above"+Idx->strNode()+"root="+Idx->root()+Idx->root()->strNode())

            // descend on oDvr to present
            const OperatorTree* oRoot=oDvr;
    oDvr=oRoot->nodeAt(pos);
    if(!oDvr)DEVABORT(Sstr+oRoot->str()+"\ndoes not have block in level"+Idx->strNode());

    Coefficients diag=oDvr->diagonal(true);
    std::vector<std::complex<double>> res(diag.data(),diag.data()+diag.size());
    return res;
}


OperatorFloorEE::OperatorFloorEE(const std::string Name, const std::string Def, const Index *IIndex, const Index *JIndex)
    :OperatorFloor(IIndex->sizeCompute(),JIndex->sizeCompute(),"FloorEE")
{
    Index::build=true;


    // Find the indices n1,l1,m1,n2,l2,m2
    construct_index(IIndex,Idx);
    construct_index(JIndex,Jdx);

    addWeights(IIndex,overlapWeights(IIndex));

    Helpers_e_e::initialize(IIndex->root(), lambda_upper_limit);

    if(_rows+_cols>100000)ABORT("bad matrix size");

    oNorm=DBL_MAX; // generic value - not suitable for controlling applications
    // use oNorm for indicating 0
    if(isZero())oNorm=0.;

    Index::build=false;
}


void OperatorFloorEE::Helpers_e_e::lobatto_quadratureNew(const BasisAbstract *B, VectorXd &quadraturePoints, VectorXd &Qweights)
{
    int N= B->integrable()->order();
    double lb = B->integrable()->lowBound();
    double ub = B->integrable()->upBound();

    quadraturePoints = VectorXd::Zero(N);
    Qweights = VectorXd::Zero(N);
    alglib::real_1d_array alpha;
    alglib::real_1d_array beta;
    alglib::real_1d_array x1;
    alglib::real_1d_array w1;
    alpha.setlength(N-1);
    beta.setlength(N-1);
    x1.setlength(N);
    w1.setlength(N);
    alglib::ae_int_t info;

    for(int i=0;i<N-1;i++)
        alpha(i)=0;

    for(int i=0;i<N-1;i++)
        beta(i)=double(i*i)/double(4*i*i-1);

    alglib::gqgenerategausslobattorec(alpha, beta, (double)2.0, (double)(-1.0) , (double)1, N, info, x1, w1);

    if (info!=1)
        cout<<"Error in Gauss-Lobatto quadrature. Error code "<<info<<endl;

    double temp = (ub-lb)/2.0;
    for(int i=0;i<N;i++){
        quadraturePoints(i)=lb+(x1[i]+1.0)*temp;
        Qweights(i)=w1[i]*temp;
    }
}

void OperatorFloorEE::axpy(const std::complex<double> & Alfa, const std::complex<double>* X, unsigned int SizX,
                           const std::complex<double> & Beta, std::complex<double>* Y, unsigned int SizY) const
{
    vector<complex<double> > x(X,X+SizX);
    vector<complex<double> > y(Y,Y+SizY);
    axpy(Alfa,x,Beta,y);
    for(size_t k=0;k<SizY;k++)Y[k]=y[k];
}

TIMERSAMPLE(axpyEE,)

void OperatorFloorEE::axpy(std::complex<double> Alfa, const vector<complex<double> > &X, std::complex<double> Beta, vector<complex<double> > &Y) const
{
    if(Idx[0]!=Jdx[0] or Idx[1]!=Jdx[1])ABORT("not diagonal in n1,n2 "+tools::str(Idx)+" - "+tools::str(Jdx));

    if(Beta!=1.){
        if(Beta==0.)for(unsigned int k=0;k<Y.size();k++)Y[k]=0.;
        else        for(unsigned int k=0;k<Y.size();k++)Y[k]*=Beta;
    }

    if((size_t)Idx[0]>=Helpers_e_e::TwoElecOverlap.size() || (size_t)Idx[1]>=Helpers_e_e::TwoElecOverlap[0].size()) return;
    if(Jdx[4]-Idx[4]!=Idx[5]-Jdx[5]) return;
    int M = Jdx[4]-Idx[4];

    //    STARTDEBUG(axpyEE);

    for(int lambda=0;lambda<=lambda_upper_limit;lambda++){

        if(Gaunt::main.coeff_isZero(Jdx[2],lambda,Idx[2],Jdx[4],M,Idx[4]) or Gaunt::main.coeff_isZero(Idx[3],lambda,Jdx[3],Idx[5],M,Jdx[5])) continue;

        double gf1,gf2;    // gaunt coefficients
        gf1 = Gaunt::main.coeff(Jdx[2],lambda,Idx[2],Jdx[4],M,Idx[4]);
        if(abs(gf1) < epsGaunt) continue;
        gf2 = Gaunt::main.coeff(Idx[3],lambda,Jdx[3],Idx[5],M,Jdx[5]);
        double gaunt_factor = gf1*gf2;

        if(abs(gaunt_factor)<epsGaunt) continue;

        // Comment this line out to disable different exapnsions at different patches
        if(Helpers_e_e::TwoElecOverlapMaxValue[Idx[0]][Idx[1]][lambda] < epsTwoElec) continue;

        gaunt_factor *= (4.0*math::pi/(2.0*double(lambda)+1.0));

        Map<VectorXcd>(Y.data(),Y.size()) += gaunt_factor*Alfa*
                Map<VectorXcd>(Helpers_e_e::TwoElecOverlap[Idx[0]][Idx[1]][lambda].data(),Helpers_e_e::TwoElecOverlap[Idx[0]][Idx[1]][lambda].size()).cwiseProduct(
                    Map<VectorXcd>(const_cast<vector<complex<double> >*>(&X)->data(),X.size()));
    }

    //    STOPDEBUG(axpyEE);
}

// Can't be called isZero, because C++ is brilliant.
static bool checkIfZero(int l1, int l2, int m1, int m2, int l1p, int l2p, int m1p, int m2p, int lambda_upper_limit, double Eps){
    if(m1p-m1!=m2-m2p) return true;
    int M = m1p-m1;
    for(int lambda=0;lambda<=lambda_upper_limit;lambda++){
        if(Gaunt::main.coeff_isZero(l1p,lambda,l1,m1p,M,m1) or Gaunt::main.coeff_isZero(l2,lambda,l2p,m2,M,m2p)) continue;

        double gf1,gf2;    // gaunt coefficients
        gf1 = Gaunt::main.coeff(l1p,lambda,l1,m1p,M,m1);
        if(abs(gf1) < Eps) continue;
        gf2 = Gaunt::main.coeff(l2,lambda,l2p,m2,M,m2p);

        double gaunt_factor=gf1*gf2*(4.0*math::pi/(2.0*double(lambda)+1.0));
        if(abs(gaunt_factor)>=Eps) return false;
    }
    return true;

}

bool OperatorFloorEE::isZero(double Eps) const
{
    if(Idx[0]!=Jdx[0] or Idx[1]!=Jdx[1]) return true;
    int unscaledElems = Helpers_e_e::TwoElecOverlap.size();
    if(Idx[0]>=unscaledElems or Idx[1]>=unscaledElems or Jdx[0]>=unscaledElems or Jdx[1]>=unscaledElems) return true; // 0 in complex scaled region
    
    return checkIfZero(Idx[2], Idx[3], Idx[4], Idx[5], Jdx[2], Jdx[3], Jdx[4], Jdx[5], lambda_upper_limit, epsGaunt);
}

void OperatorFloorEE::construct_index(const Index *I, std::vector<int> &Idx)
{// order n1,n2,m1,m2,l1,l2
    Index::build=true;
    if(not I->hasFloor()) ABORT("Not floor ?? ");

    Idx.resize(6);
    for(const Index* s=I; ; s=s->parent()){
        if(s->parent()==0) break;
        if(s->parent()->axisName()=="Rn1") Idx[0]=s->nSibling();
        else if(s->parent()->axisName()=="Rn2") Idx[1]=s->nSibling();
        else if(s->parent()->axisName()=="Eta1") Idx[2]=s->physical();
        else if(s->parent()->axisName()=="Eta2") Idx[3]=s->physical();
        else if(s->parent()->axisName()=="Phi1") Idx[4]=s->physical();
        else if(s->parent()->axisName()=="Phi2") Idx[5]=s->physical();
        else ABORT("Unknown axis name: "+s->parent()->axisName());
    }
    Index::build=false;

}

void OperatorFloorEE::constrain(UseMatrix& Mult, const Index* IIndex, const Index* JIndex){

    // ordered: m1, m2, l1
    std::vector<int> iIdx(3, INT_MIN);
    std::vector<int> jIdx(3, INT_MIN);

    auto construct = [](const Index* Idx, std::vector<int>& Result){
        for(const Index* s=Idx; s->parent()!=0; s=s->parent()){
            if(s->parent()->axisName()=="Eta1") Result[2]=s->physical();
            else if(s->parent()->axisName()=="Phi1") Result[0]=s->physical();
            else if(s->parent()->axisName()=="Phi2") Result[1]=s->physical();
        }
    };

    construct(IIndex, iIdx);
    construct(JIndex, jIdx);

    for(int i=0; i<2; i++){
        if(iIdx[i] == INT_MIN or jIdx[i] == INT_MIN) return; // Not all magnetic momenta set
    }

    if(iIdx[2] == INT_MIN or jIdx[2] == INT_MIN) return; // Not both l_1 set

    /*
     * Constrain based on Gaunt values - should not produce bloating
     */

    // Ensure Gaunts are initialised
    Helpers_e_e::initialize(IIndex->root(), lambda_upper_limit);

    for(size_t i=0; i<Mult.rows(); i++){
        for(size_t j=0; j<Mult.cols(); j++){
            if(Mult(i, j)==0.) continue;

            int l2 = std::abs(iIdx[1]) + i;
            int l2p = std::abs(jIdx[1]) + j;

            if(checkIfZero(iIdx[2], l2, iIdx[0], iIdx[1], jIdx[2], l2p, jIdx[0], jIdx[1], lambda_upper_limit, epsGaunt))
                Mult(i, j)=0.;
        }
    }
}

void OperatorFloorEE::pack(std::vector<int> &Info, std::vector<std::complex<double> > &Buf) const{
    Buf.clear();
    for(unsigned int k=0;k<6;k++)
        Buf.push_back(complex<double>(Idx[k],Jdx[k]));
    packBasic(Info,Buf);
}
OperatorFloorEE::OperatorFloorEE(const std::vector<int> &Info, const std::vector<std::complex<double> > &Buf)
    :OperatorFloor("FloorEE"){

    if(Helpers_e_e::TwoElecOverlap.size()==0) ABORT("Two Electron Overlap Matrices not initialised");

    unpackBasic(Info, Buf);
    Idx.clear();
    Jdx.clear();
    for(int i=0; i<6; i++){
        Idx.push_back(Buf[i].real());
        Jdx.push_back(Buf[i].imag());
    }
}


void OperatorFloorEE::Helpers_e_e::analyseDVR(const BasisAbstract* B, UseMatrix& pts, UseMatrix& wgs, vector<bool>& wgs_inc){

    UseMatrix pts_mat,wgs_mat;
    const BasisAbstract * B1 = B;
    BasisMat1D * b1 = new BasisMat1D("<1>", B1, B1);
    BasisMat1D * b2 = new BasisMat1D("<Q>", B1, B1);
    pts_mat = b2->useMat();
    wgs_mat = b1->useMat();
    wgs_inc.clear();  // Does the function include weight
    for(size_t i=0;i<wgs_mat.rows();i++){
        if(abs(wgs_mat(i,i).real()-1.0)<1e-14) wgs_inc.push_back(true);
        else wgs_inc.push_back(false);
    }

    // Cannot get lobatto points from quad rule??? So using own routine
    VectorXd pts_temp,wgs_temp;
    lobatto_quadratureNew(B,pts_temp,wgs_temp);
    if(abs(B->integrable()->lowBound())<1e-12){
        VectorXd temp = pts_temp.segment(1,pts_temp.size()-1);
        pts_temp = temp;
        temp = wgs_temp.segment(1,wgs_temp.size()-1);
        wgs_temp = temp;
    }

    // For some reason the points are not in ascending order, so rearrange according to pts_mat
    wgs = UseMatrix::Zero(wgs_temp.size(),1);
    pts = UseMatrix::Zero(wgs_temp.size(),1);

    for(size_t i=0;i<pts_mat.rows();i++){
        int index=-1;
        for(size_t j=0;j<pts_mat.rows();j++){
            if(abs(pts_temp(j)-pts_mat(i,i).real())<1e-9) {index=j; break;}
            if((i==0 or i==pts_mat.rows()-1) and abs(pts_temp(j)*wgs_temp(i)-pts_mat(i,i).real())<1e-9) {index=j; break;}
        }
        if(index==-1) ABORT("Could not find the point. Why??? "+tools::str(pts_mat(i,i).real()));
        wgs(i) = wgs_temp[index];
        pts(i) = pts_temp[index];
    }
}

OperatorFloorEE::~OperatorFloorEE(){
}
