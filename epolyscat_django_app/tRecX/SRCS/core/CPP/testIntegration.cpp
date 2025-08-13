// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "testIntegration.h"

TestIntegration::TestIntegration()
{

}

Eigen::MatrixXd TestIntegration::GetBasisCoeff(const Index* idx){
    const BasisDVR* b=dynamic_cast<const BasisDVR*>(idx->basis());
    //if(b->size()>3)DEVABORT("Test is only up to order 3, current:"+b->size());

    std::vector<double> dvrPoints, dvrWeights;
    b->dvrRule(dvrPoints, dvrWeights);
    Eigen::MatrixXd M(b->order(),b->order());
    Eigen::MatrixXd C(b->order(),b->size());
    Eigen::MatrixXd F(b->order(),b->size());
    Eigen::VectorXd X(b->order());
    for(int j=0;j<b->order();j++){
        X(j)=(b->upBound()-b->lowBound())*j/b->order() + b->lowBound()+(b->upBound()-b->lowBound())/(b->order()+1);
    }
    for(int i=0;i<b->order();i++){
        M(i,0)=1.;
    }
    for(int i=0;i<b->order();i++){
        for(int j=1;j<b->order();j++){
            M(i,j)=M(i,j-1)*X(i);
        }
    }
    for(int i=0;i<b->order();i++){
        for(int j=0;j<b->size();j++){
            F(i,j)=b->val(X(i))[j].real();
        }
    }
    C=M.inverse()*F;
    //std::cout<<"C for "<<nCousin(idx)<<" with a="<<idx->basis()->lowBound()<<" and b="<<idx->basis()->upBound()<<std::endl;
    //std::cout<<C<<std::endl;


    //-----Test-----//
    for(int i=0;i<b->size();i++){
        double sum=C(0,i);
            for(int k=1;k<b->order();k++){
                sum+=C(k,i)*pow(dvrPoints[0],k);
            }
            if(sum-b->val(dvrPoints[0])[i].real()>1.e-3){
                DEVABORT(Sstr+"Wrong coefficients"+sum+"val"+b->val(dvrPoints[0])[i].real());
            }
        }
    return C;
}

double TestIntegration::IntegrateElementEE(const Index *idx, int i, const Index *jdx, int j, std::string type){
    Eigen::VectorXd C1=GetBasisCoeff(idx).col(i);
    Eigen::VectorXd C2=GetBasisCoeff(jdx).col(j);
    const BasisIntegrable* ibas=idx->basis()->integrable();
    const BasisIntegrable* jbas=jdx->basis()->integrable();
    //std::cout<<"INTEGRATING FOR i = "<<i<<" and j = "<<j<<" for idx->Cousin, jdx->Cousin "<<nCousin(idx)<<" "<<nCousin(jdx)<<std::endl;
    if(type=="Left"){
        return IntegratePolynomial(C1,ibas->upBound(),ibas->lowBound(),1)*IntegratePolynomial(C2,jbas->upBound(),jbas->lowBound(),2);
    }
    else if(type=="Right"){
        return IntegratePolynomial(C1,ibas->upBound(),ibas->lowBound(),2)*IntegratePolynomial(C2,jbas->upBound(),jbas->lowBound(),1);
    }
    else if(type=="Diagonal"){
        if(idx!=jdx)DEVABORT("Element not on the diagonal");
        //return IntegrateLowerTriangle(C1,C2,idx->basis()->upBound(),idx->basis()->lowBound())+IntegrateUpperTriangle(C1,C2,idx->basis()->upBound(),idx->basis()->lowBound());
        return IntegrateLowerTriangle(C1,C2,ibas->upBound(),ibas->lowBound())+IntegrateLowerTriangle(C2,C1,ibas->upBound(),ibas->lowBound());
    }
    else {
        DEVABORT("Wrong type:"+type);
        return 0;
    }
}

double TestIntegration::IntegrateLowerTriangle(Eigen::VectorXd C1, Eigen::VectorXd C2,double UpBound,double LowBound){
    //if(C2.size()!=4)DEVABORT("generalize IntegrateLowerTriangle");

    //double A=IntegratePolynomial(C2,0,LowBound,2);
    double sum1=0.;
    double sum2=0.;
    for(int i=0;i<C1.size();i++){
        for(int j=0;j<C2.size();j++){
            sum1+=C1(i)*C2(j)*IntegrateMonomial(i+j+4,UpBound,LowBound)/(j+3);
        }
    }
    for(int i=0;i<C1.size();i++){
        for(int j=0;j<C2.size();j++){
            sum2+=C1(i)*C2(j)*IntegrateMonomial(i+1,UpBound,LowBound)*pow(LowBound,j+3)/(j+3);
        }
    }
    return sum1-sum2;;
}

double TestIntegration::IntegrateUpperTriangle(Eigen::VectorXd C1, Eigen::VectorXd C2,double UpBound,double LowBound){
    //if(C2.size()!=4)DEVABORT("generalize IntegrateLowerTriangle");

    //double A=IntegratePolynomial(C2,0,LowBound,2);
    double sum1=0.;
    double sum2=0.;
    for(int i=0;i<C1.size();i++){
        for(int j=0;j<C2.size();j++){
            sum1+=C1(i)*C2(j)*IntegrateMonomial(i+2,UpBound,LowBound)*pow(UpBound,j+2)/(j+2);
        }
    }
    for(int i=0;i<C1.size();i++){
        for(int j=0;j<C2.size();j++){
            sum2+=C1(i)*C2(j)*IntegrateMonomial(i+j+4,UpBound,LowBound)/(j+2);
        }
    }
    return sum1-sum2;;
}

double TestIntegration::IntegratePolynomial(Eigen::VectorXd C, double UpBound, double LowBound, int powerofX){
    //std::cout<<"Printing for LowBound="<<LowBound;
    double sum=0;
    for(int i=0;i<C.size();i++){
        sum+=C(i)*IntegrateMonomial(i+powerofX, UpBound, LowBound);
    }
    return sum;
}

double TestIntegration::IntegrateMonomial(int power,double UpBound, double LowBound){
    return (pow(UpBound,power+1)-pow(LowBound,power+1))/(power+1);
}
