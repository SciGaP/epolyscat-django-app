// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "schroedingerMachine.h"

#include "abort.h"
#include "mpiWrapper.h"
#include "printOutput.h"

// for testing schroma routines and classes
int main(int argc, char* argv[]) {

    MPIwrapper::Init(argc,argv);

    // read Psi's
    int nPsi,sizePsi;
    sm_SizesTrain("testData",&nPsi,&sizePsi);
    double data[nPsi*sizePsi];
    sm_DataTrain("testData",data);

    // get the training structure
    int nElemTrain(sm_NElemTrain());
    std::vector<int> sizeElemTrain(nElemTrain);
    std::vector<double> bElemTrain(nElemTrain+1);
    sm_StrucTrain(sizeElemTrain.data(),bElemTrain.data());

    // change the element boundaries
    std::vector<float> bElem(bElemTrain.size()-4);
    bElem.assign(bElem.size(),(bElemTrain.back()-bElemTrain.front())/(bElem.size()-1));
    bElem.front()=exp(bElemTrain.front());

    SchroedingerMachine m;
    double *datN=data;
    for(int n=0;n<nPsi;n++,datN+=sizePsi){
        if(true or n==0){
            Eigen::MatrixXcd psi(sizePsi/(nElemTrain)/2,nElemTrain);
            for(int j=0,ij=0;j<psi.cols();j++)
                for(int i=0;i<psi.rows();i++,ij+=2)
                    psi(i,j)=std::complex<double>(datN[ij],datN[ij+1]);
            double loss=m.LossL2(Eigen::Map<Eigen::VectorXf>(bElem.data(),bElem.size()),psi)(0,0);
            std::cout<<"new loss "<<n<<": "<<loss<<std::endl;
            if(abs(loss-1.)>1.e-5)ABORT("here");
        }
    }

    PrintOutput::message("done");
    MPIwrapper::Finalize();


}
