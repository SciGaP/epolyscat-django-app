// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "tools.h"

//#include "operator.h"
#include "rungeKutta4.h"
#include "simple.h"
#include "readInput.h"
#include "printOutput.h"
#include "basicDisc.h"
#include "parameters.h"
#include "timer.h"
#include "pulse.h"
#include "discretizationDerived.h"
#include "plot.h"
#include "coefficients.h"

using namespace std;
using namespace Eigen;
using namespace tools;

static unsigned int currentFile=0;       ///< will be set to the class static fileCount

TIMER(all,)
TIMER(operator,)
TIMER(eigen,)
TIMER(matrix,)
TIMER(io,)

int main(int argc, char* argv[]) {

    START(all);

    ReadInput::openMain("",argc,argv);

    Discretization * D = new BasicDisc(ReadInput::main);
    D->print();

    string hamDef,ovrDef,srcDef;
    ReadInput::main.read("Operator","hamiltonian",hamDef,"<1><J><1/2.dJd-Jr>+<1><1/2.dJ1q2d><rJr>","hamiltonian of Maxwell's equations");
    ReadInput::main.read("Operator","overlap",ovrDef,"<1><J><J>","overlap");
    ReadInput::main.read("Operator","source",srcDef,"|delta[1]>|cos[0]>|cos[1]>","inhomogenuity");

    PrintOutput::title("OPERATORS");
    PrintOutput::paragraph();
    PrintOutput::lineItem("Hamiltonian",hamDef);
    PrintOutput::newLine();
    PrintOutput::lineItem("Overlap",ovrDef);
    PrintOutput::newLine();
    PrintOutput::lineItem("Source",srcDef);


    Pulse::read(ReadInput::main,true,"Pulse");


    double accuracy,cutE;
    double tBeg,tEnd,tPrint,tStore;
    ReadInput::main.read("TimePropagation","begin",tBeg,"0","begin time of propagation");
    ReadInput::main.read("TimePropagation","end",  tEnd  ,"0","end time of propagation");
    ReadInput::main.read("TimePropagation","print",tPrint,tools::str((tEnd-tBeg)/50,4),"printout intervals");
    ReadInput::main.read("TimePropagation","store",tStore,tools::str(tPrint),"store intervals");
    ReadInput::main.read("TimePropagation","accuracy",accuracy,"1.e-8","accuracy control");
    ReadInput::main.read("TimePropagation","cutEnergy",cutE,str(DBL_MAX),"remove energies above");

    PrintOutput::paragraph();

    int neig;
    ReadInput::main.read("Print","eigenvalues",neig,"0","how many eigenvalues to print");

    PrintOutput::title("OUTPUT");
    PrintOutput::lineItem("Directory",ReadInput::main.output());
    PrintOutput::paragraph();


    START(operator);
    Operator SR("Overlap",ovrDef,D,D);
    Operator HR("Hamiltonian",hamDef,D,D);
    Operator HI("Source",srcDef,D,D);
    Operator HS("Ham+Source",hamDef+"+"+srcDef,D,D);
    STOP(operator);

    Plot plot(D,ReadInput::main); // set up plot, get the definitions from file

    bool eigenOnly;
    unsigned int maxEigen;
    ReadInput::main.read("Control","eigenOnly",eigenOnly,"false","flag -eigenOnly stops after computing eigenvalues",1,"eigenOnly");
    ReadInput::main.read("Control","maxEigen",maxEigen,"2000","flag -suppresses diagonalization for sytem size > maxEigen",1,"maxEigen");

    ReadInput::main.finish(); // operator parameters may require further input
    //============================================================================================

    vector<complex<double> > eval;
    vector<Coefficients*> evec;

    // plot source
    Coefficients c(D->idx());
    Wavefunction src(10.,&c);
    for(unsigned int i=0;i<=10;i++){
        src.time=Pulse::current.beginPrint()+(Pulse::current.endPrint()-Pulse::current.beginPrint())*i*0.1;
        src.setToZero();
        HI.axpy(src,src);
//        src.coefs->inverseOverlap();
        DEVABORT("re-write using inverse overlap");
        plot.plot(*src.coefs,ReadInput::main.output()+"src"+tools::str(i));
    }

    D->idx()->testInverseOverlap();

    // get the spectrum (if not too large)
    if(src.coefs->size()<=maxEigen){
        if(src.coefs->size()>2500)
            PrintOutput::warning("system size "+tools::str(src.coefs->size())
                                 +" - expect long eigenvalue calculation");
        START(eigen);
        UseMatrix::eigenMethod=UseMatrix::general;
        HR.eigen(SR,eval,evec);
        STOP(eigen);

        ofstream eig,eigPosI;
        eig.open((ReadInput::main.output()+"eig").c_str());
        unsigned int countPos=0;
        PrintOutput::message("Eigenvalues in "+ReadInput::main.output()+"eig");
        for (unsigned int n=0;n<eval.size();n++){
            eig<<setprecision(10)<<real(eval[n])<<", "<<imag(eval[n])<<endl;
            if(imag(eval[n])>1.e-12){
                if(countPos==0)eigPosI.open((ReadInput::main.output()+"eigPosI").c_str());
                countPos++;
                eigPosI<<setprecision(10)<<real(eval[n])<<", "<<imag(eval[n])<<endl;
            }
        }
        eig.close();
        if(countPos>0){
            PrintOutput::message("Eigenvalues pos. imag. ("+tools::str(countPos)+") in "+ReadInput::main.output()+"eigPosI");
            eigPosI.close();
        }
        if(eigenOnly)exit(0);
    }
    else PrintOutput::warning("system size "+tools::str(src.coefs->size())+" - eigenvalues not calculated");

    // Output: print expectation values and plot wave function
    vector<Operator*>printExp(1,&HR);
    printExp.push_back(&SR);
    RungeKutta4::Output rkOut(printExp,"",tPrint,ReadInput::main.output(),&plot);
    // Derivative: apply Hamiltonian
    RungeKutta4::Derivative rkDer(&HS);
    // propagator setup
    RungeKutta4 rk4(&rkDer,&rkOut,accuracy);

    // propagate
    Wavefunction wf(D,tBeg);
    wf.setToZero(); // make sure the wave function is empty
    rk4.propagate(&wf,tEnd);


    STOP(all);
    Timer::write(cout);
    PrintOutput::title("output directory "+ReadInput::main.output());

}
