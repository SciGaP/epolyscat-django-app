// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "finDiff.h"
#include "tools.h"
#include "useMatrix.h"
#include "printOutput.h"
#include "basisFunction.h"
#include "algebra.h"

using namespace std;

std::string FinDiff::FDgrid="exp";
std::string FinDiff::FDmatrix="exact";

FinDiff::FinDiff(double X0, double X1, unsigned int N, double R0,std::complex<double> Eta, const string Kind)
    :r0(R0),eta(Eta){
    if(N<2)ABORT("cannot use grid points N="+tools::str(N));
    for(unsigned int k=1;k<N-1;k++)grid.push_back(X0+k*(X1-X0)/double(N-1));
    weig.assign(grid.size(),1./double(grid.size()));
    construct(Kind);
}

FinDiff::FinDiff(UseMatrix Grid, UseMatrix Weight, double R0, complex<double> Eta, const string Kind):r0(R0),eta(Eta){
    for(unsigned int k=0;k<Grid.size();k++){
        grid.push_back(Grid(k).real());
        if(Weight.size()>0)weig.push_back(Weight(k).real());
        else               weig.push_back(1.);
        construct(Kind);
    }
}

void FinDiff::construct(const std::string Kind){

    xLast=0;
    scalLast=0;

    // scaling
    string sR0=tools::str(r0),sEta=tools::str(eta);
    string kind=Kind.substr(0,Kind.find('['));

    singularX.clear();
    string arg,der,der2,der3;

    argScal=0;
    derScal=0;
    scalD1=0; // only for demonstrating errors....
    scalD2=0;
    if(kind=="sudden"){
        string left ="chi["+tools::str(-DBL_MAX)+","+tools::str(-r0)+"](Q)";
        string right="chi["+tools::str(r0)+","+tools::str(DBL_MAX)+"](Q)";

        arg="Q+( "+sR0+"-Q+(Q-"+sR0+")*"+sEta+")*"+right+
                "+(-"+sR0+"-Q+(Q+"+sR0+")*"+sEta+")*"+left ;
        der="1.+("+sEta+"-1)*"+left+"+("+sEta+"-1)*"+right;
        singularX.push_back(-r0);
        singularX.push_back( r0);
    }

    else if (kind=="smooth"){
        double smoothRange=tools::string_to_double(tools::stringInBetween(Kind,"[","]"));
        string sR1=tools::str(r0+smoothRange);

        string trunc="("+sEta+"+(1-"+sEta+")*(trunc["+sR0+","+sR1+"]))";
        string left  =sEta+"*chi["+tools::str(-DBL_MAX)+",-"+sR1+"](Q)";
        string lTrans="chi[-"+sR1+",-"+sR0+"](Q)*"+trunc;
        string inner= "chi[-"+sR0+", "+sR0+"](Q)";
        string rTrans="chi[ "+sR0+", "+sR1+"](Q)*"+trunc;
        string right =sEta+"*chi["+sR1+","+tools::str(DBL_MAX)+"](Q)";

        der=left+"+"+lTrans+"+"+inner+"+"+rTrans+"+"+right;

        singularX.push_back(-r0-smoothRange);
        singularX.push_back(-r0);
        singularX.push_back( r0);
        singularX.push_back( r0+smoothRange);

        // lambda' and lambda''
        der2=     "chi[-"+sR1+",-"+sR0+"](Q)*(1-"+sEta+")*trunc["+sR0+","+sR1+",1]"
                +"+chi[ "+sR0+", "+sR1+"](Q)*(1-"+sEta+")*trunc["+sR0+","+sR1+",1]";
        der3=     "chi[-"+sR1+",-"+sR0+"](Q)*(1-"+sEta+")*trunc["+sR0+","+sR1+",2]"
                +"+chi[ "+sR0+", "+sR1+"](Q)*(1-"+sEta+")*trunc["+sR0+","+sR1+",2]";
    }

    else if (kind=="exp"){
        singularX.push_back(-abs(r0));
        singularX.push_back( abs(r0));
        string left ="chi["+tools::str(-DBL_MAX)+",-"+sR0+"](Q)";
        string right="chi["+sR0+","+tools::str(DBL_MAX)+"  ](Q)";
        string lExp="*exp(-"+tools::stringInBetween(Kind,"[","]")+"*(Q+"+sR0+"))";
        string rExp="*exp( "+tools::stringInBetween(Kind,"[","]")+"*(Q-"+sR0+"))";
        der="1.-"+left+"-"+right+"+"+sEta+"*("+left+lExp+"+"+right+rExp+")";
    }
    else ABORT("undefined scaling: "+Kind);

    // set up the functions
    vector<string> def(4);def[0]=arg;def[1]=der;def[2]=der2;def[3]=der3;
    vector<Algebra*>alg(4);
    for(unsigned int i=0;i<def.size();i++){
        if(def[i]!=""){
            alg[i]=new Algebra(def[i]);
            if(not alg[i]->isAlgebra())ABORT("malformed algebra: "+alg[i]->definition()+"\n"+Algebra::failures);
        }
    }
    argScal=alg[0];
    derScal=alg[1];
    scalD1 =alg[2];
    scalD2 =alg[3];
    if(derScal==0)ABORT("grid spacing lambda undefined");
}

// coordinate scaling transformation
std::complex<double> FinDiff::scalArg(double X){

    // numerically integrate
    vector<double> steps;
    if(xLast<X){
        steps=tools::insertInterval(xLast,X,singularX);
        for (unsigned int k=1;k<steps.size();k++)scalLast+=derScal->integral(steps[k-1],steps[k]);
    } else {
        steps=tools::insertInterval(X,xLast,singularX);
        for (unsigned int k=1;k<steps.size();k++)scalLast-=derScal->integral(steps[k-1],steps[k]);
    }

    // for debugging: check against explicit integral function
    if(argScal!=0 and abs(scalLast-argScal->val(X))>1.e-10)
        ABORT("error in numerical integration: "+tools::str(scalLast)+" "+tools::str(argScal->val(X)));

    xLast=X;
    return scalLast;
}
std::complex<double> FinDiff::scalDer(double X){
    complex<double>res=derScal->val(X);
    if(res.real()<0.)ABORT("scaling must be in upper half plane: "+derScal->definition());
    return res;
}

// recurrence coefficients for the Legendre polynomials
double a(int i) {return 0.;}
double b(int i) {return (2.*double(i)-1.)/double(i);}
double c(int i) {return   (-double(i)+1.)/double(i);}

void FinDiff::interFunc(unsigned int N, complex<double> Z, complex<double> A, complex<double> B, vector<vector<complex<double> > > &V){


    if(abs(A-B)<=abs(A+B)*1.e-12)ABORT("zero size interval: ["+tools::str(A)+","+tools::str(B)+"]");
    complex<double> q=(2.*Z-(A+B))/(B-A);
    if(abs(q)>1.+1.e-9)ABORT("argument outside interval: "+tools::str(Z)+" ["+tools::str(A)+", "+tools::str(B)+"]");

    //    // exponentially damped interpolation hypothesis
    //    double Alfa=0.;
    //    complex<double> expDamp;
    //    if(A.real()> r0){
    //        Alfa= alfa;
    //        expDamp=exp(-Alfa*(Z-A));
    //    }
    //    else if(B.real()<-r0){
    //        Alfa=-alfa;
    //        expDamp=exp( Alfa*(Z-B));
    //    }
    //    if(abs(expDamp.imag())>1.e-12)ABORT("complex exponent is not intended");

    // evaluate values, 1st and 2nd derivative of the (complex) legendred polynomials
    V.assign(3,vector<complex<double> >(0));
    if(N>0){
        V[0].push_back(1.);
        V[1].push_back(0.);
        V[2].push_back(0.);
    }
    if(N>1){
        V[0].push_back(a(1)+b(1)* q*V[0][0]);
        V[1].push_back(     b(1)*(q*V[1][0]+V[0][0]));
        V[2].push_back(0.);
    }
    if(N>2){
        for (int n=2; n<N;n++){
            V[0].push_back((a(n) +  b(n)*q)*V[0][n-1]                  + c(n)*V[0][n-2]);
            V[1].push_back((a(n) +  b(n)*q)*V[1][n-1] + b(n)*V[0][n-1] + c(n)*V[1][n-2]);
            V[2].push_back((a(n) +  b(n)*q)*V[2][n-1] + b(n)*V[1][n-1] + c(n)*V[2][n-2] + b(n)*V[1][n-1]);
        }
    }

    // inner derivatives
    for (unsigned int n=0;n<N;n++){
        V[1][n]*=2./(B-A);
        V[2][n]*=4./((B-A)*(B-A));
    }

    //    // exponential damping
    //    if(Alfa!=0.){
    //        for (unsigned int n=0;n<N;n++){
    //            V[0][n]=V[0][n]*expDamp;
    //            V[1][n]=V[1][n]*expDamp-Alfa*V[0][n];
    //            V[2][n]=V[2][n]*expDamp-Alfa*(2.*V[1][n]-Alfa*V[0][n]);
    //        }
    //    }
}

void FinDiff::matrix(string Oper, unsigned int Points, const UseMatrix &Grid, const UseMatrix &Weig,
                     double R0, complex<double> Eta, double Par, UseMatrix &Matrix, string Kind){

    Kind=Kind+"["+tools::str(Par)+"]";
    FinDiff fd(Grid,Weig,R0,Eta,Kind);
    UseMatrix dum;
    if(Oper=="d_1_d")                  fd.scheme(Points,dum,Matrix);
    else if(Oper=="d_1" or Oper=="1_d")fd.scheme(Points,Matrix,dum);
    else if(Oper=="1/sqrt(Q*Q+2)"){
        Matrix=UseMatrix::Zero(fd.grid.size(),fd.grid.size());
        for(unsigned int k=0;k<fd.grid.size();k++){
            complex<double>z=fd.scalArg(fd.grid[k]);
            Matrix(k,k)=1./sqrt(2.+z*z);
        }
    }
    else ABORT("operator not implemented: "+Oper);
    Matrix.expand();
}

/// approximate Oper to a given Order
void FinDiff::scheme(unsigned int Points, UseMatrix & Der, UseMatrix & Lap){

    cout<<"scheme for "+FDmatrix+" "+FDgrid<<endl;
    if(Points%2!=1)ABORT("must have odd number of points (even degree), is: "+tools::str(Points));

    Der=UseMatrix(grid.size(),grid.size(),Points/2,Points/2);
    Lap=UseMatrix(grid.size(),grid.size(),Points/2,Points/2);

    // extend lower and upper ends of grids
    vector<double> gg(grid);
    double x;
    for(unsigned int k=0;k<Points/2;k++){
        x=gg[0]-pow(gg[1]-gg[0],2)/(gg[2]-gg[1]);
        gg.insert(gg.begin(),1,x);
        x=gg.back()+pow(gg[gg.size()-1]-gg[gg.size()-2],2.)/(gg[gg.size()-2]-gg[gg.size()-3]);
        gg.push_back(x);
    }

    for(int k=0;k<grid.size();k++){

        // get local subgrid
        vector<double> lg;
        for(int i=k-int(Points/2);i<=k+int(Points/2);i++)
            lg.push_back(gg[Points/2+i]);

        // evaluate fitting function at scaled points
        vector<vector<complex<double> > >valder;
        UseMatrix val(Points,Points),der(Points,Points),lap(Points,Points);
        // reference interval
        complex<double> a=scalArg(lg[0]),b=scalArg(lg.back());

        for(unsigned int i=0;i<lg.size();i++){

            complex<double>sqrtLa=1.;
            if(FDmatrix=="exact"){
                // evaluate at transformed points
                interFunc(Points,scalArg(lg[i]),a,b,valder);
                sqrtLa=sqrt(scalDer(lg[i]));
            }
            else if (FDmatrix=="standard" or FDmatrix=="approx"){
                // evaluate at grid points
                interFunc(Points,lg[i],lg[0],lg.back(),valder);
            }
            else ABORT("undefined FDmatrix: "+FDmatrix);

            if(FDmatrix=="exact" or FDmatrix=="approx"){
                // for convenience, here we define val=[sqrt(la) Pi]^T, der=[sqrt(la) Pi']^T
                for(unsigned int l=0;l<Points;l++){
                    val(l,i)= sqrtLa*valder[0][l];
                    der(l,i)= sqrtLa*valder[1][l];
                    lap(l,i)=-sqrtLa*valder[2][l];
                }
            }
            else if (FDmatrix=="standard"){
                if(scalD1==0 or scalD2==0)ABORT("need to define derivatives for FDmatrix=standard with "+FDgrid);
                complex<double> la=scalDer(lg[i]),dla=scalD1->val(lg[i]),qla=scalD2->val(lg[i]);
                for(unsigned int l=0;l<Points;l++){
                    la=     scalDer(lg[l]);
                    dla=scalD1->val(lg[l]);
                    qla=scalD2->val(lg[l]);

                    val(l,i)=valder[0][l];
                    der(l,i)=valder[1][l]/la
                            -0.5*dla/pow(la,2);
                    lap(l,i)=-valder[2][l]/(la*la)
                            +valder[1][l]*   2.*    dla/pow(la,3)
                            -valder[0][l]*(1.25*pow(dla,2)/pow(la,4)
                                           -0.5*    qla/pow(la,3));
                }
            }
            else
                ABORT("not implemented: "+FDmatrix);
        }

        // get the local stencil
        val.solve(der);
        val.solve(lap);


        // scheme = central row of [val^{-1} der]^T = central col of val^{-1} der
        int ioff=max(0,k-int(Points/2))-k;
        int ilen=min(int(grid.size()),int(Points/2)+k+1)-max(0,k-int(Points/2));

        if(FDmatrix=="approx"){
            // compose
            complex<double> qlak,rlak,dlak;
            complex<double> qlai,rlai,dlai;
            int kk=Points/2;
            qlak=1./scalDer(lg[kk]);
            rlak=sqrt(qlak);
            dlak=scalD1->val(lg[kk])*pow(qlak,2);
            for(unsigned int i=0;i<ilen;i++){
                qlai=1./scalDer(lg[i]);
                rlai=sqrt(qlai);
                dlai=0.5*scalD1->val(lg[i])*pow(qlai,2);

                lap(kk,i)         =lap(kk,i)*qlak*qlai; // note: lap here carries the minus already
                lap(kk,i)         =lap(kk,i)-der(kk,i)*(-dlak*qlai+qlak*dlai);
                if(i==kk)lap(kk,i)=lap(kk,i)+pow(dlai,2);

                der(kk,i)=der(kk,i)*rlak*rlai;
            }
        }

        Der.block(k,k+ioff,1,ilen)=der.col(Points/2).transpose().block(0,Points/2+ioff,1,ilen);
        Lap.block(k,k+ioff,1,ilen)=lap.col(Points/2).transpose().block(0,Points/2+ioff,1,ilen);
    }

    Der.purge(1.e-12); // remove near-zeros
    Lap.purge(1.e-12); // remove near-zeros
}


void FinDiff::test(){

    // compute the 1d-hydrogen complex scaled spectrum

    unsigned int points=3;
    double rmax=15;
    double r0=1;
    double dx=0.25;
    complex<double> eta=exp(complex<double>(0.,0.));

    FDgrid="smooth";
    FDmatrix="approx";

    unsigned int n=2*int(rmax/dx)+1;
    //    FinDiff fd(-rmax,rmax,n,r0,eta,"sudden");
    //      FinDiff fd(-rmax,rmax,n,r0,eta,"exp[0.]");
    FinDiff fd(-rmax,rmax,n,r0,eta,"smooth[10.]");

    cout<<"functions:\n"<<fd.derScal->definition()<<"\n"<<endl;

    // verify scaling function
    ofstream plot;
    plot.open("findiffScaling");
    //    if(fd.argScal!=0)plot<<"# argScal "+fd.argScal->definition<<endl;
    plot<<"# derScal "+fd.derScal->definition()<<endl;
    for (double x=-2.*fd.r0;x<=2.*fd.r0;x+=2.*fd.r0*0.01)
    {
        plot<<x<<", "<<fd.scalArg(x).real()<<", "<<fd.scalArg(x).imag()
           <<    ", "<<fd.scalDer(x).real()<<", "<<fd.scalDer(x).imag()<<endl;
    }
    plot<<flush;
    cout<<"*** scaling functions on \"findiffScaling\""<<endl;

    UseMatrix der,lap;
    fd.scheme(points,der,lap);
    if(der.rows()<31){
        der.print("1st derivative");
        lap.print("2nd derivative");
    }
    if(FDmatrix=="approx"){
        if(not lap.isSymmetric(1.e-12))PrintOutput::warning("approx lap is not symmetric ???");
        PrintOutput::message("approx lap is  symmetric");
    }

    UseMatrix ham(lap);
    ham*=0.5;
    for(unsigned int i=0;i<fd.grid.size();i++){
//        complex<double> z=fd.scalArg(fd.grid[i]);
        //ham(i,i)+=-1./sqrt(2.+z*z);
    }
    UseMatrix eval(ham.rows(),1),evalT(ham.rows(),1);
    UseMatrix rvec,lvecT;
    UseMatrix ovr=UseMatrix::Identity(ham.rows(),ham.cols());
    ham.expand();
    //    ham.show();
    cout<<"*** computing eigenvalues, size: "<<ham.rows()<<endl;
    ham.eigen(eval,rvec,ovr);
    ham.transpose().eigen(evalT,lvecT,ovr);

    eval.block(0,0,10,1).transpose().print("eigenvalue",5);
    ofstream eig;
    eig.open(("eigFinDiff_"+FDmatrix).c_str());
    int upper=0,lower=0;
    cout<<"\neigenvalues in upper half-plane: "<<endl;
    for(unsigned int i=0;i<eval.size();i++){
        if(eval(i).imag()>1.e-12){
            upper++;
            cout<<eval(i).complex()<<endl;
        }
        if(eval(i).imag()<-1.e-12)lower++;
        eig<<eval(i).real()<<", "<<eval(i).imag()<<endl;
    }
    cout<<"count: "<<upper<<"/"<<lower<<endl;

    FDgrid="smooth";
    FDmatrix="standard";
    FinDiff fd0(-rmax,rmax,n,r0,eta,"smooth[10.]");
    UseMatrix der0,lap0;
    fd0.scheme(points,der0,lap0);


    (der0-der).show("1st derivative");
    (lap0-lap).show("2nd derivative");
}
