// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "gaunt.h"
//#include "basisAbstract.h"
//#include "operatorData.h"
#include "basisAssocLeg.h"

using namespace std;

Gaunt Gaunt::main = Gaunt();

Gaunt::Gaunt():lmax(-1),order(0){}

void Gaunt::initialize(int lmax)
{
  this->lmax = lmax;
  order = max(3,lmax+1);

  std::vector<double> _quad,_weig;
  Passoc.clear();
  for(int m=0;m<=lmax;m++){

      BasisAssocLeg assoLegM(std::string("AssociatedLegendre: ")+tools::str(m)+","+tools::str(lmax+1));

      if(m==0){
          assoLegM.quadRule(order,_quad,_weig);
          quadX = Eigen::VectorXd::Zero(_quad.size());
          quadW = Eigen::VectorXd::Zero(_quad.size());
          for(int i=0;i<quadX.rows();i++){
              quadX(i) = _quad[i];
              quadW(i) = _weig[i];
            }
        }

      UseMatrix val = assoLegM.val(_quad);
      Eigen::MatrixXd valEig = Eigen::MatrixXd::Zero(val.rows(),val.cols());
      for(size_t i=0;i<val.rows();i++)
        for(size_t j=0;j<val.cols();j++){
            complex<double > a = val(i,j).complex();
            if(abs(imag(a))>1e-14) ABORT("Complex part in values ?? "+tools::str(imag(a)));
            valEig(i,j) = real(a);
          }
      Passoc.push_back(valEig*sqrt(1.0/(2.0*math::pi)));      // normalization for m component;
    }

  // Run the test
  test();
}

bool Gaunt::coeff_isZero(int lc, int l1, int l2, int mc, int m1, int m2)
{
  return mc!=m1+m2 or abs(l1-l2)>lc or lc>l1+l2 or abs(mc)>lc or abs(m1)>l1 or abs(m2)>l2;
}

double Gaunt::coeff(int lc, int l1, int l2, int mc, int m1, int m2)
{
  
  if(std::max(lc, std::max(l1, l2)) > lmax) initialize(std::max(lc, std::max(l1, l2)));

  // returns the gaunt coefficient
  if (mc!=m1+m2)  return 0.0;
  if (abs(l1-l2)>lc or lc>l1+l2) return 0.0;
  if (abs(mc)>lc or abs(m1)>l1 or abs(m2)>l2) return 0.0;

  int f=1;
  if (mc<0 and (-mc)%2==1) f  = -1;
  if (m1<0 and (-m1)%2==1) f *= -1;
  if (m2<0 and (-m2)%2==1) f *= -1;

  return (Passoc[abs(mc)].col(lc-abs(mc)).array()
          *Passoc[abs(m1)].col(l1-abs(m1)).array()
          *Passoc[abs(m2)].col(l2-abs(m2)).array()
          *quadW.array()
          ).sum()*2.0*math::pi*(double)f;

  //  return (Passoc[abs(mc)][lc-abs(mc)].array()*Passoc[abs(m1)][l1-abs(m1)].array()*Passoc[abs(m2)][l2-abs(m1)].array()*quadw.cast<long double>().array()).sum()*2.0*myPi*(double)f;
}

void Gaunt::test(){

//    if(MPIwrapper::isMaster()==1) cout << "gaunt test ... " << flush;
    if(this->Passoc.size()==0) return;
    double eps = 1e-12;
    int errorcount = 0;
    int shortentest=0;
    for( int m=lmax;m>=0;m--){
        if(shortentest++ > 10) break;
        for(unsigned int l=m; l<=(unsigned int)lmax; l++){
            // if(l>85) break;
            double test = coeff(l,l,0,m,m,0);
            if (abs(abs(test)-sqrt(0.25/math::pi))>eps) {
                cout<<"failed test1 in gaunt.h. m="<<m<<", l="<<l<<" off by " << abs(abs(test)-sqrt(0.25/math::pi))<<" "<<abs(test)<<" "<<sqrt(0.25/math::pi)<<endl;// exit(0);
                errorcount++;
            }
            if (abs(abs(coeff(l,l,0,-m,-m,0))-sqrt(0.25/math::pi))>eps)            {
                cout<<"failed test2 in gaunt.h. m="<<m<<", l="<<l<<" off by " << abs(abs(coeff(l,l,0,-m,-m,0))-sqrt(0.25/math::pi))<<endl; //exit(0);
                errorcount++;
            }
            if (abs(abs(coeff(0,l,l,0,m,-m))-sqrt(0.25/math::pi))>eps)            {
                cout<<"failed test3 in gaunt.h. m="<<m<<", l="<<l<<" off by " << abs(abs(coeff(0,l,l,0,m,-m))-sqrt(0.25/math::pi))<<endl;//exit(0);
                errorcount++;
            }
            if (abs(abs(coeff(0,l,l,0,-m, m))-sqrt(0.25/math::pi))>eps)            {
                cout<<"failed test4 in gaunt.h. m="<<m<<", l="<<l<<" off by " << abs(abs(coeff(0,l,l,0,-m, m))-sqrt(0.25/math::pi))<<endl;//exit(0);
                errorcount++;
            }
        }
    }
    if(errorcount>0){
        cout << endl << "GauntCoeffTable::test: long double should be giving correct results until lmax/mmax > 200." << endl;
        cout << "if shit happened at l/m around 90, this is possibly because the compiler/computer/whatever has a different (smaller: like double) definition of long double than ours ... " << endl;
        ABORT("GauntCoeffTable::test fatal");
    }
    else{
    //    if(MPIwrapper::isMaster()==1) cout << "finished" << endl;
    }

}

double Gaunt::twoYintegrals(int l1, int l2, int m1, int m2, string operat)
{

  ABORT("Need to implement");

//  BasisSet::Def AssLeg(lmax+1,-1.0,2.0,"assocLegendre",true,true,true,Coordinate::fromString("Eta"),vector<int>(0),
//                       ComplexScaling(),false,vector<double>(1,m));
//  BasisSet* YLM = BasisSet::get(AssLeg);

//  UseMatrix mat;
//  OperatorData::get(operat,vector<BasisSet*>(1,YLM),vector<BasisSet*>(1,YLM),vector<UseMatrix>(1,mat));
//  return mat().real();

//  //these are used to compute the overlaps appearing in dipol field hamiltonian
//      Matrix<long double,Dynamic, 1> temp,temp2;
//      double result = 0.0;
//      if(operat == "q") result= ( ValueOfAssLeg(m1,l1,temp).array()*normForSphericHarmonic(m1,l1)*ValueOfAssLeg(m2,l2,temp2).array()*normForSphericHarmonic(m2,l2)*quadw.cast<long double>().array()*quadx.cast<long double>().array()).sum();
//      else if(operat == "q^2") result= ( ValueOfAssLeg(m1,l1,temp).array()*normForSphericHarmonic(m1,l1)*ValueOfAssLeg(m2,l2,temp2).array()*normForSphericHarmonic(m2,l2)*quadw.cast<long double>().array()*quadx.cast<long double>().array()*quadx.cast<long double>().array()).sum();
//      else if(operat == "(1-q^2)")  result= ( ValueOfAssLeg(m1,l1,temp).array()*normForSphericHarmonic(m1,l1)*ValueOfAssLeg(m2,l2,temp2).array()*normForSphericHarmonic(m2,l2)*quadw.cast<long double>().array()*
//                                              (1.0-quadx.cast<long double>().array()*quadx.cast<long double>().array())).sum();
//      else if(operat == "|") result= ( ValueOfAssLeg(m1,l1,temp).array()*normForSphericHarmonic(m1,l1)*ValueOfAssLeg(m2,l2,temp2).array()*normForSphericHarmonic(m2,l2)*quadw.cast<long double>().array()).sum();
//      else if(operat == "(1-q^2)d") result= ( ValueOfAssLeg(m1,l1,temp).array()*normForSphericHarmonic(m1,l1)*DerivOfAssLeg(m2,l2,temp2).array()*
//                                              normForSphericHarmonic(m2,l2)*quadw.cast<long double>().array()*(1.0-quadx.cast<long double>().array()*quadx.cast<long double>().array())).sum();
//      else if(operat == "sqrt(1-q^2)") result= (ValueOfAssLeg(m1,l1,temp).array()*normForSphericHarmonic(m1,l1)*ValueOfAssLeg(m2,l2,temp2).array()*
//                                                normForSphericHarmonic(m2,l2)*quadw.cast<long double>().array()*sqrt(1.0-quadx.cast<long double>().array()*quadx.cast<long double>().array())).sum();
//      else if(operat == "sqrt(1-q^2)qd") result= (ValueOfAssLeg(m1,l1,temp).array()*normForSphericHarmonic(m1,l1)*DerivOfAssLeg(m2,l2,temp2).array()*
//                                                  normForSphericHarmonic(m2,l2)*quadw.cast<long double>().array()*quadx.cast<long double>().array()*sqrt(1.0-quadx.cast<long double>().array()*quadx.cast<long double>().array())).sum();
//      else if(operat == "sqrt(1-q^2)q") result= (ValueOfAssLeg(m1,l1,temp).array()*normForSphericHarmonic(m1,l1)*ValueOfAssLeg(m2,l2,temp2).array()*
//                                                  normForSphericHarmonic(m2,l2)*quadw.cast<long double>().array()*quadx.cast<long double>().array()*sqrt(1.0-quadx.cast<long double>().array()*quadx.cast<long double>().array())).sum();
//      else if(operat == "1/sqrt(1-q^2)") result= (ValueOfAssLeg(m1,l1,temp).array()*normForSphericHarmonic(m1,l1)*ValueOfAssLeg(m2,l2,temp2).array()*
//                                                  normForSphericHarmonic(m2,l2)*quadw.cast<long double>().array()/sqrt(1.0-quadx.cast<long double>().array()*quadx.cast<long double>().array())).sum();
//      else if(operat == "q^2/sqrt(1-q^2)") result=(ValueOfAssLeg(m1,l1,temp).array()*normForSphericHarmonic(m1,l1)*ValueOfAssLeg(m2,l2,temp2).array()*
//                                                    normForSphericHarmonic(m2,l2)*quadw.cast<long double>().array()*quadx.cast<long double>().array()*quadx.cast<long double>().array()/sqrt(1.0-quadx.cast<long double>().array()*quadx.cast<long double>().array())).sum();
//      else{
//          cout << "operator in twoYintegrals in gaunt.h unknown " <<operat<< endl; throw;
//      }
//      if(result != result){
//          cout << normForSphericHarmonic(m1,l1) << " " << normForSphericHarmonic(m2,l2) << endl;
//          cout << "nan in GauntCoeffTable::twoYintegrals " << l1 << " " << l2 << " " << m1 << " " << m2 << " " << operat << endl;
//          throw;
//      }
//      return result;

}
