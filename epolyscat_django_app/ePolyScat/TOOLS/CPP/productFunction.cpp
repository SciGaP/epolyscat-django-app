// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "productFunction.h"

#include "algebra.h"
#include "readInput.h"
#include "printOutput.h"

ProductFunction::ProductFunction(std::string Coors, std::vector<std::string> Algs)
    :_coors(Coors)
{

    _algs.resize(Algs.size());
    for(int k=0;k<Algs.size();k++){
        std::vector<std::string> fact=tools::splitString(Algs[k],'.',"<({[",">)}]");
        if(fact.size()!=std::count(_coors.begin(),_coors.end(),'.')+1)
            ABORT("algebra count in"+Algs[k]+"does not match coordinates"+_coors
                  +"\nseparate factors by '.', use brackets (alg1).(alg2).alg3 etc if needed");
        for(int l=0;l<fact.size();l++){
            _algs[k].push_back(std::shared_ptr<Algebra>(new Algebra(fact[l])));
            if(not _algs[0].back()->isAlgebra())
                ABORT("malformed algebra in"+fact[k]+"\n"+Algebra::failures);
        }
    }
}

std::vector<std::complex<double> > ProductFunction::operator()(const std::vector<double> X) const{
    validArguments(X,"(ProductFunction)");
    std::vector<std::complex<double> > res(_algs.size(),1.);
    for(int k=0;k<_algs.size();k++){
        for(int l=0;l<X.size();l++)res[k]*=_algs[k][l]->val(X[l]);
    }
    return res;
}

void ProductFunction::read(ReadInput &Inp){
    if(not Inp.found("ProductFunction"))return;

    std::string name("BLANK"),coors("BLANK"),namePrev,coorsPrev;;
    std::vector<std::string> algs;
    int line(1);
    do {
        namePrev=name;
        Inp.read("ProductFunction","name",name,namePrev,"unique name for this function",line);


        coorsPrev=coors;
        algs.push_back("");
        Inp.read("ProductFunction","coordinates",coors,coorsPrev,"coordinates e.g. Phi.Eta.Rn or X.Y.Z",line);
        Inp.read("ProductFunction","factorAlgebras",algs.back(),"BLANK","product of single argument algebras, e.g. 1.Q.exp[-1](Q)",line);
        if(name=="BLANK" or coors=="BLANK" or algs.back()=="BLANK")
            ABORT("ProductFunction: must specify name,coordinates,factorAlgebras, got"+name+","+coors+","+algs.back());

        if((namePrev!="BLANK" and name!=namePrev) or Inp.endCategory("ProductFunction",line+1)){
            // save if name changes or end of input
            add(name,std::shared_ptr<VectorValuedFunction>(new ProductFunction(coors,algs)));
            algs.clear();
        }

    }while (not Inp.endCategory("ProductFunction",++line));

}
