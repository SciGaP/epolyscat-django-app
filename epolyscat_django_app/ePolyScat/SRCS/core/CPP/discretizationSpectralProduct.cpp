// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "discretizationSpectralProduct.h"

#include <string>
#include <vector>
#include "useMatrix.h"
#include "index.h"
#include "indexQuot.h"
#include "operatorDefinition.h"
#include "operatorFloor.h"
#include "operatorTree.h"
#include "operatorDiagonal.h"
#include "basicDisc.h"
#include "printOutput.h"
#include "purgeIndex.h"
#include "inverse.h"

#include "indexProd.h"
#include "operatorAbstractProduct.h"
#include "eigenSolver.h"
#include "log.h"
#include "indexConstraint.h"

#include "indexExtract.h"
#include "indexPermute.h"
#include "operatorExpandIndex.h"
#include "tree.h"
#include "eigenTools.h"

#include "parallelOperator.h"

using namespace std;

bool matchZeroStructure(const Eigen::MatrixXcd & A,const Eigen::MatrixXcd & B, double Eps=1.e-12){
    if(A.rows()!=B.rows() or A.cols()!=B.cols())return false;
    for(int k=0;k<A.size();k++)
        if((std::abs(A.data()[k])>Eps)!=(std::abs(B.data()[k])>Eps))return false;
    return true;
}


DiscretizationSpectralProduct::DiscretizationSpectralProduct(const Discretization *D, const string OpSeparable, double Emin, double Emax,
                                                             bool excludeEnergyRange)
    :DiscretizationSpectral(D,"["+tools::str(Emin,3,DBL_MAX/2.)+","+tools::str(Emax,3,DBL_MAX/2.)+"]"),
      proj(nullptr /* does not allow 0 */)
{

    name="spectralproduct_"+_selectionCriterion;

    std::vector<std::string>sep,term;
    if(std::count(OpSeparable.begin(), OpSeparable.end(), ':')==1){
        term.push_back(OpSeparable);
        if(OpSeparable.substr(0)=="+"){
            sep.push_back("+");
        }
    }
    else{
        tools::splitString(OpSeparable,"+",term,sep,"<(",">)");
    }

    // drop terms that do not match hierarchy (for sub-region calculation)
    for(int k=term.size()-1;k>=0;k--){
        if(term[k].rfind(":")==term[k].find_last_of(")]>")+1){
            std::vector<std::string> hier=tools::splitString(D->idx()->hierarchy(),'.');
            std::vector<std::string> subH=tools::splitString(term[k].substr(term[k].rfind(":")+1),'.');
            for(std::string s: subH){
                if(std::find(hier.begin(),hier.end(),s)==hier.end()){
                    PrintOutput::DEVmessage("SpectralProduct: removed "+term[k]+", does not match "+D->idx()->hierarchy());
                    term.erase(term.begin()+k);
                    sep.erase(sep.begin()+k);
                    break;
                }
            }
        }
    }

    for(auto s: sep){
        if(s!="+" and s!=" ") ABORT("Illegal separator - must be \"+\" or blank");
    }

    // by default, assume equal distribution of thresholds over terms
    double emin(Emin),emax(Emax);
    //    if(abs(emin)<DBL_MAX)emin/=term.size();
    //    if(abs(emax)<DBL_MAX)emax/=term.size();

    LOG_PUSH("FactorDiscretizations");

    // Use pointers to avoid copying/moving discretizations
    std::vector<std::unique_ptr<Discretization>> discretizations;
    std::vector<std::unique_ptr<Index>> complementIndices;
    std::vector<std::unique_ptr<OperatorTree>> ops;

    // loop through all terms
    for(size_t k=0;k<term.size();k++){
        string def=term[k].substr(0,term[k].find(":"));
        string coo=term[k].substr(term[k].find(":")+1);

        string complementCoo;


        // select axes
        std::vector<Axis> ax;
        std::vector<Axis> complementAx;
        for(size_t l=0;l<D->getAxis().size();l++){
            if(coo.find(D->getAxis()[l].name)!=string::npos)
                ax.push_back(D->getAxis()[l]);
            else{
                complementAx.push_back(D->getAxis()[l]);
                complementCoo += (complementCoo.size() == 0 ? "" : ".") + D->getAxis()[l].name;
            }
        }
        Index* xIdx=new IndexExtract(D->idx(),tools::splitString(coo,'.'),false,D->constraint);
        std::vector<std::string> compName,axName(tools::splitString(coo,'.')),allNames(tools::splitString(D->idx()->hierarchy(),'.'));
        for(string nam: allNames)
            if(     std::find(axName.begin(),axName.end(),nam)==axName.end() and
                    std::find(compName.begin(),compName.end(),nam)==compName.end())
                compName.push_back(nam);


        if(xIdx==0){
            PrintOutput::warning(Sstr+"Empty factor in DiscretizationSpectralProduct:"+coo+"not in"+D->idx()->hierarchy());
            continue;
        }
        discretizations.push_back(std::unique_ptr<Discretization>());


        if(compName.size() > 0){
            Index* idx = 0;

            // Hacky special case: Propagation in unbound subregion
            // BasicDisc does not support axis.basDef.funs=="grid"
            // (After all, it would need to know the grid points.)
            bool unboundDOF = false;
            for(auto& v: complementAx){
                if(v.name.find("nSurface")==0){
                    unboundDOF=true;
                    break;
                }
            }

            if(unboundDOF){
                auto hackIndex = [&complementAx](const Index* Idx, const std::unique_ptr<Discretization>& disc){
                    const Index* prodBottom = Idx;
                    while(prodBottom->axisName().find("spec")!=0){
                        prodBottom = prodBottom->descend();
                    }

                    std::vector<Axis> complementAxTop;
                    for(auto& v: complementAx){
                        if(v.name.find("nSurface")!=0 and v.name.find("kGrid")!=0){
                            complementAxTop.push_back(v);
                        }else{
                            break;
                        }
                    }

                    Index* complementTop = 0;
                    if(complementAxTop.size() > 0){
                        BasicDisc topDisc(complementAxTop, &IndexConstraint::main);
                        complementTop = topDisc.idx();
                        topDisc.idx() = 0;
                    }

                    bool inBottom=false;
                    std::string firstAxisNameBottom;
                    for(const Index* i=Idx; i!=0; i=i->descend()){
                        if(i->axisName().find("spec")!=0) inBottom=true;
                        if(inBottom && disc->idx()->hierarchy().find(i->axisName()) != std::string::npos){
                            firstAxisNameBottom = i->axisName();
                            break;
                        }
                    }

                    Index* bottom = disc->idx();
                    while(bottom->axisName()!=firstAxisNameBottom){
                        bottom = bottom->descend();
                    }

                    IndexQuot complementBottom(prodBottom, bottom);
                    Index* result;
                    if(complementTop != 0){
                        result = new IndexProd(complementTop, &complementBottom);
                        delete complementTop;
                    }else{
                        result = new Index(complementBottom);
                    }
                    return result;
                };


                idx = hackIndex(D->idx(), discretizations.back());
            }else{
                idx=new IndexExtract(D->idx(),compName,false,&IndexConstraint::main);
            }

            // Fallback: Works with and without unbound propagation,
            // but does not cover constraints
            // idx = new IndexQuot(D->Idx, discretizations.back()->idx());

            complementIndices.push_back(std::unique_ptr<Index>(idx));
        }else{
            complementIndices.push_back(std::unique_ptr<Index>(nullptr));
        }

        // get spectral discretizations
        ops.push_back(std::unique_ptr<OperatorTree>(
                          new OperatorTree(def, OperatorDefinition(def,coo).str(),xIdx, xIdx)));
        factors.push_back(std::unique_ptr<DiscretizationSpectral>(
                              new DiscretizationSpectral(ops.back().get(), emin, emax, excludeEnergyRange)));

        factors.back()->name = name + "_factor" + std::to_string(factors.size());
        if(factors.back()->idx()!=0)factors.back()->check(ops.back().get());

        if(factors.back()->idx() == 0 or factors.back()->idx()->sizeStored() == 0){
            discretizations.erase(discretizations.begin() + discretizations.size() - 1);
            factors.erase(factors.begin() + factors.size() - 1);
            ops.erase(ops.begin() + ops.size() - 1);
        }
    }
    for(auto& d: discretizations){
        // Drop ownership, as we need this index in eigenVectors
        if(d!=0)d->idx()=0;
    }

    LOG_POP();

    for(size_t i=0;i<factors.size(); i++){

        //----------------------------------------
        // for checking: get vector from factor subspace and apply operator
        // (this must be before modifications fot factors[i])
        Coefficients proC(ops[i]->jIndex);
        Coefficients ranC(ops[i]->iIndex);
        Coefficients opeC(ops[i]->iIndex);
        ranC.setToRandom();
        for(size_t k=0;k<ranC.size();k++)ranC.data()[k]=k;
        ranC.makeContinuous();
        factors[i]->projector()->apply(1.,ranC,0.,proC);
        ops[i]->apply(1., proC, 0., ranC);
        ops[i]->iIndex->inverseOverlap()->apply(1., ranC, 0., opeC);
        opeC.makeContinuous();
        //----------------------------------------

        // modify factors[i]: extend by identity levels
        factors[i]->parent = D;
        ParallelOperator par(factors[i]->mapFromParent());
        OperatorTree *mapFrom=new OperatorExpandIndex(dynamic_cast<const OperatorTree*>(factors[i]->mapFromParent()),D->idx(),false);
        OperatorTree *mapTo=  new OperatorExpandIndex(dynamic_cast<const OperatorTree*>(factors[i]->mapToParent()),  D->idx(),true);

        // make matching indices identical
        if(not mapTo->jdx()->treeEquivalent(mapFrom->idx()))
            DEVABORT("jdx\n"+mapTo->jdx()->str()+"\nidx\n"+mapFrom->idx()->str()+"\nnot equivalent: "+Index::failureCompatible);
        mapTo->replaceIndex(0,mapFrom->idx());

        factors[i]->idx() = const_cast<Index*>(mapFrom->idx());
        factors[i]->_mapFromParent=std::shared_ptr<OperatorTree>(mapFrom);
        factors[i]->_mapToParent=std::shared_ptr<OperatorTree>(mapTo);
        LOG_POP();


        /*
         * Spectral Values
         */
        LOG_PUSH("SpectralValues");
        OperatorDiagonal* newSpectralValues = new OperatorDiagonal(factors[i]->spectralValues->name,factors[i]->idx());
        newSpectralValues->setupAccordingToIndex();
        Coefficients cSpec(factors[i]->spectralOper()->iIndex);
        factors[i]->spectralOper()->storeInCoefficients(cSpec);
        Coefficients cSpecX(OperatorExpandIndex::coefficients(cSpec,factors[i]->idx()));

        newSpectralValues->setFromCoefficients(cSpecX);
        LOG_POP();

        // Store spectral values
        newSpectralValues->updateFunction(0., OperatorDiagonal::identityFunction);
        factors[i]->spectralValues = newSpectralValues;

        /*
         * Run a quick check as in DiscretizationSpectral::check(Op).
         * compare operator on factor subspace (previously computed c3) with operator on product subspace
         * Cannot use check(op), as this would require op \otimes ovr, we only get op \otimes Id
         */
        Coefficients cSpec1(factors[i]->idx());
        Coefficients cSpec2(factors[i]->idx());
        Coefficients cCheck(OperatorExpandIndex::coefficients(proC,factors[i]->mapFromParent()->jdx()));

        factors[i]->mapFromParent()->apply(1., cCheck, 0., cSpec1);
        factors[i]->spectralOper()->updateFunction(0., OperatorDiagonal::identityFunction);
        factors[i]->spectralOper()->apply(1., cSpec1, 0., cSpec2);
        factors[i]->mapToParent()->apply(1., cSpec2, 0., cCheck);

        double eps=1e-7;
#ifdef _DEVELOP_
        eps=1e-10;
#endif
        Coefficients cCheck2(OperatorExpandIndex::coefficients(opeC,D->idx()));
        if(cCheck.cwiseRelativeError(cCheck2).isZero(eps)){
            PrintOutput::DEVmessage(Sstr+"OK HP=U^-1 d U  for "+name+"["+i+"]");
        }else{
            PrintOutput::warning(
                        Str("HP=U^-1 d U not satisfied for")+name+"["+i+"], error ="+cCheck.norm(),50,0,
                        " This is a test on the spectral decomposition"
                        "\n Numerical errors are inevitable, but the level appears high and may compromize results"
                        "\n Possible fix: lower order in discretization, higher cutoff energy"
                        );
            if(cCheck.norm()>1.e-3)DEVABORT("serious error - cannot use decomposition - try higher or no spectral cut");
        }

        /*
         * We do not extend factors[i]->eigenVectors and factors[i]->eigenValues by identity, as these are more useful
         * in their original form (see ChannelsSubregion)
         */

    }
    if(not factors.size())idx() = 0;
    else idx()=const_cast<Index*>(D->idx());// Keep it from being deleted
}

DiscretizationSpectralProduct::Projector::Projector(const DiscretizationSpectralProduct* Parent):
    OperatorAbstract("projector", Parent->idx(), Parent->idx()),parent(Parent){
}

void DiscretizationSpectralProduct::Projector::apply(
        std::complex<double> A,
        const Coefficients& Vec,
        std::complex<double> B,
        Coefficients& Y) const{

    if(A != 1. or B != 0.) ABORT("Not implemented");

    Y = Vec;
    for(size_t i=0; i<parent->factors.size(); i++){
        Coefficients cSpec(parent->factors[i]->mapFromParent()->idx());
        parent->factors[i]->mapFromParent()->apply(1., Y, 0., cSpec);
        parent->factors[i]->mapToParent()->apply(-1., cSpec, 1., Y);
    }
}

const OperatorAbstract* DiscretizationSpectralProduct::projector() const{
    if(proj == 0){
        proj = std::unique_ptr<Projector>(new Projector(this));
    }

    return proj.get();
}

void DiscretizationSpectralProduct::checkFull(const OperatorAbstract* Projector, double cutE) const{
    LOG_PUSH("checkFull");
    /*
     * Check the projectors commute
     */
    for(size_t i=0; i<factors.size(); i++){
        for(size_t j=i+1; j<factors.size(); j++){
            PrintOutput::DEVmessage("Checking if projectors "+std::to_string(i)+" and "+std::to_string(j)+" commute");
            Coefficients c(idx());
            Coefficients c1(idx());
            Coefficients c2(idx());

            c.setToRandom();
            c1 = c;
            c2 = c;
            
            Coefficients cSpecI(factors[i]->idx());
            Coefficients cSpecJ(factors[j]->idx());
            
            factors[i]->mapFromParent()->apply(1., c, 0., cSpecI);
            factors[i]->mapToParent()->apply(1., cSpecI, -1., c1);
            factors[j]->mapFromParent()->apply(1., c1, 0., cSpecJ);
            factors[j]->mapToParent()->apply(1., cSpecJ, -1., c1);

            factors[j]->mapFromParent()->apply(1., c, 0., cSpecJ);
            factors[j]->mapToParent()->apply(1., cSpecJ, -1., c2);
            factors[i]->mapFromParent()->apply(1., c2, 0., cSpecI);
            factors[i]->mapToParent()->apply(1., cSpecI, -1., c2);

            c1 -= c2;
            if(c1.norm() / c2.norm() > 1.e-9){
                PrintOutput::warning("Projectors "+std::to_string(i)+" and "+std::to_string(j)+" do not commute");
            }else{
                PrintOutput::DEVmessage("Projectors "+std::to_string(i)+" and "+std::to_string(j)+" commute");
            }
        }
    }

    /*
     * Check eigenvectors with energy > cutE. These should all be projected out
     *
     * Note, this only works, if the operator given is in fact the projection, not the full Hamiltonian
     */
    {
        PrintOutput::DEVmessage("Checking eigenspace of projector for E > cutE");
        int total=0;
        int in_subspace=0;
        int in_complement=0;

        EigenSolver slv(cutE, DBL_MAX);
        slv.compute(Projector);
        for(auto v: slv.rightVectors()){
            total++;
            Coefficients c(Projector->jIndex);
            projector()->apply(1.,*v, 0., c);

            if(c.norm() / v->norm() < 1.e-9){
                in_complement++;
            }else{
                c -= *v;
                if(c.norm() / v->norm() < 1.e-9){
                    in_subspace++;
                }
            }
        }
        PrintOutput::DEVmessage("Checked projector for E > cutE: subspace/complement/total: "
                                +std::to_string(in_subspace)+"/"+std::to_string(in_complement)+"/"+std::to_string(total));

        if(in_subspace + in_complement != total){
            PrintOutput::warning("Projector is not diagonal with eigenvalues 0 and 1!");
        }

        if(in_subspace > 0){
            PrintOutput::warning("Projector does not project out all eigenvectors with E > cutE");
        }
    }

    /*
     * Check the eigenvectors with eigenvalues <= cutE, these should predominantly not be projected out,
     * however -- if the projection is a product -- some will. Anyway, P should be diagonal with
     * eigenvalues 0 and 1.
     *
     * Note, this only works, if the operator given is in fact the projection, not the full Hamiltonian
     */
    if(idx()->sizeStored() > 2000){
        PrintOutput::DEVwarning("Skipping test for E < cutE");
    }else{
        PrintOutput::DEVmessage("Checking eigenspace of projector for E < cutE");
        int total=0;
        int in_subspace=0;
        int in_complement=0;

        EigenSolver slv(DBL_MIN, cutE);
        slv.compute(Projector);
        for(auto v: slv.rightVectors()){
            total++;
            Coefficients c(Projector->jIndex);
            projector()->apply(1.,*v, 0., c);

            if(c.norm() / v->norm() < 1.e-9){
                in_complement++;
            }else{
                c -= *v;
                if(c.norm() / v->norm() < 1.e-9){
                    in_subspace++;
                }
            }
        }

        PrintOutput::DEVmessage("Checked projector for E < cutE: subspace/complement/total: "
                                +std::to_string(in_subspace)+"/"+std::to_string(in_complement)+"/"+std::to_string(total));

        if(in_subspace + in_complement != total){
            PrintOutput::warning("Projector is not diagonal with eigenvalues 0 and 1!");
        }
    }


    /*
     * Check that S P is a symmetric matrix. The matrix P corresponds to P^i_j, S to S_ij, so S P
     * corresponds to P_ij.
     *
     * This should hold even for ECS.
     */
    if(idx()->sizeStored() < 200){
        PrintOutput::DEVmessage("Checking orthogonality property");
        UseMatrix mat;
        OperatorAbstractProduct projector_cov{"test", { idx()->overlap(), projector() }};
        projector_cov.matrixAdd(1., mat);

        if(not mat.isSymmetric(1.e-9) and not mat.isZero(1.e-9)){
            PrintOutput::warning("Spectral projection is not symmetric");
        }else{
            PrintOutput::DEVmessage("SP is symmetric");
        }

    }else{
        PrintOutput::DEVmessage("Checking orthogonality for 1000 random vectors");
        OperatorAbstractProduct projector_cov{"test", { idx()->overlap(), projector() }};

        bool fail = false;
        for(int i=0; i<1000; i++){
            Coefficients c1(idx());
            c1.setToRandom();
            c1.makeContinuous();

            Coefficients c2(idx());
            c2.setToRandom();
            c2.makeContinuous();

            std::complex<double> m1 = projector_cov.matrixElement(c1, c2, true);
            std::complex<double> m2 = projector_cov.matrixElement(c2, c1, true);
            if(std::abs(m1 - m2) > 1.e-9){
                PrintOutput::warning("SP is not symmetric");
                fail = true;
            }
        }
        if(not fail){
            PrintOutput::DEVmessage("SP is symmetric");
        }
    }

    LOG_POP();
}
