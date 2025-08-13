// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "channelsSubregion.h"

#include "mpiWrapper.h"
#include "index.h"
#include "operatorTree.h"
#include "discretizationSpectralProduct.h"
#include "operatorDefinition.h"
#include "readInput.h"
#include "plot.h"
#include "printOutput.h"
#include "operatorDiagonal.h"
#include "asciiFile.h"
#include "readInput.h"
#include "averageOverAxes.h"
#include "pulse.h"
#include "discretizationGrid.h"
#include "parallelOperator.h"

using namespace std;

double ChannelsSubregion::tStore;
ChannelsSubregion::Config ChannelsSubregion::config;

void ChannelsSubregion::read(ReadInput& Inp){
    std::string def = "0";
    if(Units::isDefined("OptCyc")) def="0.25 OptCyc";
    Inp.read("ChannelsSubregion", "tStore", tStore, def, "Store spectra timestep");
}

ChannelsSubregion::ChannelsSubregion(const Discretization* D, const OperatorTree* Ham, std::string Region, ReadInput& Inp):
    spec("si_spec_"+Region, "", tStore), tBeginAverage(-DBL_MAX)
  #ifdef _CHANNELS_SUBREGION_GRID_
  , specGrid("si_spec_grid_"+std::to_string(Subregion), "", tStore)
  #endif
{

    /*
         * Hardcoded configuration for Helium 6D
         */
    if(D->idx()->hierarchy().find("Eta1")!=string::npos and D->idx()->hierarchy().find("Eta2")!=string::npos){
        if(D->idx()->hierarchy().find("kRn1")!=string::npos)
            config=Config({"(1/2<<Laplacian>>-2<<Coulomb>>):Phi2.Eta2.Rn2", "kRn1", {"Phi2", "Eta2", "Eigen"}, {"Phi1", "Eta1"}});
        if(D->idx()->hierarchy().find("kRn2")!=string::npos)
            config=Config({"(1/2<<Laplacian>>-2<<Coulomb>>):Phi1.Eta1.Rn1", "kRn2", {"Phi1", "Eta1", "Eigen"}, {"Phi2", "Eta2"}});
    }


    if(config.kAxis=="") ABORT("Specify configuration in ChannelsSubregion for: "+D->idx()->hierarchy());

    // Set up spectral discretization and check that at least one ionic state was found
    disc = unique_ptr<DiscretizationSpectralProduct>(
                new DiscretizationSpectralProduct(D, config.hamiltonian, -DBL_MAX, 0., false));

    if(disc->factors.size() == 0) ABORT("No channels found for"+config.str());

    for(size_t i=0; i<disc->factors[0]->eigenvalues().size(); i++){
        PrintOutput::newLine();
        PrintOutput::lineItem("Channel "+std::to_string(i), disc->factors[0]->eigenvalues()[i].real());
    }

    coeff = unique_ptr<Coefficients>(
                new Coefficients(disc->factors[0]->idx()));

    avg = unique_ptr<Coefficients>(
                new Coefficients(disc->factors[0]->idx()));

    temp = unique_ptr<Coefficients>(
                new Coefficients(disc->factors[0]->idx()));

    // Necessary for AllreduceSUM
    coeff->treeOrderStorage();

    avg->setToZero();

    // Set distribution of all children to "all" to bypass ParallelOperator in OperatorTree::apply
    const OperatorTree* mapFromParent = dynamic_cast<const OperatorTree*>(disc->factors[0]->mapFromParent());
    for(size_t i=0; i<mapFromParent->childSize(); i++)
        ParallelOperator::bcast(mapFromParent->child(i));

#ifdef _CHANNELS_SUBREGION_GRID_
    std::string phiAxis;
    std::string etaAxis;
    for(auto v: config[subregion].partialAxes){
        if(v.find("Phi") == 0) phiAxis = v;
        else if(v.find("Eta") == 0) etaAxis = v;
    }

    if(phiAxis != "" and etaAxis != ""){
        grid = std::unique_ptr<DiscretizationGrid>(new DiscretizationGrid(
                                                       disc->factors[0].get(),
        { phiAxis, etaAxis },
        { 100, 9 },
        { {0., 2*math::pi}, {-1., 1.} }
                                                   ));
        coeffGrid = std::unique_ptr<Coefficients>(new Coefficients(grid->idx()));
    }
#endif
}

void ChannelsSubregion::write(double Time){
    if(true or tStore == 0.){
        // Averaging disabled
        spec.write(Time, *coeff);

#ifdef _CHANNELS_SUBREGION_GRID_
        if(grid != 0){
            grid->mapFromParent()->apply(1., *coeff, 0., *coeffGrid);
            specGrid.write(Time, *coeffGrid);
        }
#endif

    }else{
        // Averaging enabled

        if(tBeginAverage == -DBL_MAX){
            tBeginAverage = Pulse::gettEnd();
            tEndAverage = tBeginAverage + tStore;
            tLastAverage = tBeginAverage;

            // Mark Pulse end time by a dummy entry
            // Set order==1 in spec.py to account for this
            avg->setToZero();
            spec.write(Pulse::gettEnd(), *avg);

#ifdef _CHANNELS_SUBREGION_GRID_
            if(grid != 0){
                coeffGrid->setToZero();
                specGrid.write(Pulse::gettEnd(), *coeffGrid);
            }
#endif
        }

        if((Time - tEndAverage) > -1.e-9){
            *coeff *= (tEndAverage - tLastAverage);
            *avg += *coeff;
            *coeff *= 1./(tEndAverage - tLastAverage);
            *avg *= 1./tStore;

            spec.write(tEndAverage, *avg);

#ifdef _CHANNELS_SUBREGION_GRID_
            if(grid != 0){
                grid->mapFromParent()->apply(1., *avg, 0., *coeffGrid);
                specGrid.write(tEndAverage, *coeffGrid);
            }
#endif

            tBeginAverage = tEndAverage;
            tEndAverage = tBeginAverage + tStore;
            tLastAverage = tBeginAverage;
            avg->setToZero();
        }

        *coeff *= (Time - tLastAverage);
        *avg += *coeff;
        tLastAverage = Time;
    }
}

void ChannelsSubregion::getConfigurations(const Discretization* D, ReadInput& Inp, std::vector<std::string> & configurations){
    // The function determines the channel configurations from the original Hamiltonian definition.
    std::vector<std::string> infCoorNames; // List of all infinite coordinates
    for(unsigned int i=0; i<D->continuityLevel.size(); i++){
        infCoorNames.push_back(D->getAxis()[D->continuityLevel[i]].name);
    }
    std::string Hamiltonian;
    Inp.read("Operator", "hamiltonian", Hamiltonian, ReadInput::noDefault, "Hamiltonian definition");
    Hamiltonian = OperatorDefinition::unBracket(Hamiltonian).str();

    std::vector<std::string> terms;
    std::vector<std::string> factors;
    int pos0 = 0;
    int pos1 = Hamiltonian.find("<",0);
    int pos2 = Hamiltonian.find(">",0);
    while((pos1!=string::npos)and(pos2!=string::npos)){
        terms.push_back(Hamiltonian.substr(pos1+1,pos2-pos1-1)); // without "<>" brackets
        factors.push_back(Hamiltonian.substr(pos0,pos1-pos0));
        pos0=pos2+1;
        pos1 = Hamiltonian.find("<",pos1+1);
        pos2 = Hamiltonian.find(">",pos2+1);
    }

    int dim = configurations.size();

    std::vector<std::string> PotentialFactors;
    std::vector<std::string> Potentials;

    for(unsigned int c=0; c<dim; c++){
        for(unsigned int i=0; i<terms.size()/dim; i++){
            // iterate through all terms of a certain coordinate
            if(terms[i*dim+c].find("d_1_d")!=string::npos){
                configurations[c]+=factors[i*dim]+"<"+terms[i*dim+c]+">";
            }
            // replace all "Q" by their respective coordinate
            if(terms[i*dim+c].find("Q")){
                replaceAll(terms[i*dim+c],"Q",infCoorNames[c]);
            }
            //HACK this is problem-specific and should not remain here
            if(terms[i*dim+c].find("/sqrt")){
                std::vector<std::string> Coordinates;
                for(unsigned int j=0; j<infCoorNames.size(); j++){
                    if(terms[i*dim+c].find(infCoorNames[j])!=string::npos){
                        Coordinates.push_back(infCoorNames[j]);
                    }
                }
                if((Coordinates.size()==1) and (Coordinates[0]==infCoorNames[c])){
                    if(terms[i*dim+c].find(infCoorNames[c])){
                        replaceAll(terms[i*dim+c],infCoorNames[c],"Q");
                    }
                    configurations[c]+=factors[i*dim]+"<"+terms[i*dim+c]+">";
                }
                if(Coordinates.size()>1){
                    PotentialFactors.push_back(factors[i*dim]);
                    Potentials.push_back(terms[i*dim+c]);
                }
            }
        }
    }
    std::vector<std::string> coordinatesUsed;
    for(unsigned int k=0; k<Potentials.size(); k++){
        string coor="";
        for(unsigned int l=0; l<infCoorNames.size(); l++){
            if(Potentials[k].find(infCoorNames[l])!=string::npos){
                coor+=infCoorNames[l];
            }
        }
        bool newCombination=false;
        for(unsigned int i=0; i<coordinatesUsed.size(); i++){
            if(coordinatesUsed[i]==coor){
                newCombination=true;
                break;
            }
        }
        if(newCombination==false){
            coordinatesUsed.push_back(coor);
            std::string potential;
            for(unsigned int c=0; c<dim; c++){
                for(unsigned int j=0; j<infCoorNames.size(); j++){
                    if(j!=c){ // removed coordinate
                        if(isCoMCoor(infCoorNames[j])==true){
                            // remove infCoorNames[j] from Potential
                            int pos1 = Potentials[k].find("/sqrt");
                            int pos2 = Potentials[k].find_last_of(")");
                            std::string insideSqrt = Potentials[k].substr(pos1+1, pos2-pos1-1);
                            int pos3 = insideSqrt.find_first_of("(");
                            int pos4 = insideSqrt.find_last_of(")");
                            std::string insidePow;
                            std::string outsidePow;
                            if(pos3!=string::npos and pos4!=string::npos){
                                // Is "insidePow" actually necessary?
                                insidePow = insideSqrt.substr(pos3+1, pos4-pos3-1); // axes
                                outsidePow = insideSqrt.substr(pos4+1,pos2-pos4-1); // additional constant
                            }
                            // split terms before the square root:
                            std::vector<std::string> termsBeforeSqrt;
                            std::string NewTerms;
                            std::string beforeSqrt = Potentials[k].substr(0,pos1);
                            if(beforeSqrt.find("*")==string::npos){
                                NewTerms=beforeSqrt;
                            }
                            else{
                                termsBeforeSqrt = tools::splitString(beforeSqrt,'*');
                                for(unsigned int i=0;i<termsBeforeSqrt.size();i++){
                                    // check every term for trunc[a,b](Q)
                                    if(termsBeforeSqrt[i].find("trunc")!=string::npos){
                                        int pos5 = termsBeforeSqrt[i].find("]");
                                        // Check whether a coordinate is written directly after trunc[a,b]
                                        if (termsBeforeSqrt[i].substr(pos5+1,1).find("(")!=string::npos){
                                            int pos6 = termsBeforeSqrt[i].find(")",pos5);
                                            std::string truncCoor = termsBeforeSqrt[i].substr(pos5+2,pos6-pos5-2);
                                            if(truncCoor==infCoorNames[c]){
                                                if(NewTerms.length()==0){
                                                    // leaving out the "(Q)" is not necessary
                                                    NewTerms+=termsBeforeSqrt[i].substr(0,pos5+1);
                                                    //                                                    NewTerms+=termsBeforeSqrt[i];
                                                }
                                                else{
                                                    NewTerms+="*"+termsBeforeSqrt[i].substr(0,pos5+1);
                                                    //                                                    NewTerms+="*"+termsBeforeSqrt[i];
                                                }
                                            }
                                        }
                                    }
                                    else{
                                        if(NewTerms==""){
                                            NewTerms+=termsBeforeSqrt[i];
                                        }
                                        else{
                                            NewTerms+="*"+termsBeforeSqrt[i];
                                        }
                                    }
                                }
                            }
                            // needs to be generalized for 3 axes
                            potential = NewTerms+"/sqrt("+infCoorNames[c]+"*"+infCoorNames[c]+outsidePow+")";
                            configurations[c]+=PotentialFactors[k]+"<"+potential+">";
                        }
                    }
                }
                // replace all coordinates by "Q"
                ChannelsSubregion::replaceAll(configurations[c], infCoorNames[c], "Q");
            }
        }
    }
    for(unsigned int c=0;c<configurations.size(); c++){
        cout << "configurations[c] = " << configurations[c] << endl;
    }
}

//HACK this is problem-specific and should not remain here
bool ChannelsSubregion::isCoMCoor(std::string Coordinate){
    // returns true if Coordinate is a center of mass coordinate, else: false
    // To be improved (determine it from the Hamiltonian)
    if(Coordinate=="Z") return true;
    return false;
}

// Could be moved to stringTools.cpp
void ChannelsSubregion::replaceAll(std::string& s, const std::string& a, const std::string& b){
    // The function replaces all occurrences of the substring s in a given string object
    if(s=="") return;
    int pos=s.find(a,0);
    while(pos!=string::npos){
        s.replace(pos, a.size(), b);
        pos=s.find(a,pos+1);
        pos=s.find(a);
    }
}

void ChannelsSubregion::average(const Coefficients* C, double Time){
    if(disc==0) return;
    if(Time < Pulse::gettEnd()) return;

    disc->factors[0]->mapFromParent()->apply(1.,*C,0.,*temp);
    disc->factors[0]->spectralOper()->updateFunction(Time, OperatorDiagonal::expItFunction);
    disc->factors[0]->spectralOper()->apply(1.,*temp,0.,*coeff);

    write(Time);
}

void ChannelsSubregion::parallelAverage(const Coefficients* C, double Time){
    if(disc==0) return;
    if(Time < Pulse::gettEnd()) return;

    const OperatorTree* mapFromParent = dynamic_cast<const OperatorTree*>(disc->factors[0]->mapFromParent())
            ->child(MPIwrapper::Rank());
    OperatorDiagonal* spectralValues = disc->factors[0]->spectralOper()->child(MPIwrapper::Rank());
    C = C->child(C->childSize() == 1 ? 0 : MPIwrapper::Rank());
    Coefficients* tempC = temp->child(MPIwrapper::Rank());
    Coefficients* coeffC = coeff->child(MPIwrapper::Rank());

    temp->setToZero();
    coeff->setToZero();

    mapFromParent->apply(1., *C, 0., *tempC);
    spectralValues->updateFunction(Time, OperatorDiagonal::expItFunction);
    spectralValues->apply(1., *tempC, 0., *coeffC);

    // Could be replaced by Gather
    MPIwrapper::AllreduceSUM(coeff->storageData(), coeff->size());

    write(Time);
}

