// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "visualizeAngularDistribution.h"

#include <sstream>

#include "index.h"
#include "coefficients.h"
#include "printOutput.h"

VisualizeAngularDistribution::VisualizeAngularDistribution(const Index* Root){
    const Index* level = Root;
    bool started = false;
    for(; level->descend()!=0; level=level->descend()){
        if(level->axisName().find("Phi") == std::string::npos and level->axisName().find("Eta") == std::string::npos){
            if(started) break;
        }else started = true;
    }
    for(; level!=0; level=level->nodeRight()){
        data[level] = 1.;
    }
}

VisualizeAngularDistribution::VisualizeAngularDistribution(const Coefficients* Coeff){
    const Coefficients* level = Coeff;
    for(; level->descend()!=0; level=level->descend()){
        if(level->idx()->axisName().find("Phi") == std::string::npos and
                level->idx()->axisName().find("Eta") == std::string::npos) break;
    }
    for(; level!=0; level=level->nodeRight()){
        data[level->idx()] = level->norm();
    }
}

void VisualizeAngularDistribution::print(){
    std::vector<int> x;
    std::vector<int> y;
    std::vector<double> val;

    for(auto p: data){
        for(const Index* idx = p.first; idx->parent()!=0; idx=idx->parent()){
            if(idx->parent()->axisName() == x_axis){
                x.push_back(idx->physical());
            }else if(idx->parent()->axisName() == y_axis){
                y.push_back(idx->physical());
            }
        }
        val.push_back(p.second);
    }

    if(x.size() == 0 or x.size() != y.size() or x.size() != val.size()){
        // Unexpected, but VisualizeAngularDistribution is not vital
        return;
    }

    int x_min = *std::min_element(x.begin(), x.end());
    int x_max = *std::max_element(x.begin(), x.end());
    int y_min = *std::min_element(y.begin(), y.end());
    int y_max = *std::max_element(y.begin(), y.end());

    if(x_min < 0){
        int tmp = std::max(-x_min, x_max);
        x_min = -tmp;
        x_max = tmp;
    }

    if(y_min < 0){
        int tmp = std::max(-y_min, y_max);
        y_min = -tmp;
        y_max = tmp;
    }

    PrintOutput::verbatim("\n" + y_axis + "\n");
    PrintOutput::verbatim("["+std::to_string(y_min)+","+std::to_string(y_max)+"]\n");
    for(int Y=y_max; Y>=y_min; Y--){
        std::ostringstream line;
        line<<"|";
        for(int X=x_min; X<=x_max; X++){
            int i=0;
            for(;i<val.size(); i++){
                if(x[i] == X and y[i] == Y) break;
            }
            line << ((i<val.size() && val[i] > 0.1) ? "X" : " ");
        }
        PrintOutput::verbatim(line.str() + "\n");
    }
    std::ostringstream bottom;
    bottom<<"+";
    for(int X=x_min; X<=x_max; X++) bottom<<"-";
    bottom<<"["<<x_min<<","<<x_max<<"] "<<x_axis;
    PrintOutput::verbatim(bottom.str() + "\n\n");

}

void VisualizeAngularDistribution::printAllPossible(const Index* Idx){
    std::string hier = Idx->hierarchy();
    if(hier.find("Phi") != std::string::npos and hier.find("Eta") != std::string::npos){
        PrintOutput::lineItem("Angular constraints","");
        if(hier.find("Eta.")!=std::string::npos and hier.find("Phi.")!=std::string::npos){
            VisualizeAngularDistribution(Idx)
                .withAxes("Phi", "Eta")
                .print();
        }
        if(hier.find("Eta1.")!=std::string::npos and hier.find("Phi1.")!=std::string::npos){
            VisualizeAngularDistribution(Idx)
                .withAxes("Phi1", "Eta1")
                .print();
        }
        if(hier.find("Eta2.")!=std::string::npos and hier.find("Phi2.")!=std::string::npos){
            VisualizeAngularDistribution(Idx)
                .withAxes("Phi2", "Eta2")
                .print();

        }
        if(hier.find("Phi1.")!=std::string::npos and hier.find("Phi2.")!=std::string::npos){
            VisualizeAngularDistribution(Idx)
                .withAxes("Phi1", "Phi2")
                .print();
        }
        if(hier.find("Eta1.")!=std::string::npos and hier.find("Eta2.")!=std::string::npos){
            VisualizeAngularDistribution(Idx)
                .withAxes("Eta1", "Eta2")
                .print();
        }
        
    }
}
