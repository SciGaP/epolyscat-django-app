#include "basisDyson.h"
#include "printOutput.h"

#include "asciiFile.h"
BasisDyson::BasisDyson(const VectorValuedFunction * MO, const std::vector<std::vector<double>> & DysonC, std::string Source)
    :_source(Source),_mo(MO),_dysonC(DysonC)
{
}


BasisDyson* BasisDyson::read(const VectorValuedFunction *MO, std::string DysonCoefs){
    std::vector<std::vector<double>> cols;
    if(DysonCoefs.find("Dummy")!=std::string::npos){
        std::vector<int> molist=tools::string_to_intVec("{"+tools::stringInBetween(DysonCoefs,"[","]")+"}");
        for(int i: molist){
            cols.push_back(std::vector<double>(molist.size(),0.));
            cols.back()[i]=1.;
        }
    }
    else {
        AsciiFile dysf(DysonCoefs);
        if(dysf.empty())return nullptr;
        dysf.readCols(cols);
        if(cols.size()==0)return nullptr;
        // check for match with MOs
        std::vector<std::string> com;
        dysf.readComments(com);
        //    for(auto l: com){
        //        if(l.find("#BasisMO")==0 and l.find(MO->str())==std::string::npos)
        //            ABORT("Dyson coefficients "+DysonCoefs+" are not for "+MO->str());
        //    }
    }
    BasisDyson* b=new BasisDyson(MO,cols,DysonCoefs);
    b->print();
    return b;
}

std::vector<std::complex<double>> BasisDyson::operator()(std::vector<double> X) const{
    std::vector<std::complex<double>> v=_mo->operator()(X);
    if(v.size()!=_dysonC[0].size())
        ABORT(Sstr+"number Dyson coefficients "+_dysonC[0].size()+"!="+v.size()+"number MOs, source "+_source);
    std::vector<std::complex<double>> res(_dysonC.size(),0.);
    for(size_t k=0;k<res.size(); k++){
        for(size_t l=0;l<_dysonC[k].size();l++)res[k]+=_dysonC[k][l]*v[l];
    }
    return res;
}

static void printColsTranspose(const std::vector<std::vector<double>> & Cols){
    double Eps=1.e-5;
    size_t rows=0;
    for(auto c: Cols)rows=std::max(rows,c.size());

    PrintOutput::newRow();
    PrintOutput::rowItem("No.");
    for(size_t k=0;k<rows;k++)PrintOutput::rowItem(int(k));

    for(size_t c=0;c<Cols.size();c++){
        PrintOutput::newRow();
        PrintOutput::rowItem(int(c));
        for(size_t r=0;r<Cols[c].size();r++){
            if(std::abs(Cols[c][r])<Eps)PrintOutput::rowItem(".");
            else                        PrintOutput::rowItem(Cols[c][r],2);
        }
    }
    PrintOutput::paragraph();
}


void BasisDyson::print() const{
    //    _mo->print();
    PrintOutput::title("Dyson coefficients");
    printColsTranspose(_dysonC);
}
