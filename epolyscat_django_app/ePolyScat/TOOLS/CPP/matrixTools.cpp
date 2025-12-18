// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "abort.h"
#include "stringTools.h"
#include "tools.h"
#include "str.h"

namespace tools{

/** @defgroup  StringTools
 *  @ingroup Tools
 *  \brief String-specific tools
 * @{
 */

string strMatrix(complex<double>* Mat,int Rows,int Cols,int Digits,string Text){
    ostringstream oss;
    oss<<"\n";
    // column numbering
    if(Text.length()>9)oss<<Text<<endl<<"row\\col  ";
    else oss<<setw(Text.length())<<Text<<setw(9-Text.length())<<"";

    if(Digits>0)
        for (unsigned int n=0;n<Cols;n++){
            if(n==0)oss<<setw(Digits+1)<<0;
            else oss<<setw(Digits+7)<<n;
        }
    oss<<endl;

    for (unsigned int m=0;m<Rows;m++){
        // row numbering
        oss<<setw(4)<<m;
        // real parts
        for (unsigned int n=0;n<Cols;n++){
            if(Digits>0)oss<<setprecision(Digits)<<setw(Digits+7)<<Mat[m+n*Rows].real();
            else        oss<<tools::zero(Mat[m+n*Rows]);
        }
        oss<<endl;

        // imaginary parts
        if(Digits>0){
            oss<<setw(4)<<" ";
            if(Digits>0)
                for (unsigned int n=0;n<Cols;n++)
                    oss<<setprecision(Digits)<<setw(Digits+7)<<Mat[m+n*Rows].imag();
            oss<<endl;
        }
    }
    return oss.str();
}


string strMatrixBlock(complex<double>* Mat, const vector<int>IBlock, const vector<int>JBlock, int Digits, string Text){
    int Rows=0;
    for(auto b: IBlock)Rows+=b;
    vector<complex<double> > bMat(IBlock.size()*JBlock.size(),0.);
    for(int j=0,jB0=0;j<JBlock.size();jB0+=JBlock[j++])
        for(int i=0,iB0=0;i<IBlock.size();iB0+=IBlock[i++])
            for(int jj=0;jj<JBlock[j];jj++)
                for(int ii=0;ii<IBlock[i];ii++)
                    if(abs(bMat[i+IBlock.size()*j])<abs(Mat[iB0+ii+Rows*(jB0+jj)]))
                        bMat[i+IBlock.size()*j]=Mat[iB0+ii+Rows*(jB0+jj)];

    return strMatrix(bMat.data(),IBlock.size(),JBlock.size(),Digits,Text);
}

string strMatrixBlock(complex<double>* Mat,int Rows,int Cols,int RowBlock, int ColBlock,int Digits,string Text){
    if(ColBlock==0)ColBlock=RowBlock;
    int iBlocks=(Rows-1)/RowBlock+1;
    int jBlocks=(Cols-1)/ColBlock+1;

    vector<int> iBlock,jBlock;
    iBlock.assign(iBlocks,RowBlock);iBlock.back()=min(iBlock.back(),int(Rows-RowBlock*(iBlock.size()-1)));
    jBlock.assign(jBlocks,ColBlock);jBlock.back()=min(jBlock.back(),int(Cols-ColBlock*(jBlock.size()-1)));

    return strMatrixBlock(Mat,iBlock,jBlock,Digits,Text);
}

}
