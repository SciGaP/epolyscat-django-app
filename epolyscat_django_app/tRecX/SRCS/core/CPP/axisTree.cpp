// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license

#include "axisTree.h"
#include "readInput.h"
#include "basisVector.h"


// nonsensical input with complicated hierarchy
static std::vector<std::string> inputExample=
{
    "Axis: subset,name,functions,nCoefficients,lower end, upper end,order,exactIntegral",
    "        , Vec,,2",
    "        , Vec1,,2",
    "Subspace, Orbital,Eigenbasis[0.5<<Laplacian>>-<<Coulomb>>:Complement],1,2,,1",
    "        , Orbital1,Eigenbasis[0.5<<Laplacian>>-<<Coulomb>>:Complement],1,2,,1",
    "  @sub1 , Vec2,,2",
    "  @sub2 , Vec3,,3",
    " @@sub2 , Vec6,,3",
    " @@sub3 , Vec7,,3",
    "Complement,Phi,,1",
    "        , Eta,assocLegendre{Phi},4,-1,1",
    "        , Rn,polynomial,        40, 0.,40.,20",
    "        @sub3 , Vec4,,2",
    "        @sub4 , Vec5,,3",
    "       @@sub3 , Vec8,,3",
    "       @@sub4 , Vec9,,3",
    "Nonsense,Phi,,1",
    "        , Eta,assocLegendre{Phi},4,-1,1",
    "        , Rn,polynomial,        40, 0.,40.,20"
};

std::vector<Axis> AxisTree::toVector() const {
    std::vector<Axis> ax;
    for(const AxisTree* axt=this;axt!=0;axt=axt->nodeNext()){
        if(axt->childSize()>1)DEVABORT("conversion to vector only for trivial AxisTree, is\n"+Tree::str(1));
        ax.push_back(*axt);
    }
    return ax;
}

bool nextLowerSubset(std::string Current,std::string Next){
    // true if Next subset is same or next level below Current
    return (Current=="" and Next!="" and Next.find("@")==std::string::npos)  // first level below root
            or (Current!="" and std::count(Current.begin(),Current.end(),'@')+1==std::count(Next.begin(),Next.end(),'@')); // next lower hierarchy
}

std::string AxisTree::readSubset(ReadInput &Inp, std::string AxisCategory, int Line, std::string Default){
    std::string res;
    Inp.read(AxisCategory,"subset",res,Default,"axis will be read into named subset, see UserManual.pdf for details",Line)
            .texdocu(R"tex(
                     This is used for hybrid discretizations.
                     Two subsets may appear at the same node.
                     In that case the two branches will be treated as distinct discretizations.
                     Operators and overlaps will be computed as hybrid.
                     Typical application is a 3d numerical vs. a standard tree basis with matrix
                     elements:
                     \[
                     \l \Phi_n(x,y,z)| Op | A_i(\phi)B_j(\eta)C_k(r)\r
                     \]
                     (See \nameref{eq:hybrid})
                     )tex");
 return res;
}

static int hybridCount=0;
AxisTree::AxisTree(ReadInput &Inp, int &Line, std::string PreviousSubset, std::string AxisCategory, std::string AlternateIndex):
    Axis(Inp,Line,AxisCategory,AlternateIndex)
{
    if(Line==1)hybridCount=0;
    // advance read to axis that is part of AlternateIndex
    int line(Line);

    Inp.texdocuCategoryAdd(AxisCategory,"name,lower end,upper end,nCoefficients,functions,order,subset",
  R"tex(
 A single line defines a single coordinate.\\
 Multiple lines are arranged into a tree-hierarchy.\\
 All branches of a node discretize the same subspace, but possibly with different coordinates
 and basis functions. In that case the discretization will be treated as ``hybrid''
 \\\\
 Example:
 \begin{verbatim}
 Axis: subset,name,nCoefficients,lower end, upper end,functions,order,alternate
        Neut, Neutral,1,,,CI[Neutral]
        Chan, Ion,numI,,,CI[Ion]
     @MolOrb, Orbital,2,0,1,CO2{9,10}:dense
 @Complement, Phi,3,,,,,               dense
     ,        Eta,3,-1,1, assocLegendre{Phi},,dense[order=3]
     ,        Rn,40, 0,BOX,polynomial,20,     dense[order=40]
     ,        Rn,10, BOX,Infty,polExp[1.]
 \end{verbatim}
 produces the axis hierarchy:
 \begin{verbatim}
  |....main: Neutral&Ion Approx. size (elements, order) 3 (1, 2) Hybrid:2
   0 |....Neut: Neutral Approx. size (elements, order) 1 (1, 1)  CI[...
   1 |....Chan: Ion Approx. size (elements, order) 2 (1, 2)  CI[Ion]...
   1   0 |....Chan: Orbital&Phi Approx. size (elements, order) 3 (1, 2) Hybrid:2
   1    0    0 |....@MolOrb: Orbital Approx. size (elements, order) 2 (1, 2)
   1    0    1 |....@Complement: Phi Approx. size (elements, order) 3 (1, 3)
   1    0    1    0 |....@Complement: Eta Approx. size (elements, order) 3 (1, 3)
   1    0    1    0    0 |....@Complement: Rn Approx. size ....
 \end{verbatim}
  )tex",
  "00,01,09,02,04,07,20,23,220,510");

    // get valid subset at present level by scanning all from Line present line
    _subset=PreviousSubset;
    while(Line<=line){
//        Inp.read(AxisCategory,"subset",_subset,_subset,"axis will be read into named subset, see UserManual.pdf for details",Line++)
//                .texdocu(R"tex(
//                         This is used for hybrid discretizations.
//                         Two subsets may appear at the same node.
//                         In that case the two branches will be treated as distinct discretizations.
//                         Operators and overlaps will be computed as hybrid.
//                         Typical application is a 3d numerical vs. a standard tree basis with matrix
//                         elements:
//                         \[
//                         \l \Phi_n(x,y,z)| Op | A_i(\phi)B_j(\eta)C_k(r)\r
//                         \]
//                         (See \nameref{eq:hybrid})
//                         )tex");
        _subset=AxisTree::readSubset(Inp,AxisCategory,Line++,_subset);

    }
    Line=line;

    if(basDef.size()==0 and AlternateIndex=="")DEVABORT(Str("emtpy basis, but non-empty Axis line")+Line);

    if(Line==1 and _subset!=PreviousSubset){
        if(basDef.size()>0){
            // list starts with subset, add root node
            childAdd(new AxisTree(Inp,Line,_subset,AxisCategory,AlternateIndex));
        }
        _subset="main";
    }

    // concatenate multi-part axes
    std::string oldName=name;
    // cannot have the elegant logics below: dynamical default inconsistent with multiple reads
    //    while(not Inp.endCategory(AxisCategory,Line+1) and name==Axis::readName(Inp,Line+1,name)){
    while(not Inp.endCategory(AxisCategory,Line+1)){
        std::string newName=Axis::readName(Inp,Line+1,"BLANK",AxisCategory);
        if(newName=="BLANK")newName=oldName;
        if(newName!=oldName)break;
        name=oldName;
        Axis xx(Inp,++Line,AxisCategory,AlternateIndex);
        if(xx.basDef.size()>0)appendInPlace(xx);
    }

    while(not Inp.endCategory(AxisCategory,Line+1)){
        std::string subNext;
        subNext=AxisTree::readSubset(Inp,AxisCategory,Line+1,_subset);


        if(_subset=="" or _subset=="main" or _subset==subNext or nextLowerSubset(_subset,subNext))
            // next in present subtree - descend in axis hiearchy
            childAdd(new AxisTree(Inp,++Line,_subset,AxisCategory,AlternateIndex));
        else
            // next not in present axis (sub-)tree, return to higher level
            return;
    }
    if(_subset=="main"){
        // root node was added - construct suitable name
        bases={BasisAbstract::factory("Hybrid: 2")};
        name=child(0)->name;
        for(size_t k=1;k<childSize();k++){
            name+="&"+child(k)->name;
        }
    }
    if(childSize()>1 and name.find("&")==std::string::npos){
        if(childSize()!=2)ABORT("cannot have more than two subsets following one Axis, have: \n"+Tree::str());
        AxisTree * curAx=deepCopy();
        for(size_t k=childSize();k>0;k--)childErase(k-1);
        childAdd(curAx);
        childBack()->name=curAx->child(0)->name+"&"+curAx->child(1)->name;
        childBack()->bases={BasisAbstract::factory("Hybrid: 2")};
    }
}

AxisTree::AxisTree(const std::vector<Axis> &Axes):Axis(Axes[0]),_subset(""){
    if(Axes.size()==0)ABORT("cannot construct AxisTree from empty Axis vector");
    if(Axes.size()>1)childAdd(new AxisTree(std::vector<Axis>(Axes.begin()+1,Axes.end())));
}

void AxisTree::nodeCopy(const AxisTree *Node, bool View){
    _subset=Node->_subset;
    Axis::operator=(*Node);
    if(View)DEVABORT("cannot have View of AxisTree");
}
std::string AxisTree::strNode(int Level) const{
    std::string s(Axis::str());
    if(_subset!="")s=_subset+": "+s;
    if(Level==0)return s;
    return s;
}

void AxisTree::print() const{
    std::vector<Axis> ax;
    for(const AxisTree* a=this;a!=0;a=a->nodeNext(this)){
        ax.push_back(*a);
        if(a->subset()!="")ax.back().name+="("+a->subset()+")";
    }

    Axis::print(ax);
}


static std::vector<std::string> factorAxes;
static bool isFactor(const AxisTree* Ax){
    return std::find(factorAxes.begin(),factorAxes.end(),Ax->name)!=factorAxes.end();
}
static bool isComplement(const AxisTree* Ax){return not isFactor(Ax);}

static AxisTree* newSubTree(const AxisTree* Ax,bool (*Select)(const AxisTree*)){
    AxisTree * sub;
    if(Select(Ax)){
        sub=new AxisTree();
        Ax->subTree(Select,sub);
    }
    else {
        int cnt=0;
        for(size_t k=0;k<Ax->childSize();k++){
            cnt+=Select(Ax->child(k));
        }
        if(cnt==1)
            sub=newSubTree(Ax->child(0),Select);
        else {
            sub=new AxisTree(Axis("Hybrid",cnt,0.,double(cnt-1),"automatic",cnt));
            for(size_t k=0;k<Ax->childSize();k++)sub->childAdd(newSubTree(Ax->child(k),Select));
        }
    }
    return sub;
}

AxisTree * AxisTree::factor(std::vector<std::string> Factor) const {
    factorAxes=Factor;
    AxisTree * fac=new AxisTree();
    subTree(isFactor,fac);
    return newSubTree(this,isFactor);
}
AxisTree * AxisTree::complement(std::vector<std::string> Factor) const {
    factorAxes=Factor;
    AxisTree * fac=new AxisTree();
    subTree(isComplement,fac);
    return newSubTree(this,isComplement);
}


