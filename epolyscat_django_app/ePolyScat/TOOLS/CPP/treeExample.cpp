// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "treeExample.h"

using namespace std;

TreeExample::TreeExample(unsigned int Depth, unsigned int Width,unsigned int LevelVal,const TreeExample *Parent){
    // example for a recursive constructor
    // note that we any reference to the upper tree structure at a given recursive level
    // is likely to give incorrect results, as construction of upper levels is incomplete

    parentRef()=Parent;

    if(Depth>8)ABORT("exceeds maximal depth of 8");
    if(Width>4)ABORT("exceeds maximal with of 4");

    if(Depth==0)return;

    for(unsigned int w=0;w<min(Width,Depth);w++){
        childAdd(new TreeExample(Depth-1,Width,10*LevelVal+w+1,this));
    }

    //IMPORTANT: when using tree members in recursive tree constructor
    // make sure tree above has been structurally set up (i.e. when exiting the recursive call)
    if(isLeaf()){
        vector<unsigned int> idx(index());
        val=1;
        for(unsigned int k=1;k<idx.size();k++)val+=10*val+idx[k]+1;
    }
    else val=depth();
}

void TreeExample::test(){
    unsigned int deep=5,wide=3;
    TreeExample tree(deep,wide);
    tree.valIndex();

    std::cout<<"Tree of depth="+tools::str(deep)+", width="+tools::str(wide)+"\n"+tree.str()<<std::endl;

    vector<unsigned int> perm(3,0);
    perm[0]=2;
    perm[1]=1;

    TreeExample * sub=tree.child(2)->descend()->child(1);
    TreeExample* leaf=const_cast<TreeExample*>(sub->firstLeaf());
    while(leaf!=0){
        cout<<" - "<<tools::str(leaf->index())<<" ";
        leaf=leaf->nodeRight(sub);
    }
    cout<<endl;

    TreeExample *copy=tree.deepCopy();
    std::cout<<"\nDeep copy\n"<<copy->str()<<endl;


    TreeExample view;
    tree.permute(perm,view,true);
    TreeExample* deri=view.deepCopy();

    std::cout<<"\nView resorted first "+tools::str(perm.size())+" levels\n"<<view.str()<<endl;
    std::cout<<"\n\nConvert to independent tree\n"<<deri->str()<<endl;

    leaf=tree.descend(perm.size()+1);
    TreeExample *lview=deri->descend(perm.size()+1);
    cout<<"sizes: "<<tree.descend(perm.size())->levelSize()<<" "<<view.descend(perm.size())->levelSize()<<endl;
    cout<<"\npermuted\n";
    for(;lview!=0;lview=lview->nodeRight())cout<<" "<<lview->strNode();
    cout<<endl;
    cout<<"view was deleted\noriginal\n";
    for(;leaf!=0;leaf=leaf->nodeRight())cout<<" "<<leaf->strNode();
    cout<<endl<<flush;

//    TreeExample partView;
//    tree.subTree(nodeLevel,&partView,true);
//    cout<<"partial tree\n"<<partView.str()<<endl;


    cout<<"subtree of itself "<<tree.isSubtree(&tree)<<endl;
    TreeExample * subt=tree.child(1);
    cout<<"is 2nd child subtree "<<tree.child(1)->isSubtree(&tree)<<" "<<subt->isSubtree(&tree)<<endl;

}

void TreeExample::valIndex(){
    vector<unsigned int>idx(index());
    val=0;
    for(unsigned int k=0;k<idx.size();k++)val=10*val+idx[k];
    for(unsigned int k=0;k<childSize();k++)child(k)->valIndex();
    val+=100000*(depth()+1);
}








