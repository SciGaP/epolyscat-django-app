// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#include "permuteOperatorTree.h"
#include "tools.h"
#include "str.h"
#include "timer.h"
#include "printOutput.h"

#include "operatorTreeExplorer.h"


PermuteOperatorTree::PermuteOperatorTree(
		const OperatorTree* _optree): optree(_optree), iIndex(_optree->iIndex), jIndex(_optree->jIndex), iPermutedIndex(0), jPermutedIndex(0), permutedOptree(0){

    if(iIndex->parent()!=0 or jIndex->parent()!=0) ABORT("Only root optrees are supported");

    for(const Index* tmp = iIndex; tmp->childSize()!=0; tmp=tmp->child(0)) iPermutation.push_back(tmp->depth());
    for(const Index* tmp = jIndex; tmp->childSize()!=0; tmp=tmp->child(0)) jPermutation.push_back(tmp->depth());

    permutedFloorLevel = iIndex->firstFloor()->depth();
}



/*
 * Convenience methods
 */

void PermuteOperatorTree::move(std::vector<unsigned int>& perm, unsigned int from, unsigned int to){
    if(from>perm.size()-1 or to>perm.size()-1) ABORT("Index out of bounds");

    unsigned int moving = perm[from];
    perm.erase(perm.begin()+from);
    perm.insert(perm.begin()+to, moving);
}

void PermuteOperatorTree::iMove(unsigned int from, unsigned int to){
    if(permutedOptree!=0){
        iPermutedIndex=0;
        jPermutedIndex=0;
        permutedOptree=0;
    }

    move(iPermutation, from, to);
}

void PermuteOperatorTree::jMove(unsigned int from, unsigned int to){
    if(permutedOptree!=0){
        iPermutedIndex=0;
        jPermutedIndex=0;
        permutedOptree=0;
    }

    move(jPermutation, from, to);
}

void PermuteOperatorTree::move(unsigned int from, unsigned int to){
    if(permutedOptree!=0){
        iPermutedIndex=0;
        jPermutedIndex=0;
        permutedOptree=0;
    }

    move(iPermutation, from, to);
    move(jPermutation, from, to);
}



/*
 * Helpers
 */

std::vector<unsigned int> PermuteOperatorTree::invertPermutation(std::vector<unsigned int> permutation){
    std::vector<unsigned int> res;
    for(unsigned int i=0; i<permutation.size(); i++){
        res.push_back(i);
    }

    for(unsigned int i=0; i<permutation.size(); i++){
        unsigned int index=0;
        for(;index<permutation.size();index++){
            if(permutation[index]==i){
                std::swap(permutation[i], permutation[index]);
                std::swap(res[i], res[index]);
                break;
            }
        }
    }

    return res;
}

// Copied from Coefficients::examplePermute
void PermuteOperatorTree::permuteCoeffs(
		const std::vector<unsigned int> permutation, const Index* srcIndex, const Coefficients& src, const Index* targetIndex, Coefficients& target){

    if(src.idx()!=srcIndex) ABORT("src indices don't match");
    if(target.idx()!=targetIndex) ABORT("target indices don't match");

    /*
     * Create zero floor depth indices
     */
    Index srcIndex0(*srcIndex);
    Index targetIndex0(*targetIndex);

    srcIndex0.resetFloor(srcIndex0.firstLeaf()->depth());
    targetIndex0.resetFloor(targetIndex0.firstLeaf()->depth());


    /*
     * Create zero floor depth permuted view on source coefficients
     */
    Coefficients vSrc0(&srcIndex0,const_cast<Coefficients*>(&src));
    Coefficients vSrcPermuted0;
    vSrc0.permute(permutation,vSrcPermuted0);

    /*
     * Crate zero floor depth view on target coefficients
     */
    Coefficients vTarget0(&targetIndex0,&target);

    /*
     * Copy data
     */
    vTarget0 = vSrcPermuted0;
}

/*
 * Classify the index hierarchy along a path from root to floor in optree:
 *
 * 0 = no descending
 * 1 = lhs descending
 * 2 = rhs descending
 * 3 = both sides descending
 */
std::vector<unsigned int> PermuteOperatorTree::classifyStructureForLeaf(const OperatorTree* leaf){
	std::vector<unsigned int> res;
	std::vector<unsigned int> path(leaf->index());
	const OperatorTree* optree = leaf->root();
	for(unsigned int i=0; i<path.size(); i++){
		OperatorTree* child = optree->child(path[i]);
		if     (child->iIndex!=optree->iIndex and child->jIndex==optree->jIndex) res.push_back(1);
		else if(child->iIndex==optree->iIndex and child->jIndex!=optree->jIndex) res.push_back(2);
		else if(child->iIndex!=optree->iIndex and child->jIndex!=optree->jIndex) res.push_back(3);
		else                                                                     res.push_back(0);
		optree=child;	
	}

	return res;
}

std::vector<unsigned int> PermuteOperatorTree::classifyStructure(const OperatorTree* optree){
	const OperatorTree* leaf = optree->firstLeaf();
	std::vector<unsigned int> structure=classifyStructureForLeaf(leaf);
	for(;leaf!=0;leaf=leaf->nextLeaf()){
		if(classifyStructureForLeaf(leaf)!=structure) ABORT("Structure is not equivalent along different branches");
	}
	return structure;
}


/*
 * Timer for matrix calls when shortening the optree
 * This is the main performance leak of both implementations
 */
TIMER(alignToIndex_matrix,)

#ifdef _PERMUTE_OPERATOR_TREE_NEWCODE_
/*
 ********************************* NEW CODE ************************************
 */

OperatorTree* PermuteOperatorTree::buildOptreeFromMatrix(
		const UseMatrix& mat, std::string name, const Index* iIndex, const Index* jIndex, const Index* iIndexRoot, const Index* jIndexRoot){
    OperatorTree* res = new OperatorTree(name, "TEMPORARY", iIndex, jIndex, 0);

	if(not iIndex->hasFloor()){
		for(unsigned int i=0; i<iIndex->childSize(); i++){
			OperatorTree* child = buildOptreeFromMatrix(mat, name, iIndex->child(i), jIndex, iIndexRoot, jIndexRoot);
			//Don't check if child is zero to preserve identical ordering among branches
			res->childAdd(child);
			
		}
	}else if(not jIndex->hasFloor()){
		for(unsigned int i=0; i<jIndex->childSize(); i++){
			OperatorTree* child = buildOptreeFromMatrix(mat, name, iIndex, jIndex->child(i), iIndexRoot, jIndexRoot);
			//Don't check if child is zero to preserve identical ordering among branches
			res->childAdd(child);
			
		}
	}else{
        if(iIndex->sizeStored()==1 && jIndex->sizeStored()==1){
            // Faster than using the factory
            res->oFloor = new OperatorFloorSingle(mat(iIndex->posIndex(iIndexRoot),jIndex->posIndex(jIndexRoot)));
        }else{
            UseMatrix fMat=mat.block(iIndex->posIndex(iIndexRoot), jIndex->posIndex(jIndexRoot), iIndex->sizeStored(), jIndex->sizeStored());
            res->oFloor=OperatorFloor::factory(std::vector<const UseMatrix*>(1,&fMat), name);
        }
    }

    return res;
}


OperatorTree* PermuteOperatorTree::alignToIndex(const OperatorTree *optree, const Index *iIndex, const Index *jIndex){

	//Check if floor levels match
	if((optree->childSize()==0) != (iIndex->hasFloor() and jIndex->hasFloor())){
		UseMatrix mat;
	
		if(optree->childSize()!=0){
			STARTDEBUG(alignToIndex_matrix)
			optree->matrix(mat);
			STOPDEBUG(alignToIndex_matrix)
		}else{
			optree->matrix(mat);
		}


		return buildOptreeFromMatrix(mat,optree->name,iIndex,jIndex,iIndex,jIndex);
	}

	OperatorTree* res = new OperatorTree(optree->name, optree->getDefinition(), iIndex, jIndex, 0);

	//Detect whether iIndex or jIndex or both or none are descended, while forcing that all children behave the same way
	int iIndexChanges = -1;
	int jIndexChanges = -1;
	for(unsigned int i=0; i<optree->childSize(); i++){
        int iTmp = optree->child(i)->iIndex!=optree->iIndex;
		int jTmp = optree->child(i)->jIndex!=optree->jIndex;

		if(iIndexChanges==-1){
			iIndexChanges=iTmp;
		}else{
			if(iIndexChanges!=iTmp) ABORT("lhs index mismatch");
		}

		if(jIndexChanges==-1){
			jIndexChanges=jTmp;
		}else{
			if(jIndexChanges!=jTmp) ABORT("rhs index mismatch");
		}
	}

	if(iIndexChanges==-1) iIndexChanges=0;
	if(jIndexChanges==-1) jIndexChanges=0;

	//Both indices are descended
	if(iIndexChanges && jIndexChanges){
		std::map<const Index*, OperatorTree*> iChildren;
		for(unsigned int i=0; i<optree->childSize(); i++){
			OperatorTree* iChild;
			if(iChildren.find(optree->child(i)->iIndex)==iChildren.end()){
				iChild = new OperatorTree(optree->name, optree->getDefinition(), iIndex->child(optree->child(i)->iIndex->nSibling()), jIndex, 0);
				iChildren[optree->child(i)->iIndex] = iChild;
			}else{
				iChild = iChildren[optree->child(i)->iIndex];
			}

			OperatorTree* child = alignToIndex(
					optree->child(i),
					iIndex->child(optree->child(i)->iIndex->nSibling()),
					jIndex->child(optree->child(i)->jIndex->nSibling())
			);

			//Don't check if child is zero to preserve identical ordering among branches
			iChild->childAdd(child);
			
		}

		//This ensures the children are ordered in ascending iIndex order
		for(unsigned int i=0; i<optree->iIndex->childSize(); i++){
			std::map<const Index*, OperatorTree*>::iterator it = iChildren.find(optree->iIndex->child(i));
			if(it!=iChildren.end()){
				OperatorTree* child = it->second;
				//Don't check if child is zero to preserve identical ordering among branches
				res->childAdd(child);
				
			}
		}
	}

	//No index or one index is descended
	else{
		for(unsigned int i=0; i<optree->childSize(); i++){
			const Index* iIdx = iIndex;
			const Index* jIdx = jIndex;

			if(optree->child(i)->iIndex!=optree->iIndex) iIdx=iIndex->child(optree->child(i)->iIndex->nSibling());
			if(optree->child(i)->jIndex!=optree->jIndex) jIdx=jIndex->child(optree->child(i)->jIndex->nSibling());

			OperatorTree* child = alignToIndex(optree->child(i), iIdx, jIdx);
			//Don't check if child is zero to preserve identical ordering among branches
			res->childAdd(child);
		}
	}

    return res;
}

void PermuteOperatorTree::insertFloorLevel(OperatorTree* optree){
	if(optree->oFloor!=0){
		OperatorTree* child = new OperatorTree(optree->name, optree->getDefinition(), optree->iIndex, optree->jIndex, 0);
		child->oFloor = optree->oFloor;
		optree->oFloor=0;
		optree->childAdd(child);
	}else{
		for(unsigned int i=0; i<optree->childSize(); i++) insertFloorLevel(optree->child(i));
	}
}


#else
/*
 ************************************ OLD CODE *********************************
 */

OperatorTree* PermuteOperatorTree::buildOptreeFromMatrix(
		const UseMatrix& mat, std::string name, const Index* iIndex, const Index* jIndex, const Index* iIndexRoot, const Index* jIndexRoot, std::complex<double>* timeDepFac){

	if(iIndex->hasFloor()!=jIndex->hasFloor()) ABORT("floor levels don't match");

    OperatorTree* res = new OperatorTree(name, "DUMMY-TODO", iIndex, jIndex, 0);

    if(not iIndex->hasFloor()){
        for(unsigned int i=0; i<iIndex->childSize(); i++){
            for(unsigned int j=0; j<jIndex->childSize(); j++){
                OperatorTree* child = buildOptreeFromMatrix(mat, name, iIndex->child(i), jIndex->child(j), iIndexRoot, jIndexRoot, timeDepFac);
                if(child->isZero()){
                    delete child;
                }else{
                    res->childAdd(child);
                }
            }
        }
    }else{
        if(iIndex->sizeStored()==1 && jIndex->sizeStored()==1){
            // Faster than using the factory
            res->oFloor = new OperatorFloorSingle(mat(iIndex->posIndex(iIndexRoot),jIndex->posIndex(jIndexRoot)));
        }else{
            UseMatrix fMat=mat.block(iIndex->posIndex(iIndexRoot), jIndex->posIndex(jIndexRoot), iIndex->sizeStored(), jIndex->sizeStored());
            res->oFloor=OperatorFloor::factory(std::vector<const UseMatrix*>(1,&fMat), name);
        }
        res->oFloor->setFactor(timeDepFac);
    }

    return res;
}

void PermuteOperatorTree::findIndexChildren(OperatorTree* optree, const Index* iIndex, const Index* jIndex, std::vector<OperatorTree*>& result){
	if(optree->iIndex->parent()==iIndex and optree->jIndex->parent()==jIndex){
        // Find deepest children (Goal: preserve sums)
        // Possibly no longer necessary due to findFloorChildren?
        std::vector<OperatorTree*> tmp;
        for(unsigned int i=0; i<optree->childSize(); i++){
            findIndexChildren(optree->child(i), iIndex, jIndex, tmp);
        }

        if(tmp.size()==0){
            result.push_back(optree);
        }else{
            for(unsigned int i=0; i<tmp.size(); i++) result.push_back(tmp[i]);
	    }
    }else{
		for(unsigned int i=0; i<optree->childSize(); i++){
			findIndexChildren(optree->child(i), iIndex, jIndex, result);
		}
	}
}

void PermuteOperatorTree::findFloorChildren(OperatorTree* optree, const Index* iIndex, const Index* jIndex, std::vector<OperatorTree*>& result){
    if((iIndex==0 or optree->iIndex==iIndex) and (jIndex==0 or optree->jIndex==jIndex) and optree->oFloor!=0){
        result.push_back(optree);
    }

    for(unsigned int i=0; i<optree->childSize(); i++){
        findFloorChildren(optree->child(i), iIndex, jIndex, result);
    }
}


OperatorTree* PermuteOperatorTree::alignToIndex(const OperatorTree *optree, const Index *iIndex, const Index *jIndex){
	std::vector<OperatorTree*> children;
	findIndexChildren(const_cast<OperatorTree*>(optree), optree->iIndex, optree->jIndex, children);

    OperatorTree* res = new OperatorTree(optree->name, optree->def(), iIndex, jIndex, 0);

    //Shorten optree depth
    if(iIndex->hasFloor() and jIndex->hasFloor() and children.size()!=0){	
        UseMatrix mat;
        optree->matrix(mat); 
        res->oFloor=OperatorFloor::factory(std::vector<const UseMatrix*>(1,&mat), optree->name);

        // TODO: Take care of timeDepFac. Possibly ABORT if mismatch
    }

    //Go deeper
    else if(not (iIndex->hasFloor() and jIndex->hasFloor()) and children.size()==0){
        std::vector<OperatorTree*> floors;
        findFloorChildren(const_cast<OperatorTree*>(optree), optree->iIndex, optree->jIndex, floors);
 
        for(unsigned int f=0; f<floors.size(); f++){
            UseMatrix mat;
            floors[f]->matrix(mat);

            for(unsigned int i=0; i<iIndex->childSize(); i++){
                for(unsigned int j=0; j<jIndex->childSize(); j++){


                    OperatorTree* child = buildOptreeFromMatrix(mat, optree->name, iIndex->child(i), jIndex->child(j), iIndex, jIndex, floors[f]->oFloor->factor());
                    if(child->isZero()){
                        delete child;
                    }else{
                        res->childAdd(child);
                    }
                }
            }
        }
    }

    //Continue traversal
    else{
		for(unsigned int i=0; i<children.size(); i++){
            OperatorTree* child = alignToIndex(children[i], iIndex->child(children[i]->iIndex->nSibling()), jIndex->child(children[i]->jIndex->nSibling()));
			if(child->isZero()){
				delete child;
			}else{
				res->childAdd(child);
			}
        }
    }

    return res;
}


void PermuteOperatorTree::placeInPermutedOptree0(OperatorTree* permutedOptree0, const OperatorTree* leaf, OperatorFloor* floor){
    const Index* iIndex0 = leaf->iIndex;
    const Index* jIndex0 = leaf->jIndex;

    std::vector<unsigned int> iPath(iIndex0->index());
	std::vector<unsigned int> jPath(jIndex0->index());
    

    if(iPath.size()!=jPath.size() or iPermutation.size()!=iPath.size()) ABORT("Incompatible indices");

    // First find sum indices
    std::vector<unsigned int> sumIndices;

    const OperatorTree* tmp = leaf;
    const OperatorTree* parent = leaf->parent();

    while(parent!=0){
        unsigned int sumIndex=0;
        for(unsigned int i=0; i<tmp->nSibling(); i++){
            if(parent->child(i)->iIndex==tmp->iIndex and parent->child(i)->jIndex==tmp->jIndex) sumIndex++;
        }

        sumIndices.insert(sumIndices.begin(), sumIndex);

        tmp=parent;
        parent=tmp->parent();
    }    

    // Walk down the operator tree
    const Index* iIdx = permutedOptree0->iIndex;
    const Index* jIdx = permutedOptree0->jIndex;
    OperatorTree* optree = permutedOptree0;
    for(unsigned int i=0; i<iPath.size(); i++){
        iIdx = iIdx->child(iPath[iPermutation[i]]);
        jIdx = jIdx->child(jPath[jPermutation[i]]);
        unsigned int sumIndex = sumIndices[iPermutation[i]]; //Simply use iPermutation

        bool success = false;
        unsigned int sumCounter=0;
        for(unsigned int k=0; k<optree->childSize(); k++){
            if(optree->child(k)->iIndex==iIdx and optree->child(k)->jIndex==jIdx){
                if(sumCounter==sumIndex){
                    optree=optree->child(k);
                    success=true;
                    break;
                }
                sumCounter++;
            }
        }
        if(!success){
            OperatorTree* child=0;
            for(unsigned int i=sumCounter; i<sumIndex+1; i++){
                child=new OperatorTree(optree->name, permutedOptree0->def(), iIdx, jIdx, 0);
                optree->childAdd(child);
            }
            optree=child;
        }

    }

    // Place operator floor
	if(optree->oFloor!=0) ABORT("Sum handling failed");
    optree->oFloor=floor;
}

/*
 *******************************************************************************
 */
#endif




TIMER(permute,)
TIMER(createOptree0,)
TIMER(permuteOptree0,)
TIMER(resetFloorLevel,)

void PermuteOperatorTree::permute(const Index* IPermutedIndex, const Index* JPermutedIndex){
    STARTDEBUG(permute)
	PrintOutput::message("STARTING PERMUTATION...");

	/*
	 * Preflight check: Ensure the structure of the optree is equivalent along different paths
	 * TODO: In order for new code to work, the definition of equivalent needs to be extended towards
	 *		 "Same indexed children among different siblings with same ordering"
	 */
	classifyStructure(optree);

    // TODO This is only a temporary solution, for OperatorTree::OperatorTree(Operator*) [i. e. i->i->i->...->j->j->j->...] aligned optrees
    UseMatrix fullMatrix;
    optree->matrix(fullMatrix);
    optree=buildOptreeFromMatrix(fullMatrix, optree->name, optree->iIndex, optree->jIndex, optree->iIndex, optree->jIndex, 0);

	/*
	 * First: Create necessary index structures
	 */

	PrintOutput::message("  SETTING UP INDICES");

    // TODO: invOvr get's lost!

	// Create zero floor depth indices
    Index iIndex0(*optree->iIndex);
    Index jIndex0(*optree->jIndex);
    iIndex0.resetFloor(iIndex0.firstLeaf()->depth());
    jIndex0.resetFloor(jIndex0.firstLeaf()->depth());
    iIndex0.sizeCompute();
    jIndex0.sizeCompute();

	// Create permuted zero floor depth indices
    Index iPermutedIndex0;
    Index jPermutedIndex0;
    iIndex0.permute(iPermutation, iPermutedIndex0);
    jIndex0.permute(jPermutation, jPermutedIndex0);
    iPermutedIndex0.sizeCompute();
    jPermutedIndex0.sizeCompute();

	// Create actual floor depth permuted indices
    const Index* iPermutedIndex; 
    const Index* jPermutedIndex;
    
    if(IPermutedIndex!=0){
        iPermutedIndex = IPermutedIndex;
    }else{
        Index* tmp = new Index(iPermutedIndex0);
        tmp->resetFloor(permutedFloorLevel);
        tmp->sizeCompute();
        iPermutedIndex = tmp;
    }
   
    if(JPermutedIndex!=0){
        jPermutedIndex = JPermutedIndex;
    }else if(iIndex==jIndex and iPermutation==jPermutation){
        jPermutedIndex = iPermutedIndex;
    }else{
        Index* tmp = new Index(jPermutedIndex0);
        tmp->resetFloor(permutedFloorLevel);
        tmp->sizeCompute();
        jPermutedIndex=tmp;
    }

	// Store them
    this->iPermutedIndex=iPermutedIndex;
    this->jPermutedIndex=jPermutedIndex;

	/*
	 * Second: Permute optree
	 */

#ifdef _PERMUTE_OPERATOR_TREE_NEWCODE_
    // Create zero floor depth unpermuted optree
	PrintOutput::message("  CREATING ZERO FLOOR DEPTH UNPERMUTED OPTREE");
	
	STARTDEBUG(createOptree0)
	OperatorTree* optree0 = alignToIndex(optree, &iIndex0, &jIndex0);
	STOPDEBUG(createOptree0)

	// TODO: This includes a check of whether alignToIndex has worked. Can be skipped using classifyStructureForLeaf
	// Classify the structure of the optree, 0=no descending, 1=lhs descending, 2=rhs descending
	std::vector<unsigned int> structure=classifyStructure(optree0);

	// TODO: Additional check that can be dropped in the future
	for(unsigned int i=0; i<structure.size(); i++){ if(structure[i]==3) ABORT("alignToIndex didn't work"); }

	// Find the permutation to permute optree0 with
	std::vector<unsigned int> permutation;
	std::vector<unsigned int> iIndices, jIndices;
	for(unsigned int i=0; i<structure.size(); i++){
		permutation.push_back(i);
		if(structure[i]==1) iIndices.push_back(i);
		else if(structure[i]==2) jIndices.push_back(i);
	}

	for(unsigned int i=0; i<iPermutation.size(); i++){
		permutation[iIndices[i]]=iIndices[iPermutation[i]];
	}

	for(unsigned int i=0; i<jPermutation.size(); i++){
		permutation[jIndices[i]]=jIndices[jPermutation[i]];
	}

	//Ensure operator floors are not permuted to the middle of the tree
	if(permutation.back()!=permutation.size()-1){
		insertFloorLevel(optree0);
		structure.push_back(0);
		permutation.push_back(permutation.size());
	}

	// Create zero floor depth permuted optree
	PrintOutput::message("  PERMUTING ZERO FLOOR DEPTH OPTREE");
    
	STARTDEBUG(permuteOptree0)
	OperatorTree permutedOptree0;
	optree0->permute(permutation, permutedOptree0);

	//TODO: Check why this is necessary, should be done in permute
	permutedOptree0.purge(permutation.size());


	//Merge indices
	PrintOutput::message("  MERGING PERMUTED OPTREE WITH PERMUTED INDICES");
	
	permutedOptree0.iIndex = &iPermutedIndex0;
	permutedOptree0.jIndex = &jPermutedIndex0;

	for(const OperatorTree* leaf=permutedOptree0.firstLeaf(); leaf!=0; leaf=leaf->nextLeaf()){
		std::vector<unsigned int> iPath(leaf->iIndex->index());
		std::vector<unsigned int> jPath(leaf->jIndex->index());
		std::vector<unsigned int> path(leaf->index());

		OperatorTree* optree = &permutedOptree0;
		const Index* iIdx = &iPermutedIndex0;
		const Index* jIdx = &jPermutedIndex0;

		int iCounter=0;
		int jCounter=0;
		for(unsigned int i=0; i<path.size(); i++){

			optree = optree->child(path[i]);

			if(structure[i]==1){
				iIdx = iIdx->child(iPath[iPermutation[iCounter]]);
				iCounter++;
			}else if(structure[i]==2){
				jIdx = jIdx->child(jPath[jPermutation[jCounter]]);
				jCounter++;
			}

			optree->iIndex = iIdx;
			optree->jIndex = jIdx;
		}
	}

	STOPDEBUG(permuteOptree0)

	// Reinstate floors
	PrintOutput::message("  RESETTING FLOOR LEVEL");
    
	STARTDEBUG(resetFloorLevel)
	permutedOptree = alignToIndex(&permutedOptree0, iPermutedIndex, jPermutedIndex);
	STOPDEBUG(resetFloorLevel)

    // Clean up
    //delete optree0; //TODO: causes segfault...

#else
	// Create zero floor depth unpermuted optree
	PrintOutput::message("  CREATING ZERO FLOOR DEPTH UNPERMUTED OPTREE");
	

	STARTDEBUG(createOptree0)
	OperatorTree* optree0 = alignToIndex(optree, &iIndex0, &jIndex0);
	STOPDEBUG(createOptree0)


	// Create zero floor depth permuted optree
	PrintOutput::message("  PERMUTING ZERO FLOOR DEPTH OPTREE");
	
	STARTDEBUG(permuteOptree0)
	OperatorTree* permutedOptree0 = new OperatorTree(optree->name, optree->def(), &iPermutedIndex0, &jPermutedIndex0, 0);
    for(OperatorTree* leaf=optree0->firstLeaf(); leaf!=0; leaf=leaf->nextLeaf()){
        OperatorFloor* floor=leaf->oFloor;
        leaf->oFloor=0;
        if(floor!=0) placeInPermutedOptree0(permutedOptree0, leaf, floor);
    }
	STOPDEBUG(permuteOptree0)

	// Reinstate floors
	PrintOutput::message("  RESETTING FLOOR LEVEL");

	STARTDEBUG(resetFloorLevel)
    permutedOptree = alignToIndex(permutedOptree0, iPermutedIndex, jPermutedIndex);
	STOPDEBUG(resetFloorLevel)

	// Clean up
	delete optree0;
	delete permutedOptree0;

#endif
	//Perform a check, should be cheap enough to do in any case
//	check();

	STOPDEBUG(permute)
    PrintOutput::message("...DONE");
}

void PermuteOperatorTree::check(){
    Coefficients c(jIndex);
	c.setToRandom();

    Coefficients cPermuted(jPermutedIndex);
    jPermuteCoeffs(c, cPermuted);
    Coefficients cTarget(iIndex);
    Coefficients cTargetPermuted(iPermutedIndex);

    optree->apply(1.,c,0.,cTarget);
	permutedOptree->apply(1.,cPermuted,0.,cTargetPermuted);

    Coefficients cTargetUnpermuted(iIndex);
    iUnpermuteCoeffs(cTargetPermuted, cTargetUnpermuted);

    std::vector<std::complex<double>* > p1, p2;
    cTarget.pointerToC(p1);
    cTargetUnpermuted.pointerToC(p2);
	double maxError=0.;
    for(unsigned int i=0; i<p1.size(); i++){
        if(std::abs(*p1[i]-*p2[i])>maxError) maxError = std::abs(*p1[i]-*p2[i]);
    }
	
	if(maxError>1.e-10) ABORT("WARNING! Test failed: "+std::to_string(maxError));

}

void PermuteOperatorTree::test(){
    // TODO
}
