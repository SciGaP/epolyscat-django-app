#ifndef PROBLEMTEMPLATES_H
#define PROBLEMTEMPLATES_H

#include "problemNode.h"
#include "problemName.h"
#include "problemObservable.h"

/// convenience templates for ProblemTree

/// identify Derived class
template<class P>
std::shared_ptr<P> kind(std::vector<std::shared_ptr<ProblemNode>> Path){
    for(auto p: Path){
        std::shared_ptr<P> r(std::dynamic_pointer_cast<P>(p));
        if(r)return r;
    }
    return std::shared_ptr<P>();
}

template<class P>
std::shared_ptr<P> kind(std::string Definition){
    std::shared_ptr<P> res((P(Definition)));
    if(res->definition()==Definition)return res;
    return std::shared_ptr<P>();
}

template <class P>
std::shared_ptr<P> makeNode(std::string Definition){
    std::shared_ptr<P> res(new P(Definition));
    if(res->definition()==Definition)return res;
    return std::shared_ptr<P>();
}

#endif // PROBLEMTEMPLATES_H
