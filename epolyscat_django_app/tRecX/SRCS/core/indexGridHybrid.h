#ifndef INDEXGRIDHYBRID_H
#define INDEXGRIDHYBRID_H

#include "indexDerived.h"
#include "operatorMap.h"

class DiscretizationGrid;
class IndexGridHybrid;

/// product grid that relates to several Parent indicies and contains the maps from the parents
class IndexGridHybrid: public Index
{
    const Index* _gdx;
    void extendShare(const Index* Idx, const std::vector<std::string> &Ax={},const std::vector<std::vector<double>> & Grid={});
public:
    IndexGridHybrid(){}
    IndexGridHybrid(const Index* Idx, const std::vector<std::string> &Ax, const std::vector<std::vector<double> > &Grid);
    const Index* gridIdx() const {return _gdx;}
    const IndexGridHybrid* childHybrid(size_t k) const {return dynamic_cast<const IndexGridHybrid*>(child(k));}
    // recover conversion axes and grid from IndexGridHybrid
    void axesAndGrids(std::vector<std::string> &Ax, std::vector<std::vector<double>> &Grid) const;
    void imap(const Index* Gdx, std::vector<int> &I, std::vector<double> Ipt={}, std::vector<double> Gpt={}) const;
    std::string strNode(int Level=Tree_defaultKind) const {return Index::strNode(Level);}
};

#endif // INDEXSHAREDGRID_H
