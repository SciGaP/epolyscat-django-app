#ifndef MAPGRIDHYBRID_H
#define MAPGRIDHYBRID_H

#include "operatorMap.h"
#include "qtEigenDense.h"

class IndexGridHybrid;

/// generalized map from Hybrid discretization
class MapGridHybrid: public Tree<MapGridHybrid>, public OperatorAbstract {
    std::shared_ptr<const OperatorAbstract> _map;
    Eigen::MatrixXcd _toGrid;
    std::vector<int> _imap;
    std::unique_ptr<Coefficients> _auxC;
public:
    using Tree::str; // there is also OperatorAbstract::str
    MapGridHybrid(const IndexGridHybrid* I, const Index* J);
    void apply(std::complex<double> A, const Coefficients& X, std::complex<double> B, Coefficients& Y) const;
    std::string strNode(int Level=Tree_defaultKind) const ;
};

#endif // MAPGRIDHYBRID_H
