#ifndef INDEXG_H
#define INDEXG_H


#include "discretizationDerived.h"
#include "indexNew.h"
#include "operatorAbstract.h"
#include "tree.h"
#include <memory>

class OperatorFloor;
class OperatorIdentity;

class IndexGrid: public IndexNew{
    friend class DiscretizationGrid;
public:
    /// grid index from Index (obsolescent)
    /// <br> convert selected levels to grid
    IndexGrid(const Index *I, std::vector<unsigned int> Level, std::vector<unsigned int> Point, std::vector<double> Limit);

    /// Rules:
    ///<br> axisName==Axis[k] will be converted to grid given in Grid[k] with quadrature weights Weig[k]
    ///<br> if Grid[k].size()==0, use Weigh[k].size() quadrature points, if Weigh[k].size()==0, use Basis order() points
    ///<br> BasisIntegrable is evaluated, BasisGrid interpolated for Grid points, other bases cannot be converted
    ///<br> on finite element axes, only Grid[k] points falling into given element will be converted (-> piece-wise grid)
    IndexGrid(const Index *I, const std::vector<std::string> &Ax, std::vector<std::vector<double>> Grid, std::vector<std::vector<double>> Weig={});

    /// return Grid's and Weig's for a standard plot
    ///
    /// rules for each k of Axes:
    /// <br> equidistant Points[k] in [a,b]=[Bounds[k][0],Bounds[k][1]]
    /// <br> if a>b, empty Grid[k] and Points[k] zeros in Weig[k] (interpreted as Weig[k].size() quadrature grid in DiscretizationGrid)
    static void gridWeight(std::vector<std::vector<double>> &Grid,std::vector<std::vector<double>> &Weig,
                           std::vector<std::string> Axes,std::vector<unsigned int> Points={},std::vector<std::vector<double>> Bounds={});

};


#endif // INDEXG_H
