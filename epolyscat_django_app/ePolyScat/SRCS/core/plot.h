// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef PLOT_H
#define PLOT_H

#include <memory>

#include "tools.h"
#include "useMatrix.h"
#include "plotCoefficients.h"
#include "plotKind.h"
#include "discretization.h"

class DiscretizationGrid;
class ReadInput;
class Coefficients;
class Index;
class OperatorAbstract;
class OperatorTree;


/** @defgroup Plot Plotting and analysis
 *  \brief plots (1d, 2d), spectral analysis, harmonics, Floquet etc.
 *  @{
 */

/// \brief evaluate densities  @f$\Psi^*@f$(args)@f$\Psi@f$(args) and optionally @f$\Psi^*@f$(args) (Op[k] @f$\Psi@f$)(args) to (ascii) file, gnuplot style
///
/// argument ("axis") ranges can be selected, arguments can be summed/integrated over, separate plots for some arguments fixed can be made
class Plot : public PlotCoefficients {
    friend class SpectrumPz;
    friend class SpectrumPlot;

    bool _realValues;
    const Index * discIndex;              ///< Index of Coefficients to be plotted
    std::shared_ptr<const Index> gIndex;  ///< Index of plot grid
    std::shared_ptr<const OperatorAbstract> gMap;   ///< Map discIndex->gIndex
    std::vector<const OperatorAbstract*> _ops;      ///< list of operators (=idenity if not specified)
    std::shared_ptr<const OperatorTree> _densityOp; ///< quadrature for axes that are integrated over

    std::string columnHeader; ///< column titles
    std::vector<std::vector<double> > axCols; ///< coordinate axes (in full length)
    mutable std::vector<std::vector<double> > _cols; ///< current columns, can be written to file
    unsigned int nDigits;     //!< number of digits in storage
    bool _append;             ///< append to exisiting _cols
    bool _onlyMatching;       ///< plot only if axes match

    std::shared_ptr<Index> iFull; ///< index for fullView
    std::shared_ptr<Index> iPlot; ///< index for plotView
    std::shared_ptr<Coefficients> gvec,lvec;   ///< input converted to gIndex grid (and left-hand)
    Coefficients* rightFullView,*leftFullView; ///< view of gvec/lvec with floor lower (MUST NOT BE A shared_ptr????)
    std::shared_ptr<Coefficients> rightView,leftView; ///< permutation of fullView in the sequence for plotting
    std::shared_ptr<Index> iStor;
    Coefficients *plotStore;  ///< modulus squared and summed over sum-coordinates
    unsigned int plotDepth;   ///< levels of plotStore from this onward will to go into single column

    double riemannSumWeight; //HACK to acount for summation over equidistant grids

    /// put into storage, sum over non-printing levels (returns sum of all squares)
    double sumLeftConjRight(Coefficients *LVec, Coefficients *RVec, Coefficients & Store) const;
    void generateAxes(const Index* Idx,  std::vector<std::vector<double>  > &Cols, unsigned int IAx=0) const;
    void intoColumns(const Coefficients &Store,  std::vector<std::vector<double>  > &Cols) const;

    void construct(std::vector<std::string> Axis, std::vector<std::string> Use,
                   std::vector<std::vector<double>> Grid, std::vector<std::vector<double>> Weig);
    /// true if Idx can be used for plotting w/o further conversion
    bool plotable(const Index* Idx,std::vector<std::string> Axes, std::vector<unsigned int> Points);

    const Plot & generate(const Coefficients & C, const OperatorAbstract *Op=0) const;
    static void plotJAD(const Coefficients * C, double E1, double E2, std::string fileName, unsigned int thetaNum=100);
    void append(bool Append);///< {_append=Append;} ///< tagged data will be written into 2d plot
    /// convert to plot, scale and add to columns (create new, if no columns exist)
    std::vector<std::vector<double> > sum(const Coefficients & C, std::string &File, double Scale=1., const std::string &Tag="") const;
    const Coefficients & plotCoefs() const {return *plotStore;}

public:
    virtual ~Plot();
    Plot();
    Plot(const Index *Idx, ReadInput & inp, bool WarnIfEmpty=false); //!< get plot definition from input
    Plot(const Index *Idx              /**< Index to plot */,
         std::vector<std::string> Axis /**< Axis to plot (all other will be summed over) */,
         std::vector<std::string> Use  /**< how to use axis "g"...grid,"p"...separate plots, "s"...sum*/ ,
         std::vector<unsigned int> GridPoints={} /**< number of points in equidistant grid (default: 0 = automatic)*/,
         std::vector<std::vector<double> > GridBounds={} /**< lower [0] and upper [1] grid boundaries (default=all)*/
            );
    Plot(const Index *Idx, const PlotKind *Kind /** plot definition, contains Axis, Use, GridPoints, GridBounds, see description of main constructor */);

    Plot& withPlotReal(){_realValues=true;return *this;}; ///< qualify contructor: plot (sum of) real values rather than (sum of) squares

    /// transform to plot and write to file
    void plot(const Coefficients & C, const std::string & File,
              const std::vector<std::string> & Head=std::vector<std::string>(0),
              const std::string Tag="",
              bool OverWrite=true) const;
    /// write current columns to file
    void write(const std::string & File,
               const std::vector<std::string> &Head={},
               bool OverWrite=true) const;
    void clear(){_cols.clear();} ///< clear any previous plot
    double integral() const; /// integrate density

    std::string briefName() const {return "wf";}

    void addOperator(const OperatorAbstract * Op); ///< plot density conj(Psi(coor)) [Op Psi](coor)
    unsigned int dimension() const {return axCols.size();} //!< number of plot axes
    void digits(unsigned int NDigits){nDigits=NDigits;} //!< number of digits for ASCII output
    std::string str() const; ///< plot info
    void print() const; ///< neatly print plot definition
    bool isEmpty() const {return plotStore==0;}


    bool isAppend() const {return _append;} ///< true if plot in append-mode
};

/** @} */
#endif // PLOT_H
