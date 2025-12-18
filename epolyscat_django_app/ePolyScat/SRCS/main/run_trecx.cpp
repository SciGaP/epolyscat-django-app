// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
// End of license


#include "tools.h"
#include "eigenTools.h"
#include "timer.h"
#include "readInput.h"
#include "tRecX_cache.h"
#include "printOutput.h"
#include "mpiWrapper.h"
#include "debugInfo.h"

#include "basisNdim.h"
#include "basicDisc.h"
#include "discretizationHybrid.h"
#include "discretizationSpectral.h"
#include "discretizationSpectralProduct.h"
#include "discretizationtsurffspectra.h"
#include "discretizationSurface.h"
#include "discretizationConstrained.h"

#include "pulse.h"
#include "harmonics.h"
#include "plot.h"
#include "plotSpectral.h"
#include "makePlot.h"
#include "scanEigenvalue.h"
#include "eigenSubspace.h"
#include "floquetAnalysis.h"
#include "finDiff.h"
#include "operatorDefinition.h"
#include "operatorDefinition.h"
#include "operatorFactor.h"
#include "operatorFloor3d.h"
#include "operatorFloorEE.h"
#include "parallelOperator.h"
#include "parallelLayout.h"

#include "timePropagator.h"
#include "timePropagatorOutput.h"

#include "tsurffSource.h"
#include "derivativeFlatInhomogeneous.h"
#include "derivativeLocal.h"
#include "operatormapchannelssurface.h"
#include "operatorGradient.h"
#include "operatorIdentity.h"
#include "autoconverge.h"
#include "spectrumISurff.h"

// #ifdef _USE_HACC_
// //#include "read_columbus_data.h"
// #include "discretizationHaCC.h"
// #include "mo.h"
// #include "ci.h"
// #endif

#include "odeFactory.h"

#include "tRecXchecks.h"
#include "channelsSubregion.h"

#include "densityOfStates.h"
#include "resolvent.h"

#include "initialState.h"

#include "coefficientsWriter.h"
#include "log.h"

#include "multiplyGrid.h"
#include "volkovGrid.h"
//#include "basisMat.h"
//#include "molecularorbital.h"
#include "basisOrbital.h"
#include "basisOrbitalNumerical.h"
#include "productFunction.h"

#include "toolsPrint.h"
#include "spectrumPlot.h"
#include "spectrumLinCom.h"

#ifdef _USE_UI_
#include "uiLayout.h"
#endif

#include "eigenTrace.h"

// develop only
#include "surfaceFlux.h"
#include "indexConstraint.h"
#include "indexNew.h"
#include "multipolePotential.h"
#include "randomPotential.h"
#include "timeCritical.h"
#include "farm.h"

using namespace std;
using namespace tools;

TIMER(all, )
TIMER(run, )
TIMER(setup, )
TIMER(operator, )
TIMER(matrix, )
TIMER(setupDisc, )
TIMER(setupPLot, )
TIMER(setupSpectral, )
TIMER(input, )
TIMER(setProp, )
TIMER(other, )

// auxiliary sections for main code
#include "tRecX.h"

static bool generateLinp()
{
    return ReadInput::main.flag("generate-linp", "only generate list-of-inputs file");
}

// a single "run": use single input file
void tRecX::run_trecx(std::string InputFile, int argc, char *argv[])
{
    //    Sstr+"try this A"+boost::math::sph_bessel(1,0.27557017407092*20)+Sendl;

    Threads::setup(MPIwrapper::communicator());
    START(all);
    LOG_ON();

    //============================================================================================
    // input and part of setup
    //==============================================================================================
    // only overall master outputs to screen
    if (not MPIwrapper::isMaster(MPIwrapper::worldCommunicator()))
        PrintOutput::off("Screen");

    // only master of current flock outputs to file
    if (not MPIwrapper::isMaster(MPIwrapper::communicator()))
        PrintOutput::off("Both");

    ReadInput::openMain(InputFile, argc, argv);
    debug_tools::setVerboseLevel();
    PrintOutput::subTitle("-> tRecX (v 0.9) (C) 2015-2021 Armin Scrinzi (licenced under GPLv2) <-");
    tRecX::tutorialMessage(ReadInput::main.output());
    PrintOutput::paragraph();

    // open output streams for calculation
    if(MPIwrapper::isMaster()){
        PrintOutput::set(ReadInput::main.outputTopDir()+PrintOutput::outExtension);
    }

    string title;
    ReadInput::main.read("Title","",title,"tRecX calculation","title that will be printed to output");

    PrintOutput::title("Date ["+Timer::currentDateTime()+"] "+title);
    std::string build;
#ifdef _NOMPI_
    build+=" Scalar";
#endif
#ifdef _DEVELOP_
    build+=" Develop";
#endif
    if(build!="")PrintOutput::subTitle("     Build: "+build);

    PrintOutput::paragraph();

    Timer::generalMonitor.setShow(0);
    Timer::generalMonitor.start();
    Timer::generalMonitor.monitor(0, "setup", ReadInput::main.outputTopDir() + "mon");

    Units::setDefault("au");        // general default units
    ReadInput::main.setUnits("au"); // convert input to these units

    // string col_path;
    // ReadInput::main.read("Columbus","path",col_path,"None","Columbus run path");
    // #ifdef _USE_HACC_
    // QuantumChemicalInput::read(ReadInput::main);
    // #endif
    // QcInput::read(ReadInput::main);
    // molecularOrbital M(*QcInput::getions());
    // molecularOrbital Mneut(*QcInput::getneutrals());
    // if (M != Mneut)
    //     ABORT("Ion and Neutral MOs are not same");
    
    // DensityMatrix D("ion");
    // DensityMatrix Dion("neutral");

    ProductFunction::read(ReadInput::main);
    ParallelLayout::read(ReadInput::main);

#ifdef _DEVELOP_
    //***************************************************************
    // for parallel debuging
    // compile with -O0
    // mpirun -np 4 ...
    // ps aux | grep tRecX
    // gdb tRecX procNumber
    // gdb> set DebugWait=0 to start each process
    if (ReadInput::main.flag("DEBUGparallel", "block, see main_trecx.cpp how to proceed"))
    {
        unsigned int DebugWait = 1;
        while (DebugWait)
            ;
        MPIwrapper::Barrier();
    }
#endif
    //***************************************************************

#ifdef _OPENMP
    PrintOutput::title("Open MP");
#endif

    tRecX_cache::read(ReadInput::main);
    tRecX::readConstantsAndParameters();
    tRecX::readBases(ReadInput::main);
    RandomPotential::read(ReadInput::main);
    Autoconverge autcon(ReadInput::main);

    OperatorDefinition::setup();

    // command line flags
    bool showOperators = ReadInput::main.flag("showOperators", "show operator structure and stop");
    std::string printMatrices;
    ReadInput::main.read("Flag", "printMatrices", printMatrices, "",
                         "show operator matrix, "
                         "norms[depth,digits]...norms of matrix-blocks to depth,"
                         "full[digits]...all, "
                         "block[ibl,jbl,i0,j0]...block of full matrix",
                         0, "printMatrices");

    // read possible multi-dimensional bases

    string halfPlane = "";
#ifdef _DEVELOP_
    double hpR, hpS, hpEta;
    ReadInput::main.read("_EXPERT_", "halfPlaneR", hpR, "0", "smoothing for half-plane suppression");
    ReadInput::main.read("_EXPERT_", "halfPlaneS", hpS, "5", "smoothing for half-plane suppression");
    ReadInput::main.read("_EXPERT_", "halfPlaneEta", hpEta, "0.2", "smoothing for half-plane suppression");
    if (hpR > 0.)
    {
        halfPlane = "<1><1-trunc[0," + tools::str(hpEta) + "]><1-trunc[" + tools::str(hpR - hpS) + "," + tools::str(hpR) + "]>";
    }
#endif

    OperatorFactor::readMatrix(ReadInput::main);
    ChannelsSubregion::read(ReadInput::main);
    OperatorFloorEE::read(ReadInput::main);
    OperatorFloor3d::read(ReadInput::main);
    tRecX::PlotFunctions(ReadInput::main);
    BasisOrbitalNumerical::read(ReadInput::main);
    IndexNew::read(ReadInput::main);
    SpectrumISurff iSurff(ReadInput::main);

    // ---------------------------------------------------------------------------------
    // Laser Pulse
#ifdef _USE_UI_
    UiLayout::setPage("System");
#endif
    Pulse::read(ReadInput::main, true);
#ifdef _USE_UI_
    UiLayout::setPage();
#endif

    // operators
    tRecX::readExtra(ReadInput::main);

    // Eigenproblem xor initial value problem
    ReadInput::main.exclude("Initial", "Eigen");

    string eigenSelect;
    bool eigenVectors, eigenElements;
    ReadInput::main.read("Eigen", "select", eigenSelect, "NONE",
                         "kinds: All, SmallReal[N]...first N, Nearest[N,Ereal,Eimag]...in complex plane, Rectangle[Rmin,Rmax,Imin,Imax]...in complex plane, NONE",
                         1, "eigenSelect")
        .texdocu(R"tex(
                     \begin{itemize}
                     \item[All] only meaningfull with Lapack for \nameref{docu:Eigen:method}
                     \item[SmallReal] sort by real part, return lowest N
                     \item[Nearest] locate N complex $E_i$ with smallest $|E_i-$Ereal-i Eimag$|$ (best wit Arpack, inverse iteration)
                     \item[Recangle] locate all within rectangle of complex plane
                     \end{itemize}
                     )tex");
    ReadInput::main.read("Eigen", "vectors", eigenVectors, "false", "save eigenvectors (binary format)", 1, "eigenVectors");
    ReadInput::main.read("Eigen", "elements", eigenElements, "false", "compute matrix elements of Operator expectation ", 1, "eigenElements");

    string initialKind;
    int initialN;
    ReadInput::main.read("Initial", "state", initialN, "0", "initial state (number of state in Hamiltonian)");
    initialKind = InitialState::readKind(ReadInput::main);

#ifdef _USE_UI_
    UiLayout::setPage("System");
#endif
    std::string hamDef = tRecX::readOperator("hamiltonian");
    std::string intDef = tRecX::readOperator("interaction");
    std::string iniDef = tRecX::readOperator("initial");
    std::string preconDef = tRecX::readOperator("projection");
    std::string paramDef = tRecX::readOperator("parameterTerm");
    std::string expecDef = readOccupationAppend(tRecX::readOperator("expectationValue"));
    std::string specDef = tRecX::readOperator("spectrum");

    if (iniDef == "hamiltonian")
        iniDef = hamDef;
    if (preconDef == "hamiltonian")
        preconDef = hamDef;
    OperatorDefinition::setParameters(intDef); // need gauge radius for checking discretization setup
#ifdef _USE_UI_
    UiLayout::setPage();
#endif

    // time propagation
    string propMethod;
    double applyThreshold, accuracy, cutE, fixStep, tBeg = 0., tEnd = 0., tPrint, tStore;
#ifdef _USE_UI_
    UiLayout::setPage("Propagation");
#endif
#ifdef _USE_UI_
    UiLayout::setPage();
#endif
    std::string dum1, dum2;
    IndexConstraint::readAll(ReadInput::main, dum1, dum2);

    std::vector<double> surf;
    bool surfAscii(false);
    DiscretizationSurface::read(ReadInput::main,surf,surfAscii);
    int spectrumPoints=0;
    double minEnergy,maxEnergy;
    bool kGrid(false),computeSpectrum(ReadInput::main.found("Spectrum"));

    bool specOverwrite;
    std::string plotWhat, ampFiles;
    SpectrumPlot::read(ReadInput::main, plotWhat, specOverwrite, ampFiles);

    if (TsurffSource::recompute())
    {
        TsurffSource::moveSubdirs(ReadInput::main);
    }

    DiscretizationTsurffSpectra::readFlags(ReadInput::main,spectrumPoints,minEnergy,maxEnergy,kGrid,true);

    int sampleForSpectrum=16;
    if(iSurff.isOn() and Units::convert(1,"OptCyc")<25.){
        PrintOutput::DEVwarning("HACK: introduced fine sampling for high-frequency iSurff, use -DEBUGsample to set sampling");
        sampleForSpectrum=64;
    }

    ReadInput::main.read("DEBUG","sample",sampleForSpectrum,tools::str(sampleForSpectrum),"sample for amplitude integration",1,"DEBUGsample");

    if(maxEnergy<=0.)tStore=-1.;
    else tStore=2.*math::pi/(maxEnergy*sampleForSpectrum);
    TimePropagator::read(ReadInput::main,tBeg,tEnd,tPrint,tStore,accuracy,cutE,fixStep,applyThreshold,propMethod);

    std::string gauge=tRecX::readOperator("surfaceGauge");
    if(gauge=="automatic"){
        if((hamDef+intDef).find("LaserLength")!=std::string::npos)
            gauge="length";
        else if((hamDef+intDef).find("LaserF")!=std::string::npos)
            if(spectrumPoints>0 and tBeg<tEnd)
                ABORT("Hamiltonian "+hamDef+intDef+"contains laser field rather than vector potential, must define Operator: surfaceGauge");
    }
    if(gauge=="length")Algebra::addSpecialConstant("Rg",DBL_MAX);

    if (maxEnergy <= 0.)
        tStore = -1.;
    else
        tStore = 2. * math::pi / (maxEnergy * sampleForSpectrum);
    TimePropagator::read(ReadInput::main, tBeg, tEnd, tPrint, tStore, accuracy, cutE, fixStep, applyThreshold, propMethod);

    // std::string gauge = tRecX::readOperator("surfaceGauge");
    if (gauge == "automatic")
    {
        if ((hamDef + intDef).find("LaserLength") != std::string::npos)
            gauge = "length";
        else if ((hamDef + intDef).find("LaserF") != std::string::npos)
            if (spectrumPoints > 0 and tBeg < tEnd)
                ABORT("Hamiltonian " + hamDef + intDef + "contains laser field rather than vector potential, must define Operator: surfaceGauge");
    }
    if (gauge == "length")
        Algebra::addSpecialConstant("Rg", DBL_MAX);

    TsurffSource::allReads(ReadInput::main);

    std::string region=TsurffSource::nextRegion(ReadInput::main);
    if(region.find("unbound axes:")!=string::npos)
        PrintOutput::set(ReadInput::main.output()+"outspec");
    else if (region!="")
        ReadInput::main.outSubdir("S_"+region);
    STOP(all);

    int line(1);
    AxisTree(ReadInput::main,line);//HACK
    do {
        std::shared_ptr<OperatorTree> initialOper,propOper,hamOper;
        Parallel::clear();
        auto stopped=Timer::stopAll();
        if(stopped.size()>0)PrintOutput::DEVwarning(Sstr+"force-stopped timers"+stopped);
        START(all);
        START(run);
        START(setup);
        LOG_PUSH("execute");
        MPIwrapper::setCommunicator(Threads::all());

        timeCritical::suspend();
        // update spectral region
        region = TsurffSource::nextRegion(ReadInput::main);
        if (generateLinp())
            PrintOutput::set("cout");
        else
        {
            if (region.find("unbound") != string::npos)
            {
                // plot spectra from existing ampl files
                if (region.find("unbound axes:") != string::npos)
                {
                    // plotting does not work in parallel - use single thread
                    MPIwrapper::setCommunicator(Threads::single());
                    if (MPIwrapper::isMaster(Threads::all()))
                    {
                        PrintOutput::set(ReadInput::main.outputTopDir() + "outspec");

                        SpectrumPlot::allRegions(ReadInput::main.outputTopDir(), ReadInput::main);

                        PrintOutput::paragraph();
                        PrintOutput::message("spectra calculated from amplitude file(s), for recomputing amplitudes -recompute");
                    }
                }
                Timer::generalMonitor.monitor(tEnd,"spectrum from ampl-file","",true);
                goto Terminate;
            }
            else if (region != "" and
                     (spectrumPoints > 0 or
                      string(argv[1]).find(ReadInput::inputCopy) + ReadInput::inputCopy.length() == string(argv[1]).length()))
            {
                ReadInput::main.outSubdir("S_" + region);
                PrintOutput::set(ReadInput::main.output() + PrintOutput::outExtension);
            }
        }

        if (generateLinp())
            PrintOutput::off();

        std::string info="Master host = "+tools::cropString(platformSpecific::current_host()+" ");
        if(MPIwrapper::Size(Threads::all())>1)info+="("+tools::str(MPIwrapper::Size())+" processes) ";
        PrintOutput::subTitle(info);
        PrintOutput::paragraph();

        STARTDEBUG(input);
        LOG_PUSH("setup");
        LOG_PUSH("discretization");
        STARTDEBUG(setupDisc);

        std::shared_ptr<Discretization> mainD(Discretization::factory(ReadInput::main));
        BasisOrbital::addIndex("main", mainD->idx());
        if (ReadInput::main.flag("showMainIndex", "print main Index and stop"))
        {
            PrintOutput::message("AxisTree:\n" + mainD->axisTree()->Tree::str());
            PrintOutput::message("Index:\n" + mainD->idx()->str());
            exit(0);
        }
        iSurff.checkCoordinates(mainD->idx()); // only few coordinate systems can do iSurff

        // with the index available, generat numerical orbitals
        BasisOrbitalNumerical::setup();

        STOPDEBUG(setupDisc);
        LOG_POP();
        LOG_PUSH("other");

        // ---------------------------------------------------------------------------------

        // Construct the discretization for the subregion, find all the sources feeding into this region and check if they exist
        std::vector<TsurffSource *> source_terms;
        std::shared_ptr<Discretization> runD(mainD);

        // not nice: in case of surface flux, sf will carry the index, and must not be destroyed
        std::shared_ptr<SurfaceFlux> sf;
        if (region == "")
            mainD->print();
        else
        {
            // this is better now, but still very clumsy

            // calculation in subregion with a source
            std::string path = DiscretizationSurface::surfacePath(ReadInput::main.outputTopDir(), mainD->idx(), region);
            DiscretizationSurface *surface = new DiscretizationSurface(path);

            // surface index on file may become transformed in SurfaceFlux
            // use transformed for constructing this propagation index
            sf.reset(new SurfaceFlux(path, ReadInput::main));

            // spectral discretization (will be forked along radial spectral values)
            runD.reset(new DiscretizationTsurffSpectra(sf->idx(),ReadInput::main));
            source_terms.push_back(new TsurffSource(dynamic_cast<DiscretizationTsurffSpectra*>(runD.get()),ReadInput::main,surface,sf));

            iSurff.addSource(source_terms.back());

            runD->print("", "DISCRETIZATION AND INCOMING FLUX");
            PrintOutput::newLine();
            MPIwrapper::Barrier();
            source_terms[0]->print(path);
            MPIwrapper::Barrier();
        }

        if (surf.size() > 0)
        {
            PrintOutput::paragraph();
            if (surf.size() == 1 and surf[0] > DBL_MAX / 2)
                PrintOutput::lineItem("Surfaces", "default (at complex scaling radius)");
            else
                PrintOutput::lineItem("Surfaces", Str("", "") + surf);
        }
        PrintOutput::paragraph();

        OperatorFloor3d::print();
        CoefficientsWriter::read(ReadInput::main);
        if (region != "")
            CoefficientsWriter::instance()->disable();
        // surface discretization for each unbound coordinate axis
        // NOTE: reverted to using same surface on all coordinates
        // NOTE: the logics below is messy and should be moved to a function of DiscretizationSurface
        vector<DiscretizationSurface *> discSurf;
        discSurf = DiscretizationSurface::all(ReadInput::main, runD, surf, region);

        vector<string> dipNames;
        string dipoleDef = Harmonics::dipoleDefinitions(ReadInput::main, mainD->coordinates(), hamDef + "+" + intDef, dipNames);
        PrintOutput::newRow();
        PrintOutput::title("OPERATORS");
        PrintOutput::paragraph();

        if (hamDef.find("[[eeInt6DHelium]]") != string::npos)
            PrintOutput::lineItem("Hamiltonian", Str(hamDef, "") + " (multipoleOrder=" + OperatorFloorEE::read(ReadInput::main) + ")");
        else
            PrintOutput::lineItem("Hamiltonian", hamDef);

        PrintOutput::newLine();

        std::string subDef = hamDef;
        if (mainD != runD)
        {
            subDef = OperatorDefinition(hamDef, mainD->idx()->hierarchy()).tsurffDropTerms(runD->idx(), runD->idx());
            if (subDef != "")
                PrintOutput::lineItem("Region " + region, subDef);
            else
                PrintOutput::lineItem("Region " + region, "only spectral amplitudes");
            PrintOutput::newLine();
        }
        if (intDef != "")
        {
            PrintOutput::lineItem("Interaction", intDef);
            PrintOutput::newLine();
            if (region != "")
            {
                PrintOutput::lineItem("Region " + region, OperatorDefinition(intDef, mainD->idx()->hierarchy()).tsurffDropTerms(runD->idx(), runD->idx()));
                PrintOutput::newLine();
            }
        }

        if (initialKind == "manyBody")
            iniDef = hamDef;

        string iniName = iniDef;
        if (iniDef == hamDef)
            iniName = "Hamiltonian";
        else if (iniDef == "atBegin")
            iniName = "H(t=Begin)";
        PrintOutput::lineItem("HamInitial", iniName);
        PrintOutput::newLine();

        PrintOutput::newLine();
        if (paramDef != "")
        {
            PrintOutput::lineItem("Static field", paramDef);
            PrintOutput::newLine();
        }
        if (expecDef != "")
        {
            PrintOutput::lineItem("Expectation Values", expecDef);
            PrintOutput::newLine();
        }
        if (preconDef != "")
            PrintOutput::lineItem("Projection", preconDef);
        if (dipoleDef != "")
        {
            PrintOutput::newLine();
            PrintOutput::lineItem("Dipoles", dipoleDef);
            PrintOutput::newLine();
        }
        PrintOutput::paragraph();
        VectorValuedFunction::print();

        // ranges for operator parameters and initialize parameters
        ScanEigenvalue eigenScan(ReadInput::main);
        if (eigenScan.size() > 0)
            eigenScan.print();

        FloquetAnalysis floqAna(hamDef, ReadInput::main);
        EigenTrace eigenTrace(ReadInput::main);

        EigenSubspace eigSub(ReadInput::main);

        Pulse::current.output("", hamDef + intDef);

        //    if(surf.size()>0)VolkovPhase::consistency(hamDef+intDef);

        Str strIni(initialKind, "");
        if (initialN != 0)
            strIni = strIni + " state[" + initialN + "] ";
        strIni = strIni + " for operator " + iniName;
        PrintOutput::lineItem("Initial state", strIni);
        PrintOutput::paragraph();

        string projConstraint;
        ReadInput::main.read("Project", "constraints", projConstraint, "", "constrain spectral projection");

        //        // set up possible spectral cuts
        //        SpectralCut specCut(runD.get(),ReadInput::main,hamDef);

        PrintOutput::title("OUTPUT");
        PrintOutput::lineItem("Directory", ReadInput::main.output());
        PrintOutput::paragraph();
        PrintOutput::paragraph();

        STOPDEBUG(input);
        STARTDEBUG(setupPLot);
        std::shared_ptr<Plot> plot(new Plot(runD->idx(), ReadInput::main)); // set up plot, get the definitions from file
        plot->print();

        std::shared_ptr<PlotCoefficients> plotC(plot);
        STOPDEBUG(setupPLot);
        STARTDEBUG(input);

        if (ReadInput::main.flag("showDiscretization", "print the discretization and stop"))
        {
            runD->print();
            exit(0);
        }

        // if no input specified, just print doc
        if (ReadInput::main.file() == "trecx.inp")
        {
            ReadInput::main.writeDoc();
            cout << "\nUsage: tRecX file[.inp] [-flags]" << endl;
            exit(0);
        }

        if (generateLinp())
        {
            ReadInput::main.writeDoc();
            cout << "\nList-of-inputs (linp) written to " + ReadInput::main.output() << endl;
            exit(0);
        }

        // create gradient operator, if specified (=0 else)
        OperatorAbstract *chanMap = 0, *grad = 0;
        string ChanAxis = OperatorMapChannelsSurface::axis(ReadInput::main);
        if (ChanAxis != "")
            grad = OperatorGradient::read(runD.get(), ReadInput::main);

        // constraints for initial state calculations (unmainained)
        DiscretizationConstrained::inputs(ReadInput::main);

        tRecX::read();
        tRecX::off("EigenSolver");
        tRecX::print();

        DensityOfStates dos(ReadInput::main);

        // further output
        bool saveWf;
        ReadInput::main.read("Output", "wavefunction", saveWf, "false", "save coefficients at print interval");
        if (ReadInput::main.flag("showIndex", "print Indices and stop"))
        {
            PrintOutput::message("run Index\n" + runD->idx()->str());
            if (sf)
                PrintOutput::message("Surface Index\n" + sf->idx()->str());
            for (auto s : discSurf)
            {
                PrintOutput::message(s->name + " Index\n" + s->idx()->str());
            }
            for (auto p : BasisOrbital::referenceIndex)
                PrintOutput::message(p.first + " Index\n" + p.second->str());
        }

        // all inputs should be finished above these lines - please do not remove
        ReadInput::main.finish();
        PrintOutput::paragraph();
        STOPDEBUG(input);
        //==== END OF INPUT ====================================================================

        if (floqAna.run(mainD->idx()))
            goto Terminate; // floquet analysis only
        if (eigenTrace.run(hamDef, mainD->idx(), ReadInput::main.output()))
            goto Terminate; // trace eigenvalues from guess

        // can use further cleanup... (get rid of Chandef/ChanAxis)
        std::string Chandef = "";
        if (ChanAxis != "")
            chanMap = new OperatorMapChannelsSurface(ReadInput::main, runD.get());

        LOG_POP();
        LOG_PUSH("Operators");

        STARTDEBUG(operator);
        vector<OperatorAbstract*>printOps;

        tRecX::SetOperators(mainD->idx()->hierarchy(), hamDef, intDef, iniDef, "", dipNames, expecDef, runD.get(), propOper, hamOper, initialOper, printOps);
        MPIwrapper::Barrier();


#ifdef _USE_HACC_
        // Shift the full Hamiltonian spectrum by a specified channel energy
        OperatorMapChannelsSurface *mapch = static_cast<OperatorMapChannelsSurface *>(chanMap);
        if (mapch != nullptr && mapch->energyShift() != 0.)
        {
            std::string ovl;
            if (hamOper->iIndex->isHybrid()) // hybrid basis needs special string
                ovl = "<1>(<<FullOverlap>>)+<1,0>[[haCC1]]+<0,1>[[haCC1]]";
            else
                ovl = "<<FullOverlap>>";
            ovl = tools::str(-mapch->energyShift()) + ovl;

            // HACK: As for now, it is not very clear what happens with OperatorTrees that got added/fused --> create two shift operators to be sure, don't dare to "delete" the pointers, as this may delete floors
            OperatorTree *shift = new OperatorTree("IonicShift", OperatorDefinition(ovl, hamOper->iIndex->hierarchy()), hamOper->iIndex, hamOper->iIndex);
            hamOper->add(shift);
            shift = new OperatorTree("IonicShift", OperatorDefinition(ovl, hamOper->iIndex->hierarchy()), hamOper->iIndex, hamOper->iIndex);
            propOper->add(shift);
        }
#endif

        Harmonics::addDipolesToPlot(dipoleDef, *plot, printOps);
        STOPDEBUG(operator);
        LOG_POP();
        LOG_PUSH("other");

        STARTDEBUG(other);
        tRecX::PrintShow(printMatrices, showOperators, *propOper, *initialOper, printOps);

        // make sure inverse of overlap works correctly
        if (propOper != 0)
            runD->idx()->testInverseOverlap();
        MPIwrapper::Barrier();

        if (specDef != "")
        {
            PrintOutput::DEVwarning("Operator:spectrum temporarily out of service");
            //            if(plot->dimension()>0)ABORT("wave function plot and spectral plot mutually excluded");
            //            if(specDef!="hamiltonian")ABORT("for now, only for Hamiltonian");
            //            std::shared_ptr<DiscretizationSpectral>spec(new DiscretizationSpectral(runD.get(),initialOper.get()));
            //            plotC.reset(new PlotSpectral(spec.get()));
        }

        // output dos (if set up) and stop
        dos.output(propOper.get());

        STOPDEBUG(other);
        LOG_POP();
        LOG_POP();

        if (eigenSelect != "NONE")
        {
            Timer::generalMonitor.monitor(0, "Eigensolver", ReadInput::main.outputTopDir() + "mon");
            Parameters::updateSpecial();
            tRecX::ComputeEigenvalues(eigenSelect, *plot, printOps, eigenVectors, eigenElements);
            MPIwrapper::Barrier();
            Timer::generalMonitor.monitor(0,"Eigenvalues",ReadInput::main.outputTopDir()+"mon");
            goto Terminate;
        }

        else if (eigenScan.size() > 0)
        {
            // scan eigenvalues for a range of parameters
            Parameters::updateSpecial();
            DEVABORT("re-formulate for OperatorTree");
            //        Operator op(*propOper);
            //        eigenScan.setEigenSub(D,&op,dynamic_cast<const Operator*>(printOps[1]));
            //        eigenScan.scan();
        }

        else if (tBeg < tEnd)
        {
            STARTDEBUG(setProp);
            LOG_PUSH("DerivativeFlat");

            // no eigenvalue - time-propagate
            DerivativeFlat *propDer = tRecX::SetDerivative(runD.get(), propOper.get(), hamOper.get(), initialOper.get(), applyThreshold, cutE, preconDef, projConstraint, source_terms);
            LOG_POP();

            // set up the output
            TimePropagatorOutput out;
            tRecX::SetOutput(out,tPrint,tStore,surfAscii,discSurf,plotC.get(),grad,ChanAxis,Chandef,chanMap);
            out.setInterval(tBeg,tEnd);
            for(auto p: printOps)out.addExpec(p);
            if(saveWf)out.withApplyAndWrite(new OperatorIdentity(runD->idx()),"wf_t");
            if(dipoleDef!="")out.sampleExpec(1); /// when harmonics (dipoles) are desired, take ALL store points)

            // checks and messages before start of time propagation
            std::shared_ptr<ChannelsSubregion> channels;
            if (propOper != 0)
            {
                OperatorDefinition(propOper->def()).checkVariable("non-const parameters may have penalty in time-propagation, avoid!", propOper->name);
                PrintOutput::warningList();

                LOG_PUSH("other");
                // set up channels
                if (runD != mainD and (runD->idx()->hierarchy() != "specX.Z.kX.Z"))
                {
                    //                    LOG_PUSH("ChannelsSubregion");
                    //                    runD->print();
                    //                    channels.reset(new ChannelsSubregion(runD.get(),hamOper.get(),region, ReadInput::main));
                    //                    LOG_POP();
                    //                    out.withChannelsSubregion(channels.get());
                    PrintOutput::DEVwarning("channels subregion disabled for development");
                }
                LOG_POP();
            }
            // time-propagator (for parallel code, create local derivative operator)
            LOG_PUSH("TimePropagator");

            TimePropagator *prop;
            DerivativeLocal *derLoc = new DerivativeLocal(propDer);
            OdeStep<LinSpaceMap<CoefficientsLocal>, CoefficientsLocal> *odeLoc;

            if (propOper)
                odeLoc = odeFactory<LinSpaceMap<CoefficientsLocal>, CoefficientsLocal>(propMethod, derLoc, accuracy);
            else
                odeLoc = odeFactory<LinSpaceMap<CoefficientsLocal>, CoefficientsLocal>("euler", derLoc, accuracy);

            prop = new TimePropagator(odeLoc, &out, accuracy, fixStep);

            LOG_POP();
            STOP(setup);

            // timer info for all setup
            PrintOutput::warningList();
            PrintOutput::timerWrite();
            Parallel::timer();

            if (ReadInput::main.flag("DEBUGsetup", "stop after setup"))
            {
                STOPDEBUG(setProp);
                STOP(run);
                goto Terminate;
            }

            // initial wave function
            if (runD != mainD)
                initialKind = "ZERO";
            else if (iniDef != "atBegin" and initialKind != "manyBody")
                initialKind = "Hinitial";

            // If there is a stored checkpoint, use it. Otherwise proceed normally
            Wavefunction wf;
            LOG_PUSH("InitialState");
            wf = InitialState::get(runD.get(), initialKind, initialN, tBeg, *initialOper, propDer);
            if (dynamic_cast<DerivativeFlat *>(derLoc))
                dynamic_cast<DerivativeFlat *>(derLoc)->project(*wf.coefs);
            if (propOper!=0 and expecDef.find("overInitial") != std::string::npos)
            {
                Coefficients c(wf.coefs->idx());
                wf.coefs->idx()->localOverlap()->apply(1., *wf.coefs, 0., c);
                std::shared_ptr<ProjectSubspace>proIni(new ProjectSubspace(std::vector<Coefficients *>(1, wf.coefs), std::vector<Coefficients *>(1, &c)));
                proIni->name = "overInitial";
                out.addExpec(new ProjectSubspaceToDual(proIni));
            }
            LOG_POP();

            // propagation monitor settings
            string mess = runD->name + " Tol=" + prop->str(1) + " " + Pulse::current.str(1);
            if (MPIwrapper::isMaster())
                out.timer()->monitor(wf.time, mess, ReadInput::main.outputTopDir() + "mon");
            STOPDEBUG(setProp);
            LOG_PUSH("propagation");

            if (propOper == 0)
            {
                // no actual propagation, only integration of source

                // add calculation of overlap
                out.sampleExpec(1); // put all into expec file
                out.addExpec(new OperatorTree("Overlap",OperatorDefinition::defaultOverlap(runD->idx()),
                                              runD->idx(),runD->idx()));
                out.checkGrowingNorm(false); // do not warn about norm>1

                // set turn off range
                for(TsurffSource* s: source_terms)s->setTurnOff(ReadInput::main,tEnd);

                //NOTE: we integrate by Euler, sample "sampl" times across shortest cycle in the spectrum
                prop->fixStep(2*math::pi/(sampleForSpectrum*maxEnergy));
            } else {
                for(TsurffSource* s: source_terms){
                    if(s->str().find("turn off")!=std::string::npos)
                        PrintOutput::DEVwarning("source turn off for region "+region+"- undefined behavior");
                }
            }

            // time-propagate
            if (wf.time < min(Pulse::gettEnd(), tEnd))
                prop->propagate(&wf, min(Pulse::gettEnd(), tEnd), "Start");
            LOG_POP();
            LOG_PUSH("afterwards");
            if (Pulse::gettEnd() < tEnd)
                prop->propagate(&wf, tEnd, "Stop");
            LOG_POP();

            // amplitude calculation
            if(propOper==0 or region=="Rn1" or region=="Rn2"){
                if(MPIwrapper::isMaster()){
                    ofstream ampl((ReadInput::main.output()+"ampl").c_str(),(ios_base::openmode) ios::beg|ios::binary);
                    Coefficients res(*wf.coefs);
                    iSurff.addCorrection(res);
                    const Coefficients * joinedWf=Threads::join(res);
                    if(Threads::isMaster())joinedWf->write(ampl,true,"IndexFull");
                    PrintOutput::message("Photo-electron amplitudes on "+ReadInput::main.output()+"ampl");
                }
            }
            else if (region=="" and iSurff.isOn()){
                // save H0 and wave function into iSurff (if active)
                iSurff.endPropagation(wf.time,*wf.coefs,hamOper,discSurf.back()->sharedFromParent());
            }

            if (region == "Rn1" or region == "Rn2")
            {
                // single-ionization spectra
                const Coefficients *joinedWf = Threads::join(*wf.coefs);
                if (Threads::isMaster())
                {
                    std::vector<std::string> ax = region == "Rn1" ? std::vector<std::string>(1, "kRn1") : std::vector<std::string>(1, "kRn2");
                    std::vector<std::string> use(1, "g");
                    if (region == "Rn1")
                        ax[0] = "kRn1";
                    Plot plt(joinedWf->idx(), ax, use);
                    plt.plot(*joinedWf, ReadInput::main.output() + "spec_single");
                    PrintOutput::message("single ionization spectrum on " + ReadInput::main.output() + "spec_single");
                }
            }

            PrintOutput::end();
            PrintOutput::paragraph();
            PrintOutput::subTitle("TimePropagator " + prop->info());
            delete propDer;
            delete prop;
            MPIwrapper::setCommunicator(Threads::all());
            Timer::generalMonitor.monitor(wf.time,"propagate region "+region,ReadInput::main.outputTopDir()+"mon");

        }
        else {
            PrintOutput::message(Str("tBeg =")+tBeg+">= tEnd="+tEnd+"--- no time propagation",0,true);
            Timer::generalMonitor.monitor(0,"no time propagation"+region,ReadInput::main.outputTopDir()+"mon");
        }
        if(MPIwrapper::isMaster() and dipNames.size()>0)
            Harmonics(ReadInput::main,"expec"," ",false,1.);

        LOG_POP();

        STOP(run);
        STOP(all);

        PrintOutput::warningList();
        PrintOutput::timerWrite();
        Parallel::timer();
    } while ((computeSpectrum or ReadInput::main.found("Spectrum")) and (region != TsurffSource::nextRegion(ReadInput::main) or TsurffSource::nextRegion(ReadInput::main).find("") != string::npos) and surf.size() > 0 and spectrumPoints != 0);

Terminate:
    LOG_POP();
    LOG_POP();

    auto stopped=Timer::stopAll();

    if(stopped.size()>0){
        PrintOutput::DEVmessage(Sstr+"force-stopped timers, rewrite timer results: "+stopped);
        PrintOutput::timerWrite();
        Parallel::timer();
    }


    PrintOutput::paragraph();
    PrintOutput::terminationMessages();

    Timer::generalMonitor.monitorDone(0);
    RandomPotential::clearAll(); // remove old, all re-generate of random potential by the same name
    PrintOutput::title("done - "+ReadInput::main.output());

}
