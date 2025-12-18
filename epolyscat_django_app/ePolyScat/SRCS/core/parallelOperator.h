// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef PARALLELOPERATOR_H
#define PARALLELOPERATOR_H

#include "blockView.h"

class ParallelProcess;

/// parallel operations on OperatorTree
class ParallelOperator
{
    OperatorAbstract * _op;
    BlockView _view;

    static std::map<std::string,int> _host;
public:
    static bool emptyLeafOrDummy(const OperatorTree* Op);
    static void sync(OperatorTree* Op);
    static void bcast(OperatorTree* Op);

    static constexpr int undefined=-1;
    static constexpr int none=-2;
    static constexpr int all=-3;
    static constexpr int distributed=-4;
    static int getHost(const OperatorAbstract* Op, bool Warn=true);
    static std::string strHost(const OperatorAbstract* Op);
    static void unsetHost(const OperatorTree *Op);

    ParallelOperator(const OperatorAbstract *Op);

    /// host for a given OperatorTree
    class Host{
    public:
        virtual ~Host(){}

        /// return host number for Op
        virtual int operator() (const OperatorTree* Op) const =0;
        virtual std::string kind() const=0;
    };

    /// actual hosting of OperatorTree leaf (all/none/distributed)
    class PresentHost:public Host{
    public:
        virtual std::string kind() const{return "present";}
        PresentHost(){}
        int operator() (const OperatorTree *Op) const;
        int redetermine(const OperatorTree *Op) const;
    };
    class SingleHost:public Host{
        int _host;
    public:
        virtual std::string kind() const{return "single";}
        SingleHost(int Host):_host(Host){}
        int operator() (const OperatorTree *Op) const{return _host;}
    };

    /// hosting according to parallel on parallel processes
    class ProcessHost:public Host{
        std::map<const OperatorTree*,int> _proc;
    public:
        virtual std::string kind() const{return "process";}
        ProcessHost(const std::vector<ParallelProcess *> &Proc);

        /// assign Op to one single Host
        ProcessHost(const OperatorTree* Op, int Host);

        int operator()(const OperatorTree *Op) const;
    };

    /// distribute diagonal blocks of operator (include overlap diagonality)
    class DiagonalBlockHost:public Host{
        const OperatorTree * _root;
    public:
        virtual std::string kind() const{return "block-diagonal";}
        /// node count relative to Root
        DiagonalBlockHost(const OperatorTree *Root):_root(Root){}
        int operator()(const OperatorTree *Op) const;
    };

    /// record distributon (for later recovery)
    class SavedHost:public Host{
        std::map<const OperatorAbstract*,int> savedHost;
    public:
        /// node count relative to Root
        virtual std::string kind() const{return "saved";}
        SavedHost(const ParallelOperator* Par);
        int operator()(const OperatorTree *Op) const;
    };

    static void setDistribution(const OperatorAbstract* Op, const Host& Hosts=PresentHost());
    void reDistribute(const Host& NewHosts); ///< move floor data to hosts according to NewHosts
    void setDistribution(){reDistribute(PresentHost());} ///< determine and register distribution

    void bcast();                        ///< distribute floor data to all processes
    void purge();                        ///< total purge: emtpy branches and zero-branches
    void fuse();                         /// fuse, but make sure parallel is OK
    void syncNormCost();                 ///< broadcast norm and application cost

    std::string str() const;
};

#endif // PARALLELOPERATOR_H
