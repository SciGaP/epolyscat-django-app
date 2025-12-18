// tRecX = tSurff+irECS - a universal Schroedinger solver
// Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
// 
// This program is free software; you can redistribute it and/or modify it 
// under the terms of the GNU General Public License as published by the Free Software Foundation; 
// either version 2 of the License, or (at your option) any later version.
// End of license
 
#ifndef OPERATOR_TREE_EXPLORER_H
#define OPERATOR_TREE_EXPLORER_H

#include <map>
#include <string>
#include <sstream>

#include "tools.h"
#include "operatorSVD.h"
#include "operatorTree.h"

#ifdef _OPERATOR_TREE_EXPLORER_

#include "vincurses.h"

/**
 * Simple ncurses based moving around an operator tree. Basic commands:
 *  - a: ascend
 *  - <C-w>hjkl: switch window
 *
 * In main:
 *  - m: print matrix
 *  - s: print matrix.show() (xXoO....)
 *  - i: info
 *
 * Other windows are basically hacky and changed for the particular use case (and hence are not
 * thoroughly tested, in doudbt simply comment out).
 *
 * Requires ncurses and LIBRARIES/vincurses libraries.
 */
class OperatorTreeExplorer: public ViNCurses::App{

    class MainWindow: public ViNCurses::Window{
        char state;

    protected:
        void render();
    public:
        MainWindow(): ViNCurses::Window("Main"), state('i'){}
    
        bool command(std::string command);
        void onTraverse();
    };

    class NavigatorWindow: public ViNCurses::Window{
        int selectedChild;

    protected:
        void render();
    public:
        NavigatorWindow(): ViNCurses::Window("Nav"), selectedChild(0){}

        bool command(std::string command);
        void onTraverse();
    };

    class ApplicationCostWindow: public ViNCurses::Window{
    
    protected:
        void render();
    public:
        ApplicationCostWindow(): ViNCurses::Window("Application Cost"){}
        
        bool command(std::string command);
        void onTraverse();
    };

    class SVDWindow: public ViNCurses::Window{

    protected:
        void render();
    public:
        SVDWindow(): ViNCurses::Window("SVD"){}

        bool command(std::string command);
        void onTraverse();
    };

    class PlaygroundWindow: public ViNCurses::Window{
        int hmat_counter;
    protected:
        void render();
    public:
        PlaygroundWindow(): ViNCurses::Window("Playground"){}

        bool command(std::string command);
        void onTraverse();
    };

    // Windows
    MainWindow* w_main;
    NavigatorWindow* w_nav;
    ApplicationCostWindow* w_appc;
    SVDWindow* w_sing;
    PlaygroundWindow* w_play;


    // Data
    const OperatorTree* optree;
    std::map<const OperatorAbstract*, UseMatrix*> mat;

    UseMatrix* matrix(const OperatorAbstract* op);

protected:
    void init_windows();
public:
    OperatorTreeExplorer(const OperatorTree* _optree): optree(_optree) {}
    bool command(std::string command, bool before_windows);

    // Getters
    const OperatorTree* getOptree(){ return optree; }
    UseMatrix* getMatrix(){ return matrix(optree); }

    // Traversal methods
    void descend(int nChild);
    void descend(int nI, int nJ);
    void ascend();
    void traverse(const OperatorTree* _optree);
};

#else //_OPERATOR_TREE_EXPLORER_

class OperatorTreeExplorer{
public:
    OperatorTreeExplorer(const OperatorTree* optree){ ABORT("Compiled without support for OperatorTreeExplorer"); }
    void run(){}
};

#endif

#endif //OPERATOR_TREE_EXPLORER_H
