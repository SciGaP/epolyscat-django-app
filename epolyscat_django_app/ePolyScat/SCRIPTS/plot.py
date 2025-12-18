#! /usr/bin/env python3

# tRecX = tSurff+irECS - a universal Schroedinger solver
# Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
# 
# This program is free software; you can redistribute it and/or modify it 
# under the terms of the GNU General Public License as published by the Free Software Foundation; 
# either version 2 of the License, or (at your option) any later version.
# End of license
 
"""
plot one or several data-files with column-wise storage

columns can be specified in square brackets after the file name
        in the form fileName[2,3,7-9] for plotting columns 2,3,7,8,9
        without square bracket, all columns are plotted
column 0 of each file is interpreted as x-axis by default
       for different x-axis use [3:7,8], which plot columns 7,8 against 3
column headers will be recognized and used as legend labels
"""

import sys


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.interpolate
from scipy.interpolate import griddata
import copy  as cp
from math import *

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors, ticker, cm
from matplotlib.colors import LogNorm

import tools as tool
from tools import *
from plotTools import *
from compareToFirst import Compare
import plot1d as p1
import plot2d as p2


def addPolar(CosX,Y):
    return 1



help=[
 "-label=PAR1,PAR2, ... parameter value(s) for PAR1,PAR2,...legend uses value(s) of linp name(s) containing PARi's"
,"-showColumns      ... print the columns names (if found in file)"
,"-peaks            ... mark photon peaks at n*omega-Ip-Up or (-2 Up for Helium)"
,"-eV            ...... transform x-axis as x -> x^2/2 *Rydberg"
,"-polar            ... polar plot of 1d or 2d data"
,"-compare          ... compare multiple 2d plots to first"
,"-linV             ... plot on linear value-axis (default is log)"
,"-linY             ... plot on linear y-axis (default is log)"
,"-logY             ... plot on logarithmic y-axis (default)"
,"-vrange=[vmin,vmax].. function value range for ALL plots "
,"-xrange=[xmin,xmax].. plot axis range for All plots "
,"-yrange=[ymin,ymax].. plot axis range for All plots "
,"-erange=[emin,emax].. range for plotting errors "
,"-rmax=r           ... radius in polar plot"
,"-equalAx          ... 2d with commensurate axes"
,"-lineoutX=x{,w}   ... lineout of 2d plot nearest to x or sum [x-w/2,x+w/2], similar for Y"
,"-normalize{=x0}   ... scale to maximal value = 1 (nearest to x0)"
,"-symX12           ... symmetrize 2d plots by 1 <--> 2 (special for He)"
,"-ratio            ... ratio of multiple plots to first"
,"-maxgraph         ... print location of graph's maximum"
,"-mingraph         ... print location of graph's minimum"
,"-sum              ... add files and columns single plot"
,"-plotfile=file    ... specify plot file name"
,"-batch            ... do keep window open"
,"-points           ... plot points rather than lines"
,"-surface          ... do surface rather than contour plots"
,"-figure=file.fig  ... read defaults and additional style options (will be superseeded by command line flags)"
]

(files,flags) = argsAndFlags()

#  check flags
flagsInHelp(flags,help)

(defaultFiles,defaultFlags,properties) = fileArgsAndFlags(help)
if len(files)==0: files=defaultFiles
flags=defaultFlags+flags


if len(files)<1:
    print("Usage: ")
    print("   plot.py file1[{xcol:}colrange] file2[{xcol:}colrange] ... {flags}")
    print("Example:")
    print("   plot.py Argon/0015/spec[3:22,26,30] -label='I(W/cm2),nCoefficients[2]'")
    print("")
    print("   will plot column 3 vs columns 22,26,30 of files Argon/0015/spec")
    print("   with legend labels including intensity and number of coeficients of 2nd axis")
    print("   first match in Argon/0015/linp is used to determine the value for -label")
    print("   column for x-axis defaults to 0, if colums are not specified, all columns will be plotted vs. 0")
    print("")
    print("Note: use available command shell syntax for specifying multiple files, e.g.")
    print("   plot.py Argon/00{15,14,19}/spec[1] instead of")
    print("   plot.py Argon/0015/spec[1] Argon/0014/spec[1] Argon/0019/spec[1]")
    print("")
    print("Command line flags:")
    for h in help: print("  ",h)
    sys.exit(0)

# try extract pre- and postfixes from all file names
dir,which,files=preAndPostFix(files)

# get total number of datasets to be plotted
nSets=0
for f in files:
    file=dir+f+which
    if file.find("[")==-1:
        datfil=DataFile(file)
        nSets+=len(datfil.datCols)
    else:
        nSets=nSets+len(expandRange(file[max(file.find('['),file.find(':'))+1:file.find(']')]) )

compare_interactive=0

def allPlots(toggle_compare=False):

    global compare_interactive
    if toggle_compare:
        if compare_interactive==0:
            if "-compare" in flags:
                compare_interactive=-1
            else:
                compare_interactive=1
        elif compare_interactive==1:
            compare_interactive=-1
        else:
            compare_interactive=1

    comp=None
    countAll=-1

    # loop through input files
    plot=Plot(flags)
    previous=""
    previousDir=""
    nPlot=-1
    leg=None
    #loop through files
    for info in files:
        file=dir+info+which
        pars=RunParameters(file[:file.rfind('/')])
        # get data into selected columns
        datfil=DataFile(file)

        if "-showColumns" in flags and len(datfil.colNames)>0:
            strCols="Columns:"
            for k in range(len(datfil.colNames)):
                strCols+=" "+str(k)+"="+datfil.colNames[k]
            print(strCols)


        if info==files[0]:
            # first file --- set up the figure
            if datfil.isMulti and not plot.isLineout():
                """ two axes -- separate panels """
                nr,nc=plot.layout(nSets)
                gs = gridspec.GridSpec(nr,nc,height_ratios=[1]*nr)
                fax=plt.figure(1)
                cmap = plt.cm.get_cmap("gnuplot")
            else:
                """ single axis --- all in one panel """
                gs = gridspec.GridSpec(1,1,height_ratios=[1,])
                fax=plt.figure(1)
                axn=fax.add_subplot(gs[0])

        if leg==None or dir=="" and which=="":
            if info in properties: leg=Legend(dir,which,flags,properties[info])
            else:                  leg=Legend(dir,which,flags)


        # loop through columns in file
        for col in datfil.datCols:
            nPlot+=1
            countAll+=1
            # sum all data into single plot, plot when arrived at last
            doPlot="-sum" not in flags or (info==files[-1] and col==datfil.datCols[-1])


            if datfil.isMulti and not plot.isLineout():
                 p2.plot2d(fax,datfil,flags,pars,plot,col,gs,nPlot,doPlot,info,countAll,compare_interactive)
            else:
                axn,comp=p1.plot1d(fax,axn,comp,datfil,flags,pars,plot,col,gs,nPlot,doPlot,info,leg,countAll+1==len(files)*len(datfil.datCols),compare_interactive)


    if (comp==None and (not datfil.isMulti or plot.isLineout())):
        # single plot, adjust sizes
        xMin,xMax=getRange(flags,plot.dim[0],"-xrange")
        vMin,vMax=getRange(flags,plot.dim[2],"-vrange")

        axn.set_xlim(xMin,xMax)
        axn.set_ylim(vMin,vMax)

        # mark peak positions (if so desired)
        if "-peaks" in flags: peakPositions(plot,pars,axn,flags)

    if ((not datfil.isMulti or plot.isLineout())):
        axTop=axn
        axBot=axn
        if comp!=None:
            axTop=comp.ax1
            axBot=comp.ax2

        # axis labels
        xlab=datfil.cols[0][0]
        if "-eV" in flags: xlab="eV"
        axBot.set_xlabel(xlab,fontsize=14)

        # line graph title
        if dir+which!="" and len(files)>1: longTitle=dir+"*"+which
        else:                              longTitle=dir+files[0]+which
        axTop.set_title(tool.briefTitle(longTitle,40,6))

    if comp!=None: fax=comp.fax
    return fax


def plotFile(flags,file):
    for f in flags:
        if f.find("-plotfile")==0:
            ff=f.split("=")
            if(len(ff)>1):
                 file=ff[1]
    if file.find(".png")==len(file)-4:
        file=file[:-4]
    return file




if(True):
    style=PlotStyle(flags)

    # this catches key press on the figure (do not know yet how to use)
    def releaskey(event):
        name=plotFile(flags,sys.argv[0].split(".py")[0])

        if event==None:
            fax=allPlots()
            fax.canvas.mpl_connect('key_release_event', releaskey)

        elif event.key=="u" or event.key=="enter":
            plt.clf()
            c=allPlots()
            plt.draw()

        elif event.key=="c":
            plt.clf()
            c=allPlots(toggle_compare=True)
            plt.draw()

        elif event.key!="l":
            print("unknown key: "+event.key)

        plt.legend(frameon=False)
        # get python file name (as run from the command line)
        plt.figure(1).savefig(name+".png")
        plt.figure(1).text(0.02,0.01, 'u..update, c..compare, l..lin/log, q..quit, graph '+plotFile(flags,sys.argv[0].split(".py")[0])+'.png', fontsize=8)

        if event==None: return "initial"
        else:           return event.key


    key=releaskey(None)

    if "-batch" not in flags:
        plt.show()


