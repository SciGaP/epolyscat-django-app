# tRecX = tSurff+irECS - a universal Schroedinger solver
# Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
# 
# This program is free software; you can redistribute it and/or modify it 
# under the terms of the GNU General Public License as published by the Free Software Foundation; 
# either version 2 of the License, or (at your option) any later version.
# End of license
 
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

from tools import *
from plotTools import *
from compareToFirst import *

colr=['red','green','blue','cyan','magenta','brown']
markr=['.','+','x','s','*']

def toPolar(CosTh,Value,Phi):
    for k in range(len(CosTh)):
        Value[k]*=1#sqrt(1.-CosTh[k]*CosTh[k])
    theta=np.arccos(CosTh)
    if abs(Phi-3.14)<1.e-2: theta*=-1
    elif abs(Phi)>1e-2: print("found phi=",Phi,"plot in upper half plane")
    return theta,cp.deepcopy(Value)

def plot1d(fax,axn,comp,datfil,flags,pars,plot,col,gs,nPlot,doPlot,info,leg,lastPlot,interactive_compare=0):

    style=PlotStyle(flags)
    lineColor=colr[(nPlot+1)%len(colr)]

    if plot.isLineout():
        xVal=datfil.xAxis()
        yVal=datfil.yAxis()
        fVal=datfil.column(col)
        atCoor,xVal,fVal=plot.lineout(xVal,yVal,fVal)
        print("linout at "+str(atCoor)+" on file "+datfil.name+"_lineout"+str(atCoor)[:5])
        f=open(datfil.name+"_lineout"+str(atCoor)[:5],'w')
        f.write("#  \n")
        for k in range(len(xVal)):
            f.write(str(xVal[k])+", "+str(fVal[k])+"\n")
        f.close()
    else:
        xVal=datfil.xAxis()
        fVal=datfil.column(col)

    if "-eV" in flags:
        for k in range(len(xVal)):
            fVal[k]*=xVal[k]
            xVal[k]*=xVal[k]*0.5*27.211

    """ HACK for converting phase into attoseconds """
    if "-attoseconds" in flags:
        for k in range(len(xVal)):
            xVal[k]*=24.188*110.3/(2*3.14159)

    xVal,fVal=average(flags,xVal,fVal)
    normalize(flags,xVal,fVal)

    maxPos=""
    if "-maxgraph" in flags:
        maxPos='{:02.3f}'.format(maxGraph(xVal,fVal),"value",np.max(fVal),"file",datfil.name)
        print("maxpos at",maxPos)
    minPos=""
    if "-mingraph" in flags:
        minPos='{:02.0f}'.format(minGraph(xVal,fVal),"value",np.max(fVal),"file",datfil.name)


    if "-polar" in flags:
        phi=float(datfil.colName(col).split(',')[2])
        xVal,fVal=toPolar(xVal,fVal,phi)
        if "-normalize" in flags: fVal/=np.max(fVal)

    # extend the plot sizes to present
    plot.extend(np.min(fVal),np.max(fVal),np.min(xVal),np.max(xVal))

    # get label (parameter values can be prepended by "-label=PARNAME")
    lab=str(col)+":";
    if maxPos!="": lab+=" max@="+maxPos
    if minPos!="": lab+=" min@="+minPos

    vMin,vMax=getRange(flags,fVal,"-vrange",1e-5)
    linY=np.min(fVal)<=0. or "-linY" in flags
    linY=linY and "-logY" not in flags and "-logV" not in flags
    if (interactive_compare==0 and ("-compare" in flags or "-ratio" in flags)) or interactive_compare==1:
        if nPlot==0:
            lab+=" "+leg.label(axn,col,datfil,info,pars)
            comp=Compare(linY,style)
            comp.reference(xVal,fVal,lab,lineColor)
        else:
            lab+=" "+leg.label(comp.ax2,col,datfil,info,pars)
            comp.compare(xVal,fVal,lab,colr[(nPlot+1)%len(colr)],"-ratio" in flags)
        if lastPlot:
            comp.plot()

    elif "-polar" in flags:
        axn=fax.add_subplot(gs[0],projection="polar")
        lab+=" "+leg.label(axn,col,datfil,info,pars)
        axn.plot(xVal,fVal,label=lab,color=lineColor)

    else:
        lab+=" "+leg.label(axn,col,datfil,info,pars)
        if "-points" in flags:
            if linY: axn.plot(xVal,fVal,'s',label=lab,color=lineColor,marker=markr[(nPlot+1)%len(markr)])
            else:    axn.semilogy(xVal,fVal,'s',label=lab,color=lineColor,marker=markr[(nPlot+1)%len(markr)])
        else:
            if linY:
                axn.plot(xVal,fVal,label=lab,color=lineColor)

            else:
                axn.semilogy(xVal,fVal,label=lab,color=lineColor)


    return axn,comp
