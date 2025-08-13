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

previous_compare=-1

def plot2d(fax,datfil,flags,pars,plot,col,gs,nPlot,doPlot,info,countAll,interactive_compare=0):
    global previous_compare
    global refX
    global refY
    global refValues

    if "-linY" in flags:
        print("cannot choose -linY or -logY in 2d plot, maybe you meant -linV or -logV ?")
        exit()

    if interactive_compare!=previous_compare:
        count_all=-1
        previous_compare=interactive_compare

    file=datfil.name
    xVal=datfil.xAxis()
    yVal=datfil.yAxis()
    fVal=datfil.column(col)

    xGrid,yGrid,vMatr=cp.deepcopy(xyzToAxisArray(xVal,yVal,fVal))
    if "-sum" in flags:
        if info==files[0] and col==datfil.datCols[0]:
            xSum=xGrid
            ySum=yGrid
            fSum=vMatr
        else:
            fMatr+=fSum
            if (xSum!=xGrid).any() or (ySum!=yGrid).any():
                exit("axes do no match")

    # transformations:
    # energy-axis in eV
    if "-symX12" in flags:
        tmp=vMatr
        vMatr=tmp+tmp.transpose()
        print("symmetrized 1<-->2")
    if "-eV" in flags: momentumToEnergyEv(xGrid,yGrid,vMatr)
    if "-normalize" in flags: vMatr=vMatr/np.max(abs(vMatr))


    ratio=1.e-9
    if "-linV" in flags: ratio=None
    xMin,xMax=getRange(flags,xGrid,"-xrange")
    yMin,yMax=getRange(flags,yGrid,"-yrange")
    vMin,vMax=getRange(flags,fVal,"-vrange",ratio)

    # extend the plot sizes to present
    plot.extend(np.min(fVal),np.max(fVal),np.min(xGrid),np.max(xGrid),np.min(yGrid),np.max(yGrid))

    if vMin<0. or "-linV" in flags: lev,tic=linLevels(vMin,vMax)
    else:                           lev,tic=logLevels(vMin,vMax)


    if (interactive_compare==0 and ("-compare" in flags)) or interactive_compare==1:
        if countAll==0:
            refValues=cp.deepcopy(vMatr)
            refX=xGrid
            refY=yGrid
        else:
            if len(refX)!=len(xGrid) or len(refY)!=len(yGrid):
                exit("plot dimension differ - cannot compare")

            for i in range(vMatr.shape[0]):
               for j in range(vMatr.shape[1]):
                   if refValues[i,j]!=0.: vMatr[i,j]=abs(vMatr[i,j]-refValues[i,j])/(abs(refValues[i,j])+abs(vMatr[i,j]))*2.
                   else: vMatr[i,j]=max(vMin*10.,vMax*1.e-12)

            if np.count_nonzero(vMatr)==0: print("Exact agreement")
            else:                          print("min/max error",np.min(vMatr),np.max(vMatr))
            vMin,vMax=getRange(flags,vMatr,"-vrange",ratio)
            if 0<=vMax and vMax<1.e-12:
                vMatr.fill(1.e-12)
                vMax=1.e-12
            vMin=max(vMax*1.e-5,vMin)
    else:
        vMin,vMax=getRange(flags,plot.dim[2],"-vrange=")

    if "-polar" in flags:
        axn=fax.add_subplot(gs[nPlot],projection='polar')

        # guessing coordinate meaning from file name
        if file.find("/spec")==len(file)-5:
            sys.exit("need \"Eta\" or \"Phi\" in file name for -polar, file-name is: "+file)
        theta=yGrid
        if file.find("Phi")!=-1 or (len(datfil.axisName)==2 and datfil.axisName[1].find("Eta")!=-1):
            theta=np.arccos(yGrid)
            # if integrated over Phi, correct weight for theta
            if datfil.columnTitle().find("Phi")==-1:
                for j in range(len(yGrid)): vMatr[j,:]*=sqrt(1-min(yGrid[j]*yGrid[j],1))

        # normalization and range for transformed data
        if "-normalize" in flags: fVal/=np.max(fVal)
        vMin,vMax=getRange(flags,vMatr,"-vrange",ratio)

       # forgot why we need this...
        if vMin>=0 and ("-logV" in flags or not "-linV" in flags):
            vMin=adjustLogRange(vMin,vMax)
            vMatr+=vMin

        CS,tic=plotLinOrLog(axn,theta,xGrid,np.transpose(vMatr,[1,0]),flags,vMin=vMin,vMax=vMax)
        axn.set_rmax(float(flagValue(flags,"-rmax",xMax)))
        axn.set_xticks(np.linspace(0,np.pi,7))
    else:
        if("-surface" in flags): axn=fax.add_subplot(gs[nPlot], projection='3d')
        else:                    axn=fax.add_subplot(gs[nPlot])

        if datfil.name.find("kXkY")==len(file)-4 or np.array(xGrid==yGrid).all() or "-equalAx" in flags:
            axn.set_aspect('equal', 'box')
        if doPlot: CS,tic=plotLinOrLog(axn,xGrid,yGrid,vMatr,flags,vMin=vMin,vMax=vMax)

    # set title
    title=datfil.name
    # add parameters if available
    flagval=flagValue(flags,"-label","")
    if flagval!="":
        vals=pars.itemList(flagval.split(","))
    elif pars.item("I(W/cm2)",0)!=None:
        vals=pars.itemList(["I(W/cm2)","FWHM","phiCEO"])
    else:
        vals=None

    if vals!=None:
        title+="\n"
        for n in vals: title+=n+","
        title=title[:-1]
        title+=" = "
        for n in vals:
            if vals[n]==None: title+="n/a,"
            else: title+=vals[n]+","
        title=title[:-1]
    axn.set_title(title)
    axn.set_xlabel(datfil.cols[0][0])
    axn.set_ylabel(datfil.cols[1][0])


    if doPlot:
       cbar=fax.colorbar(CS,ticks=tic,ax=axn, shrink=0.7)
       cbar.ax.set_ylabel('yield')

    if "-sum" in flags and doPlot:
        print("sum on file",file+"_sum")
        f=open(file+"_sum",'w')
        f.write("#  \n")
        for k in range(len(xVal)):
            if k>0 and xVal[k]!=xVal[k-1]: f.write("\n")
            f.write(str(xVal[k])+", "+str(yVal[k])+", "+str(fVal[k])+"\n")
        f.close()
    if "-peaks" in flags:
        HeIp=2.88738
        if doPlot: twoElectronEnergies(xGrid,pars,axn,flags,HeIp)
        lineoutAllSumEnergy(datfil,xGrid,yGrid,vMatr,HeIp)
