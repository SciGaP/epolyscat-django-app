#! /usr/bin/env python

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
from tools import *
from plotTools import *

class Compare:
    def __init__(self,linY,style):
        self.tolerance=1.e-4
        self.errRange=1.e-5  # plot in [maxErr*errRange,maxErr]
        self.plotRange=1.e-7 # plot in [max*plotRange,max]
        self.linY=linY
        self.args=argsAndFlags()[0]
        self.flags=style.flags
        self.style=style
        self.fax,(self.ax1,self.ax2)=plt.subplots(nrows=2,sharex=True,num=1)

        gs = gridspec.GridSpec(2,1,height_ratios=style.heightRatios([3,5]))
        self.ax1=plt.subplot(gs[0])
        self.ax2=plt.subplot(gs[1])

        # axis labels
        self.ax1.set_ylabel(style.label("y"),fontsize=style.labelSize())
        self.ax2.set_xlabel(style.label("x"),fontsize=style.labelSize())
        self.ax2.set_ylabel("Relative error",fontsize=style.labelSize())


        self.fmin= 1e10
        self.fmax=-1e10
        self.errmax=-1
        self.errmin=1

    def addPlot(self,xVal,fVal,lab,colr):
        self.fmin=min(self.fmin,np.min(fVal))
        self.fmax=max(self.fmax,np.max(fVal))

        if self.style.legendOff(): lab=None
        if self.linY: self.ax1.plot(xVal,fVal,label=lab,color=colr)
        else:         self.ax1.semilogy(xVal,fVal,label=lab,color=colr)

    def reference(self,xVal,fVal,lab,colr):
        self.xref=xVal
        self.fref=fVal

        fmax=max(self.fmax,np.max(fVal))
        [self.xmin,self.xmax]=[np.min(xVal),np.max(xVal)]
        self.addPlot(xVal,fVal,lab,'black')

    def compare(self,xVal,fVal,lab,colr,ratio):
        f1=self.fref
        fx=griddata(xVal,fVal,(self.xref,),method='linear')
        #plt.suptitle(dir+"*"+which,fontsize=20)
        self.addPlot(self.xref,fx,"",colr)

        self.fmin=min(self.fmin,np.min(fx))
        self.fmax=max(self.fmax,np.max(fx))
        if np.max(abs(fx-self.fref)) !=0:
            ferr=[]
            xerr=[]
            for n in range(len(fx)):
                if ratio:
                    if f1[n]!=0:
                        ferr.append(fx[n]/f1[n])
                        xerr.append(self.xref[n])
                else:
                    if (fx[n]!=0 or f1[n]!=0) and fx[n]-self.fref[n]!=0:
                        ferr.append(abs(fx[n]-self.fref[n])/abs(0.5*(fx[n]+f1[n])))
                        xerr.append(self.xref[n])

            if self.style.legendOff(): lab=None
            if ratio:
                self.ax2.plot(xerr,ferr,label=lab,color=colr)
            else:
                self.ax2.semilogy(xerr,ferr,label=lab,color=colr)

            self.errmin=min(self.errmin,np.min(ferr))
            self.errmax=max(self.errmax,np.max(ferr))

            # rms deviations
            relErr=(fx-self.fref)/(self.fref+self.tolerance*np.max(self.fref))
            errRMS=np.sqrt( np.sum(relErr*relErr) * (xVal[1]-xVal[0]) / (np.max(xVal)-np.min(xVal)) )
            print("errRMS:",errRMS,"\t",self.tolerance,"\t",np.max(self.fref))
        else: print("identical")

    def plot(self):
        # remove ticks in pupper panel
        self.fax.subplots_adjust(hspace=0.1)
        plt.setp([a.get_xticklabels() for a in self.fax.axes[:-1]], visible=False)

        xrang=getRange(self.flags,[np.min(self.xref),np.max(self.xref)],"-xrange")
        print (xrang,self.flags)
        self.ax1.set_xlim(xrang)
        self.ax2.set_xlim(xrang)

        errmin=self.errmin
        errmax=self.errmax
        print("Error range: ",errmin,errmax,self.errRange,"set range by -erange=[eMin,eMax]")
        errmin=max(errmax*self.errRange,errmin)
        errmin,errmax=getRange(self.flags,[errmin,errmax],"-erange",self.errRange)
        self.ax2.set_ylim([errmin,errmax])

        fmin=self.fmin
        fmax=self.fmax
        fmin,fmax=getRange(self.flags,[fmin,fmax],"-vrange")
        self.ax1.set_ylim([fmin,fmax])

        #finalize
        self.ax1.legend().draw_frame(False)
        if self.ax2.legend()!=None:
            self.ax2.legend().draw_frame(False)


