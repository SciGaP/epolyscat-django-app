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

def scatterAndSizes(data,mark,label,color):
    if label.find("HF")!=-1:
        lin='-'
        labsym=""
    else:
        lin='--'
        labsym=label

    plt.plot(data[0],data[1],lin,color=color)
    w=np.array(data[2])/np.max(data[2])*100
    plt.scatter(data[0],data[1],w,marker=mark,label=labsym,color=color)

def saveAndShow(savename=None):
    # get python file name (as run from the command line)
    if savename!=None:
        name=savename
    else: 
        name=sys.argv[0].split(".py")[0]+".png"
    plt.savefig(name,bbox_inches='tight')
    plt.show(block=False)
    answer = raw_input('"u" to update ')
    plt.close()

class PlotStyle:
    def __init__(self,flags):
        self.flags=flags

    def label(self,which,default=""):
        if default!="": lab=default
        else: lab=which+" (a.u.)"
        for f in self.flags:
            if f.find("-"+which+"Label=")==0:
                lab=f.split("=")[1]
                break;
        return lab

    def labelSize(self,default=14):
        siz=default;
        for f in self.flags:
            if f.find("-labelSize=")==0:
                siz=int(f.split("=")[1])
                break;
        return siz

    def legendOff(self):
        return "-legendOff" in self.flags

    def heightRatios(self,default):
        for f in self.flags:
            if f.find("-heightRatios=")==0:
                rat=f.split("=")[0][1:-1]
                rat=rat.split(",")

                rat=[int(rat[0]),int(rat[1])]
                return rat
        else:
            return default



class Plot:
    def lineoutPosWidth(self,Flag):
        if Flag.find(",")!=-1:
            pos=float(Flag[10:Flag.find(",")])
            wid=float(Flag[Flag.find(",")+1:])
        else:
            pos=float(Flag[10:])
            wid=0.
        return pos,wid

    def setLineout(self,Flags):
        self.atX=None # lineout at X-coordinate
        self.atY=None # lineout at Y-coordinate
        for f in Flags:
            if f.find("-lineoutX=")!=-1: self.atX,self.wid=self.lineoutPosWidth(f)
            if f.find("-lineoutY=")!=-1: self.atY,self.wid=self.lineoutPosWidth(f)
#        if f.find("-eV")!=1:
#            if self.atX!=None: self.atX=sqrt(2*self.atX/27.21)
#            if self.atY!=None: self.atX=sqrt(2*self.atX/27.21)
#            print("lineout",self.atX)

    """ collect parameters and modify data for figure """
    def __init__(self,Flags):
        self.setLineout(Flags)
        self.dim=[[1.e10,-1.e10],[1.e10,-1.e10],[1.e10,-1.e10],]

    def layout(self,Nsub):
        """ layout for subplots in figure (rows/cols"""
        if Nsub<3:   nr=1
        elif Nsub<9: nr=2
        else:        nr=3
        nc=int((Nsub-1)/nr+1)

        return nr,nc

    def extend(self,fmin,fmax,xmin,xmax,ymin=None,ymax=None):
        """ extend the plot size to present """
        def extend1(xmin,xmax,dim):
            dim[0]=min(dim[0],xmin)
            dim[1]=max(dim[1],xmax)
            return dim

        self.dim[0]=extend1(xmin,xmax,self.dim[0])
        self.dim[2]=extend1(fmin,fmax,self.dim[2])
        if ymin!=None: self.dim[1]=extend1(ymin,ymax,self.dim[1])

    def isLineout(self):
        return self.atX!=None or self.atY!=None


    def lineout(Plot,Xcol,Ycol,Zcol,Flags=None):
        """ exctract line-out according to Plot """
        # axes and matrix of values
        x,y,v=xyzToAxisArray(Xcol,Ycol,Zcol)

        if Flags!=None: Plot.setLineout(Flags)

        # extract row/col nearest to atX/atY
        if Plot.atX!=None:
            colMin = (np.abs(x-Plot.atX+Plot.wid/2)).argmin()
            colMax = (np.abs(x-Plot.atX-Plot.wid/2)).argmin()+1
            return x[colMin],y,np.divide(np.sum(v[:,colMin:colMax],axis=1),colMax-colMin)
        elif Plot.atY!=None:
            rowMin = (np.abs(y-Plot.atY+Plot.wid/2)).argmin()
            rowMax = (np.abs(y-Plot.atY-Plot.wid/2)).argmin()+1
            print("sum",rowMin,rowMax,Plot.atY,y[0],y[-1])
            return y[rowMin],x,np.divide(np.sum(v[rowMin:rowMax,:],axis=0),rowMax-rowMin)
        else:
            print("specify lineout point -lineoutX=XVAL or -lineoutY=YVAL")
            exit(1)

def energyPeakPositions(E0,Omega,Emin,Emax):
    ePos =[]
    nPhot=[]
    epos=E0-Omega
    nphot=-1
    while epos<Emax-Omega:
        epos=epos+Omega
        nphot+=1
        if epos>Emin:
            ePos.append(epos)
            nPhot.append(nphot)
    return ePos,nPhot


def lineoutSumEnergy(SumEnergy,Width,xGrid,yGrid,vMatr,Flags=None):
    """ create a lineout si(xGrid) along xGrid+yGrid=SumEnergy"""
    res=[]
    dGrid=[]
    for xE in np.arange(np.min(xGrid),np.max(xGrid)+0.1,0.1):
        yE=SumEnergy-xE
        if Width<xE and xE<SumEnergy-Width and yE>0:
            dGrid.append(2*xE-SumEnergy)
            """ average over square surrounding (xE,yE)"""
            iLow=findFirstLarger(xGrid,xE-Width/2)
            iUpp=findFirstLarger(xGrid,xE+Width/2)
            jLow=findFirstLarger(yGrid,yE-Width/2)
            jUpp=findFirstLarger(yGrid,yE+Width/2)
            cnt=0
            sum=0
            for i in range(iLow,iUpp):
                for j in range(jLow,jUpp):
                    cnt+=1
                    sum+=vMatr[i,j]
            if cnt>0: res.append(sum/cnt)
            else:     res.append(0.)
    return dGrid,res

def lineoutAllSumEnergy(Datfil,xGrid,yGrid,vMatr,Ip):
    pars=RunParameters(Datfil.name[:Datfil.name.rfind('/')])
    omega,up,ip=pulseParameters(pars,Ip)
    esum,nphot=energyPeakPositions(-ip-2*up,omega,0,np.max(xGrid))
    for k in range(len(esum)):
        dgrid,espec=lineoutSumEnergy(esum[k],omega*0.5,xGrid,yGrid,vMatr,Ip)
        if len(dgrid)>0:
            f=open(Datfil.name+"_linoutN"+str(nphot[k]),'w')
            f.write("# lineout at sumEnergy="+str(esum[k])+", witdh="+str(omega*0.5)+"\n")
            for l in range(len(dgrid)):
                f.write(str(dgrid[l])+", "+str(espec[l])+"\n")
            f.close()

            tspec=np.fft.rfft(espec,norm="ortho")
            f=open(Datfil.name+"_linoutT"+str(nphot[k]),'w')
            f.write("# Fourier transform at sumEnergy="+str(esum[k])+", witdh="+str(omega*0.5)+"\n")
            dt=0.5*omega/(np.max(dgrid)-np.min(dgrid))
            for l in range(len(tspec)):
                f.write(str(l*dt)+", "+str(abs(tspec[l]))+"\n")
            f.close()



def normalize(flags,x,y,atX=None):
    for f in flags:
        if f.find("-normalize")!=-1:
           try:
               x0=float(f.split("=")[1])
               print("normalized at x=",x[np.abs(x - x0).argmin()],"is",y[np.abs(x - x0).argmin()])
               y/=y[np.abs(x - x0).argmin()]
           except:
               print("normalized at maximum",np.max(np.abs(y)))
               y/=np.max(np.abs(y))
           return
    if atX!=None:
        print("normalized at x=",x[np.abs(x - atX).argmin()],"is",y[np.abs(x - atX).argmin()])
        y/=y[np.abs(x-atX).argmin()]


# average around values around x-points
def average(flags,x,y,av=None):

    for f in flags:
        if f.find("-average=")!=-1:
            av=float(f.split("=")[1])/2.
    if av==None or av==0.: return x,y

    ya=np.zeros((len(x)))
    for k in range(len(x)):
        na=0
        for l in range(len(x)):
            if x[k]-av<x[l] and x[l]<x[k]+av:
                na+=1
                ya[k]+=y[l]
        ya[k]/=na
    return x,ya


def getRange(flags,data,command,minRatio=None):

    for flag in flags:
        if flag.find(command)==0:
            return floatRange(flag)

    up= np.max(data)
    low=np.min(data)
    if low>=0 and minRatio!=None:
        low=up*minRatio
    return low,up

def adjustLogRange(vMin,vMax):
    if vMin<=0:
        print("lowest value is"+str(vMin)+" <= 0, using default vmin="+str(vMax*1e-5)+", may specify -vrange=[vmin,vmax]")
        vMin=vMax*1e-5
    return vMin


def plotLinOrLog(axn,xGrid,yGrid,vMatr,flags,vMin=None,vMax=None):
    cmap = plt.cm.get_cmap("gnuplot")

    vmin,vmax=getRange(flags,vMatr,"-vrange=")
    yMin,yMax=getRange(flags,yGrid,"-yrange=")
    xMin,xMax=getRange(flags,xGrid,"-xrange=")
    if vMin==None: vMin=vmin
    if vMax==None: vMax=vmax

    if "-linV" in flags:
        lev,tic=linLevels(vMin,vMax)
        if "-surface" in flags: CS=axn.plot_surface(xGrid,yGrid,vMatr,lev,cmap=cmap,vmin=vMin,vmax=vMax)
        else:                       CS=axn.contourf(xGrid,yGrid,vMatr,lev,cmap=cmap,vmin=vMin,vmax=vMax)
    else:
        vMin=adjustLogRange(vMin,vMax)
        lev,tic=logLevels(vMin,vMax)
        if "-surface" in flags:
            CS=axn.plot_surface(xGrid,yGrid,np.log10(vMatr),lev,cmap=cmap,vmin=log10(vMin),vmax=log10(vMax))
        else:
            CS=axn.contourf(xGrid,yGrid,vMatr,lev,cmap=cmap,vmin=vMin,vmax=vMax,locator=ticker.LogLocator())

    axn.set_ylim(yMin,yMax)
    if "-polar" in flags:
        axn.set_rlim(xMin,xMax)
    else:
        axn.set_xlim(xMin,xMax)
    return CS,tic

def xyzToAxisArray(Xcol,Ycol,Zcol):
    """ convert from 3-column to 2dim data """
    fastFirst=Xcol[0]!=Xcol[1];

    if not fastFirst:
        for k in range(len(Xcol)):
            if(Xcol[0]!=Xcol[k]):
                yAxis=Ycol[:k]
                break
        xAxis=np.array([Xcol[k] for  k in range(0,len(Xcol),len(yAxis))])

        vArray=Zcol.reshape(len(xAxis),len(yAxis)).transpose() # numpy is row-wise
    else:
        for k in range(len(Ycol)):
            if(Ycol[0]!=Ycol[k]):
                xAxis=Xcol[:k]
                break
        yAxis=np.array([Ycol[k] for  k in range(0,len(Ycol),len(xAxis))])
        vArray=Zcol.reshape(len(yAxis),len(xAxis)) # numpy is row-wise

    return xAxis,yAxis,vArray

def momentumToEnergyEv(kX,kY,sigmaK):
    au2eV=27.211386
    for i in range(len(kX)):
        for j in range(len(kY)):
            sigmaK[i,j]*=kX[i]*kY[j]/(au2eV*au2eV)
    for k in range(len(kX)): kX[k]*=kX[k]*0.5*au2eV
    for k in range(len(kY)): kY[k]*=kY[k]*0.5*au2eV

def pulseParameters(inputPars,ionizationPot):
    # check wave-length
    lam=inputPars.allItems("lambda(nm)")
    for l in lam:
        if l!=lam[0]: print("multiple wave-length, using first: ",lam)
    wavelength=float(inputPars.item("lambda(nm)",0))

    inte=inputPars.allItems("I(W/cm2)")
    intensity=float(inte[0][0])

    omega=45.5633/wavelength
    up=intensity/(4*omega*omega)/3.50944e16
    ip=ionizationPot

    au2eV=27.211386
    omega*=au2eV
    up*=au2eV
    ip*=au2eV

    return omega,up,ip

def twoElectronEnergies(eGrid,inputPars,axn,flags,ionizationPot,photonNumber=None,colr=None):
    """
    compute peak postitions
    from wavelength, Ip, and intensity
    add lines to plot
    """

    print("Ip",ionizationPot)
    omega,up,ip=pulseParameters(inputPars,ionizationPot)
    title="om="+str(omega)+", Up="+str(up)+", ip="+str(ip)

    # draw lines at n omega -ip-2*up
    emax=np.max(eGrid)
    emin=np.min(eGrid)
    epos,nphot=energyPeakPositions(-ip-2*up,omega,emin,emax)
    actualPlot=-1
    for k in range(len(epos)):
        if photonNumber==None or nphot[k] in photonNumber:
            c="b"
            w=1
            actualPlot+=1
            if colr!=None:
                c=colr[actualPlot%len(colr)]
                w=3
            axn.plot([0,epos[k]],[epos[k],0],'--',color=c,linewidth=w)
            bbox_props = dict(boxstyle="square,pad=0.", fc="white", ec="b", lw=0, alpha=0.8)
            axn.text(epos[k]/2,epos[k]/2-10,str(nphot[k]),color=c, fontsize=14, fontweight='bold',bbox=bbox_props)

    # draw diagonals
    ediag=-omega
    while ediag<emax-omega:
        ediag+=omega
        if photonNumber==None: axn.plot([ediag,emax],[0,emax-ediag],'--',color='y')

    return title


def peakPositions(pltsize,inputPars,ax1,flags):
    """
    compute peak postitions
    from wavelength, Ip, and intensity
    add lines to plot
    """
    title=inputPars.run

    # get wave-length, intensity, ip
    name=""
    # check wave-length
    lam=inputPars.allItems("lambda(nm)")
    for l in lam:
        if l!=lam[0]: print("multiple wave-length, using first: ",lam)
    wavelength=float(inputPars.item("lambda(nm)",0))
    # add up all intensities
    inte=inputPars.allItems("I(W/cm2)")
    intensity=0
#    for i in inte: intensity+=float(i)
    intensity=float(inte[0][0])

    omega=45.5633/wavelength
    up=intensity/(4*omega*omega)/3.50944e16
    ip=0.1888

    if "-eV" in flags:
        up*=27.211
        ip*=27.211
        omega*=27.211

    # draw lines at n omega - ip - up
    epos=-ip
    title="om="+str(omega)+", Up="+str(up)+", ip="+str(ip)
    epos=epos-up
    title=title+" Up subtracted"
    plt.title(title)

    nphot=0
    emax=pltsize.dim[0][1]*pltsize.dim[0][1]*0.5
    while epos<min(emax-omega,pltsize.dim[0][1]-omega):
        epos=epos+omega
        nphot=nphot+1
        if epos>0:
            kpos=epos
            if not "-eV" in flags: kpos=sqrt(2*epos)
            ax1.plot([kpos,kpos],pltsize.dim[2],color='b')
            ax1.text(kpos,pltsize.dim[2][1]*0.5,str(nphot),rotation=90)

    return name

def linLevels(vmin,vmax):
    """
    a set of linearly spaced levels in [vmin,vmax]
    tics at 5 intervals
    """



    fact=(vmax-vmin)/16
    levs=[vmin+fact,]
    for n in range(1,17): levs.append(levs[n-1]+fact)

    # spacing of ticks neares 2-digit
    spac=(vmax-vmin)/5
    tics=[0,]
    while tics[0] >vmin:
        tics.insert(0,tics[0]-spac)
    while tics[-1]<vmax:
        tics.append(tics[-1]+spac)
    return levs,tics

def logLevels(vmin,vmax):
    """
    a set of logarithmically spaced levels in [vmin,vmax]
    ticks are returned at powers of 10
    """
    lmin=int(log10(vmin)+1)
    lmax=int(log10(vmax)+1)

    tics=[]
    fact=pow(vmax/vmin,float(1./16.))
    levs=[vmin*fact,]
    for n in range(1,17): levs.append(levs[n-1]*fact)
    for l in range(lmin,lmax): tics.append(pow(10,l))
    return levs,tics

def minGraph(x,y):
    my=-np.array(y)
    return maxGraph(x,my)


def maxGraph(x,y):
    """
    locate the maximum in a graph, by quadratic interpolation
    """
    kmax=np.argmax(y)
    if kmax==0 or kmax==len(y):
        print("WARNING: no local maximum")
        return x[kmax]
    
    xMat=np.zeros((3,3))
    for i in range(3):
        xij=1
        for j in range(3):
            xMat[i,j]=xij
            xij*=x[kmax-1+i]

    c=np.linalg.solve(xMat,y[kmax-1:kmax+2])
    
    xmax=-c[1]/(2*c[2])
    if x[kmax-1]>xmax or x[kmax+1]<xmax:
        print("parabolic location of maximum failed: ",xm,"points",x[kmax-1:kmax+2])
    
    return xmax



class Legend:
    """
    for creating labels and formating the legend
    - use "-label=PARNAME" on the command line to extract values from input (linp-file) to legend
    """
    def __init__(self,Dir,Postfix,Flags,Graphproperties=dict(),Pars=None):
        self.off=False
        self.previous=""
        self.previousDir=""
        self.dir=Dir
        self.postfix=Postfix
        self.labShort=""
        self.labName=""
        self.graphproperties=Graphproperties
        for f in Flags:
            if f.find("-label=")!=-1:
                self.labName=f[7:].strip()
                if self.labName[ 0]=="'" or self.labName[ 0]=='"':self.labName=self.labName[1:]
                if self.labName[-1]=="'" or self.labName[-1]=='"':self.labName=self.labName[:-1]
                if Pars!=None:
                    for n in self.labName.split(","):
                        self.labShort+=","+Pars.short(n,0)
                if(len(self.labShort)!=0): self.labShort=self.labShort[1:]

    def label(self,ax1,col,datfil,info,pars):
        """
        create a new label for a plot, add header info as appropriate
        - if there is a dir or postfix add legend header
        - if there is a new file or info and multiple columns, add legend header
        - put into legend label the info that is not in header
        """

        file=self.dir+info+self.postfix

        if "legend" in self.graphproperties:
            return self.graphproperties["legend"]

        parVal=""
        # if input parameters are given, add values to label
        if self.labName!="":
            for l in self.labName.split(","):
                parVal+=" "+pars.item(l,0)+","
        parVal=parVal[:-1]

        # add column names
        lab=str(datfil.cols[col][0])

        # directory or post-fix are given: create header
        curDir=self.dir+"..."+self.postfix
        if curDir.rfind("]")+1==len(curDir): curDir=curDir[:curDir.rfind("[")]

        if len(datfil.datCols)==1:
            # single column, append file name
            lab+=parVal
            if self.dir+self.postfix=="":
                lab+=" "+briefTitle(file,30,0,5)
            else:
                lab+=" *="+info
            lab+=""

        else:
            # place into intermediate legend header
            if self.dir+self.postfix!="" and curDir!=self.previousDir and self.previousDir!="":
                labHead=""
                if self.labShort!="": labHead=self.labShort+": "+curDir
                else:                 labHead=curDir
                if labHead!="": ax1.plot([1,],[1,],' ',label=tool.briefTitle(labHead,30,0,5))

            # if file name changed and not in label, force into legend
            if self.previous!=info and lab.find(file)==-1:
                self.previous=info
                subtitle=briefTitle(info.split('[')[0],20,0,5)
                ax1.plot([1,],[1,],' ',label=parVal+" .."+subtitle,alpha=1)
        return lab

