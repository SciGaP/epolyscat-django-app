#!/usr/bin/env python3

# tRecX = tSurff+irECS - a universal Schroedinger solver
# Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
# 
# This program is free software; you can redistribute it and/or modify it 
# under the terms of the GNU General Public License as published by the Free Software Foundation; 
# either version 2 of the License, or (at your option) any later version.
# End of license
 
# -*- coding: utf-8 -*-

import sys
import numpy as np
import os.path

def Sout(printItem):
    print("Sout",printItem)

def briefTitle(longTitle,maxLength,frontLength,innerLength=None):
    """ reduce long title to maxLength keeping frontLength by omitting parts inside"""
    if len(longTitle)<maxLength:
        return longTitle

    if innerLength==None: innerLength=frontLength

    parts=longTitle.split("/")
    n=len(parts)-1
    briefT=parts[-1]
    while n>0 and (len(briefT)+frontLength+3<maxLength or n>len(parts)-3):
        n=n-1
        if n==len(parts)-3 and len(briefT)+frontLength+5>=maxLength:
           briefT=parts[n][:max(innerLength,maxLength-len(briefT)-2)]+"../"+briefT
        else:
            briefT=parts[n]+"/"+briefT

    if n!=0 and frontLength>0: briefT=longTitle[:frontLength]+"../"+briefT
    return briefT

def findFirstOf(string,chars):
    """ find first position in string of any from as string of characters"""
    pos=100000
    for k in range(len(chars)):
        if string.find(chars[k])!=-1: pos=min(pos,string.find(chars[k]))
    return pos

def preAndPostFix(listOfNames):
    """seeks pre- and postfix in a list of names"""
    # no or single name
    if len(listOfNames)<2: return "","",listOfNames

    prefix=listOfNames[0]
    pstfix=listOfNames[0]
    pre=len(prefix)
    pst=len(pstfix)
    for n in listOfNames:
        while pre>len(n) or (n[ :pre]!=prefix[ :pre]): pre-=1
        while pst>len(n) or (n[-pst:]!=pstfix[-pst:]): pst-=1

    var=[]
    for n in listOfNames:
        var.append(n[pre:-pst])
    return prefix[:pre],pstfix[-pst:],var

def flagsInHelp(flags,help):
    fail=""
    for f in flags:
        f0=f
        if f0.find("=")!=-1:f0=f0[:f0.find("=")]
        for h in help:
            h0=h
            ff=findFirstOf(h," [=.,{")
            if ff!=-1: h0=h[:ff]
            if f0==h0: break
        else: fail=f
    if fail!="":
        print("undocumented flag: "+fail)
        print("allowed flags:")
        for h in help: print("  ",h)
        sys.exit(0)

def properties(proplist):
    """ converts list of strings ["p:a","q:b",etc] into a dict """
    props=dict()
    for p in proplist:
        pp=p.split(":")
        props[pp[0].strip()]=pp[1].strip()
    return props


def fileArgsAndFlags(cl_admissable):
    """
    read command-line style arguments and flags from -figure=file-name.fig (extension required)
    one argument by line
    """
    onlyFile=[
    "-font=NAME..."
    ,"-labelSize=SIZ..."
    ,"-legendOff          # do not write a legend into the graphs"
    ,"-panelRatio=[5,2]    # controls size ratio if there are two panels"
    ,"-compare             # compare first to all others"
    ,"-eV                  # use eV (for energy) instead of a.u."
    ,"-erange=[1.e-3,1]    # error range"
    ,"-xrange=[0,140]      # x-range"
    ,"-labelSize=12        # font size of axis labels"
    ,"-xLabel=Energy (eV)  # label to show (x,y,e are valid)"
    ,"-yLabel=Spectrum (arb.u.)"
    ,"-xLabel ...label on x-axis"
    ]
    legal=cl_admissable+onlyFile

    (cl_fils,cl_flags)=argsAndFlags()
    for f in cl_flags:
        if f.find("-figure=")==0:
            file=flagValue(cl_flags,"-figure=")
            break
    else:
        return cl_fils,cl_flags,dict()

    if file.find(".fig")!=len(file)-4:
        print("figure files must have extension .fig, got: "+file)
        sys.exit(0)

    fig=open(file)
    line=fig.readline()
    comments=[]
    flags=[]
    files=[]
    dfils=dict()
    while len(line.split()):
        if line[0]=="#":
            if line[1]!="-": comments.append(line[:-1])
        elif line[0]=="-":
            flags.append(line[:-1].split("#")[0].strip())
        else:
            files.append(line[:-1].split("#")[0].split()[0].strip())
            dfils[files[-1]]=properties( line[:-1].split("#")[0].split()[1:])
        line=fig.readline()

    print("\n --- arguments from file "+file+" ---\n")
    for c in comments: print(c)

    for f in flags:
        if f.find("-plotfile=")==0: break
    else:
        flags.append("-plotfile="+file[:file.find(".fig")]+".png")

    flagsInHelp(flags,legal)
    return files,flags,dfils

def argsAndFlags():
    """
    separete sys.argv into flags (argumentes starting with "-") and the rest
    """
    args=[]
    flags=[]
    for item in sys.argv[1:]:
        if item[0]=="-": flags.append(item)
        else:            args.append(item)
    return args,flags

def flagValue(flags,command,default=None):

    for flag in flags:
        if flag.find(command)==0 and flag.find("=")!=-1:
            return flag.split("=")[1]
    if default==None:
        print("specify value for flag ",command)
        sys.exit(1)

    return default


def prePost(flags):
    """ resolve flags for pre-and post-fixes by -dir=PreFix and -which=PostFix """
    for flag in flags:
        if flag.find("-dir")==0:
            print("obsolete flag",flag," use \"before\" flag -b=... instead")
        if flag.find("-dir")==0:
            print("obsolete flag",flag," use \"after\" flag -a=... instead")
    dir=""
    which=""
    for flag in flags:
        if flag.find("-b=")==0: dir=flag[3:]
        if flag.find("-a=")==0: which=flag[3:]
    return dir,which

def getch():
    """ get single character from terminal """
    import termios
    import sys, tty
    def _getch():
        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setraw(fd)
            ch = sys.stdin.read(1)
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
        return ch
    return _getch

def floatRange(Range):
    if Range.find("[")==-1 or Range.find("]")==-1:
        print("not a valid float range: ",Range,"specify as in, e.g. [1.e-4,100]")
        sys.exit(1)

    rang=Range[Range.find("[")+1:Range.find("]")].split(",")
    if len(rang)!=2:
        print("not a valid float range: ",Range,"specify as in, e.g. [1.e-4,100]")
        sys.exit(1)
    return float(rang[0]),float(rang[1])


def splitOutsideBrackets(str,sep,left,right):
    """ split string at sep, except where within left-right pairs of brackets"""
    if len(left)!=len(right):
        print("number of left brackets does not match number of right brackets")
        print(" left:",left)
        print("right:",right)

    # replace separators between brackets by PLACEHOLDER
    ph="PLACEHOLDER"
    ph=ph.replace(sep,"") # modify to avoid accidental agreement

    loc=str.find(sep)
    while loc!=-1:
        for k in range(len(left)):
            # if loc past left bracket, but not past right, hide sep by replacing with placeholder
            if str[loc:].count(left[k])<str[loc:].count(right[k]):
                str=str[:loc]+str[loc:].replace(sep,ph)
                break;
        loc=str.find(sep,loc+1)

    # split and re-insert separators
    sp=str.split(sep)
    for k in range(len(sp)): sp[k]=sp[k].replace(ph,sep)

    return sp

def rangeInSquareBrackets(Range):
    beg=0
    end=-1
    if Range.find('[')==-1: return beg,end

    if Range.find('[')==-1 or Range.find(']')==-1 or Range.find(':')==-1:
        print("not a valid range: ",Range)
        sys.exit(1)

    r=Range[Range.rfind('[')+1:].strip()
    r=r[:r.rfind(']')].strip()
    if r[0:r.find(':')]!="": beg=int(r[:r.find(':')])
    if r[r.find(':')+1:]!="":end=int(r[r.find(':')+1:])

    return beg,end

def expandRange(colRange):
    """
    expand a range string
    example: "1,2,6,8-10,3" -> [1,2,6,8,9,10,3]
    """
    jobs=[]
    subr=colRange.split(",")
    for s in subr:
        rang=s.split("-")
        if len(rang)==1: jobs.append(int(rang[0]))
        else:
            end=10000
            if rang[1]!="": end=int(int(rang[1])+1)
            jobs=jobs+ list(range(int(rang[0]),end))
    return jobs

def findFirstLarger(Range,Val):
    for i in range(len(Range)):
        if Range[i]>Val: return i
    return len(Range)

class DataFile:
    """
    structured access to a data file
    specify files and columns
    examples:
        myDataFile... read all columns, plot column 0 against all other, or, if 2d cols 0,1 against the rest
        myDataFile[3,1:7,8,4]...cols 3,1,7,8,4 from myDataFile, interprete 3,1 as x- and y-axes of a 2d plots 7,8,4
                                note: file[7][2,3:4] ... last square bracket is column info
    """
    def __init__(self,FileCols):
        self.head=None
        self.name=FileCols
        if self.name.find("[")!=-1: self.name=FileCols[:FileCols.rfind("[")]
        openFile = open(self.name)
        self.kind="gnuplot"
        self.pars=RunParameters(FileCols[:FileCols.rfind("/")])

        # collect the header lines
        line=str(openFile.readline())
        self.head=[]
        while line[0]=="#":
            self.head.append(line)
            line=str(openFile.readline())

        # find how columns are separated: ',' or whitespace
        if self.firstDataLine(openFile).find(',')!=-1: self.colSep=","
        else:                                          self.colSep=""

        # detect multiple datasets separated by blank lines
        self.isMulti=False
        linePrev=self.firstDataLine(openFile)
        while linePrev!="":
            line=str(openFile.readline())
            if linePrev.strip()=="" and line.strip()!="":
                self.isMulti=True
                break
            linePrev=line

        # determine columns to be read
        line=str(self.firstDataLine(openFile))
        self.axCols=[0,]
        if self.isMulti: self.axCols=[0,1]
        if FileCols.find('[')==-1:
            # nothing selected - read all columns
            self.datCols=range(len(self.axCols),len(self.rowSplit(line)))
        else:
            # columns given --- expand ranges and get axes (if given)
            colstr=FileCols[FileCols.rfind('[')+1:FileCols.rfind(']')]
            self.datCols=expandRange(colstr.split(':')[-1])
            if 0 in self.datCols or self.isMulti and 1 in self.datCols:
                sys.exit("data columns start at 1 for 1d and at 2 for 2d")
            if colstr.find(':')!=-1:
                self.axCols=expandRange(colstr.split(':')[0])

        # axis names
        self.axisName=["x","y","z"]
        nAxNam=0
        if len(self.head)>0:
            self.axisName=self.head[-1][1:].split(":")[0].split(",")

            # if "=", first 1 or 2 columns are for axes
            if self.head[-1].find("=")!=-1: nAxNam=len(self.axCols)
            else:                           nAxNam=0

        # number of columns
        nCols=len(self.rowSplit(self.firstDataLine(openFile)))
        
        # interprete header line as column names (various formats)
        self.colNames=self.axisName
        if len(self.head)>0:
            head=self.head[-1][1:]
            if head.find("=")!=-1:
                # multicolumn header with names after =
                self.colNames=self.colNames+head.split("=")[-1].split()
            elif head.find("sum[")!=-1:
                # names after sum
                head=head[head.find("sum[")+4:]
                self.colNames=self.colNames+head.split("]")[-1].split()
            elif(len(head.split())==nCols):
                # try blank-separated
                self.colNames=head.split()
            elif(len(head.split(','))==nCols):
                # try column-separated
                self.colNames=head.split(',')


        # number of column headers does not match number of columns - name by column number
        if len(self.colNames)!=nCols:
            self.colNames=list(range(len(self.rowSplit(line))))

        # strip brackets from column names
        for i in range(len(self.colNames)):
            n=str(self.colNames[i]).strip()
            if n.find("(") ==0:        n=n[1:]
            if n.rfind(")")==len(n)-1: n=n[:-1]
            self.colNames[i]=n

        # set up a dictionary for axis and data columns
        self.cols={}
        for c in self.axCols:  self.cols[c]=[self.colNames[c],[]]
        for c in self.datCols: self.cols[c]=[self.colNames[c],[]]

        # get axes and data columns
        self.datCols=list(self.datCols)
        self.datCols.extend(self.axCols)
        line=self.firstDataLine(openFile)
        while line!="":
            if line.strip()!="":
                data=self.rowSplit(line)
                for c in self.datCols:
                    self.cols[c][1].append(float(data[c]))
            line=openFile.readline()
        self.datCols=self.datCols[:-len(self.axCols)]

        # done, close file
        openFile.close()

    def xAxis(self):
        return np.array(self.cols[self.axCols[0]][1])

    def yAxis(self):
        return np.array(self.cols[self.axCols[1]][1])

    def column(self,col):
        return np.array(self.cols[col][1])

    def str(self):
        s=self.name+"\n"
        for key in self.cols:
            s+=str(key)+":"+str(self.cols[key][0])+"["+str(len(self.cols[key][1]))+"] "
        return s

    def colName(self,col):
        if not col in self.cols:
            print("Column",col,"not read",self.cols)
        return self.cols[col][0]


    def rowSplit(self,line):
        line=str(line)
        if self.colSep=="": return line.strip().split()
        else: return line.strip().split(self.colSep)

    def firstDataLine(self,file):
        file.seek(0)
        for n in range(len(self.head)): line=file.readline()
        return file.readline()

    def ranges(self,col):
        args,flags=argsAndFlags()

        ratio=1.e-5
        if "-linY" in flags: ratio=None
        vMin,vMax=getRange(flags,self.column(col),"-vrange",ratio)
        xMin,xMax=getRange(flags,self.xAxis(),"-xrange")
        if self.isMulti:
           yMin,yMax=getRange(flags,self.yAxis(),"-yrange")
        else:
           yMin,yMax=vMin,vMax
        return vMin,vMax,xMin,xMax,yMin,yMax

    def columnTitle(self):
        h=self.head[-1].split("=")[0]
        h=h[h.find(":"):]
        if h.find("(")==-1: return ""
        return h[h.find("("):].strip()

class RunParameters:
    """
    access to input parameters from a given tRecX "linp" file
    (the linp-file echos all inputs as actually read)
    """

    def __init__(self,RunDir):
        self.run=RunDir
        if os.path.exists(RunDir+'/linp-extract'):
            self.linp=open(RunDir+'/linp-extract','r')
        elif os.path.exists(RunDir+'/linp'):
            self.linp=open(RunDir+'/linp','r')
        else:
            self.linp=None;

    def hasParameters(self):
        return self.linp!=None

    def item(self,Name,I):
        """ string value of input item """
        if self.allItems(Name)==None: return None
        if len(self.allItems(Name)[0])==0: return None
        return self.allItems(Name)[0][I]

    def itemList(self,List):
        res={}
        for l in List:
            res[l]=self.item(l,0)
        return res

    def short(self,Name,I):
        """ short name of input item """
        return self.allItems(Name)[1][I]


    def floatItem(self,Name,I):
        """ float value of input item """
        return float(self.allItems(Name)[0][I])

    def allItems(self,Name):
        """ all input items matching Name """
        if self.linp==None:
            print("no linp-file found, cannot display paramters")
            return
        self.linp.seek(0) # position to beginning of file
        items=[]
        short=[]
        line=str(self.linp.readline())
        while line!="":
            if line.find(Name)!=-1:
                val = line.split("=")
                val[1]=val[1][:val[1].find("\inputValue{")]
                if len(val)>1: items.append(val[1].rstrip()) # remove possible trailing whitspace
                else:          items.append("")
                if len(val)>2: short.append(val[2].rstrip()) # remove possible trailing whitspace
                else:          short.append(Name)
            line=str(self.linp.readline())
        return items,short


