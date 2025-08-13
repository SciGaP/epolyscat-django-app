# This Python file uses the following encoding: utf-8

# if __name__ == "__main__":
#     pass

import os
import sys
import shutil
import copy
import string

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import tools
import datetime

def addLinp(list,root,job):
    """
    add a linp-file to list of linps
    """

    if not os.path.isfile(root+'/'+job+'/inpc'): return list

    try:
        inp=open(root+'/'+job+"/linp",'r').readlines()
    except:
        print("linp-file missing in",root+'/'+job)
        return list

    if len(list)==0:
        list=dict()
        list["Root"]=[]
        list["Run"]=[]
    list["Root"].append(root)
    list["Run"].append(job)

    for lin in inp:
        vals=tools.splitOutsideBrackets(lin,"=","[({<","])}>")
        valExt=valFromLine(lin)
        if vals[0] not in list: list[vals[0]]=[] # new input item
        # supplement with n/a for previous runs
        while len(list[vals[0]])<len(list["Run"])-1: list[vals[0]].append("n/a")
        list[vals[0]].append(valExt)

    return list

def diffLinp(list,chooseAll=False):
    """
    show only list entries that differ
    """

    printList=[list["Run"],]
    printList[0].insert(0,"Run")

    for key in sorted(list.keys()):
        if(key=="Run"): continue # Run is first line
        val0=None
        for val in list[key]:
            if val0==None: val0=val
            if val0!=val or chooseAll:
                printList.append(list[key])
                printList[-1].insert(0,key)
                break

    longList=dict() # list long entries at end of table
    printLen = [0,]*len(printList[0])
    for l in printList:
        for k in range(len(l)):
            if k!=0 and len(l[k])>20:
                # replace with footnote entry
                if l[k] not in longList: longList[l[k]]="[-"+str(len(longList))+"-]"
                l[k]=longList[l[k]]
            # get the column width
            printLen[k]=max(len(l[k]),printLen[k])

    totalLen=0
    for n in printLen: totalLen+=n+2

    #HACK: need to disentangle generation of list and pring
    if chooseAll: return printList

    print("-"*totalLen)
    for l in printList:
        for k in range(len(l)):
            fstr='{0:'+str(printLen[k])+'} '
            print(fstr.format(l[k]),end='')
        if l[0]=="Run": print("\n","-"*totalLen)
        else:           print("")

    if len(longList)>0: print("."*totalLen)
    printLongVals(longList)

    return printList

def stalled(clock,term):
    if term=="+": return clock,term
    if clock=="--.--":return clock,"?"
    now=datetime.datetime.now()
    mins=int(clock.split(':')[0])*60+int(clock.split(':')[1])
    if now.hour*60+now.minute<mins+10: return clock,term
    return clock,"?"

def addToExtract(root,printList,add,flags):
    """
    determine all variables that are not in linp-extract, optionally append to linp-extract
    """

    if not os.path.exists(root):
        exit("no such directory: "+root)

    if os.path.isfile(root+"/linp-extract"):
        extractLines=open(root+"/linp-extract",'r').readlines()
    else:
        extractLines=[]

    addPar=[]
    for f in flags:
        if f.find("-addToExtract=")==0:
            addPar=f.split("=")[1].split(",")

    newPar=[]
    for l in printList:
        if l[0]=="Run" or l[0]=="Root" or l[0].find("---")!=-1: continue
        short=shortName(l[0])
        for x in extractLines:
            if x.find(l[0])==0: break
        else:
            if addPar==[] or short in addPar:
                # add all new or thosespecified on input
                newPar.append(short)
                extractLines.append(l[0]+"= ="+shortName(l[0])+"\n")
    if newPar!=[]:
        if "-addToExtract" in flags or add or addPar!=[]:
            for l in extractLines:
                if l[0]=="#": extractLines.remove(l)
            open(root+"/linp-extract",'w').writelines(extractLines)
        else:
            print(" --- add new parameters by lRuns.py with -diff -addToExtract="+",".join(newPar)+" ---")

def defaultLinpExtract(root,job,flags):
    """ create a default for lin-extract """
    list=dict()
    list=addLinp(list,root,job)
    difl=diffLinp(list,True)
    addToExtract(root,difl,True,flags)

    # rewrite linp-extract file
    select=['Run',
    'Absorption:upper[1]',
    'Axis:nCoefficients[1]',
    'Axis:nCoefficients[2]',
    'Axis:nCoefficients[3]',
    'Operator:hamiltonian[1]'
    ]

    lines=open(root+"/linp-extract",'r').readlines()
    addl=[]
    for s in select:
        for l in lines:
            if l.find(s+"=")==0: addl.append(l)
    for l in lines:
        addl.append("#"+l)

    f=open(root+"/linp-extract",'w')
    f.writelines(addl)

    f.close()

    print("\n*************************************************************")
    print("...GENERATED DEFAULT EXTRACTION LIST (linp-extract) ...")
    print("\n         --- re-run to see results ----")
    print("\nInstructions for extending the list of parameters")
    print("- run with option lRuns.py ... -diff -addToExtract")
    print("- edit the generated file as needed")
    print("- repeat this with more files, new differences can be added")
    print("*************************************************************")
    exit(0)

def getMonitor(rundir):
    mon=dict()
    try:
        f=open(rundir+"/mon",'r')
        for l in f.readlines():
            if l.find("Stat: ")==0:
                l=l[:l.find('(')]+l[l.find(')')+1:]
                its=l.split()
                mon["date"]=its[1]
                mon["clock"]=its[2]
                mon["param"]=its[3]
                mon["info"]=" ".join(its[4:])
                break
    except:
        pass

    return mon




def inputExtraction(Root):
    """
    get entry names from "linp-extract" file
    get column header names where these are defined by an extra =colName
    """
    entr=[]
    head=[]
    rang=[]
    outf=[]
    if os.path.isfile(Root+"/linp-extract"):
        lines=open(Root+"/linp-extract",'r').readlines()
        for l in lines:

            if l.find("DO NOT EDIT")!=-1:
                print("remove leftover line from file:\n",Root+"/linp-extract")
                print ('"'+l+'"')
                exit(1)

            if l[0]=='#': continue

            if l[0]=='@':
                # substrings from outf
                outf.append(l[l.find(':')+1:-1]+l[:l.find(':')])
                continue

            info=tools.splitOutsideBrackets(l,"=","[({<","])}>")
            entr.append(info[0])
            rang.append([0,-1])
            if len(info)==3:
                rang[-1][0],rang[-1][1]=tools.rangeInSquareBrackets(info[2])
                head.append(info[2][:info[2].find('[')].strip())
            else:
                head.append('... ')
    return entr,head,rang,outf

# extract according to linp-extract (generate linp-extract if not present)
def extractInput(root,job,extractOutf,longVals,extract,extractRange,flags):

    if not os.path.isfile(root+"/linp-extract"):
        # create initial extraction file (needs editing)
        defaultLinpExtract(root,job,flags)

    # extract from input file
    inp=open(root+'/'+job+"/linp",'r').readlines()

    extInp,extExtra=extractOutput(root,job,extractOutf)

    k=0
    for name in extract:
        for val in inp:
            if name==val.split('=')[0]:
                vals=tools.splitOutsideBrackets(val,"=","{[(<","}])>")

                val1="=".join(vals[1:])
                if val1.strip().find("e+308")==-1: valExt=val1.strip()
                else: valExt="Infty"

                valExt=valFromLine(val)
                if len(valExt)>15:
                   if valExt not in longVals:
                       longVals[valExt]="[-"+str(len(longVals))+"-]"
                   valExt=longVals[valExt]

                # we can define a sub-range of the input string
                if extractRange[k][1]==-1: extInp+="\t"+valExt[extractRange[k][0]:]
                else:                      extInp+="\t"+valExt[extractRange[k][0]:extractRange[k][-1]]
                break;
        else:
            extInp+="\tn/a"
        k+=1
    return extInp,extExtra,longVals

def extractOutput(root,job,extractOutf):
    """
    extract information from outf,out,outana files in root/job/
    """

    mon=getMonitor(root+"/"+job)

    runDir=os.listdir(root+"/"+job)

    if not "outf" in runDir:
        return " "+job+":       --- no output file ---",""

    outlines=open(root+'/'+job+"/outf",'r').readlines()

    # check for regular termination
    term=" "

    if term==" " and "err" in runDir:
        errlines=open(root+'/'+job+"/out",'r').readlines()
        for lin in errlines:
            if lin.find("Killed")!=-1:
                term="E"

    if term==" " and "out" in runDir:
         errlines=open(root+'/'+job+"/out",'r').readlines()
         for lin in errlines:
             if lin.find("Killed")!=-1: term="E"

    # new style spectra
    if term==" " and "outspec" in runDir:
        term="a"
        speclines=open(root+'/'+job+"/outspec",'r').readlines()
        for lin in speclines[-5:]:
            if lin.find(" done ")!=-1: term='+'

    # old style spectra
    if term==" " and "anaout" in runDir:
        term="a"
        analines=open(root+'/'+job+"/anaout",'r').readlines()
        for lin in analines[-5:]:
            if lin.find(" done ")!=-1: term='+'

    # main run ended
    if term==" " and "outf" in runDir:
        for lin in outlines[-5:]:
            if lin.find(" done ")!=-1: term='*'

    # count regions in progress
    if term==" " or term=="*" or term=="a":
        subD=0
        for f in runDir:
            if f.find("S_")==0:
                 subD+=1
        if subD==1: term="."
        if subD==2: term=":"
        if subD>2:  term="|"

    # find master host name
    host="(no info)"
#<<<<<<< HEAD
    for k in range(len(outlines)):
        if outlines[k].find("Master host = ")!=-1: host=outlines[k].strip().split()[-1];
    host=(10*" "+host.strip()+" ")[-10:]
#=======
#    for lin in outlines:
#        if lin.find("Master host = ")!=-1: host=lin.strip().split()[3];
#    host=(10*" "+host.strip()+" ")[-11:]
#>>>>>>> 1a4de7662a2ed804d1ded690451359cd2e04cb9e

    try:
        clock=mon["clock"]+" "
        info =mon["info"]
        if clock=="":
            if len(info)>8: clock=inf[:7]+'..'
            else: clock=info
            if clock=="": clock="  -.-    "

            clock,term=stalled(clock,term)
    except:
        info=""
        clock="--n/a-- "
        term="-"

    # find time-propagation header
    proc="?"
    for k in range(len(outlines)):
        if outlines[k].find("running parallel with")!=-1: proc=outlines[k].strip().split()[4];
        if outlines[k].find("CPU")!=-1 and outlines[k].find("(%)")!=-1: break
    else:
        if info!="":
            return term+job+": "+host+clock+"--- "+info[:44].ljust(44),""
        else:
            return term+job+": "+host+clock+"--- computation not started ---".center(48),""

    # find progress printout
    last=" --- --- --- --- ---"
    kLast=len(outlines)-1
    for l in range(k+1,len(outlines)):
        # accept/reject - finished
        if outlines[l].find("accept/reject")!=-1: kLast=l;break
    while not isRunLine(outlines[kLast]): kLast-=1
    else: last=outlines[kLast]

    # compose info into single line:
    s=term+job+": "+host+clock                  # status (' '/'*') + job number + host
    s+=(last.split()[0]+" "*10)[:9]       # current running time
    s+=proc.rjust(5)                      # number of processe, right justified
    s+=' '+(last.split()[2].strip()+" "*20)[:7]  # time-propagation time
    s+=' '+(last.split()[3].strip()+" "*20)[:10] # standard norm
    if len(last.split())>5:
        s+=' '+(last.split()[5].strip()+" "*20)[:14] # 0th expectation value (may not always be there?)
    else:
        s+="     ----     "

    extra=parseFile(root+'/'+job,extractOutf)

    return s,extra

def printLongVals(longList):
    linebreak=150
    indent=5
    for key, value in sorted(longList.items()):
        if len(key)<linebreak:
            print(" "*indent,value+" = "+key)
        else:
            print(" "*indent,value+" = "+key[:linebreak]+"\n"+" "*(indent+10
            )+key[linebreak:])


def shortName(linpEntry):
    """
    construct a short name from linpEntry
    """
    s=linpEntry
    line=s[s.rfind("[")+1:s.rfind("]")]

    # remove all brackets
    for b in '{}[]()<>/': s=s.replace(b,'')

    short=s.strip()[:3]
    short+=s[s.find(":")+1:s.find(":")+2].upper()
    short+=s[s.find(":")+2:s.find(":")+4]
    short+=line
    return short

def delayParameters(files):
    """ extract parameters related to RABITT delays from spec files"""
    phim=[]
    dlay=[]
    alig=[]
    yiel=[]
    for f in files:
        lines=open(f,'r').readlines()
        for l in lines:
            if l.find("maxPhi=")!=-1:
                part=l.split(',')
                for p in part:
                    if p.find("delay=") !=-1:  dlay.append(float(p.split("=")[1]))
                    if p.find("maxPhi=")!=-1: phim.append(float(p.split("=")[1]))
                    if p.find("align=") !=-1:  alig.append(float(p.split("=")[1])/3.14159127*180)
                    if p.find("yield=") !=-1:  yiel.append(float(p.split("=")[1]))
    return alig,dlay,yiel,phim



def parseFile(jobDir,extractOutf):
    '''parse a jobDir/file and return strings between markers'''

    extra=''
    for codeOut in extractOutf:
        file=jobDir+"/"+codeOut[codeOut.find('@')+1:]
        if not os.path.exists(file):
            extra+='\t??'
            continue
        outlines=open(file,'r').readlines()
        for out in outlines:
            beg=out.find(codeOut.split('=')[0])
            end=out.find(codeOut.split('=')[1],beg+len(codeOut.split('=')[0]))
            if beg!=-1 and end!=-1:
                extra+='\t'+out[beg+len(codeOut.split('=')[0]):end]
                if len(codeOut.split('='))>2 and codeOut.split('=')[2]:
                    rang=rangeInSquareBrackets(codeOut.split('=')[2])
                    extra=extra[rang[0]:rang[1]]
                break
        else:
            extra+='\t??'

    return extra


def isRunLine(l):
    """
    True if line contains run-number
    """
    # last line criteria: all numbers, at least 4
    if len(l.split())<4: return False
    for s in l.split():
        try: float(s)
        except ValueError: return False
    return True

def valFromLine(lin):
    vals=tools.splitOutsideBrackets(lin,"=","[({<","])}>")
    val1="=".join(vals[1:])
    if val1.find("\\")!=-1: val1=val1[:val1.find("\\")]
    if val1.find("(mutable)")!=-1:val1=val1[:val1.find("(mutable)")]
    if val1.strip().find("e+308")==-1: valExt=val1.strip()
    else: valExt="Infty"
    return valExt


