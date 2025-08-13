# This Python file uses the following encoding: utf-8

# if __name__ == "__main__":
#     pass

import os
import sys
import shutil
import copy

def argsAndFlags():
    args=[]
    flags=[]
    for item in sys.argv[1:]:
        if item[0]!="-": args.append(item)
        else: flags.append(item)
    return args,flags

def copySubstitute(jobdir,infile,replaceList):
    """substitute values form replace list in copy of infile to jobdir/inpc"""
    lines=open(infile+".inp",'r').readlines()
    incf=open(jobdir+'/inpc','w')
    for l in lines:
        if l.find("#define ")==0:
            t=l.split()
            if t[1] in replaceList:
                t[2]=replaceList[t[1]]
                l=" ".join(t)+"\n"
        incf.write(l)
    incf.flush()

def newRunDirectory(exedir,jobdir,infile,inpars={}):
    """ create a new run-directory by counting up numbers"""
    if not os.path.exists(jobdir+"/inpc"):
        if infile.rfind(".inp")!=-1 and infile.rfind(".inp")==len(infile)-4:
            infile=infile[:infile.rfind(".inp")]

        # create overall directory if needed
        if not os.path.exists(infile+"/"): os.mkdir(infile)

        # make a job directory
        for i in range(10000):
            jobdir=exedir+'/'+infile+'/'+str(i).zfill(4)
            if not os.path.exists(jobdir): break
        os.mkdir(jobdir)

        copySubstitute(jobdir,infile,inpars)
    return jobdir

def parameterValues(inflags):
    """extract lists of input parameters form -with=Nam1,val1,val2,..:Nam2,val3,val3,...:..."""
    lists=[]
    for f in inflags:
        if f.find("-with=")==0:
            ranges=f.split("=")[1].split(':')
            for r in ranges:
                lists.append(r.split(','))
    return lists

def nextInLists(inpars,lists):
    """increment inpars by ranges in lists, create first if empty"""
    if len(lists)==0:
        return {}

    if len(inpars)==0:
        newpars={}
        for l in lists:
            newpars[l[0]]=l[1]
        return newpars

    newpars=copy.deepcopy(inpars)
    indx=lists[0].index(inpars[lists[0][0]])+1
    if indx<len(lists[0]):
        # next in present list
        newpar=lists[0][indx]
    elif len(lists)>1:
        # this list exhausted, reset to first, try next
        newpar=lists[0][1]
        newpars=nextInLists(inpars,lists[1:])
    else:
        return {}
    newpars[lists[0][0]]=newpar

    if len(newpars)==len(inpars):
        # update successful
        return newpars
    else:
        return {}
