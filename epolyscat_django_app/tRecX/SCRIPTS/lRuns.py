#! /usr/bin/env python3

# tRecX = tSurff+irECS - a universal Schroedinger solver
# Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
# 
# This program is free software; you can redistribute it and/or modify it 
# under the terms of the GNU General Public License as published by the Free Software Foundation; 
# either version 2 of the License, or (at your option) any later version.
# End of license
 
"""
prints a summary of input for jobs from the same input file
output is controlled by a file that is generated at first run
"""
import os
import sys
import shutil
import string

from tools import argsAndFlags,expandRange
import runTools as rt
from runTools import extractInput,inputExtraction,printLongVals,addLinp,addToExtract,diffLinp,defaultLinpExtract


if __name__ == "__main__":

    args,flags=argsAndFlags();

    if len(args)<1:
        print("\n Input overview for runs in a base directory")
        print(" -------------------------------------------")
        print(" usage: lRuns base_dir select")
        print("        base_dir....  w/o the /0000,/0001, etc.")
        print("        select......  format 3,5,8,10-13,8 (w/o blanks)")
        print("                      42- all >=42, default: 0-9999")
        print(" flags: -C   ...print completed (default: all)")
        print("        -?   ...show only stalled runs")
        print("        -r   ...show only running or stalled")
        print("        -diff...all entries that differ in the selected runs")
        print("        -tar ...create base_dir.tar.gz to contain human-relevant outputs of selected runs")
        print("")
        print(" Example for linp-extract:")
        print("     @outf:'Rg='=>>=Rg (NOTE: blanks between the '=' will not be stripped)")
        print("     Axis:functions[2]=assocLegendre{Phi.L-M<5}=AssLeg[17:23]")
        print("     Axis:lower end[4]=20.=Box")
        print(" plots three columns with headers 'AssLeg' 'Box' 'Rg'")
        print(" values characters 17-23 of 'assocLegendre{Phi...}', here '.L-M<5', of Axis:functions[2]")
        print(" and all after first '=' and before second '=', i.e. 20. of Axis:lower end[4] ")
        print(" the run's outf will be scanned for \"Rg\" and \">>\", string in between will be printed")
        exit()

    # input file
    infile=args[0]

    # what to display
    display="all"
    if "-C"   in flags: display="completed"
    elif "-?" in flags: display="stalled"
    elif "-r" in flags: display="running"

    # full path to base directory
    if infile[0]=='/':
        root=infile
    else:
        root=os.getcwd()+"/"+infile


    # which values to extract form input and column headers
    extract,extractName,extractRange,extractOutf=inputExtraction(root);

    #header line
    head=head=" Run:"
    head+="    Master  MoniTime Running #Proc  Time    NormTotal <Expectation0>"
    for n in extractName: head+="\t"+n
    for o in extractOutf:
        head+="\t"
        if o.split("=")>2: head+=o.split("=")[-1].rstrip()
        print("outf: ",o)

    # output header line
    if "-diff" not in flags: print(" --- data selected according to "+root+"/linp-extract ---")
    headCnt=0

    # generate list of jobs for display
    jobs=range(0,9999)
    if len(args)>1:
        jobs=expandRange(args[1])

    if not os.path.isfile(root+"/linp-extract"):
        # create initial extraction file (needs editing)
        defaultLinpExtract(root,str(jobs[0]).zfill(4),flags)

    if "-diff" in flags:
        list=dict()
        for job in jobs: list=addLinp(list,root,str(job).zfill(4))
        if len(list)==0: exit("no input directories for "+root)
        addToExtract(root,diffLinp(list),False,flags)
        exit()


    tar_include=['inpc','outf','spec*','harmonics_expec','S_*/outf','Laser','linp','outspec']
    lruns_list=[]
    longVals=dict()
    noProp=[]
    for i in jobs:
        job=str(i).zfill(4)
        if os.path.isfile(root+'/'+job+'/inpc'):
            if headCnt%30==0:
                print("\033[1m"+head+"\033[0m") # the escapes make bold-face
                if "-tar" in flags: lruns_list.append(head)
                headCnt+=1
            if not os.path.isfile(root+'/'+job+"/linp"):
                print("linp-file missing - generate dummy in "+job)
                os.system('touch '+root+"/"+job+"/linp")
            line,extra,longVals=extractInput(root,job,extractOutf,longVals,extract,extractRange,flags)

            if (display=="all" and line.find("no propagation")==-1) \
            or (display=="completed" and (line[0]=="*" or line[0]=="+" or line[0]==".")) \
            or (display=="running" and (line[0]==" " or line[0]=="?")) \
            or (display=="stalled" and line[0]=="?"):
                print(line+extra)
                if "-tar" in flags:
                    lruns_list.append(line+extra)
                    for f in tar_include:
                        os.system('tar -rf '+root+'.tar '+root+'/'+job+'/'+f+ '&> tar_messages')
                headCnt+=1
            if line.split()[3].find("---")!=-1 and line[0]!="." and line[0]!="+" and line[0]!="*":
                noProp.append(i)

    printLongVals(longVals)
    if len(noProp)>0: print("W/o propagation",len(noProp)," last few:",noProp[-min(len(noProp),20):])

    if "-tar" in flags:
        lruns_list.extend(longVals)
        open("lruns_list",'w').writelines(s + '\n' for s in lruns_list)
        os.system('tar -rf '+root+'.tar lruns_list &> tar_messages')
        os.system('gzip --force '+root+'.tar')
        print(root+'.tar.gz contains '+",".join(tar_include)+' for all runs')
