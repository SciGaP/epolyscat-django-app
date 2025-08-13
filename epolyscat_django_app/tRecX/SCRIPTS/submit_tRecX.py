#! /usr/bin/env python

# tRecX = tSurff+irECS - a universal Schroedinger solver
# Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
# 
# This program is free software; you can redistribute it and/or modify it 
# under the terms of the GNU General Public License as published by the Free Software Foundation; 
# either version 2 of the License, or (at your option) any later version.
# End of license
 
import os
import sys
import shutil
import copy

import submit_tools as st

def resume(base,runFlags):
    global useTim

    # name for the queueing system
    name=base.split('/')[-2]
    numb=base.split('/')[-1]
    name=name[0:min(5,len(name))]+numb[-min(3,len(numb)):]

    # KCS (SLURM)
    if 'kcs-login'==os.environ['HOSTNAME'] or os.environ['HOSTNAME']=='ASC':
        maxMem=360000
        maxTim=72
        # guessing memory on ASC
        if os.environ['HOSTNAME']=='ASC':
            maxMem=128000
            maxTim="10"

        if useTim=="": useTim="72"

        if int(nproc>31):
            mem=maxMem
        elif useMem=="":
            mem=min(8000*int(nproc),maxMem)
        elif useMem=="all":
            mem=maxMem
        else:
            mem=useMem

        commands=\
            "#!/bin/bash                                    \n"+\
            "#SBATCH -o "+base+"/out                        \n"+\
            "#SBATCH -e "+base+"/err                        \n"+\
            "#SBATCH -D "+exedir+"                          \n"+\
            "#SBATCH -J "+name+"                            \n"+\
            "#SBATCH --export=NONE                          \n"+\
            "#SBATCH --get-user-env                         \n"+\
            "#SBATCH --mail-type=end                        \n"+\
            "#SBATCH --mail-user=armin.scrinzi@lmu.de       \n"+\
            "#SBATCH --time="+useTim+":00:00                \n"+\
            "#SBATCH --mem="+str(mem)+"                     \n"+\
            "#SBATCH --nodes="+str(nnode)+"                 \n"+\
            "#SBATCH --ntasks-per-node="+str(ntask)+"       \n"


        # designate compute nodes
        if useNodes!="": commands=commands+"#SBATCH --nodelist="+useNodes+"                   \n"

        if os.environ['HOSTNAME']!='ASC':
            commands+=\
            "module unload intel                            \n"+\
	    "module unload mpi.intel                        \n"+\
	    "module load intel/18.0                         \n"+\
	    "module load mpi.intel/2018                     \n"+\
            ""
        if nproc<1: commands+="mpirun -genv I_MPI_DEVICE shm "
        else:       commands+="mpiexec -np "+str(nproc)+" "
        commands+=exedir+"/tRecX "+base+"/inpc  "+runFlags+"\n"
        #commands+=exedir+"/a.out "+base+"/inpc  "+runFlags+"\n"

        print commands
        open(base+'/script','w').write(commands)
        os.system('sbatch '+base+'/script')


    # CAPP (SLURM)
    elif 'capp1'==os.environ['HOSTNAME'] or 'sulamith'==os.environ['HOSTNAME']:
        if useMem=="":
            if queue=="capp_medium": mem=8000*int(nproc)
            else              : mem=6000*int(nproc)
        elif useMem=="all":
            mem=128000
        else:
            mem=useMem

        print mem
        commands=\
            "#!/bin/bash                                    \n"+\
            "#SBATCH -p "+queue+"                           \n"+\
            "#SBATCH -o "+base+"/out                        \n"+\
            "#SBATCH -e "+base+"/err                        \n"+\
            "#SBATCH -D "+exedir+"                          \n"+\
            "#SBATCH -J "+name+"                            \n"+\
            "#SBATCH --get-user-env                         \n"+\
            "#SBATCH --mail-type=end                        \n"+\
            "#SBATCH --mail-user=armin.scrinzi@lmu.de       \n"+\
            "#SBATCH --export=NONE                          \n"+\
            "#SBATCH --time=40:00:00                        \n"+\
            "#SBATCH --ntasks="+nproc+"                     \n"

        # designate compute nodes
        if useNodes!="": commands=commands+"#SBATCH --nodelist="+useNodes+"                   \n"

        if  int(nproc)<17: commands+="#SBATCH --mem="+str(mem)+"\n"

        # SLURM is broken with mpirun - the following is a workaround
        #if nproc!="1": commands+="mpirun "+exedir+"/tRecX "+base+"/inpc  "+runFlags+"\n"
        #else:          commands+=exedir+"/tRecX "+base+"/inpc "+runFlags+"\n"
        commands+="mpirun -genv I_MPI_DEVICE shm "+exedir+"/tRecX "+base+"/inpc  "+runFlags+"\n"

        open(base+'/script','w').write(commands)
        os.system('sbatch '+base+'/script')

    # LRZ linux cluster (SGE)
    elif 'LRZlinux'==os.environ['HOST']:
        commands=\
            "#!/bin/bash                 \n"+\
            "#$ -M armin.scrinzi@lmu.de  \n"+\
            "#$ -S /bin/bash             \n"+\
            "#$ -N "+name+"              \n"+\
            "#$ -o "+base+"out           \n"+\
            "#$ -e "+base+"err           \n"+\
            "#$ -l mf=1000M,h_rt=10:00:0 \n"+\
            "#$ -l march=x86_64          \n"+\
            ". /etc/profile              \n"+\
            "cd "+os.getcwd()+"          \n"+\
            exedir+"/tsurff "+base+"/inpc -noScreen \n"
        open(base+'/script','w').write(commands)
        os.system('qsub '+base+'/script')
        
    else:
        exit('undefined system: '+os.environ['HOST'])+' or '+os.environ['HOSTNAME']


# separate args from possible flags
inargs,inflags=st.argsAndFlags()

if len(inargs)<1:
    print "\n   usage: submit_tRecX input nproc queue   [defaults: nproc=1, queue=medium]"
    print "            -nodes=capp[28-30,32]  to select given nodes "
    print "            -mem=all | 1234          to set memory (default = 8000/proc, all=128GB) "
    print "            -time=24                 run-time (hours), default on kcs=48"
    print "            -with=NAM1,v1,v2,...:NAM2,v3,v4,...   multiple submit loop, substitute for #define NAMi val"
    exit()

# input file
infile=inargs[0]

# number of tasks
nproc=1
if len(inargs)>1:
    nproc=int(inargs[1])

nnode=1
ntask=nproc
if nproc>32:
    nnode=(nproc-1)/32+1
    ntask=min((nproc-1)/nnode+1,32)

# run time
useTim=""
for f in inflags:
    if f.find("-time=")!=-1:
        useTim=f.split('=')[1]

# queue
queue="capp_medium"
if len(inargs)>2:
    queue="capp_"+inargs[2]

useNodes=""
for f in inflags:
    if f.find("-nodes=")!=-1:
        useNodes=f.split('=')[1]

useMem=""
for f in inflags:
    if f.find("-mem=")!=-1:
        useMem=f.split('=')[1]

runFlags=""
for f in inflags:
    if f.find("-mem=")==-1 and f.find("-nodes=")==-1 and f.find("-time=")==-1 and f.find("-with=")==-1:
       runFlags+=" "+f

# exedir
exedir=os.getcwd()

# job directory
jobdir=infile
if jobdir.rfind('/inpc')!=-1 and infile.rfind('/inpc')==len(infile)-5:
    jobdir=infile[:infile.rfind('/inpc')]

inpars=st.nextInLists({},st.parameterValues(inflags))
if inpars!={}:
    while len(inpars)>0:
        nextdir=st.newRunDirectory(exedir,jobdir,infile,inpars)
        print "Submitted to",nextdir,"with parameters",inpars
        resume(nextdir,runFlags)
        inpars=st.nextInLists(inpars,st.parameterValues(inflags))
else:
    jobdir=st.newRunDirectory(exedir,jobdir,infile)
    # "resume", i.e. run from job directory
    resume(jobdir,runFlags)

