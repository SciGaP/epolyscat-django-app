#! /usr/bin/env python

# tRecX = tSurff+irECS - a universal Schroedinger solver
# Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
# 
# This program is free software; you can redistribute it and/or modify it 
# under the terms of the GNU General Public License as published by the Free Software Foundation; 
# either version 2 of the License, or (at your option) any later version.
# End of license
 
if len(sys.argv)<4:
    print "bin spectral values to produce spectral density, usage:"
    print "\npolar.py file bin emax"
    sys.exit(0)

# get python file name (as run from the command line)
name=sys.argv[0].split(".py")[0]

# base directory to use
file=sys.argv[1]

# binning
DeltaE=float(sys.argv[2])

# plot range
Emax=2.
if len(sys.argv)==4: Emax=float(sys.argv[3])

Thresh=1.e-4

#======================================================================
# end of input
#======================================================================

#===========================================================
# set up plot panels
#==========================================================

# plot colors for easier comparison
colr=['red','green','blue','cyan','magenta','brown']



# read lines for run
f = open(file)
lines = f.readlines()
f.close()

eb=np.array([-0.5,-0.125,-0.0555,0])
sb=np.zeros(4)
ec=np.array([DeltaE*(i+0.5) for i in range(int(Emax/DeltaE))])
sc=np.zeros(len(ec))

for l in lines:
    if l[0]!="#":
        vals=l.split()
        e=float(vals[0])
        p=float(vals[2])**2+float(vals[3])**2
        if e<-0.499: sb[0]+=p
        elif e<-0.149: sb[1]+=p
        elif e<-0.05:  sb[2]+=p
        elif e<0.0:  sb[3]+=p
        elif int((e+0.03125)/DeltaE)<len(sc): sc[int((e+0.03125)/DeltaE)]+=p

# remove values below threshold
rm=[]
for k in range(len(sc)):
    if sc[k]<Thresh: rm.append(k)
print rm
sc=np.delete(sc,rm)
ec=np.delete(ec,rm)

#===========================================================
# finalize the plots
#===========================================================
#plt.yscale('log')
#plt.plot(eb,sb,ms=8,marker='s',mfc='b',linestyle='none')
#plt.plot(ec,sc,ms=8,marker='*',mfc='r',linestyle='none')
plt.axis('off')
plt.bar(eb,np.log10(sb/Thresh),width=0.05)
plt.bar(ec,np.log10(sc/Thresh),width=0.95*DeltaE,color='r')


axes = plt.gca()
axes.set_xlim([-0.7,Emax])


plt.savefig(name+".png")
plt.show(block=False)
answer = raw_input("plot on "+name+".png\n<return> to finish")

