# tRecX = tSurff+irECS - a universal Schroedinger solver
# Copyright (c) 2015 - 2021 by Armin Scrinzi (armin.scrinzi@lmu.de)
# 
# This program is free software; you can redistribute it and/or modify it 
# under the terms of the GNU General Public License as published by the Free Software Foundation; 
# either version 2 of the License, or (at your option) any later version.
# End of license
 
import sys
print sys.argv[1]

import os
import os.path
import string

curr = os.getcwd()         #save the current directory
path = sys.argv[3]
path = curr +"/"+ path         # change to absolute path

print "copying "+sys.argv[1]+" to "+path
os.system("cp -r "+sys.argv[1]+" "+path)
os.system("chmod -R u=xrw "+path)

os.chdir(path+"/WORK")

# create separate integral files for overlap, kinetic energy, potenitial energy and electron-electron repulsion integrals in the WORK directory.
# File names created - Ovelap, Kinetic_Energy, Potential_Energy and Electron_Repulsion

# Remove the files if they already exist
if(os.path.isfile("Overlap")):
    p=os.system("rm Overlap")
if(os.path.isfile("Kinetic_Energy")):
    p=os.system("rm Kinetic_Energy")
if(os.path.isfile("Potential Energy")):
    p=os.system("rm Potential Energy")
if(os.path.isfile("Electron_Repulsion")):
    p=os.system("rm Electron_Repulsion")

# create the required files
os.system("touch Overlap")
os.system("touch Kinetic_Energy")
os.system("touch Potential Energy")
os.system("touch Electron_Repulsion")

# Extract the required information from integrals file
if(os.path.isfile("ints")):
    p=os.system("rm ints")

# Required execution followed by the input file names
p=os.system("(echo moints ; echo moints) | $COLUMBUS/iwfmt.x > ints 2>&1")
#p=os.system("(echo moints ; echo moints) | $COLUMBUS/iwfmt.x >& ints")

# The ints file created has the info of all integrals along with long lot of extra record values, which is not neccesary for us.

f = open("ints","r")
l = f.readlines()

for i in range(len(l)):
    l[i]=l[i].split()

# an extra function for convinience
def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

# get to the starting of Overlap integrals
for i in range(len(l)):
    if(len(l[i])==3): 
        if(isFloat(l[i][0]) and isFloat(l[i][1]) and isFloat(l[i][2])):
            ov_start = i
            break

# writing the overlap matrix elements
w = open("Overlap","w")
for i in range(ov_start,len(l)):
    if(len(l[i])==3):
        print>>w, '{0:20}  {1:4}  {2:4}'.format(l[i][0], l[i][1], l[i][2])
    else:
        v_start=i+2
        break
w.close()

# writing the Potential Energy matrix elements
w = open("Potential_Energy","w")
for i in range(v_start,len(l)):
    if(len(l[i])==3):
        print>>w, '{0:20}  {1:4}  {2:4}'.format(l[i][0], l[i][1], l[i][2]) 
    else:
        kin_start=i+2
        break
w.close()

# writing the Kinetic Energy matrix elements
w = open("Kinetic_Energy","w")
for i in range(kin_start,len(l)):
    if(len(l[i])==3):
        print>>w, '{0:20}  {1:4}  {2:4}'.format(l[i][0], l[i][1], l[i][2])
    else:
        break
w.close()

# Look for starting of 2 electron integrals
for i in range(len(l)):
    if(l[i]==['filename', '2']):
        vee_start = i+2
        break

#writing electron repulsion integrals
w = open("Electron_Repulsion","w")
for i in range(vee_start,len(l)):
    if(len(l[i])==5):
        print>>w, '{0:20}  {1:4}  {2:4}  {3:4}  {4:4}'.format(l[i][0], l[i][1], l[i][2], l[i][3], l[i][4]) 
    #else:
     #   break
w.close()

# Extracting the determinant information into binary file slater.red and eivectors.red

if(os.path.isfile("eivectors.red")):
    p=os.system("rm eivectors.red")
if(os.path.isfile("slaterfile.red")):
    p=os.system("rm slaterfile.red")
if(os.path.isfile("consolidatefile")):
    p=os.system("rm consolidatefile")

if(not(os.path.isfile("cipc_out"))):
    os.system("echo 4 | $COLUMBUS/cipc.x > cipc_out 2>&1")
    #os.system("echo 4 | $COLUMBUS/cipc.x >& cipc_out")

if(os.path.isfile("consolidatefile")):
    os.system("rm consolidatefile")

tol = float(sys.argv[2])
print("Determinants: Tolerance used "+str(tol))
os.system("$COLUMBUS/civecconsolidate -t "+str(tol)+" eivectors.combined slaterfile eivectors.red slaterfile.red consolidatefile < civecconsolidate.in > civecconsolidate.out")

os.chdir(path+"/MOCOEFS")
if(os.path.isfile("mocoef_mc.sp")):
    f = open("mocoef_mc.sp","r")
    l = f.readlines()
    for i in range(len(l)):
        temp = ''
        for j in range(len(l[i])):
            if l[i][j] == 'D':
                temp+='E'
            else:
                temp+=l[i][j]
        l[i] = temp
            

    f.close()
    f = open("mocoef_mc.sp","w")
    for i in range(len(l)):
        f.write(l[i])
