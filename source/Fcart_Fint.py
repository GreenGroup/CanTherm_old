#!/usr/bin/env python

"""
Copyright (c) 2002-2009 William H. Green and the CanTherm Team

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""


#take as arguments the "file name" which contains the force constant matrix and geometry.
#Also assumes that inertia.dat file is available in the present directory. From this file
#all the information for inertia is obtained.

#Once the inertia file is read the cartesian coordinates for displacement corresponding to
#hindered rotors is generated. Also external translation and rotational vectors are generated.
#Using all these vectors a projector matrix P is generated.

#the force constant matrix not containing the internal roations, external rotations and translations
#is generated Fc = (I-P)*Fc*(I-P)

#this force matrix is converted to mass-weighted-cartesian coordinates and the frequencies are 
#calculated

from numpy import *
from scipy import *
from scipy.linalg import *
import sys
import readGeomFc
file = open(sys.argv[1],'r')
inertia = open('inertia.dat','r')

(geom,Mass,Fc)=readGeomFc.readGeomFc(file)

#make the full cartesian force constant matrix
for i in range(0,3*Mass.size):
   for j in range(i,3*Mass.size):
       Fc[i,j] = Fc[j,i]

(pivots,rotAtoms,numRotors) = readGeomFc.readInertia(inertia)

intRotMatrix=matrix(array(zeros((3*Mass.size,numRotors),dtype=float)))

#form cartesian vectors for all rotors
for i in range(numRotors):
   e12=matrix('0 0 0');
   e21=matrix('0 0 0');
   e12=geom[pivots[2*i]-1,:]-geom[pivots[2*i+1]-1,:]
   e12=e12/linalg.norm(e12)
   e21=-e12
   atoms1 = rotAtoms[2*i]
   atoms2 = rotAtoms[2*i+1]
   for j in atoms1:
     e31 = geom[j-1,:]-geom[pivots[2*i]-1,:]
     #e31 = e31/linalg.norm(e31)
     intRotMatrix[3*(j-1):3*j,i]=transpose(cross(e31,e12))
   for j in atoms2:
     e42=geom[j-1,:]-geom[pivots[2*i+1]-1,:]
     intRotMatrix[3*(j-1):3*j,i]=transpose(cross(e42,e21))

#make all the modes of unit length
for i in range(numRotors):
    intRotMatrix[:,i]=intRotMatrix[:,i]/linalg.norm(intRotMatrix[:,i])

#make the int Rotors Orthonormal
intRotMatrix = matrix(orth(intRotMatrix))


#make translation and rotation unit vectors
tranrot = matrix(zeros((3*Mass.size,6),dtype=float))
for i in range(Mass.size):
   tranrot[3*i,0] = 1.0
   tranrot[3*i+1,1] = 1.0
   tranrot[3*i+2,2] = 1.0

   tranrot[3*i:3*i+3,3] = transpose(matrix([0, -geom[i,2], geom[i,1]]))
   tranrot[3*i:3*i+3,4] = transpose(matrix([geom[i,2], 0, -geom[i,0]]))
   tranrot[3*i:3*i+3,5] = transpose(matrix([-geom[i,1], geom[i,0], 0]))
   
tranrot = matrix(orth(tranrot))
P = tranrot*transpose(tranrot)
#P = intRotMatrix*transpose(intRotMatrix)
#Fc = (P)*Fc*(P)

inttranrot = matrix(zeros((3*Mass.size,6+numRotors),dtype=float))
inttranrot[:,0:numRotors]=intRotMatrix
inttranrot[:,numRotors:numRotors+6]=tranrot

inttranrot = matrix(orth(inttranrot))

#P = inttranrot*transpose(inttranrot)
I = matrix(eye(3*Mass.size,3*Mass.size))


Fc = (I-P)*Fc*(I-P)
#Fc = P*Fc*P

Tcmc = mat(zeros((3*Mass.size,3*Mass.size),dtype=float))
for i in range(Mass.size):
   for j in range(3):
      Tcmc[(i)*3+j,(i)*3+j]=1.0/sqrt(Mass[i])



Fc = Tcmc*(Fc*Tcmc)

[l,v]=linalg.eigh(Fc)

v = Tcmc*v

for i in range(3*Mass.size):
   v[:,i]=v[:,i]/linalg.norm(v[:,i])

num = Mass.size
l = sort(l)
for i in range(len(l)/3):
 for j in range(3):
  #print sqrt(l[3*i+j])*337.0/6.5463e-02,'\t',
  print sqrt(l[3*i+j] * (627.5095*4180/6.023e23)  * (1.88972e10**2) * (1/1.67e-27) )/2/math.pi/3e10,'\t',
 print 

print


