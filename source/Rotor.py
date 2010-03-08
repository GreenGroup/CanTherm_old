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


from numpy import *
import pdb

class Rotor:

    pivotAtom=0
    atomsList=[]
    pivot2=0
    moments = []
    level = 0
    r = mat('0.0 0.0 0.0')
    symm = 1.0

    def __init__(self,atomsList,pivot2, level,symm,Mass):
        self.pivotAtom = atomsList[0]
        self.atomsList = atomsList
        self.pivot2 = pivot2
        self.level = level
        self.parent = 0
        self.symm = symm
        self.nonRotorList=[]
        for j in range(len(Mass)):
          if (self.atomsList.count(j+1) == 0):
             self.nonRotorList.append(j+1)
        
    def getAxes(self,geom, Mass):
        z = -(geom[self.pivot2-1,:]-geom[self.pivotAtom-1,:])
        z = z/linalg.norm(z)

        cm = matrix('0.0 0.0 0.0')

        M = 0.0
        for i in self.atomsList:
            cm = cm+Mass[i-1]*geom[i-1,:]
            M = M + Mass[i-1]
        cm = cm/M

        xtemp = (cm - geom[self.pivotAtom-1,:])
        xtemp = xtemp/linalg.norm(xtemp)

        diff = xtemp-z
        different = False

        for i in range(3):
            if not( -1e-10<(xtemp[0,i]-z[0,i])<1e-10):
                different = True
                break

        if (different):
            x = xtemp - (z*transpose(xtemp))*z
            x = x/linalg.norm(x)
            y = matrix(cross(z,x))
        else :
            xtemp = z + mat(' 0.0 0.0 1.0')
            x = xtemp - (z*transpose(xtemp))*z
            x = x/linalg.norm(x)
            y = matrix(cross(z,x))

        self.dircos = matrix(zeros((3,3),dtype=float))
        self.dircos[0,:] = x
        self.dircos[1,:] = y
        self.dircos[2,:] = z
        
    def getMoments(self, geom, Mass):
        geomTemp = geom.copy();

        r = geom[self.pivotAtom-1,:]

        self.r = r.copy()
 
        #translate so that pivot atom is at the origin
        for i in range(Mass.size):
            geomTemp[i,:] = geom[i,:] - geom[self.pivotAtom-1,:]

        #now rotate so that the axes are parallel to rotor axes
        for i in range(Mass.size):
            geomTemp[i,:]= transpose(self.dircos*transpose(geomTemp[i,:]))


        A = 0.0 #moment of inertia about the z axis
        B = 0.0 #xz cross product of ineria
        C = 0.0 #yz cross product of inertia
        Ux = 0.0 #first order x moment
        x = geomTemp[:,0]
        y = geomTemp[:,1]
        z = geomTemp[:,2]

        Uy = 0.0


        for k in self.atomsList:
            i = k-1
            A = A + Mass[i]*(x[i]**2+y[i]**2) 
            B = B + Mass[i]*x[i]*z[i]
            C = C +  Mass[i]*y[i]*z[i]
            Ux = Ux + Mass[i]*x[i]
            Uy = Uy + Mass[i]*y[i]

        #print Uy
        #print self.atomsList, Mass
        #print geom
        #exit()                
        self.moments=[A,B,C,Ux]
        #print self.moments

