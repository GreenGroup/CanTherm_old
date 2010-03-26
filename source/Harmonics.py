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
import readGeomFc
import math
import Gnuplot,Gnuplot.funcutils

class Harmonics:
    numFit = 0
    Kcos = []
    Ksin = []
    A = 0.0
    Barrier_Height = 0.0

    def __init__(self,numFit, Kcos, Ksin):
        self.numFit = numFit
        self.Kcos = Kcos
        self.Ksin = Ksin

    def getPotential(self, angle):
#assuming the angle is in degrees
        rad = angle*math.pi/180
        pot = self.A
        for i in range(self.numFit):
           pot = pot + self.Ksin[i]*math.sin((i+1)*rad) + self.Kcos[i]*math.cos((i+1)*rad)
        return pot

    def fitPotential(self,file):
        read = open(file,'r')
        lines = read.readlines()
        angles = []
        potentials = []
        pot0 = 0.0
        nfit = len(lines)-4 #first three lines are comments and last time is repeat
        potgiven = []
        for i in range(3,len(lines)-1):
           tokens = lines[i].split()

           potentials.append(float(tokens[1]))
           if i == 3:
              pot0 = potentials[0]
           potentials[i-3] = (potentials[i-3]-pot0)*627.5095
           potgiven.append([float(i-3)*360/nfit, potentials[i-3]])

        #now fit the potentials
#        Y = transpose(matrix(potentials[:nfit]))
#        X = matrix(zeros((nfit,11),dtype=float))
# MRH 28Jan2010: comment out previous two lines, added next 2 lines
  #      Y = matrix(zeros((nfit+1,1),dtype=float))
  #      X = matrix(zeros((nfit+1,11),dtype=float))
# END MRH 28Jan2010
# MRH 2Feb2010: V(phi=0) ~= 0 necessarily
#     CFG suggests eliminating the A term so MRH is implementing this to test
        Y = matrix(zeros((nfit+1,1),dtype=float))
        X = matrix(zeros((nfit+1,10),dtype=float))
        for i in range(nfit):
           #MRH 28Jan2010: next line added
           Y[i,0]=potentials[i]
           angle = float(i)*2*math.pi/nfit
           #MRH commented out next line 2Feb2010
           #X[i,0]=1.0
           for j in range(5):
              #MRH commented out next 2 lines 2Feb2010
              #X[i,j+1] = math.cos((j+1)*angle)
              #X[i,6+j] = math.sin((j+1)*angle)
              X[i,j] = math.cos((j+1)*angle)
              X[i,j+5] = math.sin((j+1)*angle)

        #MRH 28Jan2010: Adding code to enforce dV/dphi = 0 @ phi=0
        Y[len(Y)-1] = 0
        for j in range(5):
           #MRH changed j+6 to j+5 on 2Feb2010
           X[len(Y)-1,j+5] = j+1
        #END MRH 28Jan2010

        XtX = transpose(X)*X
        XtY = transpose(X)*Y
        b = linalg.inv(XtX)*XtY
        
        for i in range(5):
           self.Kcos.append(0.0)
           #self.Kcos[i] = float(b[i+1])
           self.Kcos[i] = float(b[i])
           self.Ksin.append(0.0)
           #self.Ksin[i] = float(b[i+6])
           self.Ksin[i] = float(b[i+5])
        #self.A = float(b[0])
        self.A = - sum(self.Kcos[:])
        self.numFit = 5
        self.Barrier_Height = max(potentials)

        #print self.Kcos
        #print self.Ksin

        #print "Potential-Read Potential-Fit"
        pot = []
        for i in range(3*nfit):
           angle = i*360/3/nfit
           pot.append([angle,self.getPotential(angle)])
        #   print '%14.2f'%potentials[i]+'%14.3f'%pot[i]
        #print
        #g=Gnuplot.Gnuplot()
        #g('set data style linespoints')
        #g.plot(potgiven,pot)
        #raw_input('Please press enter to continue ...\n')
 
        
        return
