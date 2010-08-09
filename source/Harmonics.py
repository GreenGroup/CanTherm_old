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
import Rotor
import geomUtility
import Gnuplot,Gnuplot.funcutils

class Harmonics:
    numFit = 0
    Kcos = []
    Ksin = []
    A = 0.0
    Barrier_Height = 0.0

    #inertia variables for variable moment of inertia case
    InumFit=0
    Icos = []
    Isin = []
    B = 0.0

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
#        g=Gnuplot.Gnuplot(persist=1)
        #g('set data style linespoints')
#        plot1 = Gnuplot.PlotItems.Data(potgiven, with_="points 3", title=None)
# 	plot2 = Gnuplot.PlotItems.Data(pot, with_="lines", title="fit" )
# 	g.plot(plot1, plot2)
        #raw_input('Please press enter to continue ...\n')
 
        
        return

    #V0 is analogous to pot0 in fitPotential above
    def fitMM4Potential(self,file,V0,dihedralMinimum, variableInertia, rotor, Ki):
	#read the potentials from the file
	#example block below:
#Methanol                                                    0   6 0  0 0  0 10.0
#0   1      0.0000000         4    0         0    0    0    0    0         0    0
#    1    2
#    1    3    1    4    1    5    2    6
#   0.73340   0.03788  -0.00626 C  1(  1)
#  -0.69331  -0.07708  -0.03183 O  6(  2)
#   1.18040  -0.92717   0.31879 H  5(  3)
#   1.04571   0.83498   0.69975 H  5(  4)
#   1.10634   0.29039  -1.02310 H  5(  5)
#  -1.06298   0.57455   0.57983 H 21(  6)
#TORSION(1)=   0.0, TORSION(2)=   0.0, ENERGY=      1.8643 KCAL/MOLE
#
#   0.73340   0.03788  -0.00626 C  1(  1)
#  -0.69331  -0.07708  -0.03183 O  6(  2)
#   1.18040  -0.92717   0.31879 H  5(  3)
#   1.04571   0.83498   0.69975 H  5(  4)
#   1.10634   0.29039  -1.02310 H  5(  5)
#  -1.06298   0.57455   0.57983 H 21(  6)
#TORSION(1)=   0.1, TORSION(2)=   0.0, ENERGY=      1.8643 KCAL/MOLE
#
#   0.73340   0.03788  -0.00626 C  1(  1)
#  -0.69331  -0.07708  -0.03183 O  6(  2)
#   1.18040  -0.92717   0.31879 H  5(  3)
#   1.04571   0.83498   0.69975 H  5(  4)
#   1.10634   0.29038  -1.02310 H  5(  5)
#  -1.06298   0.57455   0.57983 H 21(  6)
#TORSION(1)=   0.2, TORSION(2)=   0.0, ENERGY=      1.8643 KCAL/MOLE

	read = open(file,'r')
	nfit = 1
	potentials = [0.0]
	potgiven = [[0.0,0.0]]#initialize with value at mimimum (reset to angle of zero)
	geomList = []
	MassList = []
	for line in read:
	    if line[39:40]==")":#look for geometry lines; cf. readGeomFc.readMM4Geom()
		MassList.append(readGeomFc.getMassByAtomicSymbol(line[30:33].strip()))
		xc = float(line[0:10])
		yc = float(line[10:20])
		zc = float(line[20:30])
		geomList.append([xc,yc,zc])
	    if(line.startswith('TORSION(1)=')):#torsional energy line terminates the geometry section
		#1: read in the torsion information
		tokens = line.split()
		nfit=nfit+1
		potVal = float(tokens[-2])-V0
		angleDeg = float(line[11:17])-dihedralMinimum
		potentials.append(potVal)
		potgiven.append([angleDeg, potVal])
		#2: process the geometry information (i.e. calculate inertia) and store the result (in variableInertia cases)
		if variableInertia:
		    #convert to the matrix format used by CanTherm
		    geom = matrix(geomList)
		    Mass = matrix(MassList).transpose()
		    inertVal = geomUtility.calculateD32forIndividualRotor(geom,Mass,rotor)
		    inertgiven.append([angleDeg, 1/inertVal])
		#3: reset geomList and MassList for reading the next geometry
		geomList = []
		MassList = []

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
           angle = potgiven[i][0]*math.pi/180.0 #note conversion to radians
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
        #pot = []
        #for i in range(3*nfit):
        #   angle = i*360/3/nfit
        #   pot.append([angle,self.getPotential(angle)])
        #   print '%14.2f'%potentials[i]+'%14.3f'%pot[i]
        #print
#        g=Gnuplot.Gnuplot(persist=1)
        #g('set data style linespoints')
#        plot1 = Gnuplot.PlotItems.Data(potgiven, with_="points 3", title=None)
# 	plot2 = Gnuplot.PlotItems.Data(pot, with_="lines", title="fit" )
# 	g.plot(plot1, plot2)
        #raw_input('Please press enter to continue ...\n')

	#fit 1/moment of inertia; unlike above, we don't need (or want) to enforce d/d(phi) = 0 at phi = 0
	if variableInertia:
	    Y = matrix(zeros((nfit,1),dtype=float))
	    X = matrix(zeros((nfit,11),dtype=float))
	    for i in range(nfit):
	       Y[i,0]=inertgiven[i][1]
	       angle = inertgiven[i][0]*math.pi/180.0 #note conversion to radians
	       for j in range(5):
		  X[i,j] = math.cos((j+1)*angle)
		  X[i,j+5] = math.sin((j+1)*angle)
	       X[i,10] = 1.0 #contribution from the constant term (in this sense, it differs from approach with potential, above)


	    XtX = transpose(X)*X
	    XtY = transpose(X)*Y
	    b = linalg.inv(XtX)*XtY

	    for i in range(5):
	       self.Icos.append(0.0)
	       self.Icos[i] = float(b[i])
	       self.Isin.append(0.0)
	       self.Isin[i] = float(b[i+5])
	    #self.B = 1/Ki - sum(self.Icos[:])#ensure it passes through 1/Ki at phi=0; note that this approach (as in the similar assignment of self.A in the method above) does not necessarily ensure that we have the least squares error given this constraint; if the fit is good, the sum of cosines should already be close to 1/Ki and this value should be close to zero
	    self.B = float(b[10]) #unlike above, we will not try to make sure the fit passes through the value at phi=0; instead, we use the constant term as one of the fitted parameters
	    self.InumFit = 5


        return
