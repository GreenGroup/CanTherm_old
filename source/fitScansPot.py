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


#the file passed to the program is the gaussian output file from which it reads
#the geometry of the molecule
#assumes that inertia.dat is present in the current directory and reads 
#the inertia data

import sys
sys.path.append('/home/sandeeps/program/CanTherm')
import readGeomFc
import os
from numpy import *
import geomUtility
import Gnuplot, Gnuplot.funcutils

#inputFiles = ['1220.log','0120.log','1000.log','1120.log','2101.log','2120.log']
inputFiles = ['022.log','112.log','012.log','120.log']

#inputFiles = ['12211.log','02211.log','22211.log','11011.log','11022.log']
#inputFiles = ['111122.log','212211.log', '102022.log', '002022.log',   '212122.log',   '022022.log']
harmonics = file('harmonics.dat_ring','w')
energy_base =  readGeomFc.readHFEnergy(inputFiles[0][:-4]+'_rot0_6.log')
numRotors = 4
harmonics.write(str(len(inputFiles))+'\n')
for files in inputFiles:
    y = matrix(zeros((13,1),dtype=float))
    x = matrix(zeros((13,11),dtype=float))
    energy = readGeomFc.readHFEnergy(files[:-4]+'_rot0_6.log')
#    harmonics.write(files+'\n')
    (geom,Mass) = readGeomFc.readGeom(open(files,'r'))
    inertia = open('inertia.dat','r')
    (rotors) = readGeomFc.readGeneralInertia(inertia,Mass)
    if (files == inputFiles[0]):
        K = geomUtility.calculateD32(geom,Mass,rotors)
        detD = 1.0
        for i in range(numRotors):
            print K[i],
            detD = detD*K[i]
        harmonics.write(str(float(detD))+'\n')
    harmonics.write(str((energy-energy_base)*627.5095)+'\n')
    if (energy < energy_base):
        print files, " has lower energy"
        exit()
    for i in range(numRotors):
        potgiven = []
        for j in range(13):
            angle = (-60.0 + j*120.0/12.0)*2*math.pi/360.0
            fname = files[:-4]+'_rot'+str(i)+'_'+str(j)+'.log'
            y[j] = (readGeomFc.readHFEnergy(fname)-energy)*627.5095
            x[j,0] = 1.0
            for k in range(5):
                x[j,k+1] = cos((k+1)*angle)
                x[j,k+6] = sin((k+1)*angle)

            potgiven.append([angle,float(y[j])])

        XtX = transpose(x)*x
        XtY = transpose(x)*y
        b = linalg.inv(XtX)*XtY

        print 'rotor', i, files

        pot = []
        for j in range(21):
            angle = (-60.0+j*120/20.0)*2*math.pi/360.0
            v = b[0]
            for k in range(5):
                v = v + b[k+1]*cos((k+1)*angle)+b[k+6]*sin((k+1)*angle)
            pot.append([angle,v])

        g = Gnuplot.Gnuplot()
        harmonics.write(str(float(b[0]))+'\n')
        for k in range(5):
            harmonics.write(str(float(b[k+1]))+'\t'+str( float(b[k+6]))+'\n')
        harmonics.write('\n')
        #g('set data style linespoints')
        #g.plot(potgiven,pot)

        #raw_input("Enter...")
                          

    harmonics.write('\n')
