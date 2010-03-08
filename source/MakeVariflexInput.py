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

import sys
import readGeomFc
import pdb
import math
from numpy import *
from scipy import *

class CanTherm:
 CalcType = ''
 ReacType = ''
 Temp = []
 MoleculeList = []
 Entropy = []
 Thermal = []
 Cp = []
 scale = 0.0
 #CBSQB3 E for H, N, O, C, P
 atomE = {'H':-0.499818 , 'N':-54.520543 , 'O':-74.987624 , 'C':-37.785385 , 'P':-340.817186}
 #expt H contains H + TC + SOC (spin orbital correction)
 atomH = {'H':50.62 , 'N':111.49 , 'O':58.163 , 'C':169.8147 }
 #BAC for C-H C-C C=C C.TB.C  O-H  C-O C=O
 bondC = [-0.11, -0.3, -0.08, -0.64, 0.02, 0.33, 0.55]



def main():
  data = CanTherm()
  inputFile = open(sys.argv[1],'r') 
  readGeomFc.readInputFile(inputFile,data)
  
  mol = data.MoleculeList[0]
  symm = mol.extSymm
  Mass = mol.Mass
  geom = mol.geom

  print "NearProlateP"
  print  "NElecStatesD"
  print  "1"
  print  "ElecStatesL"
  print mol.nelec," , 0.0D0"
  print "NModesD"
  print len(mol.Freq), " ,0.99"
  print  "ModesL"
  for i in range(len(mol.Freq)/3+1):
      for j in range(3):
          if 3*i+j <len(mol.Freq):
              print '%10.3f'%mol.Freq[3*i+j],
      print
  print
  if (len(mol.Harmonics) != 0):
    print "HindRotorD"
    print len(mol.Harmonics)
    for i in range(len(mol.Harmonics)):
       harmonic = mol.Harmonics[i]
       rot = mol.rotors[i+1]
       symm = symm * mol.rotors[i+1].symm
       print "CartesianP"
       print "NMAtomsD"
       print len(rot.atomsList)," ", len(rot.nonRotorList)   
       print "PositionL"
       for j in rot.atomsList:
         print float(Mass[j-1]), "%12.3f"%geom[j-1,0], "%12.3f"%geom[j-1,1], "%12.3f"%geom[j-1,2]
       print float(Mass[rot.pivot2-1]), "%12.3f"%geom[rot.pivot2-1,0], "%12.3f"%geom[rot.pivot2-1,1], "%12.3f"%geom[rot.pivot2-1,2]
       for j in rot.nonRotorList:
         if (j != rot.pivot2):
           print float(Mass[j-1]), "%12.3f"%geom[j-1,0], "%12.3f"%geom[j-1,1], "%12.3f"%geom[j-1,2]

       print "HindTypeD"
       print "10"
       print "NCOSNSIND"
       print "5 5"
       print "HindParL"
       for i in range(5):
          print harmonic.Kcos[i]*349.8, i+1
       for i in range(5):
          print harmonic.Ksin[i]*349.8, i+1

       print "0.0"
    print
  print "SigRotD"
  print symm 
  print "CartesianP"
  print "NAtomsD"
  print len(mol.Mass)
  print "PositionL"
  for i in range(len(Mass)):
    print float(Mass[i]), "%12.3f"%geom[i,0], "%12.3f"%geom[i,1], "%12.3f"%geom[i,2]

if __name__ == "__main__":
   main()



