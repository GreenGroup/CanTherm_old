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
sys.path.append('/home/sandeeps/site-packages')
sys.path.append('/home/sandeeps/site-packages/Numeric')
import readGeomFc
import Eckart
import pdb
import math
from numpy import *
from scipy import *
import FitCp 

class CanTherm:
 CalcType = ''
 ReacType = ''
 TunnelType = ''
 numProds = []
 Temp = []
 MoleculeList = []
 Entropy = []
 Thermal = []
 Cp = []
 Parition = []
 scale = 0.0
 #CBSQB3 E for H, N, O, C, P
 atomEcbsqb3 = {'H':-0.499818 , 'N':-54.520543 , 'O':-74.987624 , 'C':-37.785385 , 'P':-340.817186}
 #G3 E for H, N, O, C, P
 atomEg3 = {'H':-0.5010030, 'N':-54.564343, 'O':-75.030991, 'C':-37.827717, 'P':-341.116432}
 #Klip QCI(dz,tz)+ MP2(tz,qz) E for H, N, O, C, P
 #atomEKlip_1 = {'H':-0.49991705, 'O':-74.99507456, 'C':-37.78778408,}
 atomEKlip_1 = {'H':-0.50003976, 'O':-75.00915718, 'C':-37.79249556,} 
#Klip QCI(tz,qz) E for H, N, O, C, P
#atomEKlip_2 = {'H':-0.50003976, 'O':-75.00692740, 'C':-37.79044862,}
 atomEKlip_2 = {'H':-0.50003976, 'O':-75.00692746, 'C':-37.79044863,}
#Klip CCSD(T)(tz,qz) E for H, N, O, C, P
 atomEKlip_2_cc = {'H':-0.50003976, 'O':-75.00681155, 'C':-37.79029443,}

 #expt H contains H + TC + SOC (spin orbital correction)
 atomH = {'H':50.62 , 'N':111.49 , 'O':58.163 , 'C':169.8147 }
 #BAC for C-H    C-C   C=C    C.TB.C  O-H   C-O   C=O   N.TB.N O=O   H-H  C.TB.N
 bondC = [-0.11, -0.3, -0.08, -0.64,  0.02, 0.33, 0.55, -2.0,  -0.2, 1.1, -0.89]



def main():
  data = CanTherm()
  inputFile = open(sys.argv[1],'r') 
  oFile = open('cantherm.out','w')
  readGeomFc.readInputFile(inputFile,data)

  data.Entropy=len(data.MoleculeList)*len(data.Temp)*[0.0]
  data.Cp=len(data.MoleculeList)*len(data.Temp)*[0.0]
  data.Thermal=len(data.MoleculeList)*len(data.Temp)*[0.0]
  data.Partition=len(data.MoleculeList)*len(data.Temp)*[1.0]
  Entropy = data.Entropy
  Cp = data.Cp
  Thermal = data.Thermal
  Partition = data.Partition

  for i in range(len(data.MoleculeList)):
     molecule = data.MoleculeList[i]
     oFile.write('Molecule '+str(i+1)+':\n')
     oFile.write('-----------\n\n')
     molecule.printData(oFile)

     oFile.write('\nThermodynamic Data\n')

     Temp = data.Temp
     #translation
     (ent,cp,dh,q) = molecule.getTranslationThermo(oFile,data.Temp) 
     for j in range(len(Temp)):
         Entropy[i*len(Temp)+j]=Entropy[i*len(Temp)+j]+ent[j]
         Cp[i*len(Temp)+j]=Cp[i*len(Temp)+j]+cp[j]
         Thermal[i*len(Temp)+j]=Thermal[i*len(Temp)+j]+dh[j]
         Partition[i*len(Temp)+j]=Partition[i*len(Temp)+j]*q[j]
  
     #vibrational
     (ent,cp,dh,q) = molecule.getVibrationalThermo(oFile,data.Temp,data.scale) 
     for j in range(len(Temp)):
         Entropy[i*len(Temp)+j]=Entropy[i*len(Temp)+j]+ent[j]
         Cp[i*len(Temp)+j]=Cp[i*len(Temp)+j]+cp[j]
         Thermal[i*len(Temp)+j]=Thermal[i*len(Temp)+j]+dh[j]
         Partition[i*len(Temp)+j] = Partition[i*len(Temp)+j]*q[j]
         #print '%12.2f'%float(ent[j]),
     #print '\n'

     #Internal rotational
     if molecule.numRotors != 0:
      (ent,cp,dh,q) = molecule.getIntRotationalThermo_Q(oFile,data.Temp) 
      for j in range(len(Temp)):
         Entropy[i*len(Temp)+j]=Entropy[i*len(Temp)+j]+ent[j]
         Cp[i*len(Temp)+j]=Cp[i*len(Temp)+j]+cp[j]
         Thermal[i*len(Temp)+j]=Thermal[i*len(Temp)+j]+dh[j]
         Partition[i*len(Temp)+j] = Partition[i*len(Temp)+j]*q[j]
         #print '%12.2f'%float(ent[j]),
     #print '\n'

     #External rotational
     (ent,cp,dh,q) = molecule.getExtRotationalThermo(oFile,data.Temp) 
     for j in range(len(Temp)):
         Entropy[i*len(Temp)+j]=Entropy[i*len(Temp)+j]+ent[j]
         Cp[i*len(Temp)+j]=Cp[i*len(Temp)+j]+cp[j]
         Thermal[i*len(Temp)+j]=Thermal[i*len(Temp)+j]+dh[j]
         Partition[i*len(Temp)+j] = Partition[i*len(Temp)+j]*q[j]

     for j in range(len(Temp)):
         Entropy[i*len(Temp)+j]=Entropy[i*len(Temp)+j]+1.985*math.log(molecule.nelec)
         Partition[i*len(Temp)+j] = Partition[i*len(Temp)+j] * molecule.nelec

     #print Enthalpy

     H = molecule.Energy
     if not molecule.Etype == 'mm4':#for the MM4 case, the value passed in should be in kcal/mol and should not require unit adjustments or atomization energy information
	 atoms = readGeomFc.getAtoms(molecule.Mass)
	 atomsH = 0.0
	 if molecule.Etype == 'cbsqb3':
	    atomE = data.atomEcbsqb3
	 if molecule.Etype == 'g3':
	    atomE = data.atomEg3
	 if molecule.Etype == 'klip_1':
	    atomE = data.atomEKlip_1
	 if molecule.Etype == 'klip_2':
	    atomE = data.atomEKlip_2
	 if molecule.Etype == 'klip_2_cc':
	    atomE = data.atomEKlip_2_cc
	 for atom in atoms:
	     H -= atomE[atom]
	     atomsH += data.atomH[atom]
	 H = H*627.5095+atomsH

     if molecule.Etype == 'cbsqb3':
       b = 0
       for bonds in molecule.bonds:
         H += bonds*data.bondC[b]
         b += 1


     #MRH 30Jan2010
     #This E0 will be used in the Eckart tunneling calculation
     molecule.E0 = H

     H += Thermal[i*len(Temp)+0]

     print '%12.2f'%H + '%12.2f'%Entropy[i*len(Temp)+0]
#     print '%12.2f'%float(H*4.187) + '%12.2f'%float(Entropy[i*len(Temp)+0]*4.187)
     for c in range(1,len(Temp)):
        print '%12.2f'%Cp[i*len(Temp)+c],
     print '\n'

     #for c in range(len(Temp)):
        #print '%12.2e'%Partition[i*len(Temp)+c],
     #print
     if data.fitcp:
       FitCp.FitHeatCapacity(Temp,Cp, molecule.linearity, molecule.Freq, molecule.numRotors, molecule.R, oFile)

  if len(data.MoleculeList) == 1:
     return

  #fit the rate coefficient
  A = matrix(zeros((len(Temp),3),dtype=float))
  A1 = matrix(zeros((len(Temp),2),dtype=float))
  y = matrix(zeros((len(Temp),1),dtype=float))

  rate = [0.0]*len(Temp)
  for j in range(len(Temp)):

    if (data.ReacType == 'Unimol'):
      #rate[j] = (1.381e-23*Temp[j]/6.626e-34)*math.exp((Entropy[len(Temp)+j]-Entropy[j])/1.985)*math.exp(-(data.MoleculeList[1].Energy - data.MoleculeList[0].Energy)*627.5095*1.0e3/1.985/Temp[j])
      kbT_h = (1.381e-23*Temp[j]/6.626e-34)
      G_TS = Thermal[len(Temp)+j]*1e3+data.MoleculeList[1].Energy*627.5095*1e3-Temp[j]*Entropy[len(Temp)+j]
      G_react = Thermal[j]*1e3+data.MoleculeList[0].Energy*627.5095*1e3-Temp[j]*Entropy[j]
      #qTS_qA = Partition[len(Temp)+j]/Partition[j]
      #exp_e_RT = math.exp(-(data.MoleculeList[1].Energy-data.MoleculeList[0].Energy)*627.5095*1e3/1.985/Temp[j])
      #print kbT_h * exp_S_R * exp_H_RT
      #print qTS_qA * exp_e_RT
      rate[j] = kbT_h * math.exp(-(G_TS-G_react)/1.985/Temp[j])

      #Tunneling:
      #Wigner - need imaginary frequency + temperature
      if data.TunnelType == 'Wigner' :
          kappa = 1.0 + 1.0/24.0 * (1.44*data.MoleculeList[1].imagFreq/Temp[j])**2
          print kappa
          rate[j] *= 1.0 + 1.0/24.0 * (1.44*data.MoleculeList[1].imagFreq/Temp[j])**2
      #Symmetric Eckart - need imaginary frequency + temperature + (E0_TS - E_react)
      elif data.TunnelType == 'sEckart' :
          delV1 = data.MoleculeList[1].E0-data.MoleculeList[0].E0
          delV1 = delV1 / 6.022e23 * 4184
          alpha1 = 2*math.pi*delV1/6.626e-34/3.00e10/(-1*data.MoleculeList[1].imagFreq)
          alpha2 = alpha1
          rate[j] *= Eckart.computeTunnelingCorrection(delV1,Temp[j],alpha1,alpha2)
      #Symmetric Eckart - need imaginary frequency + temperature + (E0_TS - E_react) + (E0_TS - E_prod)
      elif data.TunnelType == 'aEckart' :
          delV1 = data.MoleculeList[1].E0-data.MoleculeList[0].E0
          delV1 = delV1 / 6.022e23 * 4184
          alpha1 = 2*math.pi*delV1/6.626e-34/3.00e10/(-1*data.MoleculeList[1].imagFreq)
          if data.numProds == 1 :
              delV2 = data.MoleculeList[1].E0-data.MoleculeList[2].E0
          else :
              delV2 = data.MoleculeList[1].E0-data.MoleculeList[2].E0-data.MoleculeList[3].E0
          alpha2 = 2*math.pi*delV2/6.022e23/6.626e-34/3.00e10*4184/(-1*data.MoleculeList[1].imagFreq)
          rate[j] *= Eckart.computeTunnelingCorrection(delV1,Temp[j],alpha1,alpha2)

    elif (data.ReacType == 'Bimol'):
      #rate[j] = (1.381e-23*Temp[j]/6.626e-34)*(82.05746*Temp[j]/1.0)*math.exp((Entropy[2*len(Temp)+j]-Entropy[len(Temp)+j]-Entropy[j])/1.985)*math.exp(-(data.MoleculeList[2].Energy - data.MoleculeList[0].Energy - data.MoleculeList[1].Energy)*627.5095*1.0e3/1.985/Temp[j])
      
      kbT_h = (1.381e-23*Temp[j]/6.626e-34)
      exp_S_R = math.exp((Entropy[2*len(Temp)+j]-Entropy[len(Temp)+j]-Entropy[j])/1.985)
      exp_H_RT = math.exp(-(Thermal[2*len(Temp)+j]-Thermal[len(Temp)+j]-Thermal[j])*1e3/1.985/Temp[j])
      rate[j] = kbT_h * exp_S_R * exp_H_RT

      #wigner correction
      #rate[j] *= 1.0 + 1.0/24.0 * (1.44*data.MoleculeList[2].imagFreq/Temp[j])**2
      #Tunneling:
      #Wigner - need imaginary frequency + temperature
      if data.TunnelType == 'Wigner' :
          rate[j] *= 1.0 + 1.0/24.0 * (1.44*data.MoleculeList[2].imagFreq/Temp[j])**2
      #Symmetric Eckart - need imaginary frequency + temperature + (E0_TS - E_react)
      elif data.TunnelType == 'sEckart' :
          delV1 = data.MoleculeList[2].E0-data.MoleculeList[0].E0-data.MoleculeList[1].E0
          delV1 = delV1 / 6.022e23 * 4184
          alpha1 = 2*math.pi*delV1/6.626e-34/3.00e10/(-1*data.MoleculeList[2].imagFreq)
          alpha2 = alpha1
          rate[j] *= Eckart.computeTunnelingCorrection(delV1,Temp[j],alpha1,alpha2)
      #Symmetric Eckart - need imaginary frequency + temperature + (E0_TS - E_react) + (E0_TS - E_prod)
      elif data.TunnelType == 'aEckart' :
          delV1 = data.MoleculeList[2].E0-data.MoleculeList[0].E0-data.MoleculeList[1].E0
          delV1 = delV1 / 6.022e23 * 4184
          alpha1 = 2*math.pi*delV1/6.626e-34/3.00e10/(-1*data.MoleculeList[2].imagFreq)
          if data.numProds == 1 :
              delV2 = data.MoleculeList[2].E0-data.MoleculeList[3].E0
          else :
              delV2 = data.MoleculeList[2].E0-data.MoleculeList[3].E0-data.MoleculeList[4].E0
          alpha2 = 2*math.pi*delV2/6.022e23/6.626e-34/3.00e10*4184/(-1*data.MoleculeList[2].imagFreq)
          rate[j] *= Eckart.computeTunnelingCorrection(delV1,Temp[j],alpha1,alpha2)

    A[j,:] = mat([1.0, math.log(Temp[j]), -1.0/1.985/Temp[j]])
    y[j] = log(rate[j])
    A1[j,:] = mat([1.0, -1.0/1.985/Temp[j]])
  b = linalg.inv(transpose(A)*A)*(transpose(A)*y)
  b1 = linalg.inv(transpose(A1)*A1)*(transpose(A1)*y)
  oFile.write('\n\nRate Data: Modified Arrhenius format\n')
  oFile.write('r = A*(T/1000)^n*exp(-Ea/R/T)'+'%12.2e'%(exp(b[0])*1000.0**float(b[1]))+'%6.2f'%b[1]+'%12.2f'%(b[2]/1.0e3)+'\n')
  oFile.write('r = A*T^n*exp(-Ea/R/T)'+'%12.2e'%(exp(b[0]))+'%6.2f'%b[1]+'%12.2f'%(b[2]/1.0e3)+'\n')
  oFile.write('%12s'%'Temperature'+'%12s'%'Rate'+'%12s\n'%'Fit Rate')
  for j in range(len(Temp)):
      fitrate = exp(b[0])*Temp[j]**float(b[1])*exp(-b[2]/1.985/Temp[j])
      oFile.write('%12.2f'%Temp[j]+'%12.2e'%rate[j]+'%12.2e\n'%fitrate)
  oFile.write('\n\n')
  oFile.write('Rate Data: Arrhenius format\n')
  oFile.write('r = A*exp(-Ea/R/T)'+'%12.2e'%(exp(b1[0]))+'%12.2f'%(b1[1]/1.0e3)+'\n')
  oFile.write('%12s'%'Temperature'+'%12s'%'Rate'+'%12s\n'%'Fit Rate')
  for j in range(len(Temp)) :
      fitrate2 = exp(b1[0])*exp(-b1[1]/1.985/Temp[j])
      oFile.write('%12.2f'%Temp[j]+'%12.2e'%rate[j]+'%12.2e\n'%fitrate2)
  oFile.write('\n\n')
  oFile.close()
  

if __name__ == "__main__":
   main()

