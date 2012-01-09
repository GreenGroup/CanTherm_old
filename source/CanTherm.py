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

 #units below are Hartree
 SOC = {'H':0.0, 'N':0.0, 'O': -0.000355, 'C': -0.000135, 'P': 0.0} #spin orbit correction (SOC) in Hartrees, values taken from note 22 of http://jcp.aip.org/resource/1/jcpsa6/v109/i24/p10570_s1 and converted to hartree (values in millihartree are also available (with fewer significant figures) in http://jcp.aip.org/resource/1/jcpsa6/v106/i3/p1063_s1)
 # CBS-QB3 and G3 methods include SOC for atoms, so no correction is necessary; gmagoon has confirmed with Gaussian, Inc. that the atom energies calculated by Gaussian using these methods include an SOC correction (though it is not explicitly stated in the output)
 #CBSQB3 E for H, N, O, C, P
 atomEcbsqb3 = {'H':-0.499818 , 'N':-54.520543 , 'O':-74.987624 , 'C':-37.785385 , 'P':-340.817186}
  #CBSQB3ultrafine E for H, N, O, C, P (P, N values taken from "regular" CBS-QB3 value as a first approximation; H unchanged, I believe)
 atomEcbsqb3uf = {'H':-0.499818 , 'N':-54.520543 , 'O':-74.987619 , 'C':-37.785376 , 'P':-340.817186}
 #G3 E for H, N, O, C, P
 atomEg3 = {'H':-0.5010030, 'N':-54.564343, 'O':-75.030991, 'C':-37.827717, 'P':-341.116432}
 #the values below do not include spin-orbit correction, so this is added below
 #Klip QCI(dz,tz)+ MP2(tz,qz) E for H, N, O, C, P
 #atomEKlip_1 = {'H':-0.49991705, 'O':-74.99507456, 'C':-37.78778408,}
 atomEKlip_1 = {'H':-0.50003976+SOC['H'], 'O':-75.00915718+SOC['O'], 'C':-37.79249556+SOC['C'],}
#Klip QCI(tz,qz) E for H, N, O, C, P
#atomEKlip_2 = {'H':-0.50003976, 'O':-75.00692740, 'C':-37.79044862,}
 atomEKlip_2 = {'H':-0.50003976+SOC['H'], 'O':-75.00692746+SOC['O'], 'C':-37.79044863+SOC['C'],}
#Klip CCSD(T)(tz,qz) E for H, N, O, C, P
 atomEKlip_2_cc = {'H':-0.50003976+SOC['H'], 'O':-75.00681155+SOC['O'], 'C':-37.79029443+SOC['H'],}

#units below are kcal/mol;
 atomH0 = {'H': 51.63 , 'N': 112.53 ,'O': 58.99 ,'C': 169.98 } #expt Hf at 0K (see Gaussian thermo whitepaper: http://www.gaussian.com/g_whitepap/thermo.htm); note: these values are relatively old and some improvement may be possible by using newer values, particularly for carbon; however, care should be taken to ensure that they are compatible with the BAC values (if BACs are used)
 atomTC = {'H': 1.01 , 'N': 1.04, 'O': 1.04 ,'C': 0.25 }#thermal contribution Hss(298K)-Hss(0K) reported by Gaussian thermo whitepaper; this will be subtracted from the corresponding value in atomH0 to produce an energy used in calculating Hf298
 atomH={'H': atomH0['H'] - atomTC['H'], 'N': atomH0['N']-atomTC['N'], 'O': atomH0['O']-atomTC['O'], 'C': atomH0['C']-atomTC['C'] }
 #old approach: included SOC in experimental values; this produced an incorrect double-counting of SOC when using methods like CBS-QB3 and G3 that include the SOC correction in their reported atom energies (atomEcbsqb3 and atomEg3)
 #expt H contains H + TC + SOC (spin orbital correction)
 #atomH = {'H':50.62 , 'N':111.49 , 'O':58.163 , 'C':169.8147 }
 
 #BAC for C-H    C-C   C=C    C.TB.C  O-H   C-O   C=O   N.TB.N O=O   H-H  C.TB.N
 bondC = [-0.11, -0.3, -0.08, -0.64,  0.02, 0.33, 0.55, -2.0,  -0.2, 1.1, -0.89]



def main():
  data = CanTherm()
  inputFile = open(sys.argv[1],'r')
  #determine the output file name by removing the contents of the input filename following the first period and appending the .canout suffix
  periodPosition = sys.argv[1].rfind(".")
  outfilename = sys.argv[1][0:periodPosition]+".canout"
  oFile = open(outfilename,'w')
  readGeomFc.readInputFile(inputFile,data)

  data.Entropy=len(data.MoleculeList)*len(data.Temp)*[0.0]
  data.Cp=len(data.MoleculeList)*len(data.Temp)*[0.0]
  data.Thermal=len(data.MoleculeList)*len(data.Temp)*[0.0]
  data.Partition=len(data.MoleculeList)*len(data.Temp)*[1.0]
  Entropy = data.Entropy
  Cp = data.Cp
  Thermal = data.Thermal
  Partition = data.Partition
  Entropy298=0.0
  Thermal298=0.0

  for i in range(len(data.MoleculeList)):
     molecule = data.MoleculeList[i]
     oFile.write('Molecule '+str(i+1)+':\n')
     oFile.write('-----------\n\n')
     molecule.printData(oFile)

     oFile.write('\nThermodynamic Data\n')

     Temp = data.Temp
     temp298 = [298.15] #use 298.15 K as well; this will be used below for printing Hf298, S298
     #translation
     (ent,cp,dh,q) = molecule.getTranslationThermo(oFile,data.Temp) 
     for j in range(len(Temp)):
         Entropy[i*len(Temp)+j]=Entropy[i*len(Temp)+j]+ent[j]
         Cp[i*len(Temp)+j]=Cp[i*len(Temp)+j]+cp[j]
         Thermal[i*len(Temp)+j]=Thermal[i*len(Temp)+j]+dh[j]
         Partition[i*len(Temp)+j]=Partition[i*len(Temp)+j]*q[j]
     (ent298,cp298,dh298,q298) = molecule.getTranslationThermo(oFile,temp298)
     Entropy298=Entropy298+ent298[0]
     Thermal298=Thermal298+dh298[0]
  
     #vibrational
     (ent,cp,dh,q) = molecule.getVibrationalThermo(oFile,data.Temp,data.scale) 
     for j in range(len(Temp)):
         Entropy[i*len(Temp)+j]=Entropy[i*len(Temp)+j]+ent[j]
         Cp[i*len(Temp)+j]=Cp[i*len(Temp)+j]+cp[j]
         Thermal[i*len(Temp)+j]=Thermal[i*len(Temp)+j]+dh[j]
         Partition[i*len(Temp)+j] = Partition[i*len(Temp)+j]*q[j]
         #print '%12.2f'%float(ent[j]),
     #print '\n'
     (ent298,cp298,dh298,q298) = molecule.getVibrationalThermo(oFile,temp298,data.scale)
     Entropy298=Entropy298+ent298[0]
     Thermal298=Thermal298+dh298[0]

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
      (ent298,cp298,dh298,q298) = molecule.getIntRotationalThermo_Q(oFile,temp298)
      Entropy298=Entropy298+ent298[0]
      Thermal298=Thermal298+dh298[0]

     #External rotational
     (ent,cp,dh,q) = molecule.getExtRotationalThermo(oFile,data.Temp) 
     for j in range(len(Temp)):
         Entropy[i*len(Temp)+j]=Entropy[i*len(Temp)+j]+ent[j]
         Cp[i*len(Temp)+j]=Cp[i*len(Temp)+j]+cp[j]
         Thermal[i*len(Temp)+j]=Thermal[i*len(Temp)+j]+dh[j]
         Partition[i*len(Temp)+j] = Partition[i*len(Temp)+j]*q[j]
     (ent298,cp298,dh298,q298) = molecule.getExtRotationalThermo(oFile,temp298)
     Entropy298=Entropy298+ent298[0]
     Thermal298=Thermal298+dh298[0]

     for j in range(len(Temp)):
         Entropy[i*len(Temp)+j]=Entropy[i*len(Temp)+j]+1.9872*math.log(molecule.nelec)
         Partition[i*len(Temp)+j] = Partition[i*len(Temp)+j] * molecule.nelec
     Entropy298=Entropy298+1.9872*math.log(molecule.nelec)

     #print Enthalpy

     H = molecule.Energy
     if not molecule.Etype == 'mm4':#for the MM4 case, the value passed in should be in kcal/mol and should not require unit adjustments or atomization energy information
	 atoms = readGeomFc.getAtoms(molecule.Mass)
	 atomsH = 0.0
	 #atomsH0 = 0.0
	 if molecule.Etype == 'cbsqb3':
	    atomE = data.atomEcbsqb3
	 if molecule.Etype == 'cbsqb3uf':
	    atomE = data.atomEcbsqb3uf
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
	#     atomsH0 += data.atomH0[atom]
	 H = H*627.5095+atomsH
	# H0 = H*627.5095+atomsH0

     if (molecule.Etype == 'cbsqb3' or molecule.Etype == 'cbsqb3uf') :
       b = 0
       for bonds in molecule.bonds:
         H += bonds*data.bondC[b]
         b += 1


     #MRH 30Jan2010
     #This E0 will be used in the Eckart tunneling calculation
     molecule.E0 = H

     H298 = H + Thermal298
     H += Thermal[i*len(Temp)+0]

     print '%12.2f'%H298 + '%12.2f'%Entropy298,
#     print '%12.2f'%float(H*4.187) + '%12.2f'%float(Entropy[i*len(Temp)+0]*4.187)
     for c in range(0,len(Temp)):
        print '%12.2f'%Cp[i*len(Temp)+c],
     print '\n'

     oFile.write("Hf298 S298 Cps:\n")
     oFile.write(str(H298)+" "+ str(Entropy298))
#     print '%12.2f'%float(H*4.187) + '%12.2f'%float(Entropy[i*len(Temp)+0]*4.187)
     for c in range(0,len(Temp)):
        oFile.write(" "+str(Cp[i*len(Temp)+c]))
#     print '\n'

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
      #rate[j] = (1.381e-23*Temp[j]/6.626e-34)*math.exp((Entropy[len(Temp)+j]-Entropy[j])/1.9872)*math.exp(-(data.MoleculeList[1].Energy - data.MoleculeList[0].Energy)*627.5095*1.0e3/1.9872/Temp[j])
      kbT_h = (1.381e-23*Temp[j]/6.626e-34)
      G_TS = Thermal[len(Temp)+j]*1e3+data.MoleculeList[1].Energy*627.5095*1e3-Temp[j]*Entropy[len(Temp)+j]
      G_react = Thermal[j]*1e3+data.MoleculeList[0].Energy*627.5095*1e3-Temp[j]*Entropy[j]
      #qTS_qA = Partition[len(Temp)+j]/Partition[j]
      #exp_e_RT = math.exp(-(data.MoleculeList[1].Energy-data.MoleculeList[0].Energy)*627.5095*1e3/1.9872/Temp[j])
      #print kbT_h * exp_S_R * exp_H_RT
      #print qTS_qA * exp_e_RT
      rate[j] = kbT_h * math.exp(-(G_TS-G_react)/1.9872/Temp[j])

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
          delV2 = delV2 / 6.022e23 * 4184
          alpha2 = 2*math.pi*delV2/6.626e-34/3.00e10/(-1*data.MoleculeList[1].imagFreq)
          if (alpha1<alpha2):
              rate[j] *= Eckart.computeTunnelingCorrection(delV1,Temp[j],alpha1,alpha2)
          else:
              # reverse the direction: pass delV2 and alpha1 in place of alpha2
              rate[j] *= Eckart.computeTunnelingCorrection(delV2,Temp[j],alpha2,alpha1)

    elif (data.ReacType == 'Bimol'):
      #rate[j] = (1.381e-23*Temp[j]/6.626e-34)*(82.05746*Temp[j]/1.0)*math.exp((Entropy[2*len(Temp)+j]-Entropy[len(Temp)+j]-Entropy[j])/1.9872)*math.exp(-(data.MoleculeList[2].Energy - data.MoleculeList[0].Energy - data.MoleculeList[1].Energy)*627.5095*1.0e3/1.9872/Temp[j])
      
      kbT_hC = (1.381e-23*Temp[j]/6.626e-34)*(82.05746*Temp[j]/1.0)
      G_TS = Thermal[2*len(Temp)+j]*1e3+data.MoleculeList[2].Energy*627.5095*1e3-Temp[j]*Entropy[2*len(Temp)+j]
      G_react1 = Thermal[len(Temp)+j]*1e3+data.MoleculeList[1].Energy*627.5095*1e3-Temp[j]*Entropy[len(Temp)+j]
      G_react2 = Thermal[j]*1e3+data.MoleculeList[0].Energy*627.5095*1e3-Temp[j]*Entropy[j]
      rate[j] = kbT_hC * math.exp(-(G_TS-G_react1-G_react2)/1.9872/Temp[j])

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
          delV2 = delV2 / 6.022e23 * 4184
          alpha2 = 2*math.pi*delV2/6.626e-34/3.00e10/(-1*data.MoleculeList[2].imagFreq)
          if (alpha1<alpha2):
              rate[j] *= Eckart.computeTunnelingCorrection(delV1,Temp[j],alpha1,alpha2)
          else:
              # reverse the direction: pass delV2 and alpha1 in place of alpha2
              rate[j] *= Eckart.computeTunnelingCorrection(delV2,Temp[j],alpha2,alpha1)

    A[j,:] = mat([1.0, math.log(Temp[j]), -1.0/1.9872/Temp[j]])
    y[j] = log(rate[j])
    A1[j,:] = mat([1.0, -1.0/1.9872/Temp[j]])
  b = linalg.inv(transpose(A)*A)*(transpose(A)*y)
  b1 = linalg.inv(transpose(A1)*A1)*(transpose(A1)*y)
  oFile.write('\n\nRate Data: Modified Arrhenius format\n')
  oFile.write('r = A*(T/1000)^n*exp(-Ea/R/T)'+'%12.2e'%(exp(b[0])*1000.0**float(b[1]))+'%6.2f'%b[1]+'%12.2f'%(b[2]/1.0e3)+'\n')
  oFile.write('r = A*T^n*exp(-Ea/R/T)'+'%12.2e'%(exp(b[0]))+'%6.2f'%b[1]+'%12.2f'%(b[2]/1.0e3)+'\n')
  oFile.write('%12s'%'Temperature'+'%12s'%'Rate'+'%12s\n'%'Fit Rate')
  for j in range(len(Temp)):
      fitrate = exp(b[0])*Temp[j]**float(b[1])*exp(-b[2]/1.9872/Temp[j])
      oFile.write('%12.2f'%Temp[j]+'%12.2e'%rate[j]+'%12.2e\n'%fitrate)
  oFile.write('\n\n')
  oFile.write('Rate Data: Arrhenius format\n')
  oFile.write('r = A*exp(-Ea/R/T)'+'%12.2e'%(exp(b1[0]))+'%12.2f'%(b1[1]/1.0e3)+'\n')
  oFile.write('%12s'%'Temperature'+'%12s'%'Rate'+'%12s\n'%'Fit Rate')
  for j in range(len(Temp)) :
      fitrate2 = exp(b1[0])*exp(-b1[1]/1.9872/Temp[j])
      oFile.write('%12.2f'%Temp[j]+'%12.2e'%rate[j]+'%12.2e\n'%fitrate2)
  oFile.write('\n\n')
  oFile.close()
  

if __name__ == "__main__":
   main()

