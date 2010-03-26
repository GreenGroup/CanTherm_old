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
from Harmonics import *
from scipy import *
import scipy.linalg 
#import random
import geomUtility
import Gnuplot
import os

class Molecule:
    R = 1.985
    kb = 1.38065e-23
    h = 6.626e-34
    amu = 1.6605e-27
#has other attributes
#geom, Mass, Fc, linearity, Energy, Etype, extSymm, nelec, rotors, potentialFile, Iext, ebase, bonds
    def __init__(self,file, isTS):
        self.Freq = []
        self.Harmonics = []
        self.hindFreq = []
        self.bonds = []
        self.Etype = ''
        self.TS = isTS
        self.E0 = []

#read dummy line "MOL #"
        line = readGeomFc.readMeaningfulLine(file)
        #self.linearity = line.split()[0].upper

#read linearity
        line = readGeomFc.readMeaningfulLine(file)
        if line.split()[0].upper() == 'NONLINEAR' :
          self.linearity = 'Nonlinear'
        elif line.split()[0].upper() == 'LINEAR' :
          self.linearity = 'Linear'
        elif line.split()[0].upper() == 'ATOM' :
          self.linearity = 'Atom'
        else :
          print 'Linearity keyword not recognized ' + line
          exit()
        #linearity = line.split()[0].upper

# read geometry
        line = readGeomFc.readMeaningfulLine(file)
        tokens = line.split()
        if tokens[0].upper() != 'GEOM' :
           print 'Geom keyword not found in the input file' + line
           exit()

        if len(tokens) == 1:
           #geometry is given following the GEOM keyword
           line = readGeomFc.readMeaningfulLine(file)
           numAtoms = int(line.split()[0])
           self.geom = matrix(array(zeros((numAtoms,3),dtype=float)))
           self.Mass = matrix(array(zeros((numAtoms,1),dtype=float)))
           for i in range(numAtoms):
               line = readGeomFc.readMeaningfulLine(file)
               tokens = line.split()
               self.geom[i,0]=double(tokens[3])
               self.geom[i,1]=double(tokens[4])
               self.geom[i,2]=double(tokens[5])
               atomicNum = int(tokens[1])
               if (atomicNum == 6): self.Mass[i]=12.0
               if (atomicNum == 8): self.Mass[i]=15.99491
               if (atomicNum == 1): self.Mass[i]=1.00783
               if (atomicNum == 7): self.Mass[i]=14.0031
               if (atomicNum == 17): self.Mass[i]=34.96885
               if (atomicNum == 16): self.Mass[i]=31.97207
               if (atomicNum == 9): self.Mass[i]=18.99840

        #read geometry from the file
        elif tokens[1].upper() == 'FILE' :
           print "reading Geometry from the file: ",tokens[2]
           geomFile = open(tokens[2],'r')
           (self.geom,self.Mass) = readGeomFc.readGeom(geomFile); 
	   #print self.geom
        else:
           print 'Either give geometry or give keyword File followed by the file containing geometry data'
           exit()

        self.calculateMomInertia()

# read force constant or frequency data
        line = readGeomFc.readMeaningfulLine(file)
	tokens = line.split()
        if tokens[0].upper() == 'FORCEC' and tokens[1].upper() == 'FILE':
            #MRH added following line on 2/Dec/2009
            if self.linearity != 'Atom':
                fcfile = open(tokens[2],'r')
                self.Fc = readGeomFc.readFc(fcfile)

                for i in range(0,3*self.Mass.size):
                    for j in range(i,3*self.Mass.size):
                        self.Fc[i,j] = self.Fc[j,i]

        elif tokens[0].upper() == "FREQ" :
             line = readGeomFc.readMeaningfulLine(file)
             numFreq = int(line.split()[0])
             i = 0
             while (i < numFreq):
                 line = readGeomFc.readMeaningfulLine(file)
                 tokens = line.split()
                 i = i+ len(tokens)
                 for j in tokens:
                     self.Freq.append(float(j))
 
             if len(self.Freq) > numFreq:
                 print 'More frequencies than ', numFreq, ' are specified'

        else:
            print 'Frequency information cannot be read, check input file again'
            exit()

#read energy 
        line = readGeomFc.readMeaningfulLine(file)
        tokens = line.split()
        if (tokens[0].upper() != 'ENERGY'):
            print 'Energy information not given'
            exit()
        if tokens[1].upper() == 'FILE' :
            print 'Reading energy from file: ', tokens[2]
            energyFile = open(tokens[2],'r')
            if (tokens[3].upper() == 'CBS-QB3'):
               self.Etype = 'cbsqb3'
            elif (tokens[3].upper() == 'G3'):
               self.Etype = 'g3'
            elif (tokens[3].upper() == 'KLIP_1'):
               self.Etype = 'klip_1'
            elif (tokens[3].upper() == 'KLIP_2'):
               self.Etype = 'klip_2'
            elif (tokens[3].upper() == 'KLIP_2_CC'):
               self.Etype = 'klip_2_cc'
            self.Energy = readGeomFc.readEnergy(energyFile, self.Etype)
            print self.Energy, self.Etype
        elif (len(tokens) == 3):
            self.Energy = float(tokens[1])
            if (tokens[2].upper() == 'CBS-QB3'):
               self.Etype = 'cbsqb3'
            elif (tokens[2].upper() == 'G3'):
               self.Etype = 'g3'
            elif (tokens[2].upper() == 'Klip_1'):
               self.Etype = 'klip_1'
            elif (tokens[2].upper() == 'Klip_2'):
               self.Etype = 'klip_2'
            elif (tokens[2].upper() == 'KLIP_2_CC'):
               self.Etype = 'klip_2_cc'
            print self.Etype.upper(),' Energy: ',self.Energy
        else :
            print 'Cannot read the Energy'
            exit()

#read external symmetry
        line = readGeomFc.readMeaningfulLine(file)
        if (line.split()[0].upper() != 'EXTSYM'):
           print 'Extsym keyword required'
           exit()
        self.extSymm = int(line.split()[1])

#read electronic degeneracy
        line = readGeomFc.readMeaningfulLine(file)
        if (line.split()[0].upper() != 'NELEC'):
           print 'Nelec keyword required'
           exit()
        self.nelec = int(line.split()[1])

#read rotor information     
        line = readGeomFc.readMeaningfulLine(file)
        if (line.split()[0].upper() != 'ROTORS'):
           print 'Rotors keyword required'
           exit()
        self.numRotors = int(line.split()[1])
        if self.numRotors == 0:
           #Line added by MRH o 2/Dec/2009
           #If no rotors exist and molecule is not a single atom, still need to:
           #   (1) Check if frequencies are already computed
           #   (2) Read in the BAC information
           if self.linearity != 'Atom':
               self.rotors = []
               if (len(self.Freq) == 0):
               #calculate frequencies from force constant
                   self.getFreqFromFc()
               line = readGeomFc.readMeaningfulLine(file)
               tokens = line.split()
               for bond in tokens:
                   self.bonds.append(float(bond))
               return
           # If the molecule is a single atom, no rotors exist, no frequencies need
           #    calculating, and no BAC information needs to be read
           else:
               return

        rotorFile = line.split()[2]
        inertiaFile = open(rotorFile,'r')
        #print self.Mass
        (self.rotors)= readGeomFc.readGeneralInertia(inertiaFile,self.Mass)
        if len(self.rotors)-1 != self.numRotors :
            print "The number of rotors specified in file, ",rotorFile,' is different than, ',self.numRotors

        if (len(self.Freq) == 0):
            #calculate frequencies from force constant
            #Added by MRH on 2/Dec/2009
            if self.linearity != 'Atom':
                self.getFreqFromFc()


#read potential information for rotors
        line = readGeomFc.readMeaningfulLine(file)
        tokens = line.split()
        if tokens[0].upper() != 'POTENTIAL': 
            print 'No information for potential given'
            exit()

        if tokens[1].upper() == 'SEPARABLE':
            if tokens[2].upper() == 'FILES':
               line = readGeomFc.readMeaningfulLine(file)
               tokens = line.split()
               if len(tokens) != self.numRotors :
                   print 'give a separate potential file for each rotor'
               for files in tokens:
                   Kcos=[]
                   Ksin =[]
                   harmonic = Harmonics(5,Kcos,Ksin)
                   harmonic.fitPotential(files)
                   self.Harmonics.append(harmonic)

            elif tokens[2].upper() == 'HARMONIC':
               for i in range(self.numRotors):
                  line = readGeomFc.readMeaningfulLine(file)
                  numFit = int(line.split()[0])
                  Ksin = []
                  Kcos = []
                  for i in range(numFit):
                      line = readGeomFc.readMeaningfulLine(file)
                      tokens = line.split()
                      Kcos.append(tokens[0])
                      Ksin.append(tokens[0])
                  harmonic = Harmonics(numFit,Kcos,Ksin)
                  self.Harmonics.append(harmonic)

        elif tokens[1].upper() == 'NONSEPARABLE':
             line = readGeomFc.readMeaningfulLine(file)
             self.potentialFile = line.split()[0]
# read the energy base
        #line = readGeomFc.readMeaningfulLine(file)
        #tokens = line.split()
        #if (tokens[0].upper() != 'ENERGYBASE'):
        #    print 'Keyword EnergyBase required'
        #    exit()
        #self.ebase=float(tokens[1])
# read the bonds
        line = readGeomFc.readMeaningfulLine(file)
        tokens = line.split()
        for bond in tokens:
          self.bonds.append(float(bond))
        #read the random seed number
        '''line = readGeomFc.readMeaningfulLine(file)
        if line.split()[0].upper() != 'RANDOMSEED':
            print 'Keyword RandomSeed required'
            exit()
        self.mcSeed = float(line.split()[1])

        #read the # mciter
        line = readGeomFc.readMeaningfulLine(file)
        if line.split()[0].upper() != 'ITER':
            print 'Keyword ITER required'
            exit()
        self.mcIter = int(line.split()[1])'''

#*************************************************************************


    def getFreqFromFc(self):
        #Following added by MRH on 2/Dec/2009
        if self.linearity == 'Atom':
            return
        
        Fc = self.Fc.copy()
        rotors = self.rotors
        numRotors = self.numRotors
        geom = self.geom
        Mass = self.Mass

        if numRotors >0:
           intRotMatrix=matrix(array(zeros((3*Mass.size,numRotors),dtype=float)))

        #inttranrot = matrix(zeros((3*Mass.size,6+numRotors),dtype=float))

        #form cartesian vectors for all rotors
        for i in range(numRotors):
           rotor = rotors[i+1]
           e12=matrix('0 0 0');
           e21=matrix('0 0 0');
           e12=geom[rotor.pivotAtom-1,:]-geom[rotor.pivot2-1,:]
           e12=e12/linalg.norm(e12)
           e21=-e12
           atoms1 = rotor.atomsList
           for j in atoms1:
              e31 = geom[j-1,:]-geom[rotor.pivotAtom-1,:]
              intRotMatrix[3*(j-1):3*j,i]=transpose(cross(e31,e12))

        #make all the modes of unit length
        for i in range(numRotors):
              intRotMatrix[:,i]=intRotMatrix[:,i]/linalg.norm(intRotMatrix[:,i])

        #make the int Rotors Orthonormal
        if numRotors >0 :
              intRotMatrix = matrix(scipy.linalg.orth(intRotMatrix))


        #make translation and rotation unit vectors
        #if self.linearity == 'Nonlinear' :
        if self.linearity != 'Atom' :
            tranrot = matrix(zeros((3*Mass.size,6),dtype=float))
            for i in range(Mass.size):
                tranrot[3*i,0] = 1.0
                tranrot[3*i+1,1] = 1.0
                tranrot[3*i+2,2] = 1.0
                tranrot[3*i:3*i+3,3] = transpose(matrix([0, -geom[i,2], geom[i,1]]))
                tranrot[3*i:3*i+3,4] = transpose(matrix([geom[i,2], 0, -geom[i,0]]))
                tranrot[3*i:3*i+3,5] = transpose(matrix([-geom[i,1], geom[i,0], 0]))
        #elif self.linearity == 'Linear' :
        #    zeroArray = matrix(zeros((Mass.size,1),dtype=float))
        #    #If the vector of z-coordinates is zero, the molecule is not aligned along the z-axis
        #    # If this is so, re-align molecule along the z-axis
        #    if (geom[:,2] == zeroArray).all() :
        #        if (geom[:,0] == zeroArray).all() :
        #            geom[:,2] = geom[:,1]
        #            geom[:,1] = zeroArray
        #        else :
        #            geom[:,2] = geom[:,0]
        #            geom[:,0] = zeroArray
        #    tranrot = matrix(zeros((3*Mass.size,5),dtype=float))
        #    for i in range(Mass.size) :
        #        tranrot[3*i,0] = 1.0
        #        tranrot[3*i+1,1] = 1.0
        #        tranrot[3*i+2,2] = 1.0
        #        tranrot[3*i:3*i+3,3] = transpose(matrix([0, geom[i,2], 0]))
        #        tranrot[3*i:3*i+3,4] = transpose(matrix([geom[i,2], 0, 0]))
        tranrot = matrix(scipy.linalg.orth(tranrot))

        if self.linearity == 'Nonlinear':
            inttranrot = matrix(zeros((3*Mass.size,6+numRotors),dtype=float))
        elif self.linearity == 'Linear':
            inttranrot = matrix(zeros((3*Mass.size,5+numRotors),dtype=float))
        #elif self.linearity == 'Atom':
        #    inttranrot = matrix(zeros((3*Mass.size,3+numRotors),dtype=float))

        if numRotors >0 :
           inttranrot[:,0:numRotors]=intRotMatrix
        
        if self.linearity == 'Nonlinear':
            inttranrot[:,numRotors:numRotors+6]=tranrot
        elif self.linearity == 'Linear':
            inttranrot[:,numRotors:numRotors+5]=tranrot
        #elif self.linearity == 'Atom':
        #    inttranrot[:,numRotors:numRotors+3]=tranrot

        inttranrot = matrix(scipy.linalg.orth(inttranrot))

        P = inttranrot*transpose(inttranrot)
        I = matrix(eye(3*Mass.size,3*Mass.size))

        Fc = (I-P)*Fc*(I-P)

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

        if self.TS:
           self.imagFreq = -1*sqrt(-1*l[0] * (627.5095*4180/6.023e23)  * (1.88972e10**2) * (1/1.67e-27)) /2/math.pi/3e10
        
        if self.linearity == 'Nonlinear':
            l = l[6+numRotors+int(self.TS):]
            for i in range(3*Mass.size-6-numRotors - int(self.TS)):
                self.Freq.append(sqrt(l[i] * (627.5095*4180/6.023e23)  * (1.88972e10**2) * (1/1.67e-27)) /2/math.pi/3e10)
        elif self.linearity == 'Linear':
            l = l[5+numRotors+int(self.TS):]
            for i in range(3*Mass.size-5-numRotors - int(self.TS)):
                self.Freq.append(sqrt(l[i] * (627.5095*4180/6.023e23)  * (1.88972e10**2) * (1/1.67e-27)) /2/math.pi/3e10)
        #elif self.linearity == 'Atom':
        #    l = l[3+numRotors+int(self.TS):]

        Fc = self.Fc.copy()
        Fc = P*Fc*P
        Fc = Tcmc*Fc*Tcmc
        [l,v]=linalg.eigh(Fc)
        l = sort(l)
        l = l[-numRotors:]
        for i in range(numRotors):
            self.hindFreq.append(sqrt(l[i] * (627.5095*4180/6.023e23)  * (1.88972e10**2) * (1/1.67e-27)) /2/math.pi/3e10)

        #for i in range(len(l)/3+1):
        #    for j in range(3):
        #      if 3*i+j <len(l):
        #        print '%10.3f'%(sqrt(l[3*i+j] * (627.5095*4180/6.023e23)  * (1.88972e10**2) * (1/1.67e-27) )/2/math.pi/3e10),
        #    print 
 
        #print



#*************************************************************************

    def printData(self,oFile):
       geom = self.geom
       Mass = self.Mass
       oFile.write('Geometry:\n')
       oFile.write('%10s'%'Mass(amu)'+'%10s'%'X(ang)'+'%10s'%'Y(ang)'+'%10s'%'Z(ang)'+'\n')
       for i in range(Mass.size):
          oFile.write('%10.3f'%float(Mass[i])+'%10.4f'%float(geom[i,0])+'%10.4f'%float(geom[i,1])+'%10.4f'%float(geom[i,2])+'\n')

       oFile.write('\nFrequencies (cm-1):\n')
       Freq = self.Freq
       if self.TS:
	oFile.write('Imaginary Frequency: '+str(self.imagFreq)+'\n')
       for i in range(len(Freq)/3+1):
         for j in range(3):
            if 3*i+j <len(Freq):
              oFile.write('%10.3f'%Freq[3*i+j])
         oFile.write('\n') 
       oFile.write('\nExternal Symmetry = '+str(self.extSymm)+'\n')
       oFile.write('Principal Moments of Inertia = '+str(self.Iext[0])+'  '+str(self.Iext[1])+'  '+str(self.Iext[2])+'\n')
       oFile.write('Electronic Degeneracy = '+str(self.nelec)+'\n')

       if self.numRotors == 0:
          return       
       oFile.write('\nFitted Harmonics V(p) = sum (A_i cos(i*p) + B_i sin(i*p)) :\n' )

       k = 1
       K = geomUtility.calculateD32(self.geom,self.Mass,self.rotors)
       for harmonic in self.Harmonics:
           oFile.write('Harmonic '+str(k)+'\n')
           oFile.write('Moment of Inertia: '+ str(float(K[k-1]))+'\n')
           oFile.write('Symmetry '+str(self.rotors[k].symm)+'\n')
           oFile.write('BarrierHeight '+str(harmonic.A)+'\n')
           oFile.write('%12s'%'A_i'+'%12s'%'B_i'+'\n')
           for j in range(5):
              oFile.write('%12.3e'%harmonic.Kcos[j]+'%12.3e'%harmonic.Ksin[j]+'\n')
           oFile.write('\n')
           k = k+1

#*************************************************************************

    def getTranslationThermo(self,oFile,Temp):
        ent = []
        cp = []
        dH = []
        qtrans = []
        R = self.R
        kb = self.kb
        h = self.h
        amu = self.amu
        oFile.write('Translational Contributions\n')
        oFile.write('%12s'%'Temperature')
        for T in Temp:
           oFile.write('%12.2f'%T)
        oFile.write('\n%12s'%'Entropy')
        i = 0
        for T in Temp:
           ent.append(R*math.log( (2.0*math.pi*self.Mass.sum()*amu*kb*T/h**2)**(1.5) * (kb*T*math.e**(2.5)/1.013e5) ))
           oFile.write('%12.2f'%ent[i])
           i = i+1

        oFile.write('\n%12s'%'Cp')
        i = 0
        for T in Temp:
           cp.append(5.0/2*R)
           oFile.write('%12.2f'%cp[i])
           i = i+1

        oFile.write('\n%12s'%'dH')
        i = 0
        for T in Temp:
           dH.append(5.0/2*R*T/1000.0)
           oFile.write('%12.2f'%dH[i])
           i = i+1
        #oFile.write('\n')

        #MRH
        oFile.write('\n%12s'%'qtrans')
        i = 0
        for T in Temp:
           qtrans.append((2.0*math.pi*self.Mass.sum()*amu*kb*T/h**2)**(1.5) * (kb*T/1.013e5))
           oFile.write('%12.2e'%qtrans[i])
           i = i+1
        oFile.write('\n')

        return ent, cp, dH, qtrans

#*************************************************************************

    def getVibrationalThermo(self,oFile,Temp,scale):
        ent = []
        cp = []
        dH = []
        parti = []
        R = self.R
        kb = self.kb
        h = self.h
        amu = self.amu
        oFile.write('\nVibrational Contributions\n')
        Freq=[] 
        for freq in self.Freq:
           Freq.append(freq*scale)

        oFile.write('Entropy:\n')
        oFile.write('%12s'%'Temperature')

        for T in Temp:
           oFile.write('%12.2f'%T)
           ent.append(0.0)
           parti.append(1.0)
        oFile.write('\n')

        #get vibrational contribution to entropy
        j = 0
        for freq in Freq:
           i = 0
           f = 'Freq: '+'%2.0f'%(j+1)
           #oFile.write('%12s'%f)
           for T in Temp:
               s = -R*math.log( 1.0 - math.exp(-h*freq*3.0e10/kb/T))
               s = s + 6.023e23*(h*freq*3.0e10/T) / ( math.exp(h*freq*3.0e10/kb/T) - 1.0 )/4.18     
               #oFile.write('%12.2f'%s)
               ent[i] = ent[i] +s
               parti[i] = parti[i] * 1.0/(1.0 - math.exp(-h*freq*3.0e10/kb/T))
               i = i+1
               #print T, freq, s
           #oFile.write('\n')
           j = j+1
        oFile.write('%12s'%'Entropy')
        for i in range(len(Temp)):
            oFile.write('%12.2f'%ent[i])
        oFile.write('\n')

        
        #get vibrational contribution to cp
        #oFile.write('\nCp:\n')
        #oFile.write('%12s'%'Temperature')
        for T in Temp:
           #oFile.write('%12.2f'%T)
           cp.append(0.0)
        #oFile.write('\n')

        j = 0
        for freq in Freq:
           i = 0
           f = 'Freq: '+'%2.0f'%(j+1)
           #oFile.write('%12s'%f)
           for T in Temp:
               c = R * (h*freq*3.0e10/kb/T)**2 * math.exp(h*freq*3.0e10/kb/T)/( 1.0 - math.exp(h*freq*3.0e10/kb/T))**2
               #oFile.write('%12.2f'%c)
               cp[i] = cp[i] +c
               i = i+1
           #oFile.write('\n')
           j = j+1
        oFile.write('%12s'%'Cp')
        for i in range(len(Temp)):
            oFile.write('%12.2f'%cp[i])
        oFile.write('\n')


        #get vibrational conbribution to thermal correction
        #oFile.write('\ndH\n')
        #oFile.write('%12s'%'Temperature')
        for T in Temp:
           #oFile.write('%12.2f'%T)
           dH.append(0.0)
        #oFile.write('\n')

        j = 0
        for freq in Freq:
           i = 0
           f = 'Freq: '+'%2.0f'%(j+1)
           #oFile.write('%12s'%f)
           for T in Temp:
               h1 = 6.023e23*(h*freq*3.0e10)/( math.exp(h*freq*3.0e10/kb/T) - 1.0 )/4180.0
               #oFile.write('%12.2f'%h1)
               dH[i] = dH[i] +h1
               i = i+1
           #oFile.write('\n')
           j = j+1
        oFile.write('%12s'%'dH')
        for i in range(len(Temp)):
            oFile.write('%12.2f'%dH[i])
        oFile.write('\n')

        oFile.write('%12s'%'q_vib')
        for i in range(len(Temp)):
            oFile.write('%12.2e'%parti[i])
        oFile.write('\n')


        return ent, cp, dH, parti

#*************************************************************************

    def getIntRotationalThermo_PG(self,oFile,Temp):
         ent = []
         cp = []
         dH = []
         R = self.R
         kb = self.kb
         h = self.h
         amu = self.amu
         seed = 500
         numIter = 100000

         oFile.write('\n\nInternal Rotational Contributions\n')
         oFile.write('%12s'%'Temperature')
    
         #iterate over temperatures
         for T in Temp:
             oFile.write('%12.2f'%T)

         sigma = 1.0
         for rotor in self.rotors:
             sigma = sigma*rotor.symm

         K = geomUtility.calculateD32(self.geom,self.Mass,self.rotors)
         print K
         p = 1.0
         a = 1.0
         for T in Temp:
             #print 'Calculating rotational entropy for T: ',T
             Sq = 0.0
             Scl = 0.0
             S = 0.0
             Hq = 0.0
             Hcl = 0.0
             H = 0.0
             cpcl = 0.0
             cpq = 0.0
             Cp = 0.0
             for l in range(self.numRotors):                   
                sum = 0.0
                vsumexpv = 0.0
                v2sumexpv = 0.0
                minpot = 5.0
                for i in range(numIter):
                   ang = random.rand()
                   pot = self.Harmonics[l].getPotential(ang*360)
                   if (pot < minpot):
                      minpot = pot

                for i in range(100):
                   ang = i*360.0/100
                   pot = self.Harmonics[l].getPotential(ang)-minpot
		   print pot
                   exit()

#                   fi = sqrt(K[l])*exp(-pot*1.0e3/R/T) 
                   fi = exp(-pot*1.0e3/R/T) 
                   sum = sum + fi
                   vsumexpv = vsumexpv + pot*1.0e3 * fi
                   v2sumexpv = v2sumexpv + pot**2*1.0e6 *fi
                   average = sum/(i+1)
                   parti =  (2.0*math.pi*kb*T*amu*1e-20/h**2)**(0.5) *(2*math.pi) * average/self.rotors[l+1].symm

                a = a*average
                S = S + R*math.log(parti) + R/2 + vsumexpv/sum/T
		H = H + R*T/2.0 + vsumexpv/sum #reference to the minimum of the well
                Cp = Cp + R/2.0 + (v2sumexpv*sum-vsumexpv**2)/sum**2/R/T**2

             sumfreq = 0.0
             for k in range(len(self.hindFreq)):
                harm = self.Harmonics[k]
                ddv = 0.0
                for l in range(5):
                  ddv = ddv-1*harm.Kcos[l]*(l+1)**2
                freq = 1.0/2.0/pi*sqrt(ddv*4180.0/K[k]/1.0e-20/amu/6.023e23)/3.0e10
                #print freq, K[l], pi, ddv, K

                Sq = Sq -R*math.log( 1.0 - math.exp(-h*freq*3.0e10/kb/T)) + 6.023e23*(h*freq*3.0e10/T) / ( math.exp(h*freq*3.0e10/kb/T) - 1.0 )/4.18
                Scl = Scl + R + R*math.log(kb*T/h/freq/3.0e10)
                Hq = Hq + 6.023e23*(h*freq*3.0e10)/( math.exp(h*freq*3.0e10/kb/T) - 1.0 )/4.18
                Hcl = Hcl + R*T
                cpq = cpq + R * (h*freq*3.0e10/kb/T)**2 * math.exp(h*freq*3.0e10/kb/T)/( 1.0 - math.exp(h*freq*3.0e10/kb/T))**2
                cpcl = cpcl + R
                sumfreq = sumfreq +freq

             H = H + Hq - Hcl
             S = S + Sq - Scl
             Cp = Cp + cpq - cpcl

             ent.append(S)
             dH.append(H/1e3)
             cp.append(Cp)

         oFile.write('\n%12s'%'Entropy')
         for i in range(len(Temp)):
            oFile.write('%12.2f'%ent[i])

         oFile.write('\n%12s'%'Cp')
         for i in range(len(Temp)):
            oFile.write('%12.2f'%cp[i])

         oFile.write('\n%12s'%'dH')
         for i in range(len(Temp)):
            oFile.write('%12.2f'%dH[i])
         #oFile.write('\n')

         return ent, cp, dH

#**************************************************************************

#*************************************************************************

    def getIntRotationalThermo_Q(self,oFile,Temp):
         ent = [0.0]*len(Temp)
         cp = [0.0]*len(Temp)
         dH = [0.0]*len(Temp)
         parti = [1.0]*len(Temp)
         R = self.R
         kb = self.kb
         h = self.h
         amu = self.amu

         oFile.write('\n\nInternal Rotational Contributions\n')
         oFile.write('%12s'%'Temperature')
    
         #iterate over temperatures
         for T in Temp:
             oFile.write('%12.2f'%T)

         sigma = 1.0
         for rotor in self.rotors:
             sigma = sigma*rotor.symm

         K = geomUtility.calculateD32(self.geom,self.Mass,self.rotors)
         for irot in range(len(self.rotors)-1):
	        harm = self.Harmonics[irot]
                ddv = 0.0
                for l in range(5):
                  ddv = ddv-1*harm.Kcos[l]*(l+1)**2
                freq = 1.0/2.0/pi*sqrt(ddv*4180.0/K[irot]/1.0e-20/amu/6.023e23)/3.0e10
                #print '%12.2f'%float(freq),
         #print

         #calculate the energy levels for the hindered rotors
         E = self.calculateElevels()
         #pdb.set_trace()
         for iT in range(len(Temp)):
            T = Temp[iT]
            #print T,
            for irot in range(len(self.rotors)-1):
               sum = 0.0
               vsum = 0.0
               v2sum = 0.0

               for e in E[irot]:
                 e = e - E[irot][0]
                 sum = sum + exp(-e*1.0e3/R/T)
                 vsum = vsum + e*1e3*exp(-e*1.0e3/R/T)
                 v2sum = v2sum + e**2*1e6*exp(-e*1.0e3/R/T)

               ent[iT] = ent[iT] + R*math.log(sum)+vsum/sum/T-R*log(self.rotors[irot+1].symm)
               dH[iT] = dH[iT] + vsum/sum/1.0e3
               cp[iT] = cp[iT] + (v2sum*sum-vsum**2)/sum**2/R/T**2
               parti[iT] = parti[iT] *sum/self.rotors[irot+1].symm
               
               #print (v2sum*sum-vsum**2)/sum**2/R/T**2,

               #print R*math.log(sum)+vsum/sum/T-R*log(self.rotors[irot+1].symm),
            #print 
         oFile.write('\n%12s'%'Entropy')
         for i in range(len(Temp)):
            oFile.write('%12.2f'%ent[i])

         oFile.write('\n%12s'%'Cp')
         for i in range(len(Temp)):
            oFile.write('%12.2f'%cp[i])

         oFile.write('\n%12s'%'dH')
         for i in range(len(Temp)):
            oFile.write('%12.2f'%dH[i])

         oFile.write('\n%12s'%'q_int')
         for i in range(len(Temp)):
            oFile.write('%12.2e'%parti[i])
         oFile.write('\n')


         return ent, cp , dH, parti

#**************************************************************************

    def calculateElevels(self):
         R = self.R

         kb = self.kb
         h = self.h
         amu = self.amu
         K = geomUtility.calculateD32(self.geom,self.Mass,self.rotors)
         E=[]
         #let us take k = -500, 500 
         m = 200
         #print K
         for irot in range(len(self.Harmonics)):
          H = mat(zeros((2*m+1,2*m+1),dtype=complex))
          kcos = self.Harmonics[irot].Kcos
          ksin = self.Harmonics[irot].Ksin

          for k in range(0,2*m+1):
             H[k,k] = 6.023e23*h**2*(k-m)**2/8.0/math.pi**2/K[irot]/amu/1e-20/4180 + self.Harmonics[irot].A

             for n in range(1,6):
                if k-n >= 0:
                  H[k,k-n] = kcos[n-1]/2 + ksin[n-1]/2j
                if k+n < 2*m+1:
                  H[k,k+n] = kcos[n-1]/2 - ksin[n-1]/2j
          (l,v) = linalg.eigh(H)
          #pdb.set_trace()
          E.append(l)
         return E

#**************************************************************************

    def getExtRotationalThermo(self,oFile,Temp):
       kb = self.kb
       h = self.h
       amu = self.amu
       R = self.R
 
       S = []
       ent = []
       cp = []
       dH = []        
       q = []
       q_rot = []

       oFile.write('\n\nExternal Rotational Contributions\n')
       oFile.write('%12s'%'Temperature')
    
         #iterate over temperatures
       for T in Temp:
           oFile.write('%12.2f'%T)

       if self.linearity == 'Nonlinear' :
         for T in Temp:
           S = log(math.pi**0.5*exp(1.5)/self.extSymm)
           q = math.pi**0.5/self.extSymm
           for j in range(3):
             S = S + log( (8*math.pi**2*self.Iext[j]*kb*T*amu*1e-20/h**2)**0.5)
             q = q * (8*math.pi**2*self.Iext[j]*kb*T*amu*1e-20/h**2)**0.5
           ent.append(S*R)
           cp.append(3.0*R/2.0)
           dH.append(3.0*R*T/2.0/1.0e3)
           q_rot.append(q)
       elif self.linearity == 'Linear' :
         for T in Temp:
           I = self.Iext[2]
           S = log(8*math.pi**2*I*kb*T*amu*1e-20*exp(1)/self.extSymm/h**2)
           q = 8*math.pi**2*I*kb*T*amu*1e-20 / h**2 / self.extSymm
           ent.append(S*R)
           cp.append(R)
           dH.append(R*T/1.0e3)
           q_rot.append(q)
       else :
         for T in Temp:
           ent.append(0.0)
           cp.append(0.0)
           dH.append(0.0)
           q_rot.append(1.0)

       oFile.write('\n%12s'%'Entropy')
       for i in range(len(Temp)):
          oFile.write('%12.2f'%ent[i])

       oFile.write('\n%12s'%'Cp')
       for i in range(len(Temp)):
          oFile.write('%12.2f'%cp[i])

       oFile.write('\n%12s'%'dH')
       for i in range(len(Temp)):
          oFile.write('%12.2f'%dH[i])

       oFile.write('\n%12s'%'q_rot')
       for i in range(len(Temp)):
          oFile.write('%12.2f'%q_rot[i])
       oFile.write('\n\n\n\n')
       return ent, cp, dH, q_rot


#**************************************************************************
    def calculateMomInertia(self):
       geom = self.geom
       Mass = self.Mass
       #change coordinates to have cm
       cm = matrix('0.0 0.0 0.0')
    
       for i in range(Mass.size):
         cm = cm+Mass[i]*geom[i,:]
        
       cm = cm/sum(Mass)
        
       for i in range(Mass.size):
         geom[i,:]=geom[i,:]-cm



#calculate moments of inertia
       I = matrix(zeros((3,3),dtype=double))
       x = array(geom[:,0])
       y = array(geom[:,1])
       z = array(geom[:,2])
       I[0,0]=sum(array(Mass)*(y*y+z*z))
       I[1,1]=sum(array(Mass)*(x*x+z*z))
       I[2,2]=sum(array(Mass)*(x*x+y*y))
       I[0,1]=I[1,0]=-sum(array(Mass)*x*y)
       I[0,2]=I[2,0]=-sum(array(Mass)*x*z)
       I[1,2]=I[2,1]=-sum(array(Mass)*z*y)

#rotate coordinate axes to be parallel to principal axes
       (l,v)=linalg.eigh(I)
       self.Iext = l; 


'''
#***************************************************************************
    def getRotationalThermo_multi_simple(self,oFile,Temp):
         ent = []
         cp = []
         dH = []
         R = self.R
         kb = self.kb
         h = self.h
         amu = self.amu
         seed = self.mcSeed
         numIter = self.mcIter
         
         integralFile = open('integrals.out','w')
         
         atoms = readGeomFc.getAtoms(self.Mass)

         oFile.write('Rotational Contributions\n')
         oFile.write('%12s'%'Temperature')
    
         #iterate over temperatures
         for T in Temp:
             oFile.write('%12.2f'%T)

         oFile.write('\n%12s'%'Entropy')
         sigma = 1.0
         for rotor in self.rotors:
             sigma = sigma*rotor.symm
       
         K = geomUtility.calculateI23(self.geom,self.Mass,self.rotors)
         (D,v) = linalg.eigh(K)
         detD = linalg.det(K)

         stepSize = 360.0
         diheds = self.numRotors*[0.0]
         presentStep = self.numRotors*[0]
         ent = []
         currentDi = self.numRotors*[0]
         for T in Temp:
             sum = 0.0
             vsumexpv = 0.0
             print 'Calculating rotational entropy for T: ',T
             Sq = 0.0
             Scl = 0.0
             S = 0.0
             E0 = 0.0
             acceptedN = 1.0
             for i in range(numIter):
                   if i/10000*10000 - i == 0:
                      print "performed ",str(i), " iterations." 
                   n = self.numRotors
                   pot = 0.0
                      
                   for j in range(n):
                      presentStep[j] = (random.random()-0.0)*stepSize                
                      pot = pot + self.Harmonics[j].getPotential(presentStep[j])
                   
                   
                   #newGeom = geomUtility.rotateDihedrals(self.geom,self.rotors,currentDi,self.Mass)
                   #K = geomUtility.calculateI23(newGeom,self.Mass,self.rotors)
                   #detD = linalg.det(K)

                   if detD < 0.0:
                      print detD, diag(K)
                      exit()
 
                   acceptedN=acceptedN+1
                   fi = sqrt(detD)*exp(-pot*1e3/R/T) 
                   sum = sum + fi
                   vsumexpv = vsumexpv + pot*1e3 * fi
                   
                   average = sum/(acceptedN)
                   parti =  (2.0*math.pi*kb*T*amu*1e-20/h**2)**(n*0.5) *(2*math.pi)**n * average/sigma
                   ent.append(R*math.log(parti) + n*R/2 + vsumexpv/sum/T)
                   integralFile.write('%12.5f'%ent[int(acceptedN)-2]+'%12.2e'%(average/sqrt(detD))+' '+str(i-i/n*n)+'\n') 
                   integralFile.flush()

             g=Gnuplot.Gnuplot()
             g('set data style linespoints')
             g.plot(ent)
             raw_input('Please press enter to continue ...\n')

             S = R*math.log(parti) + n*R/2 + vsumexpv/sum/T

             for freq in self.hindFreq:
                Sq = Sq -R*math.log( 1.0 - math.exp(-h*freq*3.0e10/kb/T)) + 6.023e23*(h*freq*3.0e10/T) / ( math.exp(h*freq*3.0e10/kb/T) - 1.0 )/4.18
                Scl = Scl + R + R*math.log(kb*T/h/freq/3.0e10)
             S = S + Sq - Scl
             oFile.write('%12.2f'%S)
             print S, R*math.log(parti), vsumexpv/sum/T, parti, average, acceptedN


         oFile.write('\n')
#**************************************************************************


    def getRotationalThermo_single_metro(self,oFile,Temp):
         ent = []
         cp = []
         dH = []
         R = self.R
         kb = self.kb
         h = self.h
         amu = self.amu
         seed = self.mcSeed
         numIter = self.mcIter
         random.seed(seed)
         dihedrals = self.numRotors*numIter*[0]
         for i in range(self.numRotors):
              for j in range(numIter):
                   dihedrals[numIter*i+j] = 360*random.random()
         #print dihedrals
         oFile.write('Rotational Contributions\n')
         oFile.write('%12s'%'Temperature')
    
         #iterate over temperatures
         for T in Temp:
             oFile.write('%12.2f'%T)

         oFile.write('\n%12s'%'Entropy')
         sigma = 1.0
         for rotor in self.rotors:
             sigma = sigma*rotor.symm

         K = geomUtility.calculateI23(self.geom,self.Mass,self.rotors)
         (detD,v) = linalg.eigh(K)
         print (detD)
         p = 1.0
         a = 1.0
         for T in Temp:
             print 'Calculating rotational entropy for T: ',T
             Sq = 0.0
             Scl = 0.0
             S = 0.0
             for l in range(self.numRotors):                   
                sum = 0.0
                vsumexpv = 0.0
                ent=[]
                pot_prev = 0.0
                for i in range(numIter):

                   if detD[l] < 0.0:
                      print detD, diag(K)
                      exit()
                   
                   pot = self.Harmonics[l].getPotential(random.rand()*360)
                   #print pot, 

                   fi = sqrt(detD[l])*exp(-pot*1e3/R/T) 
                   sum = sum + fi
                   vsumexpv = vsumexpv + pot*1e3 * fi

                   average = sum/(i+1)
                   parti =  (2.0*math.pi*kb*T*amu*1e-20/h**2)**(0.5) *(2*math.pi) * average/self.rotors[l+1].symm
                   ent.append(R*math.log(parti) + R/2 + vsumexpv/sum/T)
                a = a*average
                S = S + R*math.log(parti) + R/2 + vsumexpv/sum/T
                print S, R*math.log(parti) + R/2 + vsumexpv/sum/T

                g=Gnuplot.Gnuplot()
                g('set data style linespoints')
                g.plot(ent)
                raw_input('Please press enter to continue ...\n')


             for freq in self.hindFreq:
                Sq = Sq -R*math.log( 1.0 - math.exp(-h*freq*3.0e10/kb/T)) + 6.023e23*(h*freq*3.0e10/T) / ( math.exp(h*freq*3.0e10/kb/T) - 1.0 )/4.18
                Scl = Scl + R + R*math.log(kb*T/h/freq/3.0e10)
             S = S + Sq - Scl
             oFile.write('%12.2f'%S)
             print S
         oFile.write('\n')
        
#**************************************************************************
'''           
