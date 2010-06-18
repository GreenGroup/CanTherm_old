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


#takes the file as an argument and then reads the geometry and lower triangular part
#of the force constant matrix and the Atomic mass vector

import os 
from numpy import *
from Rotor import *
import pdb
import CanTherm
from Molecule import *
import re,shutil,getopt

def readInputFile(file,data):

#read calculation type
# Will either be 'THERMO' or 'REAC'
    line = readMeaningfulLine(file)
    if line.split()[1].upper() == 'THERMO' :
        data.CalcType = 'Thermo'
    elif line.split()[1].upper() == 'REAC' :
        data.CalcType = 'Reac'
    else:
        print 'First line of the input file is neither Reac or Thermo'
        exit()

#read reaction & tunneling type, if reaction calculation specified
    line = readMeaningfulLine(file)
    if data.CalcType == 'Reac':

        #read reaction type: unimolecular or bimolecular
        if line.split()[1].upper() == 'UNIMOL' :
            data.ReacType = 'Unimol'
        elif line.split()[1].upper() == 'BIMOL' :
            data.ReacType = 'Bimol'
        else :
            print 'ReactionType is not recognized.  Reaction field has options "Unimol" and "Bimol"'
            exit()
        line = readMeaningfulLine(file)

        #read tunneling type
        if line.split()[1].upper() == 'WIGNER' :
            data.TunnelType = 'Wigner'
        elif line.split()[1].upper() == 'SYMMECKART' :
            data.TunnelType = 'sEckart'
        elif line.split()[1].upper() == 'ASYMMECKART' :
            data.TunnelType = 'aEckart'
            data.numProds = int(line.split()[2])
        elif line.split()[1].upper() == 'NONE' :
            data.TunnelType = 'none'
        else :
            print 'TunnelType is not recognized.  Tunneling field has options "Wigner", "SymmEckart", and "AsymmEckart"'
            exit()
        line = readMeaningfulLine(file)
    
#read Temperature range
    tokens = line.split()
    if tokens[0].upper() == 'TLIST:':
        data.fitcp = False
        #line = readMeaningfulLine(file);
        numTemp = int(tokens[1])
        i = 0
        while i < numTemp:
            #line = readMeaningfulLine(file)
            #tokens = line.split()
            #i = i+ len(tokens)
            #for j in tokens:
                #data.Temp.append(float(j))
            data.Temp.append(float(tokens[i+2]))
            i = i+1
        if len(data.Temp) > numTemp:
            print 'More Temperatures than ', numTemp, ' are specified'

    elif tokens[0].upper() == 'TRANGE:':
        data.fitcp = False
        #line = readMeaningfulLine(file);
        #tokens = line.split()         
        T0 = float(tokens[1])
        dT = float(tokens[2])
        numTemp = int(tokens[3])
        for i in range(numTemp):
            data.Temp.append(T0+i*dT)
            
    elif line.split()[0].upper() == 'NASA':
        data.fitcp = True
        Tref = float('298')
        data.Temp.append(Tref)
        T0 = float('300')
        dT = float('50')
        for i in range(34):
            data.Temp.append(float(T0 + i*dT))
        T0 = float('2000')
        dT = float('250')
        for i in range(17):
            data.Temp.append(float(T0 + i*dT))


    else:
        print 'Temperaure information not given: Either use keyword Tlist or Trange'
        exit()


#read scaling factor
    line = readMeaningfulLine(file)
    if (line.split()[0].upper() != 'SCALE:'):
       print 'Give a scaling factor'
       exit()
    data.scale = float(line.split()[1])

# read molecular data
    if (data.CalcType == 'Thermo'):
       numMol = 1
    elif data.CalcType == 'Reac' and data.ReacType == 'Unimol':
       if data.TunnelType != 'aEckart' :
           numMol = 2
       else :
           numMol = 2 + data.numProds
       #numMol = 2
    elif data.CalcType == 'Reac' and data.ReacType == 'Bimol':
       if data.TunnelType != 'aEckart' :
           numMol = 3
       else :
           numMol = 3 + data.numProds
       #numMol = 3

    for i in range(numMol):
        if (data.ReacType=='Unimol' and i==1) or (data.ReacType=='Bimol' and i==2):
          #Molecule is a transition state
          molecule = Molecule(file, True)
        else:
          molecule = Molecule(file, False)	   
        data.MoleculeList.append(molecule)
    return 

def readMeaningfulLine(file):
    readMore = True
    while (readMore):
        line = file.readline()
        index = line.find('!')
        line = line[:index]
        if (len(line.split()) != 0):
            return line
    
def readInertia(file):
    lines = file.readlines()
    numAtoms = int(lines[0])
    numRotors = int(lines[numAtoms+1])
    pivots = []
    rotAtoms = []
    for i in range(numRotors):
       tokens = lines[numAtoms+2+3*i].split()
       pivots.append(int(tokens[0]))
       pivots.append(int(tokens[1]))
       atoms1=[]
       atoms2=[]
       atom1tokens = lines[numAtoms+2+3*i+1].split()
       atom2tokens = lines[numAtoms+2+3*i+2].split()
       for j in range(int(tokens[2])):
           atoms1.append(int(atom1tokens[j]))
       for j in range(int(tokens[3])):
           atoms2.append(int(atom2tokens[j]))
       rotAtoms.append(atoms1)
       rotAtoms.append(atoms2)
    return pivots, rotAtoms,numRotors


def readGeneralInertia(file,Mass):
    lines = file.readlines()
    rotors = []

    prevLevel=1
    i = 0
    parentLevel = 0 #gives the index of parent in rotors

    for line in lines:
        tokens = line.split()
        if (len(tokens) == 0):
            continue
        level = int(tokens[0][1])

        if (level!=1):

            if (rotors[parentLevel].level == level -2):
                parentLevel = i-1
            elif (rotors[parentLevel].level >= level):
                jumps = rotors[parentLevel].level - level +1
                for k in range(jumps):
                    parentLevel = rotors[parentLevel].parent
        symm = int(tokens[1])
        pivot2 = int(tokens[2])
        atomsList=[]
        for atoms in tokens[3:]:
            atomsList.append(int(atoms))
        rotor = Rotor(atomsList,pivot2,level,symm,Mass)
        rotor.parent = parentLevel
        #print rotor.symm, rotor.pivot2, rotor.pivotAtom, rotor.atomsList
        rotors.append(rotor) 
        i = i+1
    return rotors


def readGeom(file):
    lines = file.readlines()
    lines.reverse()
#read geometries

    i = lines.index("                          Input orientation:                          \n")
    geomLines = []
    stillRead = True
    k = i-5;
    while (stillRead):
       geomLines.append(lines[k])
       k = k-1
       if (lines[k].startswith(" ------------")):
          stillRead = False

    geom = matrix(array(zeros((len(geomLines),3),dtype=float)))
    Mass = matrix(array(zeros((len(geomLines),1),dtype=float)))

    for j in range(len(geomLines)):
        line = geomLines[j]
        tokens = line.split()
        geom[j,0]=double(tokens[3])
        geom[j,1]=double(tokens[4])
        geom[j,2]=double(tokens[5])
        Mass[j] = getMassByAtomicNumber(int(tokens[1]))
    return geom, Mass

    # Example:
#	      FINAL ATOMIC COORDINATE
#           ATOM          X           Y           Z      TYPE
#         C(    1)    -3.21470    -0.22058     0.00000   (  1)
#         H(    2)    -3.30991    -0.87175     0.89724   (  5)
#         H(    3)    -3.30991    -0.87174    -0.89724   (  5)
#         H(    4)    -4.08456     0.47380     0.00000   (  5)
#         C(    5)    -1.88672     0.54893     0.00000   (  1)
#         H(    6)    -1.84759     1.21197    -0.89488   (  5)
#         H(    7)    -1.84759     1.21197     0.89488   (  5)
#         C(    8)    -0.66560    -0.38447     0.00000   (  1)
#         H(    9)    -0.70910    -1.04707    -0.89471   (  5)
#         H(   10)    -0.70910    -1.04707     0.89471   (  5)
#         C(   11)     0.66560     0.38447     0.00000   (  1)
#         H(   12)     0.70910     1.04707     0.89471   (  5)
#         H(   13)     0.70910     1.04707    -0.89471   (  5)
#         C(   14)     1.88672    -0.54893     0.00000   (  1)
#         H(   15)     1.84759    -1.21197    -0.89488   (  5)
#         H(   16)     1.84759    -1.21197     0.89488   (  5)
#         C(   17)     3.21470     0.22058     0.00000   (  1)
#         H(   18)     3.30991     0.87174     0.89724   (  5)
#         H(   19)     4.08456    -0.47380     0.00000   (  5)
#         H(   20)     3.30991     0.87175    -0.89724   (  5)
def readMM4Geom(file):
    geomList = []
    MassList = []
    line = file.readline()
    while (not line[0:29] == '      FINAL ATOMIC COORDINATE'):#get to the beginning of the geometry section
	line = file.readline()
    line = file.readline() #header line (ATOM X Y Z TYPE)
    line = file.readline()
    while len(line.split()) > 0:
        MassList.append(getMassByAtomicSymbol(line[0:10].strip()))
	xc = float(line[17:29])
	yc = float(line[29:41])
	zc = float(line[41:53])
	geomList.append([xc,yc,zc])
	line = file.readline()
    #convert to the matrix format used by CanTherm
    geom = matrix(geomList)
    Mass = matrix(MassList).transpose()

    return geom, Mass

    # Example FORCE.MAT (water):
#Current cartesian coordinates              R   N=           9
#  1.38227456E-08  1.24339856E-01  1.78448645E-15 -1.45071721E+00 -9.86824393E-01
#  2.86852231E-08  1.45071697E+00 -9.86824393E-01 -2.86852870E-08
#Cartesian Force Constants                  R   N=          45
#  6.13473058E-01 -1.28387285E-07  4.70688939E-01 -1.21303083E-08 -1.82546103E-14
# -2.99988514E-08 -3.06736559E-01 -1.92514360E-01  6.06514705E-09  3.22985142E-01
# -2.34942228E-01 -2.35344484E-01  4.64554839E-09  2.13728428E-01  2.07648292E-01
#  6.06512751E-09  3.80659948E-09  5.67157634E-15 -6.38641140E-09 -4.22606394E-09
#  5.11430362E-15 -3.06736469E-01  1.92514524E-01  6.06514394E-09 -1.62485503E-02
#  2.12137811E-02  3.21284443E-10  3.22984964E-01  2.34942198E-01 -2.35344678E-01
# -4.64555550E-09 -2.12138072E-02  2.76963785E-02  4.19465712E-10 -2.13728413E-01
#  2.07648337E-01 -6.86992934E-08 -7.11914581E-08  2.99988585E-08  1.00143808E-07
#  7.79862503E-08 -5.05238509E-15 -8.85736817E-09  1.20357928E-08 -2.62903779E-08
#Atomic Number List                         R   N=           3
#    8    1    1
#Internal Coordinate List                   R   N=           0
#    1    2    0    0    1.0000
#    1    3    0    0    1.0000
#    2    1    3    0    1.0000
def readMM4Fc(file):
    #size the matrix using the number of coordinates:
    line = file.readline()
    threen = int(line.split()[5])
    Fc = matrix(array(zeros((threen,threen),dtype=float)))
    #read force constants
    line = file.readline()
    while (not line.startswith('Cartesian Force Constants')):#get to the beginning of the force constants section
	line = file.readline()
    #n = int(line.split()[5])
    #read the force constants section
    atomID1=0#row number - 1
    atomID2=0#column number - 1
    line = file.readline()
    while (not line.startswith('Atomic Number List')):
	split = line.split()#split line should have 5 elements
	for i in range(0,len(split)):
	    Fc[atomID1,atomID2]=float(split[i])
	    if(atomID2==atomID1):#reset atomID2 and increment atomID1 once we get to the end of the values for the row
		atomID2 = 0
		atomID1 = atomID1 + 1
	    else:#otherwise, increment atomID2 by 1
		atomID2 = atomID2 + 1
	line = file.readline()#read the next line

    #now Fc should be lower triangular (+ diagonal)
    return  Fc


def readGeomFc(file):
    lines = file.readlines()
    lines.reverse()
#read geometries

    i = lines.index("                          Input orientation:                          \n")
    geomLines = []
    stillRead = True
    k = i-5;
    while (stillRead):
       geomLines.append(lines[k])
       k = k-1
       if (lines[k].startswith(" ------------")):
          stillRead = False

    geom = matrix(array(zeros((len(geomLines),3),dtype=float)))
    Mass = matrix(array(zeros((len(geomLines),1),dtype=float)))

    for j in range(len(geomLines)):
        line = geomLines[j]
        tokens = line.split()
        geom[j,0] = double(tokens[3])
        geom[j,1]=double(tokens[4])
        geom[j,2]=double(tokens[5])
        if (int(tokens[1])== 6): Mass[j]=12.0
        if (int(tokens[1])== 8): Mass[j]=15.99491
        if (int(tokens[1])== 1): Mass[j]=1.00783
        if (int(tokens[1])== 7): Mass[j]=14.0031
        if (int(tokens[1])== 17): Mass[j]=34.96885
        if (int(tokens[1])== 16): Mass[j]=31.97207
        if (int(tokens[1])== 9): Mass[j]=18.99840

#read force constants
    Fc = matrix(array(zeros((len(geomLines)*3,len(geomLines)*3),dtype=float)))
    i = lines.index(" Force constants in Cartesian coordinates: \n")
    fclines = []
    stillRead = True
    k = i-1;
    while (stillRead):
       fclines.append(lines[k])
       k = k-1
       if (lines[k].startswith(" Force constants in internal coordinates:")):
          stillRead = False

    numRepeats = len(geomLines)*3/5 + 1
    j = 0
    while j in range(numRepeats):
       i = 5*j
       while i in range(5*j,len(geomLines)*3):
          line = fclines[(j*(len(geomLines)*3 + 1) - j*(j-1)*5/2)+i+1-5*j]
          tokens = line.split()
          k = 0
          while k in range(0,min(i+1-5*j,5)):
              Fc[i,5*j+k]=float(tokens[k+1].replace('D','E'))
              k = k+1
          i= i+1
       j = j+1

    return geom, Mass, Fc

def readFc(file):
    lines = file.readlines()
    lines.reverse()
#read geometries

    i = lines.index("                          Input orientation:                          \n")
    geomLines = []
    stillRead = True
    k = i-5;
    while (stillRead):
       geomLines.append(lines[k])
       k = k-1
       if (lines[k].startswith(" ------------")):
          stillRead = False

    geom = matrix(array(zeros((len(geomLines),3),dtype=float)))
    Mass = matrix(array(zeros((len(geomLines),1),dtype=float)))

    for j in range(len(geomLines)):
        line = geomLines[j]
        tokens = line.split()
        geom[j,0]=double(tokens[3])
        geom[j,1]=double(tokens[4])
        geom[j,2]=double(tokens[5])
        Mass[j] = getMassByAtomicNumber(int(tokens[1]))

#read force constants
    Fc = matrix(array(zeros((len(geomLines)*3,len(geomLines)*3),dtype=float)))
    i = lines.index(" Force constants in Cartesian coordinates: \n")
    fclines = []
    stillRead = True
    k = i-1;
    while (stillRead):
       fclines.append(lines[k])
       k = k-1
       if (lines[k].startswith(" Force constants in internal coordinates:")):
          stillRead = False

    numRepeats = len(geomLines)*3/5 + 1
    j = 0
    while j in range(numRepeats):
       i = 5*j
       while i in range(5*j,len(geomLines)*3):
          line = fclines[(j*(len(geomLines)*3 + 1) - j*(j-1)*5/2)+i+1-5*j]
          tokens = line.split()
          k = 0
          while k in range(0,min(i+1-5*j,5)):
              Fc[i,5*j+k]=float(tokens[k+1].replace('D','E'))
              k = k+1
          i= i+1
       j = j+1

    return  Fc

def printNormalModes(l,v,num,Mass):
   for i in range(num):
    for k in range(3):
      print ('%18.3f '%(sqrt(l[i*3+k])*337.0/6.5463e-02)),  
    print

    for j in range(Mass.size):
      for m in range(3):
         for k in range(3):
           print ('%6.2f'%(v[j*3+k,i*3+m])),
         print ' ',
      print
    print  

def readEnergy(file, string):
   com = file.read()
   if string == 'cbsqb3':
     Energy=re.search('CBS-QB3 \(0 K\)= '+' \s*([\-0-9.]+)',com).group(1)
   if string == 'g3':
     Energy=re.search('G3\(0 K\)= '+' \s*([\-0-9.]+)',com).group(1)
   if string == 'klip_1':
     Energy=re.search('QCI= '+' \s*([\-0-9.]+)',com).group(1)
   if string == 'klip_2':
      Energy=re.search('QCI= '+' \s*([\-0-9.]+)',com).group(1)
   if string == 'klip_2_cc':
      Energy=re.search('CC= '+' \s*([\-0-9.]+)',com).group(1)
 
   return float(Energy)

def readHFEnergy(fileName):
 file = open(fileName,'r')
 Result = ""

 line = file.readline()

 startReading = False
 while line!="":
                
                if (line[:9] == r" 1\1\GINC"):
                    startReading = True
                if (line[-4:-1]== r"\\@"):
                    Result = Result+line[1:-1]
                    break
                if (startReading):
                    Result= Result+line[1:-1]
                    #pdb.set_trace()
                line = file.readline()

 hf_5 = getNum(Result,r"\HF=")
 file.close() 
 return hf_5

def getNum(Result,id):
    if (len(Result) <= 5):
        return 0
    fields = Result.split(id)
    resStr = fields[1]
    numstr = ""
    for i in range(len(resStr)):
        if resStr[i] == "\\" :
            break
        numstr = numstr + resStr[i]
    return float(numstr)

def getAtoms(Mass):
    atoms = len(Mass)*['']
    j = 0
    for m in Mass:
        if (int(m)== 12): atoms[j]='C'
        if (int(m)== 15): atoms[j]='O'
        if (int(m)== 1): atoms[j]='H'
        if (int(m)== 14): atoms[j]='N'
        if (int(m)== 34): atoms[j]='Cl'
        if (int(m)== 31): atoms[j]='S'
        if (int(m)== 18): atoms[j]='F'
        j = j+1
    return atoms

def getMassByAtomicNumber(atomicNum) :
    if (atomicNum == 6): Mass=12.0
    if (atomicNum == 8): Mass=15.99491
    if (atomicNum == 1): Mass=1.00783
    if (atomicNum == 7): Mass=14.0031
    if (atomicNum == 17): Mass=34.96885
    if (atomicNum == 16): Mass=31.97207
    if (atomicNum == 9): Mass=18.99840
    return Mass

def getMassByAtomicSymbol(atomicSym) :
    if (atomicSym == 'C'): Mass=12.0
    if (atomicSym == 'O'): Mass=15.99491
    if (atomicSym == 'H'): Mass=1.00783
    if (atomicSym == 'N'): Mass=14.0031
    if (atomicSym == 'Cl'): Mass=34.96885
    if (atomicSym == 'S'): Mass=31.97207
    if (atomicSym == 'F'): Mass=18.99840
    return Mass
