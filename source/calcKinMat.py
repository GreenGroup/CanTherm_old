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

import readGeomFc
import sys
import geomUtility
from numpy import *
from scipy import *
from scipy.linalg import *
import pdb

def main():
    logFile = open(sys.argv[1],'r')
    file = logFile
    inertiaFile = open('inertia.dat','r')    
    (geom,Mass)=readGeomFc.readGeom(logFile)
    (rotors)= readGeomFc.readGeneralInertia(inertiaFile)
    #D=calculateI23(geom,Mass,rotors)
    D=calculateD(geom,Mass,rotors)
    (l,v)= linalg.eigh(D)
    print l, linalg.det(D)


def calculateD(geom,Mass,rotors):
    #D=geomUtility.calculateD(geom,Mass,rotors)
    #print linalg.det(D)
    #exit()
    numRotors = len(rotors)
#    print geom
#    exit()
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
    (l,v)=linalg.eig(I)
    prinAxes = transpose(v);

    for i in range(Mass.size):
        geom[i,:]=transpose(prinAxes*transpose(geom[i,:]))

#again calculate the moments of inertia and confirm that they are diagonal
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
#    print I
    
    K=matrix(zeros((6+numRotors-1,6+numRotors-1),dtype=double))
    
    M = eye(3,3)
    M = matrix(M*sum(Mass))
    U1 = matrix('0.0; 0.0; 0.0')
        
    for i in range(3):
        U1[i]= sum(array(Mass)*array(geom[:,i]))


    U=matrix(zeros((3,3),dtype=double))
    U[0,1]=-1*U1[2]
    U[1,0]=U1[2]
    U[0,2]=U1[1]
    U[2,0]=-1*U1[1]
    U[1,2]=-1*U1[0]
    U[2,1]=U1[0]
    
    
    K[0:3,0:3]=M
    K[3:6,0:3]=U
    K[0:3,3:6]=transpose(U)
    K[3:6,3:6]=I
    
    
    irotor = 0
    for rotor in rotors[1:]:
        rotor.getAxes(geom,Mass)
        #pdb.set_trace()  
        rotor.getMoments(geom,Mass)

        A = rotor.moments[0]
        B = rotor.moments[1]
        C = rotor.moments[2]
        Ux = rotor.moments[3]
        dircos = rotor.dircos
        r = transpose(rotor.r)
        level = rotor.level
        
    #    print dircos[2,:]
        K[6+irotor,0:3] = Ux*dircos[1,:]
        K[0:3,6+irotor] = transpose(Ux*dircos[1,:])
        
        beta_i1 = mat('0.0;0.0;0.0')
        for i in range(3):
            beta_i1[i] = dircos[2,i]*A - dircos[0,i]*B - dircos[1,i]*C 

        beta_i1[0] = beta_i1[0]+ (dircos[1,2]*r[1] - dircos[1,1]*r[2])*Ux
        beta_i1[1] = beta_i1[1]+ (dircos[1,0]*r[2] - dircos[1,2]*r[0])*Ux
        beta_i1[2] = beta_i1[2]+ (dircos[1,1]*r[0] - dircos[1,0]*r[1])*Ux
        
        K[6+irotor,3:6] = transpose(beta_i1)
        K[3:6,6+irotor] = beta_i1
        
        K[6+irotor,6+irotor] = A
        
        
        
        if (level >= 3):
            pivot2 = rotor.pivot2
            numAncestors = level-2
            
            presentParent = rotor
            #print numAncestors
            for i in range(numAncestors):
                parentNum = presentParent.parent
                
                presentParent = rotors[parentNum]
                dircos_ipar = dircos*transpose(presentParent.dircos)
                #pdb.set_trace()
                #ripar = presentParent.dircos*(transpose(geom[rotor.pivotAtom-1,:]-geom[presentParent.pivotAtom-1,:]))
                ripar = transpose((geom[rotor.pivotAtom-1,:]-geom[presentParent.pivotAtom-1,:])*presentParent.dircos)
                beta_z = dircos_ipar[2,2]*A - dircos_ipar[0,2]*B - dircos_ipar[1,2]*C 
                beta_z = beta_z + (dircos_ipar[1,1]*ripar[0] - dircos[1,0]*ripar[1])*Ux
                K[6+irotor,6+parentNum-1]=K[6+parentNum-1,6+irotor] = beta_z
                
        irotor = irotor+1


        #print S
        #print K
    S = K[3:,3:] - K[3:,0:3]*linalg.inv(K[0:3,0:3])*K[0:3,3:]

    D = S[3:,3:] - S[3:,0:3]*linalg.inv(S[0:3,0:3])*S[0:3,3:]
    
    '''for i in range(numRotors-1):
      for j in range(numRotors-1):
        print '%8.3f'%D[i,j],
      print
    print

    for i in range(6+numRotors-1):
      for j in range(6+numRotors-1):
        print '%8.3f'%K[i,j],
      print
    print

    for i in range(3+numRotors-1):
      for j in range(3+numRotors-1):
        print '%8.3f'%S[i,j],
      print
    print
    #print linalg.det(D)
    #exit()'''
    
    return D



def calculateI23(geom,Mass,rotors):
    numRotors = len(rotors)
    
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
    (l,v)=linalg.eig(I)
    prinAxes = transpose(v);

    for i in range(Mass.size):
        geom[i,:]=transpose(prinAxes*transpose(geom[i,:]))

#again calculate the moments of inertia and confirm that they are diagonal
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
#    print I
    
    K=matrix(zeros((6+numRotors-1,6+numRotors-1),dtype=double))
    
    M = eye(3,3)
    M = matrix(M*sum(Mass))
    U1 = matrix('0.0; 0.0; 0.0')
        
    for i in range(3):
        U1[i]= sum(array(Mass)*array(geom[:,i]))


    U=matrix(zeros((3,3),dtype=double))
    U[0,1]=-1*U1[2]
    U[1,0]=U1[2]
    U[0,2]=U1[1]
    U[2,0]=-1*U1[1]
    U[1,2]=-1*U1[0]
    U[2,1]=U1[0]
    
    
    
    K[0:3,0:3]=M
    K[3:6,0:3]=U
    K[0:3,3:6]=transpose(U)
    K[3:6,3:6]=I
    
    
    irotor = 0
    beta_i = []
    redMom = matrix(zeros((numRotors-1,numRotors-1),dtype=double))
    Amm = matrix(zeros((numRotors-1,numRotors-1),dtype=double))
    Im0 = mat(zeros((numRotors-1,1),dtype=double))

    for rotor in rotors[1:]:
        rotor.getAxes(geom,Mass)
        rotor.getMoments(geom,Mass)
    
        A = rotor.moments[0]
        B = rotor.moments[1]
        C = rotor.moments[2]
        Ux = rotor.moments[3]
        dircos = rotor.dircos
        r = transpose(rotor.r)
                
        beta_i1 = mat('0.0;0.0;0.0')
        for i in range(3):
            beta_i1[i] = dircos[2,i]*A - dircos[0,i]*B - dircos[1,i]*C 

        beta_i1[0] = beta_i1[0]+ (dircos[1,2]*r[1] - dircos[1,1]*r[2])*Ux
        beta_i1[1] = beta_i1[1]+ (dircos[1,0]*r[2] - dircos[1,2]*r[0])*Ux
        beta_i1[2] = beta_i1[2]+ (dircos[1,1]*r[0] - dircos[1,0]*r[1])*Ux
        
        beta_i.append(beta_i1)


        Im0[irotor] = A
        for i in range(3):
           Im0[irotor] = Im0[irotor] - (dircos[1,i]*Ux)**2/sum(Mass) - beta_i1[i]**2/I[i,i]
        irotor = irotor +1

    for k in range(len(rotors)-1):
        for l in range(len(rotors)-1):
           rotk = rotors[k+1]
           rotl = rotors[l+1]
           for i in range(3):
               Amm[k,l] = Amm[k,l] + rotk.dircos[1,i]*rotl.dircos[1,i]*rotk.moments[3]*rotl.moments[3]/sum(Mass) \
                            + beta_i[k][i]*beta_i[l][i]/I[i,i]
           
    for i in range(len(rotors)-1):
       redMom[i,i] = Im0[i]
       #for j in range(len(rotors)-1):
       #   if j!=i :
             #redMom[i,i] = redMom[i,i] - (1.0/2.0)*Amm[i,j]**2/Im0[j]
    
    
    #S = K[3:,3:] - K[3:,0:3]*linalg.inv(K[0:3,0:3])*K[0:3,3:]

    #D = S[3:,3:] - S[3:,0:3]*linalg.inv(S[0:3,0:3])*S[0:3,3:]
    
    return redMom
if __name__ == "__main__":
   main()

