
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

#Computes the Eckart tunneling correction
#   Inputs are:
#       delV1 - energy difference between TS and reactants [=] Joules
#       T - absolute temperature [=] Kelvin
#       alpha1 - dimensionless energy difference between TS and reactants
#       alpha2 - dimensionless energy difference between TS and reactants (if symmetric)
#                  or between TS and products (if asymmetric)
#   Output is kappa(T), the dimensionless tunneling correction factor

import os 
from numpy import *
from scipy import *
import CanTherm

def computeTunnelingCorrection(delV1,T,alpha1,alpha2):
    #print delV1, ' ', alpha1, ' ', alpha2, ' ', T
    k = 1.381e-23
# The following lines of code were written by MRH for debugging purposes
# The file table1.out contains a table of kappa(T) for a range of dimensionless
#    variables alpha1, alpha2, and ustar and was compared against Table1 of:
#    "Tunneling corrections for unsymmetrical eckart potential energy barriers", 
#    H.S. Johnston and J. Heicklen, J. Phys. Chem., v,66, (1962), p.532-533
# MRH believes he found an error in the data reported in the paper:
#    alpha1=0.5; alpha2=0.5,1,2,4; ustar=all
#    alpha1=1; alpha2=1,2; ustar=all
# Problem is that division in calculating kappa_E should be done element-by-element
#    MRH believes this was not done in Johnston paper, rather a least-squares fitting
#    was performed to solve for "a" where "aX=Y" by calling "Y/X"
#
#    oFile = open('table1.out','w')
#    alpha1vec = [0.5, 1, 2, 4, 8, 12, 16, 20]
#    alpha2vec = [0.5, 1, 2, 4, 8, 12, 16, 20]
#    ustarvec = [2, 3, 4, 5, 6, 8, 10, 12, 16]
#    for i in range(len(alpha1vec)):
#        alpha1 = alpha1vec[i]
#        for j in range(i,len(alpha2vec)):
#            alpha2 = alpha2vec[j]
#            for m in range(len(ustarvec)):
#                T = 2*math.pi*delV1/alpha1/k/ustarvec[m]
#                integral = integrate.quad(f_of_E,0,50,args=(delV1,k*T,alpha1,alpha2))[0]
#                kappa_T = integral*math.exp(delV1/k/T)
#                oFile.write(str(kappa_T) + '\t')
#            oFile.write('\n')
#        oFile.write('\n')
#    oFile.close()
# MRH changed limits of integration from 0 --> infinity to 0 --> 50
#    Problem: Some calculations were underestimating the integral, most likely due to large
#       step sizes when the important portion of the curve is close to zero.
#    MRH computed the value of E_kT* such that f_of_E(E_kT*) < 0.001*max(f_of_E)
#       for all alpha1, alpha2, ustar combinations present in Johnston paper.  This value was
#       ~35, so MRH used decided on 50 as the upper bound of the integral
    #integral = integrate.quad(f_of_E,0,100,args=(delV1,k*T,alpha1,alpha2))[0]
    #kappa_T = integral * math.exp(delV1/k/T)
    
    x = array(zeros((10000,1),dtype=float))
    f_x = array(zeros((10000,1),dtype=float))
    for i in range(10000):
        x[i] = 1000.0*i/9999
        f_x[i] = f_of_E(x[i],delV1,k*T,alpha1,alpha2)
    max = f_x.max(0)
    lowerlimit = max/1000.0
    vector_of_Es = (f_x>lowerlimit).nonzero()
    maxE = x[vector_of_Es[0][-1]]
    minE = x[vector_of_Es[0][0]]
    #print str(minE) + ' ' + str(maxE)
    integral = integrate.quad(f_of_E,minE,maxE,args=(delV1,k*T,alpha1,alpha2))[0]
    kappa_T = integral * math.exp(delV1/k/T)
    print kappa_T
    return kappa_T

def f_of_E(E_kt,delV1,kT,alpha1,alpha2):
    radicand = alpha1*alpha2-4*math.pi*math.pi/16
    if radicand < 0 :
        twopid = 2*math.sqrt(-1*radicand)
    else :
        twopid = 2*math.sqrt(radicand)
    nondimE = E_kt*kT/delV1
    twopia = 2*math.sqrt(alpha1*nondimE)/(1/math.sqrt(alpha1)+1/math.sqrt(alpha2))
    radicand2 = (nondimE-1)*alpha1+alpha2
    if radicand2 < 0:
        twopib = 2*math.sqrt(-1*radicand2)
    else:
        twopib = 2*math.sqrt(radicand2)
    twopib = twopib/(1/math.sqrt(alpha1)+1/math.sqrt(alpha2))
# python cannot handle computing the value of cosh(700)
#    To be safe, MRH checks if any of the cosh arguments are greater than 200
# If all cosh() arguments less than 200, compute kappa_E as normal
    if (twopia < 200) & (twopib < 200) & (twopid < 200) :
        kappa_E = 1 - (math.cosh(twopia-twopib)+math.cosh(twopid)) / (math.cosh(twopia+twopib)+math.cosh(twopid))
# If not, need to be smarter about which terms to evaluate
    else :
# If at least one of the following expressions is greater than 5, we can eliminate most of the exponential terms
#    after writing out the definition of cosh() and dividing all terms by exp(twopid)
        if (twopia-twopib-twopid > 10) | (twopib-twopia-twopid > 10) | (twopia+twopib-twopid > 10) :
            kappa_E = 1 - exp(-2*twopia) - exp(-2*twopib) - exp(-twopia-twopib+twopid) - exp(-twopia-twopib-twopid)
# If all of the arguments are less than 5, then evaluate the kappa_E expression normally, except use the expanded
#    definition - expanding the cosh argument and dividing all terms by exp(twopid)
        else :
            numerator = math.exp(twopia-twopib-twopid)+exp(-twopia+twopib-twopid)+1+exp(-2*twopid)
            denominator = math.exp(twopia+twopib-twopid)+exp(-twopia-twopib-twopid)+1+exp(-2*twopid)
            kappa_E = 1 - numerator/denominator
    integrand = math.exp(-E_kt)*kappa_E
    return integrand
