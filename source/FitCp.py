from numpy import *
import pdb
import readGeomFc
import math
import Gnuplot,Gnuplot.funcutils
from scipy import *
from scipy.optimize import fsolve
from scipy.optimize import fminbound

def FitHeatCapacity(Temp, Cp, linearity, Freq, numRotors, R,oFile):

       B = fitWilhoitB(Temp,Cp, linearity, Freq, numRotors, R, 500)
       Cp_0, Cp_inf, B,a,b,c,d = getWilhoit(Temp,Cp, linearity, Freq, numRotors, R, B)
      
       oFile.write('Wilhoit Parameters: \n')
       oFile.write('Cp_0 = ' + str(Cp_0) + '\n')
       oFile.write('Cp_inf = ' + str(Cp_inf) + '\n')
       oFile.write('B = ' + str(B) +'\n')
       oFile.write('a= ' + str(a) + '\n')
       oFile.write('b= ' + str(b) + '\n')
       oFile.write('c= ' + str(c) + '\n')
       oFile.write('d= ' + str(d) + '\n')

       T_inter = 1000.
       T_inter = fitNASAT(T_inter, Temp, Cp, R)

       z1 = getNASA(T_inter, Temp, Cp, R)
       z2 = getShomate(Temp, Cp, R)

       plot_results(Temp,Cp, R, Cp_0, Cp_inf, B, a, b, c, d, T_inter, z1, z2)       
       
#-------------- Wilhoit Section
#--------------- given a Wilhoit B, find the other parameters
def getWilhoit(Temp, Cp, linearity, Freq, numRotors, R, B):
     
    if linearity == 'Linear':
        Cp_0 = R*float(1.5 + 1.0 + 1.0);
    elif linearity == 'Nonlinear':
        Cp_0 = R*float(1.5 + 1.5 + 1.0);
    else:    
        print 'linearity not found'

    Cp_inf = R*float(len(Freq) + numRotors);
    
    Y = matrix(zeros((len(Cp),1),dtype=float))
    X = matrix(zeros((len(Cp),4),dtype=float))
    reduced_cp = []
    y = []
#    B = float('500')
    for i in range(len(Temp)):
        reduced_cp.append((float(Cp[i]) - Cp_0)/(Cp_inf - Cp_0))
        y.append(Temp[i]/(Temp[i]+B))

        Y[i,0] = (reduced_cp[i] - y[i]**2)
        for j in range(4):
            X[i,j] = (y[i]**3 - y[i]**2) * y[i]**j
    
    XtX = transpose(X)*X
    XtY = transpose(X)*Y
    z = linalg.inv(XtX)*XtY
    
    B = B
    a = float(z[0])
    b = float(z[1])
    c = float(z[2])
    d = float(z[3])

    return Cp_0, Cp_inf, B, a, b, c, d

#------------- find the Wilhoit B
def fitWilhoitB(Temp, Cp, linearity, Freq, numRotors, R, B0):
    if linearity == 'Linear':
        Cp_0 = R*float(1.5 + 1.0 + 1.0);
    elif linearity == 'Nonlinear':
        Cp_0 = R*float(1.5 + 1.5 + 1.0);
    else:    
        print 'linearity not found'
           
    Cp_inf = R*float(len(Freq) + numRotors)
    
    def fitme(B, Temp, Cp, Cp_0, Cp_inf, R):
            
        Y = matrix(zeros((len(Cp),1),dtype=float))
        X = matrix(zeros((len(Cp),4),dtype=float))
        
        reduced_cp = []
        y = []
        for i in range(len(Temp)):
            reduced_cp.append((float(Cp[i]) - Cp_0)/(Cp_inf - Cp_0))
            y.append(Temp[i]/(Temp[i]+B))

            Y[i,0] = (reduced_cp[i] - y[i]**2)
            for j in range(4):
                X[i,j] = (y[i]**3 - y[i]**2) * y[i]**j
    
        XtX = transpose(X)*X
        XtY = transpose(X)*Y
        z = linalg.inv(XtX)*XtY
    
        B = B
        a = float(z[0])
        b = float(z[1])
        c = float(z[2])
        d = float(z[3])
        intermediate = []
        Cp_Wilhoit = []
        minimize_me = []
        for i in range(len(Temp)):
#        for i in range(30):       
               intermediate.append( (y[i]-1.0) * ( a + b * y[i] + c * y[i]**2 + d * y[i]**3) )
               Cp_Wilhoit.append( Cp_0 + (Cp_inf - Cp_0) * y[i]**2 * (1 + intermediate[i]))
               minimize_me.append( (Cp[i] - Cp_Wilhoit[i]) )
        minimize_me = sum( minimize_me[:])**2 
        return minimize_me       
                
 #   result = fsolve(fitme, B0, args=(Temp, Cp, Cp_0, Cp_inf, R))
    result = fminbound(fitme, 300.0, 3000.0, args=(Temp, Cp, Cp_0, Cp_inf, R))
    print 'Wilhoit T = ', result

    return result


#---------- End Wilhoit Section --------------------

#-----------Begin NASA Section

def getNASA(T_inter, Temp, Cp, R):

    Y = matrix(zeros((len(Temp)+3,1),dtype=float))
    X = matrix(zeros((len(Temp)+3,10),dtype=float))
    
    for i in range(len(Temp)):
        Y[i,0] = float(Cp[i]/R)
        
        if Temp[i] <= T_inter:
            for j in range(5):
                X[i,j] = float(Temp[i])**j
        elif Temp[i] >= T_inter:
            for j in range(5):
                X[i,j+5] = float(Temp[i])**j
        else:
            print 'WTF?'  
 
         # continuity
    for j in range(5):
        X[len(X)-3,j] = float(T_inter)**j        
        X[len(X)-3,j+5] = -float(T_inter)**j     

        # first derivative is zero
    for j in range(1,5):
        X[len(X)-2,j] = j*float(T_inter)**(j-1)        
        X[len(X)-2,j+5] = -j*float(T_inter)**(j-1)
 
        # second derivative is zero
    for j in range(2,5):
        X[len(X)-1,j] = j*(j-1)*float(T_inter)**(j-2)        
        X[len(X)-1,j+5] = -j*(j-1)*float(T_inter)**(j-2)

    XtX = transpose(X)*X
    XtY = transpose(X)*Y
    z = linalg.inv(XtX)*XtY

    return z


#------------- find the NASA T_inter
def fitNASAT(T_inter0, Temp, Cp, R):
    """Find the optimal intermediate temperature for the NASA polynomials."""
    def fitme(T_inter, Temp, Cp, R):
        z=getNASA(T_inter, Temp, Cp, R)
        Cp_NASA = []
        minimize_me = []
        for i in range(len(Temp)):
            T = Temp[i]
            T2 = T*T
            if T <= float(T_inter):
               Cp_NASA.append( R*(z[0] + z[1] * T +  z[2] * T2 + z[3] * T*T2 +  z[4] * T2*T2 ) )
               minimize_me.append( (Cp[i] - Cp_NASA[i]) )
            elif T > float(T_inter):
               Cp_NASA.append( R*(z[5] + z[6] * T +  z[7] * T2 + z[8] * T*T2 +  z[9] * T2*T2 ) )
               minimize_me.append( (Cp[i] - Cp_NASA[i]) ) 
            else:
                print 'WTF?' 
        minimize_me = sum( minimize_me[:])**2 
            
        return minimize_me       
                
    result = fminbound(fitme, 300., 3000., args=(Temp, Cp, R) ) 
    print "NASA T= ", result
    return result

#-----------Begin SHOMATE Section

def getShomate(Temp, Cp, R):

    Y = matrix(zeros((len(Temp),1),dtype=float))
    X = matrix(zeros((len(Temp),5),dtype=float))
    
    for i in range(len(Temp)):
        Y[i,0] = float(Cp[i]/R)
        for j in range(4):
            X[i,j] = float(Temp[i])**j
        X[i,4] = 1.0 / float(Temp[i])**2    
          
    XtX = transpose(X)*X
    XtY = transpose(X)*Y
    z = linalg.inv(XtX)*XtY

    return z





#########################################################
def plot_results(Temp, Cp, R, Cp_0, Cp_inf, B, a, b, c, d, T_inter, z1, z2):
       
    intermediate = []
    Cp_Wilhoit = []
    y = []
    cpgiven = []
    cpwillie = []
    willie_error = []
    for i in range(len(Temp)):
        y.append( Temp[i]/(Temp[i] + B) )
        intermediate.append( (y[i]-1.0) * ( a + b * y[i] + c * y[i]**2 + d * y[i]**3) )
        Cp_Wilhoit.append( Cp_0 + (Cp_inf - Cp_0) * y[i]**2 * (1 + intermediate[i]))
        cpgiven.append([Temp[i], Cp[i]])
        cpwillie.append([Temp[i], Cp_Wilhoit[i]])
        willie_error.append([Temp[i], (Cp[i] - Cp_Wilhoit[i])/Cp[i]*100 ])

    Cp_NASA = []
    cpnasa = []
    nasa_error = []
    for i in range(len(Temp)):
        if Temp[i] <= float(T_inter):
           Cp_NASA.append( R*(z1[0] + z1[1] * Temp[i] +  z1[2] * Temp[i]**2 + z1[3] * Temp[i]**3 +  z1[4] * Temp[i]**4 ) )
           cpnasa.append([Temp[i], Cp_NASA[i]])
           nasa_error.append([Temp[i], (Cp[i] - Cp_NASA[i])/Cp[i]*100 ])
        elif Temp[i] > float(T_inter):
           Cp_NASA.append( R*(z1[5] + z1[6] * Temp[i] +  z1[7] * Temp[i]**2 + z1[8] * Temp[i]**3 +  z1[9] * Temp[i]**4 ) )
           cpnasa.append([Temp[i], Cp_NASA[i]])
           nasa_error.append([Temp[i], (Cp[i] - Cp_NASA[i])/Cp[i]*100 ])
        else:
           print 'WTF?' 
      
    Cp_Shomate = []
    cpshom = []
    shom_error = []
    for i in range(len(Temp)):
        Cp_Shomate.append( R*(z2[0] + z2[1] * Temp[i] + z2[2] * Temp[i]**2 + z2[3] * Temp[i]**3 + z2[4]/Temp[i]**2 ) )
        cpshom.append([ Temp[i], Cp_Shomate[i] ])
        shom_error.append([Temp[i], (Cp[i] - Cp_Shomate[i])/Cp[i]*100 ])
   
    g=Gnuplot.Gnuplot(persist = 1)
    g('set multiplot')
    g('set size 1,.5')
    plot1 = Gnuplot.PlotItems.Data(cpgiven, with_="points 3", title=None,)
    plot2 = Gnuplot.PlotItems.Data(cpwillie, with_="lines", title="Wilhoit" )
    plot3 = Gnuplot.PlotItems.Data(cpnasa, with_="lines", title="NASA" )
    plot4 = Gnuplot.PlotItems.Data(cpshom, with_="lines", title="Shomate" )
    g('set origin 0.0,0.5')
    g.plot(plot1,plot2, plot3, plot4)
    g('set origin 0.0, 0.0')
    plot12 = Gnuplot.PlotItems.Data(willie_error, with_="linespoints", title="Wilhoit")
    plot13 = Gnuplot.PlotItems.Data(nasa_error, with_="linespoints", title="NASA")
    plot14 = Gnuplot.PlotItems.Data(shom_error, with_="lines", title="Shomate")
    g.plot(plot12, plot13)
#    g.plot(plot12, plot13, plot14)
    g('unset multiplot')
 
