
# Rocket Project
# Noah Alex and Robert Nievar
# December 14th, 2019
# Dynamics; Dr. Bruce Gale

### PROBLEMS 6, 7, AND 8 ###

#########  Libraries  ##########
import math
import time
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import axes
from matplotlib import figure
import pandas as pd
from scipy.integrate import quad

##########  Functions  ###########
def thrustcurve():
    ## DATA FROM THE PLOT GIVEN IN THE PDF ###
    times = [
    0,0.009243735,0.087357361,0.471362851,0.969394183,1.48137112,1.986603125,
    2.475519681,2.926609918,2.995517621,3.147096337,3.2153
    ]
    times = np.asarray(times) #CHANGING TO NUMPY ARRAY
    thrusts = [
    0,1311.883715,1544.447211,1496.457918,1555.521663,1640.425797,
    1555.521663,1393.096365,1286.043327,1153.1499,182.2895916,7.89
    ]
    thrusts = np.asarray(thrusts) #CHANGING TO NUMPY ARRAY

    def estimate_coef(x, y): 
        # number of observations/points 
        n = np.size(x) 
        # mean of x and y vector 
        m_x, m_y = np.mean(x), np.mean(y) 
        # calculating cross-deviation and deviation about x 
        SS_xy = np.sum(y*x) - n*m_y*m_x 
        SS_xx = np.sum(x*x) - n*m_x*m_x 
        # calculating regression coefficients 
        b_1 = SS_xy / SS_xx 
        b_0 = m_y - b_1*m_x 
        return(b_0, b_1) 

    def integrand(var, b_0, b_1):
        return (b_1*var + b_0)
    thrustarray = []
    
    for i in range(len(thrusts)):
        try:
            y = [thrusts[i], thrusts[i+1]]
            y = np.asarray(y)

            x = [times[i], times[i+1]]
            x = np.asarray(x)
            b_0, b_1 = estimate_coef(x, y)
            ans, err = quad(integrand, x[0], x[1], args = (b_0,b_1))
            thrustarray.append(ans)
        except:
            break
    return thrustarray, times

def problem6_3(thrustarray, times):

    def Yrocket(thrustarray, times):
        g = 32.2 # gravitational acceleration
        Fg = 50.1
        mass = Fg/g
        massloss = (5.7/32.2) / (11)
        Cd = .3
        A = np.pi*(6/2)**2 #area of rocket
        p = 0.000002379 #density
        launchtheta = (85.0 * np.pi) / 180.0
        sin = np.sin(launchtheta)
        cosine = np.cos(launchtheta)
        Rdrogue = (18/12.0)/2.0 #radius of drogue chute
        Vdrogue = 83.1 # terminal veloctiy with drogue chute
        Adrogue = (math.pi)*Rdrogue**2
        Vmain = 11.2
        Amain = math.pi * ((144/12.0)/2.0)**2
        Cddrogue = 2 * (Fg)/ ( p * Vdrogue**2 * Adrogue) 
        Cdmain= 2 * (Fg)/ ( p * Vmain**2 * Amain)
        Ft = 310
        Fnet, h0, a, hf, Fd, t, V0, V,Fddrogue, Fdmain, Fsum  = (False,)*11
        Vvector, FDvector, Avector, Tvector, Hvector = ([] for i in range(5))
        h = 0

        for i in range(len(thrustarray)):
            mass = mass - massloss
            Fd =  Cd*A*((p*(V**2)))*.5
            t = times[i]
            Fnet = (thrustarray[i]*sin) - (mass*32.2) - Fd
            a = Fnet/mass
            Avector.append(a)
            V = a*t
            Vvector.append(V)
            h += .5 * a * (t**2)
            Hvector.append(h)
            Tvector.append(t)
        while V > 0:
            Tvector.append(t)
            Fd =  Cd*A*((p*(V**2)))*.5
            Fnet = - Fg - Fd
            a = Fnet/mass
            Avector.append(a)
            V += a*.01
            Vvector.append(V)
            h += V*.01
            Hvector.append(h)
            t +=  .01
        while h > 600: #Drogue chute deployed
            Tvector.append(t)
            Fsum = Fd - Fg
            a = Fsum / mass
            h += V*.01
            Hvector.append(h)
            Avector.append(a)
            V += (a*.01)
            Vvector.append(V)

            Fd = .5 * p * (V**2) * Adrogue * Cddrogue
            FDvector.append(Fd)
            t += .01
        while h > 0: #main chute
            Tvector.append(t)
            V += (a*.01)
            Vvector.append(V)
            
            # drogue chute opens
            Fd = .5 * p * (V**2) * Amain * Cdmain
            FDvector.append(Fd)
            Fsum = Fd - Fg
            a = Fsum / mass
            Avector.append(a)
            h += V*.01
            Hvector.append(h)
            t += .01

        #generate subplots
        fig, (plt1, plt2, plt3) = plt.subplots(1,3, figsize = (12,6), dpi = 80)
        plt.autoscale()
        #rotate and align labels
        fig.autofmt_xdate()
        fig.suptitle('Problem 6: Revision of 3 (Y)', fontsize = 15)

        #plot the subplots
        plt1.plot(Tvector, Hvector, color='red',label = 'Height')
        plt1.set_title('Height')
        plt1.set(xlabel = 'Time (s)', ylabel = 'Altitude (ft)')
        plt2.plot(Tvector,Vvector, color='black',label = 'Velocity')
        plt2.set_title('Velocity')
        plt2.set(xlabel = 'Time (s)', ylabel = 'Velocity (ft/s)')
        plt3.plot(Tvector, Avector, color='olive', label='Acceleration')
        plt3.set_title('Acceleration')
        plt3.set(xlabel = 'Time (s)', ylabel = 'Acceleration(ft/s^2)')
        
        #plt.tight_layout()

        plt.show()
        return Tvector

    def Xrocket(Tyvector, thrustarray, times):
        g = 32.2 # gravitational acceleration
        Fg = 50.1
        mass = Fg/g
        massloss = (5.7/32.2) / (11)
        Cd = .3
        A = np.pi*(6/2)**2 #area of rocket
        p = 0.000002379 #density
        launchtheta = (85.0 * np.pi) / 180.0
        sin = np.sin(launchtheta)
        cosine = np.cos(launchtheta)
        Rdrogue = (18/12.0)/2.0 #radius of drogue chute
        Vdrogue = 83.1 # terminal veloctiy with drogue chute
        Adrogue = (math.pi)*Rdrogue**2
        Vmain = 11.2
        Amain = math.pi * ((144/12.0)/2.0)**2
        Cddrogue = 2 * (Fg)/ ( p * Vdrogue**2 * Adrogue) 
        Cdmain= 2 * (Fg)/ ( p * Vmain**2 * Amain)
        Ft = 310
        A = np.pi*(6/2)**2 #area of rocket
        p = 0.00002379 #density
        Fnet, h0, a, hf, Fd, t, V0, V, Fddrogue, Fdmain, Fsum, x  = (False,)*12

        Vvector, FDvector, Avector, Tvector, Xvector = ([] for i in range(5))

        for i in range(len(thrustarray)):
            t = times[i]
            mass = mass - massloss
            Tvector.append(t)
            Fd =  Cd*A*((p*(V**2)))*.5 
            Fx = (thrustarray[i] * cosine) - Fd
            a = Fx / mass
            Avector.append(a)
            V = a*t
            Vvector.append(V)
            x += V*t
            Xvector.append(x)

        while t < Tyvector[-1]:

            Fd =  Cd*A*((p*(V**2)))*.5 
            Tvector.append(t)
            Fx = Fd
            a = Fx / mass
            Avector.append(a)
            V = a*t
            Vvector.append(V)
            x += V*t
            Xvector.append(x)
            t += .01

        #generate subplots
        fig, (plt1, plt2, plt3) = plt.subplots(1,3, figsize = (12,6), dpi = 80)
        plt.autoscale()

        #rotate and align labels
        fig.autofmt_xdate()
        fig.suptitle('Problem 6: Revision of 3 (X)', fontsize = 15)

        #plot the subplots
        plt1.plot(Xvector,Tvector, color='red',label='X travel')
        plt1.set_title('X travel')
        plt1.set(xlabel = 'X-Travel (Ft)', ylabel = 'Time (s)')
        plt2.plot(Tvector,Vvector, color='black',label = 'Velocity')
        plt2.set_title('Velocity')
        plt2.set(xlabel = 'Time (s)', ylabel = 'Velocity (ft/s)')
        plt3.plot(Tvector, Avector, color='olive', label='Acceleration')
        plt3.set_title('Acceleration')
        plt3.set(xlabel = 'Time (s)', ylabel = 'Acceleration (ft/s^2)')
        
        #plt.tight_layout()

        plt.show()  

    Tyvector = Yrocket(thrustarray, times)
    Xrocket(Tyvector, thrustarray, times)

def problem6_4(thrustarray,times):
    def Yrocket(thrustarray, times):
        # variable/array initialization
        g = 32.2 # gravitational acceleration
        Fg = 50.1
        mass = Fg/g
        massloss = (5.7/32.2) / (3.25 * 100)
       
        Cd = .3
        A = np.pi*(6/2)**2 #area of rocket
        p = 0.00002379 #density
        launchtheta = (85.0 * np.pi) / 180.0
        sin = np.sin(launchtheta)
        Fnet, h0, a, hf, Fd, t, V0, V,Fddrogue, Fdmain, Fsum  = (False,)*11
        Vw = 7.33*sin
        Vvector, FDvector, Avector, Tvector, Hvector = ([] for i in range(5))
        Rdrogue = (18/12.0)/2.0 #radius of drogue chute
        Vdrogue = 83.1 # terminal veloctiy with drogue chute
        Adrogue = (math.pi)*Rdrogue**2
        Vmain = 11.2
        Amain = math.pi * ((144/12.0)/2.0)**2
        Cddrogue = 2 * (Fg)/ ( p * Vdrogue**2 * Adrogue) 
        Cdmain= 2 * (Fg)/ ( p * Vmain**2 * Amain)
        h = 0
        for i in range(len(thrustarray)):
            mass = mass - massloss
            Tvector.append(t)
            Fd =  Cd*A*((p*(V**2)))*.5
            
            t = times[i]
            Fnet = (thrustarray[i]*sin) - (mass*32.2) - Fd
            a = Fnet/mass
            Avector.append(a)
            
            V = (a*t) + Vw
            Vvector.append(V)
            h += .5 * a * t**2
            Hvector.append(h)
            
        V += a*.01 + Vw
        while V > 0:
            Tvector.append(t)
            Fd =  Cd*A*((p*(V**2)))*.5

            Fnet = - (mass*32.2) - Fd

            a = Fnet/mass
            Avector.append(a)
            V += a*.01 
            Vvector.append(V)
            h += V*.01
            Hvector.append(h)
            t +=  .01

        V += (a*.01) + Vw
        while h > 600: #Drogue chute deployed

            Tvector.append(t)
            Fd = .5 * p * (V**2) * Adrogue * Cddrogue
            
            FDvector.append(Fd)
            Fsum = Fd - Fg

            a = Fsum / mass

            Hvector.append(h)
            Avector.append(a)
            V += (a*.01)
            h += V*.01
            Vvector.append(V)
            t += .01
            
        V += (a*.01) + Vw   
        while h > 0: #main chute
            Tvector.append(t)
            V += (a*.01)
            Vvector.append(V)
            # drogue chute opens
            Fd = .5 * p * (V**2) * Amain * Cdmain
            FDvector.append(Fd)
            Fsum = Fd - Fg
            a = Fsum / mass
            Avector.append(a)
            h += V*.01
            Hvector.append(h)
            t += .01

        #generate subplots
        fig, (plt1, plt2, plt3) = plt.subplots(1,3, figsize = (12,6), dpi = 80)
        plt.autoscale()

        #rotate and align labels
        fig.autofmt_xdate()
        fig.suptitle('Problem 6: Revision of 4 (Y)', fontsize = 15)

        #plot the subplots
        plt1.plot(Tvector, Hvector, color='red',label = 'Height')
        plt1.set_title('Height')
        plt1.set(xlabel = 'Time (s)', ylabel = 'Altitude (ft)')
        plt2.plot(Tvector,Vvector, color='black',label = 'Velocity')
        plt2.set_title('Velocity')
        plt2.set(xlabel = 'Time (s)', ylabel = 'Velocity (ft/s)')
        plt3.plot(Tvector, Avector, color='olive', label='Acceleration')
        plt3.set_title('Acceleration')
        plt3.set(xlabel = 'Time (s)', ylabel = 'Acceleration(ft/s^2)')
        
        #plt.tight_layout()

        plt.show()
        return Tvector

    def Xrocket(Tyvector, thrustarray,times):
        # variable/array initialization
        g = 32.2 # gravitational acceleration
        Fg = 50.1
        mass = Fg/g
        Ft = 310
        Cd = .3
        massloss = (5.7/32.2) / (3.25 * 100)
        launchtheta = (85.0 * np.pi) / 180.0
        cosine = np.cos(launchtheta)

        A = np.pi*(6/2)**2 #area of rocket
        p = 0.00002379 #density
        Fnet, h0, a, hf, Fd, t, V0, V, Fddrogue, Fdmain, Fsum  = (False,)*11
        Vw = 7.33*cosine

        Vvector, FDvector, Avector, Tvector, Xvector = ([] for i in range(5))
        x = 0
        Rdrogue = (18/12.0)/2.0 #radius of drogue chute
        Vdrogue = 83.1 # terminal veloctiy with drogue chute
        Adrogue = (math.pi)*Rdrogue**2
        Vmain = 11.2
        Amain = math.pi * ((144/12.0)/2.0)**2
        Cddrogue = 2 * (Fg)/ ( p * Vdrogue**2 * Adrogue) 
        Cdmain= 2 * (Fg)/ ( p * Vmain**2 * Amain)
        
        for i in range(len(thrustarray)):
            mass = mass - massloss
            Tvector.append(t)
            Fd =  Cd*A*((p*(V**2)))*.5 
            Fx = (thrustarray[i] * cosine) - Fd
            a = Fx / mass
            Avector.append(a)
            V += a*.01 + Vw/100
            Vvector.append(V)
            x += V*.01
            Xvector.append(x)
            t = times[i]

        V += (a*.01) + Vw
        while t < Tyvector[-1]:
            Tvector.append(t)
            Fd =  Cd*A*((p*(V**2)))*.5 
            Fx = Fd
            a = Fx / mass
            Avector.append(a)
            
            V += (a*.01)
            Vvector.append(V)
            x += V*.01
            Xvector.append(x)
            t += .01

        #generate subplots
        fig, (plt1, plt2, plt3) = plt.subplots(1,3, figsize = (12,6), dpi = 80)
        plt.autoscale()

        #rotate and align labels
        fig.autofmt_xdate()
        fig.suptitle('Problem 6: Revision of 4 (Y)', fontsize = 15)

        #plot the subplots
        plt1.plot(Xvector,Tvector, color='red',label='X travel')
        plt1.set_title('X travel')
        plt1.set(xlabel = 'X-Travel (Ft)', ylabel = 'Time (s)')
        plt2.plot(Tvector,Vvector, color='black',label = 'Velocity')
        plt2.set_title('Velocity')
        plt2.set(xlabel = 'Time (s)', ylabel = 'Velocity (ft/s)')
        plt3.plot(Tvector, Avector, color='olive', label='Acceleration')
        plt3.set_title('Acceleration')
        plt3.set(xlabel = 'Time (s)', ylabel = 'Acceleration (ft/s^2)')
        
        #plt.tight_layout()

        plt.show()  
    Tyvector = Yrocket(thrustarray, times)
    Xrocket(Tyvector, thrustarray,times)

def problem7():
    def Yrocket(thrustarray, times):
        # variable/array initialization
        g = 32.2 # gravitational acceleration
        Fg = 50.1
        mass = Fg/g
        massloss = (5.7/32.2) / (3.25 * 100)
        
        Cd = .3
        A = np.pi*(6/2)**2 #area of rocket
        p = 0.00002379 #density
        launchtheta = (85.0 * np.pi) / 180.0
        sin = np.sin(launchtheta)
        Fnet, h0, a, hf, Fd, t, V0, V,Fddrogue, Fdmain, Fsum  = (False,)*11
        Vw = 7.33*sin
        Vvector, FDvector, Avector, Tvector, Hvector = ([] for i in range(5))
        Rdrogue = (18/12.0)/2.0 #radius of drogue chute
        Vdrogue = 83.1 # terminal veloctiy with drogue chute
        Adrogue = (math.pi)*Rdrogue**2
        Vmain = 11.2
        Amain = math.pi * ((144/12.0)/2.0)**2
        Cddrogue = 2 * (Fg)/ ( p * Vdrogue**2 * Adrogue) 
        Cdmain= 2 * (Fg)/ ( p * Vmain**2 * Amain)
        h = 0

        for i in range(len(thrustarray)):
            mass = mass - massloss
            Tvector.append(t)
            Fd =  Cd*A*((p*(V**2)))*.5
            
            t = times[i]
            Fnet = (thrustarray[i]*sin) - (mass*32.2) - Fd
            a = Fnet/mass
            Avector.append(a)
            
            V = (a*t) + Vw
            Vvector.append(V)
            h += .5 * a * t**2
            Hvector.append(h)
            
        V += a*.01 + Vw
        while V > 0:
            Tvector.append(t)
            Fd =  Cd*A*((p*(V**2)))*.5
            Fnet = - (mass*32.2) - Fd
            a = Fnet/mass
            Avector.append(a)
            V += a*.01 
            Vvector.append(V)
            h += V*.01
            Hvector.append(h)
            t +=  .01

        V += (a*.01) + Vw
        while h > 600: #Drogue chute deployed

            Tvector.append(t)
            Fd = .5 * p * (V**2) * Adrogue * Cddrogue
            
            FDvector.append(Fd)
            Fsum = Fd - Fg
            a = Fsum / mass
            Hvector.append(h)
            Avector.append(a)
            V += (a*.01)
            h += V*.01
            Vvector.append(V)
            t += .01
            
        V += (a*.01) + Vw
            
        while h > 0: #main chute
            Tvector.append(t)
            V += (a*.01)
            Vvector.append(V)
            # drogue chute opens
            Fd = .5 * p * (V**2) * Amain * Cdmain
            FDvector.append(Fd)
            Fsum = Fd - Fg
            a = Fsum / mass
            Avector.append(a)
            h += V*.01
            Hvector.append(h)
            
            t += .01

        #generate subplots
        fig, (plt1, plt2, plt3) = plt.subplots(1,3, figsize = (12,6), dpi = 80)
        plt.autoscale()

        #rotate and align labels
        fig.autofmt_xdate()
        fig.suptitle('Problem 7 Y-Direction', fontsize = 15)

        #plot the subplots
        plt1.plot(Tvector, Hvector, color='red',label = 'Height')
        plt1.set_title('Altitude')
        plt1.set(xlabel = 'Time (s)', ylabel = 'Altitude (ft)')
        plt2.plot(Tvector,Vvector, color='black',label = 'Velocity')
        plt2.set_title('Velocity')
        plt2.set(xlabel = 'Time (s)', ylabel = 'Velocity (ft/s)')
        plt3.plot(Tvector, Avector, color='olive', label='Acceleration')
        plt3.set_title('Acceleration')
        plt3.set(xlabel = 'Time (s)', ylabel = 'Acceleration(ft/s^2)')
        
        #plt.tight_layout()

        plt.show()
        return Tvector

    def Xrocket(Tyvector, thrustarray,times):

        # variable/array initialization
        g = 32.2 # gravitational acceleration
        Fg = 50.1
        mass = Fg/g
        Ft = 310
        Cd = .3

        massloss = (5.7/32.2) / (3.25 * 100)
        launchtheta = (85.0 * np.pi) / 180.0
        cosine = np.cos(launchtheta)

        A = np.pi*(6/2)**2 #area of rocket
        p = 0.00002379 #density
        Fnet, h0, a, hf, Fd, t, V0, V, Fddrogue, Fdmain, Fsum  = (False,)*11
        Vw = 7.33*cosine

        Vvector, FDvector, Avector, Tvector, Xvector = ([] for i in range(5))
        x = 0
        Rdrogue = (18/12.0)/2.0 #radius of drogue chute
        Vdrogue = 83.1 # terminal veloctiy with drogue chute
        Adrogue = (math.pi)*Rdrogue**2
        Vmain = 11.2
        Amain = math.pi * ((144/12.0)/2.0)**2
        Cddrogue = 2 * (Fg)/ ( p * Vdrogue**2 * Adrogue) 
        Cdmain= 2 * (Fg)/ ( p * Vmain**2 * Amain)
        
        for i in range(len(thrustarray)):
            mass = mass - massloss
            Tvector.append(t)
            Fd =  Cd*A*((p*(V**2)))*.5 
            Fx = (thrustarray[i] * cosine) - Fd
            a = Fx / mass
            Avector.append(a)
            V += a*.01 + Vw/100
            Vvector.append(V)
            x += V*.01
            Xvector.append(x)
            t = times[i]

        V += (a*.01) + Vw
        while t < Tyvector[-1]:
            Tvector.append(t)
            Fd =  Cd*A*((p*(V**2)))*.5 
            Fx = Fd
            a = Fx / mass
            Avector.append(a)
            
            V += (a*.01)
            Vvector.append(V)
            x += V*.01
            Xvector.append(x)
            t += .01

        #generate subplots
        fig, (plt1, plt2, plt3) = plt.subplots(1,3, figsize = (12,6), dpi = 80)
        plt.autoscale()

        #rotate and align labels
        fig.autofmt_xdate()
        fig.suptitle('Problem 7 X-Direction', fontsize = 15)

        #plot the subplots
        plt1.plot(Xvector,Tvector, color='red',label='X travel')
        plt1.set_title('X travel')
        plt1.set(xlabel = 'X-Travel (Ft)', ylabel = 'Time (s)')
        plt2.plot(Tvector,Vvector, color='black',label = 'Velocity')
        plt2.set_title('Velocity')
        plt2.set(xlabel = 'Time (s)', ylabel = 'Velocity (ft/s)')
        plt3.plot(Tvector, Avector, color='olive', label='Acceleration')
        plt3.set_title('Acceleration')
        plt3.set(xlabel = 'Time (s)', ylabel = 'Acceleration (ft/s^2)')
        
        #plt.tight_layout()
        plt.show()  
    Tyvector = Yrocket(thrustarray, times)
    Xrocket(Tyvector, thrustarray,times)

def problem8():

    def measured():
        data = pd.read_csv('/Users/noahalex/Downloads/Rocket Experiment Data for Project (1).csv')
        Ytime = data['Time Axial A']
        Ya = data['Axial A']
        Xtime = data['Time Lateral A']
        Xa = data['Lateral A']
        VyVector, YVector, VxVector, XVector = ([] for i in range(4))
        X, Y, Vy, Vx = (False,)*4
        # 1 G  = 9.81 m/s^2
        for i in range(len(Ya)):
            try:
                ay = Ya[i]
                #convert to m/s^2
                ayy = ay*9.81

                t = Ytime[i+1] - Ytime[i]
                Vy += ayy*t
                VyVector.append(Vy)
                Y += Vy*t
                YVector.append(Y)
            except:
                break
        for i in range(len(Xa)):
            try:
                ax = Xa[i]

                #convert G's to m/s^2
                axx = ax * 9.81
                t = Xtime[i+1] - Xtime[i] 
                Vx += axx*t
                VxVector.append(Vx)
                X += Vx*t
                XVector.append(X)
            except:
                break

        XVector.append(XVector[-1])
        YVector.append(YVector[-1])

        fig, (plt1, plt2) = plt.subplots(1,2, figsize = (11,6), dpi = 80)
        plt.autoscale()
        fig.suptitle('Measured Data Problem 8', fontsize = 15)
        plt1.plot(Ytime, YVector, color = 'red', label = 'Height')
        plt1.set_title('Height')
        plt1.set(xlabel = 'Time (s)', ylabel = 'Altitude (ft)')
        plt2.plot(Xtime, XVector, color = 'green', label = 'X travel')
        plt2.set_title('X travel')
        plt2.set(xlabel = 'Time (s)', ylabel = 'X-Travel (ft)')
        plt.show()

    def barometer():
        data = pd.read_csv('/Users/noahalex/Downloads/Rocket Experiment Data for Project (1).csv')

        BaroTime = data['Time Baro']
        BaroAlt = data['Altitude Baro']
        VxVector, XVector = ([] for i in range(2))

        plt.plot(BaroTime, BaroAlt)
        plt.suptitle('Barometer Problem 8', fontsize = 15)
        plt.xlabel('Time (s)')
        plt.ylabel('Altitude (ft)')
        plt.show()

    def stratologgers():
        data = pd.read_csv('/Users/noahalex/Downloads/Stratalogger.csv')
        PrimeTime = data['Time 1']
        PrimeAlt = data['Altitude 1']
        SecondTime = data['Time 2']
        SecondAlt = data['Altitude 2']
        fig, (plt1, plt2) = plt.subplots(1,2, figsize = (12,6), dpi = 80)
        fig.suptitle('Stratalogger Data Problem 8', fontsize = 15)
        plt1.plot(PrimeTime, PrimeAlt, color='red',label = '')
        plt1.set_title('Primary Strata')
        plt1.set(xlabel = 'Time (s)', ylabel = 'Altitude (ft)')
        plt2.plot(SecondTime, SecondAlt, color='black',label = '')
        plt2.set_title('Secondary Strata')
        plt2.set(xlabel = 'Time (s)', ylabel = 'Altitude (ft)')
        plt.show()
    measured()
    barometer()
    stratologgers()

########## Script ###########
thrustarray, times = thrustcurve()

problem6_3(thrustarray, times)
problem6_4(thrustarray, times)

problem8()

