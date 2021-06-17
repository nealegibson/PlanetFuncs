#!/usr/bin/python

import numpy as np
import pylab

from MyFuncs import PlotEccOrbit, RVCurve, EccLightCurve, PlotEccOrbit_aRs

#import Cnst
# import PlanetOrbit
# from EccLightCurveModel import EccLightCurve,PlotEccOrbit
# from RadialVelocityModel import RadialVelocity

#test ecc light curve model
#define system parameters
P = 4.353 #period / days
T_zero = 4.0 #central transit time / days
R_star = 1. #solar radii
p = 0.1 #ratio
c1, c2 = 0.2,0.2 #ld parameters
M_star = 1 #solar mass
M_planet = 1. #jupiter mass
i = 88. #radians
T_dur = P*12. #hours either side of centre
T_step = 0.1 #time step in hours
e = 0.3
w = 30
SecDepth = 0.001
RV_offset = 10

#set transit parameters and plot light curve
#note T_zero is ***central transit time*** for both models
tr_par = [T_zero,P,M_star,M_planet,R_star,p,i,c1,c2,e,w,1.0,0.,SecDepth]
rv_par = [T_zero,P,M_star,M_planet*np.sin(i*np.pi/180.),RV_offset,w,e]

#define times to run models
t_start = T_zero-P
t_end = T_zero+P
t = np.arange(t_start,t_end,0.001)

#make transit model
pylab.figure(1,(8,8))
pylab.subplot(211)
EccLightCurve(tr_par,t,plot=True)

#make rv model
pylab.subplot(212)
RVCurve(rv_par,t,plot=True)

#plot orbit in 3D
pylab.figure(2)
t = np.arange(t_start,t_end+0.05,0.05)
PlotEccOrbit(tr_par,t)

pylab.figure(3)
tr_aRs_par = [T_zero,P,9.1,p,0.2,c1,c2,e,w,1.0,0.,SecDepth]
PlotEccOrbit_aRs(tr_aRs_par,t)

raw_input()
