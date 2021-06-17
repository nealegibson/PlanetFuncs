#!/usr/bin/python

import numpy as np
import pylab
import time as pytime

import PlanetFuncs as MF # edit to see if works on PlanetFuncs
#from LightCurveModel import LightCurveQuad,LightCurveNonlin

#HD189733 properties
#ephemeris from Agol et al. (2010)
P = 2.21857567 #period / days
T_zero_Run1 = 2454279.436714 #central transit time / days
R_star = 0.755 #solar radii
p_ACS = 0.1572 #radius ratio from ACS data (white light curve)
p_NIC = 0.155 #radius ratio from NICMOS data (approx from white noise)
p_input = 0.155 #initial radius guess
M_star = 0.825 #solar masses
M_planet = 1.15 #jupiter mass
inc = 85.68 #degrees
c1,c2,c3,c4 = 0.3,0.2,0.0,0.0

g1,g2,g3,g4 = 0,c1,0,c2

par_quad = [T_zero_Run1,P,M_star,M_planet,R_star,p_input,inc,c1,c2,1.00,0.0]
par_nonlin = [T_zero_Run1,P,M_star,M_planet,R_star,p_input,inc,g1,g2,g3,g4,1.00,0.0]

time = np.linspace(2454279.2,2454279.7,200)

#pylab.plot(time,LightCurveQuad(time,par_quad),'r-')
#pylab.plot(time,LightCurveNonlin(time,par_nonlin),'b--')

print "TEST:"
pylab.figure(1)
z = np.arange(0.,1.2,0.01)
#z = np.linspace(0.,1.2,100)
f = MF.FluxNonlin(z,0.1,0.1,0.1,0.1,0.1)
pylab.plot(z,f,'.--')

pylab.figure(2)
z = np.linspace(0., 1.2, 1000)
cns = np.vstack((np.zeros(4), np.eye(4)))
for coef in cns:
  f = MF.FluxNonlin(z,0.1,*coef)
  pylab.plot(z, f)

z = np.linspace(0., 3., 5000)

pylab.figure(3)
t0 = pytime.time()
ld = [0.76895072,-0.57463144,0.78942074,-0.34030351]
for i,p in enumerate(np.arange(0.01,1.2,0.05)):
#  print i
  pylab.plot(z, MF.FluxNonlin(z,p,*ld))
print pytime.time() - t0, 'secs'

pylab.figure(4)
t0 = pytime.time()
for i,p in enumerate(np.arange(0.01,1.2,0.05)):
#  print i, p
  pylab.plot(z, MF.FluxQuad(z,p,0.1,0.1))
print pytime.time() - t0, 'secs'

z = np.linspace(0., 1.8, 1000)

pylab.figure(5)
t0 = pytime.time()
for i,p in enumerate(np.arange(0.01,0.4,0.01)):
#  print i
  pylab.plot(z, MF.FluxNonlin(z,p,*cns[0]))
print "nonlin model:",pytime.time() - t0, 'secs'

# K = np.mat(np.random.rand(1000000).reshape(1000,1000))
# t0 = pytime.time()
# for i,p in enumerate(np.arange(0.01,0.4,0.01)):
# #  print i
#   a = K.I
# print "matrix inverse:", pytime.time() - t0, 'secs'

pylab.figure(6)
t0 = pytime.time()
for i,p in enumerate(np.arange(0.01,0.4,0.01)):
#  print i, p
  pylab.plot(z, MF.FluxQuad(z,p,0.1,0.1))
print "quadratic model:",pytime.time() - t0, 'secs'

raw_input()