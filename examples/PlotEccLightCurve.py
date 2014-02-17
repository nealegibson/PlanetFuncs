#!/usr/bin/python

import numpy as np
import pylab

import MyFunctions as MF

#light curve parameters
pars = [2454000.0,1.,1.,1.,1.,.1,88.,0.2,0.3,0.8,90.,1.,0.,0.01]
pars_aRs = [2454000.0,1.,10.,0.1,0.2,0.2,0.3,0.8,90.,1.,0.,0.01]
wn = 0.0003

#create the data set and plot (ie training data)
time = np.linspace(2453999.9,2454000.1,50000)
flux = MF.EccLightCurve(pars,time)
pylab.figure(1)
pylab.plot(time,flux,'r-')
pylab.plot(time,flux+np.random.normal(0,wn,time.size),'.')

#test with the aRs model
flux = MF.EccLightCurve_aRs(pars_aRs,time)
pylab.figure(2)
pylab.plot(time,flux,'r-')
pylab.plot(time,flux+np.random.normal(0,wn,time.size),'.')

#test for memory leaks
# flux = np.zeros(time.size)
# for i in range(100000):
#   print i
#   pars_aRs[4] = np.random.rand()/10. + 0.2
#   flux = MF.EccLightCurve_aRs(pars_aRs,time)

raw_input()
