#!/usr/bin/python

import numpy as np
import pylab

import MyFunctions as MF

#light curve parameters
lc_pars = [2454000.0,1.,11.,.1,0.6,0.2,0.3,1.,0.]
wn = 0.0003

#create the data set (ie training data)
time = np.arange(2453999.5,2454000.5,0.001)
flux = MF.Transit_aRs(lc_pars,time)

pylab.plot(time,flux,'r-')
pylab.plot(time,flux+np.random.normal(0,wn,time.size),'.')

#test for memory leaks
# time = np.linspace(2453999.5,2454000.5,50000)
# flux = np.zeros(time.size)
# for i in range(100000):
#   print i
#   lc_pars[4] = np.random.rand()/10. + 0.2
#   flux = MF.Transit_aRs(lc_pars,time)
#  flux = np.concatenate([flux,MF.Transit_aRs(lc_pars,time)])
#print flux.size
raw_input()
