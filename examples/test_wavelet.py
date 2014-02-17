#!/usr/bin/python

import numpy as np
import pylab

import MyFunctions as MF

#light curve parameters
lc_pars = [2454000.0,1.,11.,.1,0.6,0.2,0.3,1.,0.]
wn = 0.0003 #white noise
w_pars = [wn,0.000,0.1]

#create the data set (ie training data)
time = np.linspace(2453999.8,2454000.2,1000)
flux = MF.Transit_aRs(lc_pars,time) + np.random.normal(0,wn,time.size)

pylab.plot(time,MF.Transit_aRs(lc_pars,time),'r-')
pylab.plot(time,flux+np.random.normal(0,wn,time.size),'.')

for i in range(100):
  x = np.zeros(500.) + np.random.normal(0,1,500)
  print MF.LogLikelihood_iid(x,1.), MF.WaveletLogLikelihood(x,[1.,0.000,1.])

print MF.LogLikelihood_iid_mf(lc_pars+[wn,],MF.Transit_aRs,time,flux), MF.WaveletLogLikelihood_mf(lc_pars+w_pars,MF.Transit_aRs,time,flux)

raw_input()
