#!/usr/bin/python

import numpy as np
import pylab

import MyFuncs as MF
import DataAnalysis as DA

mean_func = MF.Transit_aRs

#define the likelihood function to maximise
#pointless strings to test that only p needs to be the first variable passed
def LogLikelihood(p,pointless_string,func,x,y,sig=1.,dummy='pointless string'):
  r = y - func(p[:-1],x)
  return MF.LogLikelihood_iid(r,p[-1],sig=1.)
#or function to minimise
def NegLogLikelihood(p,func,x,y,sig=1.):
  r = y - func(p[:-1],x)
  return -MF.LogLikelihood_iid(r,p[-1],sig=1.)

#define function parameters (transit)
n_time = 250
n_time_pred = 1000
P = 4. #period / days
T0 = 0.
a_Rstar = 12. #solar radii
p_in = 0.1 #radius ratio
b = 0.25 #degrees
c1,c2 = 0.2,0.2
par_tr = [T0,P,a_Rstar,p_in,b,c1,c2,1.0,0.0]
gp_tr = [T0+0.001,P,a_Rstar-3.,p_in+0.005,b-0.1,c1,c2,1.0,0.0]
fixed_par_tr = [0,1,0,0,1,1,1,0,0]
err_par_tr = [0.0001,0,0.1,0.0005,0.0,0,0,0.00005,0.]
time = np.linspace(-0.2,0.2,n_time)
time_pred = np.linspace(-0.3,0.3,n_time_pred)

#create the transit function
flux = mean_func(par_tr,time) + DA.WhiteNoise(n_time) * 0.001

#create the guess params and fixed params
gp = gp_tr + [0.001,]
fixed_par = fixed_par_tr + [0,]

#call the optimiser to maximise the log likelihood
fp = MF.Optimise(LogLikelihood,gp,("another pointless string",mean_func,time,flux,'randstring'),fixed=fixed_par,method='NM')
#or minimise the neglog likelihood
#fp = ML.Optimise(NegLogLikelihood,gp,(mean_func,time,flux),type='min',fixed=fixed_par)

#create a plot
pylab.plot(time_pred,mean_func(gp[:-1],time_pred),'g-')
pylab.plot(time_pred,mean_func(fp[:-1],time_pred),'r-')
pylab.errorbar(time,flux,yerr=fp[-1],color='k',linestyle='.')
pylab.axhline(0.98,color='r')
pylab.errorbar(time,flux-mean_func(fp[:-1],time)+0.98,yerr=fp[-1],color='k',linestyle='.')

raw_input()