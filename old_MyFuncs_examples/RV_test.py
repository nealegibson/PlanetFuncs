#!/usr/bin/python

import numpy as np

from MyFuncs import RVCurve
import RadialVelFit as RV

t = np.arange(0,4.1,0.001)
# print t
# print RV.RadVel(t,1,2,3,4)

T0,P,Mstar,MPsini,V0,w,e = 0,2.,1.,1.,100,np.pi*3,0.3

par = (T0,P,Mstar,MPsini,0,0,e)
RVCurve(par,t,plot=True)

par = (T0,P,Mstar,MPsini,10,90.,e)
RVCurve(par,t,plot=True)

par = (T0,P,Mstar,MPsini,20,180.,e)
RVCurve(par,t,plot=True)

par = (T0,P,Mstar,MPsini,200,270.,e)
RVCurve(par,t,plot=True)

raw_input()
