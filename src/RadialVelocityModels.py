
import numpy as np
import pylab

from . import Cnst
from .RadialVelocity import RadVel

###############################################################################

def RVCurve(par,t,plot=False):
  """Radial velocity function."""
  
  #convert parameters
  T0,P,Mstar,MPsini,V0,w,e = par
  MPsini *= Cnst.MJup
  Mstar *= Cnst.MSun
  w *= np.pi / 180.
  
  #calculate K
  K = np.power(2*np.pi*Cnst.G/(P*24*3600),1./3.) \
             *(MPsini/np.power((Mstar+MPsini)**2,1./3.))*(1./(np.sqrt(1-e**2)))
  
  #make w lie in range 0-2pi
  if w >= 2*np.pi:
    w -= 2*np.pi #make f lie in range 0-2pi
  elif w < 0:
    w += 2*np.pi #make f lie in range 0-2pi
  
  #true anomaly of central transit time
  f = 1.*np.pi/2. + w
  if f >= 2*np.pi:
    f -= 2*np.pi #make f lie in range 0-2pi
  elif f < 0:
    f += 2*np.pi #make f lie in range 0-2pi
  
  if f < np.pi:
    E = np.arccos( (np.cos(f) + e) / (e*np.cos(f)+1.) )
    M_tr = E - e*np.sin(E)
    T_peri = T0 + M_tr * P/(2*np.pi)
  
  if f >= np.pi:
    #f = np.pi - f #correct for acos calc
    E = np.arccos( (np.cos(f) + e) / (e*np.cos(f)+1.) )
    M_tr = E - e*np.sin(E)
    #M_tr = 2*np.pi - M_tr
    T_peri = T0 - M_tr * P/(2*np.pi)
  
  #true anomaly of secondary transit time
  f = 1.*np.pi/2. - w
  if f >= 2*np.pi:
    f -= 2*np.pi #make f lie in range 0-2pi
  elif f < 0:
    f += 2*np.pi #make f lie in range 0-2pi
  
  if f < np.pi:
    E = np.arccos( (np.cos(f) + e) / (e*np.cos(f)+1.) )
    M_sec = E - e*np.sin(E)
    T_sec = T_peri + M_sec * P/(2*np.pi)
  
  if f >= np.pi:
    E = np.arccos( (np.cos(f) + e) / (e*np.cos(f)+1.) )
    M_sec = E - e*np.sin(E)
    T_sec = T_peri - M_sec * P/(2*np.pi)
  
  #calculate mean anomaly
  M = (2*np.pi/P) * (t - T_peri)
  
  #calculate Radial Velocity using C-function
  RV = RadVel(M, e, P, w+np.pi, K) + V0 #pass w+np.pi as w is arg_peri of planet
  
  if plot:
    pylab.plot(t,RV)
    pylab.axvline(x=T0,color='r',linestyle='-') #mark the central transit time
    pylab.axvline(x=T_peri,color='g',linestyle='-') #mark the pericentre passage
    pylab.axvline(x=T_sec,color='r',linestyle='-.') #mark the secondary transit time
    pylab.axhline(y=V0,color='g',linestyle='--') #mark the secondary transit time
    pylab.xlabel('Time / days')
    pylab.ylabel('Radial Velocity (m/s)')
    pylab.xlim(t.min(),t.max())
  
  return RV

###############################################################################
