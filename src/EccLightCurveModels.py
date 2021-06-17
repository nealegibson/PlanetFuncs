
import numpy as np
import pylab
from mpl_toolkits.mplot3d import Axes3D

from .TransitFlux import FluxQuad
from . import Cnst
from . import PlanetOrbit

"""
These have not been checked/tested in a long time - just including so can be modified
in the future.

Also would need to optimised, e.g. remove plotting options etc.

"""

###############################################################################

def EccLightCurve(par,t,plot=False):
  """Lightcurve function for eccentric orbits."""

  #read in parameters
  T0,P,Mstar,Mplanet,Rstar,p,i,c1,c2,e,w,foot,Tgrad,Sec_depth = par
  
  #convert units
  Rstar *= Cnst.RSun
  Mstar *= Cnst.MSun
  Mplanet *= Cnst.MJup
  i *= np.pi / 180.
  w *= np.pi / 180.

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
    E = np.arccos( (np.cos(f) + e) / (e*np.cos(f)+1.) )
    M_tr = E - e*np.sin(E)
    T_peri = T0 - M_tr * P/(2*np.pi)
  
  #calculate mean anomaly
  M = (2*np.pi/P) * (t - T_peri)
    
  #semi-major axis
  a = np.power( ( ((P*24.*60.*60./(2*np.pi))**2 * Cnst.G * (Mstar+Mplanet) ) ) , (1./3.) )
  
  #calc normalised separation z using PlanetOrbit functions
  norm = PlanetOrbit.get_norm(M,a/Rstar,e,w,i)
  y_coord = PlanetOrbit.get_y(M,a/Rstar,e,w,i) #get y coord to separate primary and secondary transits
  
  #print "a/R =",a/Rstar
  
  #calculate flux
  f = np.ones(norm.size)
  
  #primary transit
  index = np.where(y_coord<0) #ie planet closer than star
  f[index] = FluxQuad(norm[index],p,c1,c2) * (foot + (t[index] - T0) * 24. * Tgrad) #time in hours for foot and Tgrad!
  
  #secondary transit
  index = np.where(y_coord>0) #ie star closer than planet
  f[index] = FluxQuad(norm[index],p,0,0) * (foot + (t[index] - T0) * 24. * Tgrad) #time in hours for foot and Tgrad!
  f[index] = (f[index]-1.) * Sec_depth/p**2 + 1 #scale the transit depth to secondary depth
  f[index] *= (foot + (t[index] - T0)*24. * Tgrad) 
  
  if plot:
    pylab.plot(t,f,'r-')
    pylab.xlabel('Time / days')
    pylab.ylabel('Relative Flux')
    pylab.ylim(f.min()-(f.max()-f.min())/2.,f.max()+(f.max()-f.min())/2.)
    pylab.xlim(t.min(),t.max())

  #return flux
  return f

###############################################################################
def EccLightCurve_aRs(par,t,plot=False):
 """Lightcurve function for eccentric orbits with aRs parameterisation."""

 #read in parameters
 T0,P,a_Rstar,p,b,c1,c2,e,w,foot,Tgrad,Sec_depth = par

 #ensure b and p >= 0
 if b<0.: b=-b
 if p<0.: p=-p

 w *= np.pi / 180.
 i = np.arccos(b/a_Rstar)

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
   E = np.arccos( (np.cos(f) + e) / (e*np.cos(f)+1.) )
   M_tr = E - e*np.sin(E)
   T_peri = T0 - M_tr * P/(2*np.pi)

 #calculate mean anomaly
 M = (2*np.pi/P) * (t - T_peri)

 #calc normalised separation z using PlanetOrbit functions
 norm = PlanetOrbit.get_norm(M,a_Rstar,e,w,i)
 y_coord = PlanetOrbit.get_y(M,a_Rstar,e,w,i) #get y coord to separate primary and secondary transits

 #print "a/R =",a/Rstar

 #calculate flux
 f = np.ones(norm.size)

 #primary transit
 index = np.where(y_coord<0) #ie planet closer than star
 f[index] = FluxQuad(norm[index],p,c1,c2) * (foot + (t[index] - T0) * 24. * Tgrad) #time in hours for foot and Tgrad!

 #secondary transit
 index = np.where(y_coord>0) #ie star closer than planet
 f[index] = FluxQuad(norm[index],p,0,0) * (foot + (t[index] - T0) * 24. * Tgrad) #time in hours for foot and Tgrad!
 f[index] = (f[index]-1.) * Sec_depth/p**2 + 1 #scale the transit depth to secondary depth
 f[index] *= (foot + (t[index] - T0)*24. * Tgrad) 

 if plot:
   pylab.plot(t,f,'r-')
   pylab.xlabel('Time / days')
   pylab.ylabel('Relative Flux')
   pylab.ylim(f.min()-(f.max()-f.min())/2.,f.max()+(f.max()-f.min())/2.)
   pylab.xlim(t.min(),t.max())

 return f

###############################################################################

def PlotEccOrbit(par,t):
  """Function to plot planet orbit in 3D"""
  
  #read in parameters
  T0,P,Mstar,Mplanet,Rstar,p,i,c1,c2,e,w,foot,Tgrad,Sec_depth = par

  #convert units
  Rstar *= Cnst.RSun
  Mstar *= Cnst.MSun
  Mplanet *= Cnst.MJup
  i *= np.pi / 180.
  w *= np.pi / 180.

  #semi-major axis
  a = np.power( ( ((P*24.*60.*60./(2*np.pi))**2 * Cnst.G * (Mstar+Mplanet) ) ) , (1./3.) )
  
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

  #calculate mean anomaly
  M = (2*np.pi/P) * (t - T_peri)
  
  #get coords
  x = PlanetOrbit.get_x(M,a/Rstar,e,w)
  y = PlanetOrbit.get_y(M,a/Rstar,e,w,i)
  z = PlanetOrbit.get_z(M,a/Rstar,e,w,i)
  
  #make plot
  ax = Axes3D(pylab.gcf())
  ax.plot(x, y, z, c='k')
  ax.scatter(x, y, z, c='r', s=50)
  ax.scatter([0],[0],[0],c='y', s=500) #plot star position
  ax.scatter([x[0],],[y[0],],[z[0],],c='g', s=100) #plot initial planet position
  ax.scatter([x[1],],[y[1],],[z[1],],c='y', s=100) #plot initial planet position
  ax.set_xlabel('X')
  ax.set_ylabel('Y')
  ax.set_zlabel('Z')
  range = abs(np.array([x,y,z])).max()
  ax.set_xlim3d(-range,range)
  ax.set_ylim3d(-range,range)
  ax.set_zlim3d(-range,range)

###############################################################################

def PlotEccOrbit_aRs(par,t):
  """Function to plot planet orbit in 3D"""
  
  #read in parameters
  T0,P,a_Rstar,p,b,c1,c2,e,w,foot,Tgrad,Sec_depth = par

  #ensure b and p >= 0
  if b<0.: b=-b
  if p<0.: p=-p

  w *= np.pi / 180.
  i = np.arccos(b/a_Rstar)
  
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

  #calculate mean anomaly
  M = (2*np.pi/P) * (t - T_peri)
  
  #get coords
  x = PlanetOrbit.get_x(M,a_Rstar,e,w)
  y = PlanetOrbit.get_y(M,a_Rstar,e,w,i)
  z = PlanetOrbit.get_z(M,a_Rstar,e,w,i)
  
  #make plot
  ax = Axes3D(pylab.gcf())
  ax.plot(x, y, z, c='k')
  ax.scatter(x, y, z, c='r', s=50)
  ax.scatter([0],[0],[0],c='y', s=500) #plot star position
  ax.scatter([x[0],],[y[0],],[z[0],],c='g', s=100) #plot initial planet position
  ax.scatter([x[1],],[y[1],],[z[1],],c='y', s=100) #plot initial planet position
  ax.set_xlabel('X')
  ax.set_ylabel('Y')
  ax.set_zlabel('Z')
  range = abs(np.array([x,y,z])).max()
  ax.set_xlim3d(-range,range)
  ax.set_ylim3d(-range,range)
  ax.set_zlim3d(-range,range)

###############################################################################
