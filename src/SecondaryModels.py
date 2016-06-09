
import numpy as np

from TransitFlux import FluxQuad,FluxNonlin
import Cnst

#####################################################################

def Sec_from_Transit_aRs(par,t):
  """Lightcurve function for cicular orbits, from aRs etc..."""

  T0,P,a_Rstar,p,b,foot,Tgrad,depth = par
  
  #ensure b and p >= 0
  if b<0.: b=-b
  if p<0.: p=-p
  
  #calculate phase angle
  theta = (2*np.pi/P) * (t - T0)
  
  #normalised separation z
  z = np.sqrt( (a_Rstar*np.sin(theta))**2 + (b*np.cos(theta))**2 );

  #calculate flux - use rescaled transit model with no limb darkening
  f = (FluxQuad(z,p,0.,0.) - 1.) / p**2. * depth + 1
  f *= (foot + (t - T0) * 24. * Tgrad) #time in hours for foot and Tgrad!
  
  #return flux
  return f

#####################################################################

def Sec_from_Transit_aRs_sq(par,t):
  """Lightcurve function for cicular orbits, from aRs etc..."""

  T0,P,a_Rstar,p,b,foot,Tgrad,Tgrad2,depth = par
  
  #ensure b and p >= 0
  if b<0.: b=-b
  if p<0.: p=-p
  
  #calculate phase angle
  theta = (2*np.pi/P) * (t - T0)
  
  #normalised separation z
  z = np.sqrt( (a_Rstar*np.sin(theta))**2 + (b*np.cos(theta))**2 );

  #calculate flux - use rescaled transit model with no limb darkening
  f = (FluxQuad(z,p,0.,0.) - 1.) / p**2. * depth + 1
  #f *= (foot + (t - T0) * 24. * Tgrad) #time in hours for foot and Tgrad!
  f *= (foot + (t - T0) * 24. * Tgrad + ((t - T0) * 24.)**2 * Tgrad2)
  #return flux
  return f

#####################################################################
