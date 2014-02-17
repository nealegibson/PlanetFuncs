
import numpy as np

from TransitFlux import FluxQuad,FluxNonlin
import Cnst

def Transit(pars,t):
  """Lightcurve function for cicular orbits, from 'normal' physical parameters of star"""

  #assign parameters
  T0,P,Mstar,Mplanet,Rstar,p,i,c1,c2,foot,Tgrad = pars
  
  #ensure p is positive
  if p<0: p=-p
  
  #convert units
  Rstar *= Cnst.RSun
  Mstar *= Cnst.MSun
  Mplanet *= Cnst.MJup
  i *= np.pi / 180.
  
  #calculate phase angle
  theta = (2*np.pi/P) * (t - T0)
  
  #semi-major axis
  a = np.power( ( ((P*24.*60.*60./(2*np.pi))**2 * Cnst.G * (Mstar+Mplanet) ) ) , (1./3.) )
  
  #normalised separation z
  z = (a/Rstar) * np.sqrt(np.sin(theta)**2 + (np.cos(i)**2)*(np.cos(theta)**2))

  #calculate flux
  f = FluxQuad(z,p,c1,c2) * (foot + (t - T0) * 24. * Tgrad) #time in hours for foot and Tgrad!
  
  #return flux
  return f

def Transit_Norm(pars,t):
  """Simply returns the normalisation of the transit function, ignoring other params"""

  #assign parameters
  T0,P,Mstar,Mplanet,Rstar,p,i,c1,c2,foot,Tgrad = pars
  
  #return baseline flux
  return np.ones(t.size)*(foot + (t - T0) * 24. * Tgrad) 

def TransitNL(pars,t):
  """Lightcurve function for cicular orbits with non-linear limb darkening law."""
  
  T0,P,Mstar,Mplanet,Rstar,p,i,c1,c2,c3,c4,foot,Tgrad = pars
  
  #ensure p is positive
  if p<0: p=-p

  #convert units
  Rstar *= Cnst.RSun
  Mstar *= Cnst.MSun
  Mplanet *= Cnst.MJup
  i *= np.pi / 180.
  
  #calculate phase angle
  theta = (2*np.pi/P) * (t - T0)
  
  #semi-major axis
  a = np.power( ( ((P*24.*60.*60./(2*np.pi))**2 * Cnst.G * (Mstar+Mplanet) ) ) , (1./3.) )
  
  #normalised separation z
  z = (a/Rstar) * np.sqrt(np.sin(theta)**2 + (np.cos(i)**2)*(np.cos(theta)**2))

  #calculate flux
  f = FluxNonlin(z,p,c1,c2,c3,c4) * (foot + (t - T0) * 24. * Tgrad) #time in hours for foot and Tgrad!
    
  #return flux
  return f

def TransitNL_Norm(pars,t):
  """Lightcurve function for cicular orbits with non-linear limb darkening law."""
  
  #assign parameters
  T0,P,Mstar,Mplanet,Rstar,p,i,c1,c2,c3,c4,foot,Tgrad = pars
      
  #return baseline flux
  return np.ones(t.size)*(foot + (t - T0) * 24. * Tgrad) 

def Transit_aRs(par,t):
  """Lightcurve function for cicular orbits, from aRs etc..."""

  T0,P,a_Rstar,p,b,c1,c2,foot,Tgrad = par
  
  #ensure b and p >= 0
  if b<0.: b=-b
  if p<0.: p=-p
  
  #calculate phase angle
  theta = (2*np.pi/P) * (t - T0)
  
  #normalised separation z
  z = np.sqrt( (a_Rstar*np.sin(theta))**2 + (b*np.cos(theta))**2 );

  #calculate flux
  f = FluxQuad(z,p,c1,c2) * (foot + (t - T0) * 24. * Tgrad) #time in hours for foot and Tgrad!
    
  #return flux
  return f

def Transit_aRs_norm(par,t):
  """Lightcurve function for cicular orbits, from aRs etc..."""

  T0,P,a_Rstar,p,b,c1,c2,foot,Tgrad = par
  
  #simply return the normalisation from the same par vector
  return foot + (t - T0) * 24. * Tgrad

def Transit_aRs_sq(par,t):
  """Lightcurve function for cicular orbits, from aRs etc..."""

  T0,P,a_Rstar,p,b,c1,c2,foot,Tgrad,Tgrad2 = par
  
  #ensure b and p >= 0
  if b<0.: b=-b
  if p<0.: p=-p
  
  #calculate phase angle
  theta = (2*np.pi/P) * (t - T0)
  
  #normalised separation z
  z = np.sqrt( (a_Rstar*np.sin(theta))**2 + (b*np.cos(theta))**2 );

  #calculate flux
  f = FluxQuad(z,p,c1,c2) * (foot + (t - T0) * 24. * Tgrad + ((t - T0) * 24.)**2 * Tgrad2) #time in hours for foot and Tgrad!

  #return flux
  return f

def Transit_aRs_sq_norm(par,t):
  """Lightcurve function for cicular orbits with non-linear limb darkening law."""
  
  #assign parameters
  T0,P,a_Rstar,p,b,c1,c2,foot,Tgrad,Tgrad2 = par

  #return baseline flux
  return np.ones(t.size) * (foot + (t - T0) * 24. * Tgrad + ((t - T0) * 24.)**2 * Tgrad2)

def TransitNL_aRs(par,t):
  """Lightcurve function for cicular orbits, from aRs etc..."""

  T0,P,a_Rstar,p,b,c1,c2,c3,c4,foot,Tgrad = par
  
  #ensure b and p >= 0
  if b<0.: b=-b
  if p<0.: p=-p
  
  #calculate phase angle
  theta = (2*np.pi/P) * (t - T0)
  
  #normalised separation z
  z = np.sqrt( (a_Rstar*np.sin(theta))**2 + (b*np.cos(theta))**2 );

  #calculate flux
  f = FluxNonlin(z,p,c1,c2,c3,c4) * (foot + (t - T0) * 24. * Tgrad) #time in hours for foot and Tgrad!
  
  #return flux
  return f

def TransitNL_aRs_sq(par,t):
  """Lightcurve function for cicular orbits, from aRs etc..."""

  T0,P,a_Rstar,p,b,c1,c2,c3,c4,foot,Tgrad,Tgrad2 = par
  
  #ensure b and p >= 0
  if b<0.: b=-b
  if p<0.: p=-p
  
  #calculate phase angle
  theta = (2*np.pi/P) * (t - T0)
  
  #normalised separation z
  z = np.sqrt( (a_Rstar*np.sin(theta))**2 + (b*np.cos(theta))**2 );

  #calculate flux
  f = FluxNonlin(z,p,c1,c2,c3,c4) * (foot + (t - T0) * 24. * Tgrad + (t - T0) * 24. * Tgrad2) #time in hours for foot and Tgrad!
  
  #return flux
  return f


