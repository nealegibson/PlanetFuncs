
import numpy as np
from .TransitFlux import FluxQuad,FluxNonlin,FluxQuad_ctypes,FluxNonlin_ctypes

def transit_no_ld(par,t):
  """
  Simple transit model for cicular orbit and no limb darkening. This is computed purely
  in python.

  par is a parameter vector:
    par = [T0,P,aRs,rho,b,foot,Tgrad]
    where:
    T0 - central transit time
    P - orbital period
    aRs - semi-major axis in units of stellar radius (system scale) or a/R_star
    rho - planet-to-star radius ratio - the parameter of interest for trasmission spectroscopy
    b - impact parameter
    foot - out-of-transit flux
    Tgrad - linear gradient of basis function
  t is time (in same units as T0 and P)

  """
  
  #get the parameters from vector
  T0,P,a_Rstar,p,b,foot,Tgrad = par

  #calculate phase angle
  theta = (2*np.pi/P) * (t - T0)

  #normalised separation z
  z = np.sqrt( (a_Rstar*np.sin(theta))**2 + (b*np.cos(theta))**2 );

  #calculate flux
  f = np.ones(t.size) #out of transit
  f[z<=1-p] = 1-p**2 #fully in transit
  ind = (z > np.abs(1-p)) * (z<=1+p) #and finally for ingress/egress
  k0 = np.arccos((p**2+z[ind]**2-1)/2/p/z[ind]) # kappa0
  k1 = np.arccos((1-p**2+z[ind]**2)/2/z[ind]) # kappa1
  f[ind] = 1-(k0*p**2 + k1 - np.sqrt( 0.25*(4*z[ind]**2-(1+z[ind]**2-p**2)**2) ))/np.pi

  #modify according to linear basis function
  f = f * (foot + (t - T0) * Tgrad)

  #return flux
  return f

def transit_quad(par,t):
  """
  Lightcurve function for cicular orbits with quadratic limb darkening.
  """
  
  #get the parameters from vector
  T0,P,a_Rstar,p,b,c1,c2,foot,Tgrad = par
  
  #calculate phase angle
  theta = (2*np.pi/P) * (t - T0)
  
  #normalised separation z
  z = np.sqrt( (a_Rstar*np.sin(theta))**2 + (b*np.cos(theta))**2 );

  #calculate flux
  f = FluxQuad(z,p,c1,c2) * (foot + (t - T0) * 24. * Tgrad) #time in hours for foot and Tgrad!
    
  #return flux
  return f

def transit_quad_ctypes(par,t):
  """
  Lightcurve function for cicular orbits with quadratic limb darkening.
  This version uses the ctypes wrapper rather than python C module. Included for testing.
  """

  #get the parameters from vector
  T0,P,a_Rstar,p,b,c1,c2,foot,Tgrad = par
    
  #calculate phase angle
  theta = (2*np.pi/P) * (t - T0)
  
  #normalised separation z
  z = np.sqrt( (a_Rstar*np.sin(theta))**2 + (b*np.cos(theta))**2 );

  #calculate flux
  f = FluxQuad_ctypes(z,p,c1,c2) * (foot + (t - T0) * 24. * Tgrad) #time in hours for foot and Tgrad!
    
  #return flux
  return f

def transit_quad_norm(par,t):
  """Lightcurve function for cicular orbits, from aRs etc..."""

  #get the parameters from vector
  T0,P,a_Rstar,p,b,c1,c2,foot,Tgrad = par
  
  #simply return the normalisation from the same par vector
  return foot + (t - T0) * 24. * Tgrad

def transit_quad_sq(par,t):
  """Lightcurve function for cicular orbits, from aRs etc..."""

  #get the parameters from vector
  T0,P,a_Rstar,p,b,c1,c2,foot,Tgrad,Tgrad2 = par
    
  #calculate phase angle
  theta = (2*np.pi/P) * (t - T0)
  
  #normalised separation z
  z = np.sqrt( (a_Rstar*np.sin(theta))**2 + (b*np.cos(theta))**2 );

  #calculate flux
  f = FluxQuad(z,p,c1,c2) * (foot + (t - T0) * 24. * Tgrad + ((t - T0) * 24.)**2 * Tgrad2) #time in hours for foot and Tgrad!

  #return flux
  return f

def transit_quad_sq_norm(par,t):
  
  #get the parameters from vector
  T0,P,a_Rstar,p,b,c1,c2,foot,Tgrad,Tgrad2 = par

  #return baseline flux
  return np.ones(t.size) * (foot + (t - T0) * 24. * Tgrad + ((t - T0) * 24.)**2 * Tgrad2)

def transit_nonlin(par,t):
  """Lightcurve function for cicular orbits, from aRs etc..."""

  #get the parameters from vector
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

def transit_nonlin_sq(par,t):
  """Lightcurve function for cicular orbits, from aRs etc..."""
  
  #get the parameters from vector
  T0,P,a_Rstar,p,b,c1,c2,c3,c4,foot,Tgrad,Tgrad2 = par
      
  #calculate phase angle
  theta = (2*np.pi/P) * (t - T0)
  
  #normalised separation z
  z = np.sqrt( (a_Rstar*np.sin(theta))**2 + (b*np.cos(theta))**2 );

  #calculate flux
  f = FluxNonlin(z,p,c1,c2,c3,c4) * (foot + (t - T0) * 24. * Tgrad + ((t - T0) * 24.)**2 * Tgrad2) #time in hours for foot and Tgrad!

  if np.any(np.isnan(f)):
    print ("some returned fluxes are nan!")
  
  #return flux
  return f

def transit_nonlin_sq_norm(par,t):
  """Return the normalisation part only..."""
  
  #get the parameters from vector
  T0,P,a_Rstar,p,b,c1,c2,c3,c4,foot,Tgrad,Tgrad2 = par
  
  #return norm  return (foot + (t - T0) * 24. * Tgrad + (t - T0) * 24. * Tgrad2)
  return (foot + (t - T0) * 24. * Tgrad + ((t - T0) * 24.)**2 * Tgrad2)

