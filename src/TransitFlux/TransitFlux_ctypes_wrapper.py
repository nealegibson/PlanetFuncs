
from __future__ import print_function

import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
import os
import glob

#read in c function
#dir = os.path.dirname(__file__) #needs to be defined relative to current file
#careful with import - different systems have different orders of glob files?
lib_file = glob.glob('{}/TransitFlux_ctypes*so'.format(os.path.dirname(__file__)))[0] #find the file
lib = ctypes.cdll.LoadLibrary(lib_file)

#define the function and argtypes and return type
flux_quad_c = lib.ctypes_flux_quad_np
flux_quad_c.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double]
flux_quad_c.restype = ctypes.c_int

#define wrapper
def FluxQuad_ctypes(z,p,c1,c2):
  
  #define empty array to get flux
  f = np.empty(z.shape,dtype=z.dtype)
  
  #call the function - output written to f
  flux_quad_c(z,f,z.size,p,c1,c2)
  
  #and return flux
  return f

#define the function and argtypes and return type
flux_nonlin_c = lib.ctypes_flux_nonlin_np
flux_nonlin_c.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS", ndim=1),ndpointer(ctypes.c_double, flags="C_CONTIGUOUS", ndim=1),ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double]
flux_nonlin_c.restype = ctypes.c_int

#define wrapper
def FluxNonlin_ctypes(z,p,cg1,g2,g3,g4):
  
  #define empty array to get flux
  f = np.empty(z.shape,dtype=z.dtype)
  
  #call the function - output written to f
  flux_nonlin_c(z,f,z.size,p,g1,g2,g3,g4)
  
  #and return flux
  return f

