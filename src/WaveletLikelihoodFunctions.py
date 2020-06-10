"""
Docstring...
"""

import numpy as np
from . import WaveletLikelihood as WL

###############################################################################
#'main' log likelihood functions here
def WaveletLogLikelihood(r,p):
  """
  logP function from wavelet transform of Carter & Winn 2009
  r - residuals
  p - (sig_w,sig_r,gamma) parameter vector
  """
  
  logP = WL.WaveletLikelihood(r,p[0],p[1],p[2])
  
  return logP

def WaveletLogLikelihood_mf(p,func,x,y):
  """
  logP function from wavelet transform of Carter & Winn 2009
  first residuals are created from the mean function
  
  p - mean function parameter + wavelet parameters (sig_w,sig_r,gamma)
  func - mean function
  x - mean function input
  y - measurements
  
  """
  
  r = y - func(p[:-3],x)
  
  return WaveletLogLikelihood(r,p[-3:])

###############################################################################
