
import numpy as np

def AddLBM(func,X):
  """
  Wrapper function to add multiplicative linear basis model to any function.
  Parameters are added at the end, and X array 'stored' in function to enable
  standard optimisers to work, and to be included in GP functions as mean function etc.
  
  func - function called as func(par,x)
  X - 2D array [N x K] where N is number of data points and K is number of basis models
  
  """
  
  #just return the function if X is None, or shape is None (ie empty array)
  if X is None: return func
  else: assert X.ndim == 2, "X should be 2D"
  n = X.shape[1] #get number of parameters required for LBM
  if n == 0: return func

  def Model(par,x): #define model with standard syntax
    return np.dot(X,par[-n:]) * func(par[:-n],x)
  
  return Model
