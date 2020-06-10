
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
    return (1.+np.dot(X,par[-n:])) * func(par[:-n],x)
  
  def Model_LBMOnly(par):
    return (1.+np.dot(X,par[-n:]))

  def ModelOnly(par):
    return func(par[:-n],x)
    
  Model.X = X #attach X to the model
  Model.sys = Model_LBMOnly #attach systematics only model to function
  Model.model = ModelOnly
  
  return Model

# def OnlyLBM(par,X):
#   "Returns only the equivalent LBM from AddLBM function"
#   
#   #just return the function if X is None, or shape is None (ie empty array)
#   if X is None: return 1
#   else: assert X.ndim == 2, "X should be 2D"
#   n = X.shape[1] #get number of parameters required for LBM
#   if n == 0: return 1
# 
#   def Model(par): #define model with standard syntax
#     return (1.+np.dot(X,par[-n:]))
#   
#   return Model

def NormX(X):
  """
  Normalise a basis matrix X, ignoring a column of ones
  X - 2D array [N x K] where N is number of data points and K is number of basis models
  """

  X_copy = np.copy(X)
  for j in range(X.shape[1]):
    if not np.allclose(X_copy[:,j],1.): #ignore if array of ones is included
      X_copy[:,j] -= X_copy[:,j].mean()
      X_copy[:,j] /= X_copy[:,j].std()
  
  return X_copy
