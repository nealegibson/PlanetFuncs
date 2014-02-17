
import numpy as np

def CallFunc_FixedPar(func,var_par,func_args,fixed=None,fixed_par=None,**kwargs):
  """
  Wrapper to allow functions with multiple parameters to have fixed parameters
  when using fmin-type optimisation functions.
  
  func - function to be called
  var_par - value of variable parameters
  func_args - arguments to be passed to the function
  fixed_par = array of 1s or 0s, 1 = fixed parameter - must equal total no of parameters
  fixed_par_val - array of fixed parameter values
  
  var_par and fixed_par_val combined *must* be equal to the no of function parameters,
    and size of fixed_par.
  """
  
  #if no fixed parameters passed - just assign var_par to par and call function
  if fixed==None:
    par = var_par[:]
  #otherwise construct the parameter vector from var_par and fixed_par_val
  else:
    fixed = np.array(fixed) #ensure fixed is a np array
    par = np.empty(fixed.size) #create empty pars array
    #assign parameters to normal param vector
    par[np.where(fixed==True)] = fixed_par
    par[np.where(fixed!=True)] = var_par
  
  #now call the funciton as normal
  return func(par,func_args,**kwargs)

