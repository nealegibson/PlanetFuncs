import numpy as np

#wrapper to do maximum likelihood with fmin/scipy routines

def NelderMead(ErrFunc, params0, function_args, maxiter=10000, maxfun=10000, verbose=True):

  if verbose:
    print ("Running Nelder-Mead simplex algorithm... ")
    t0 = time.clock()
  params = fmin(ErrFunc, params0, args=function_args, maxiter=maxiter, maxfun=maxfun)
  if verbose:
    print ("(Time: %f secs)" % (time.clock()-t0))
    print ("Maximum likelihood hyperparameters: ", params,"\n")
  return params

def NegLogLikelihood(*args,**kwargs):
  """
  negative log likelihood function from mean functions and hyperparameters
  needed for scipy fmin functions to minimise.
  The 'fixed' parameter allows parameters to be constant - bit messy - maybe a better way...

  """
  
  return - LogLikelihoodGP(*args,**kwargs)

def LogLikelihoodGP(p,X,f,KernelFunction,HPFixed=None,HPFixed_params=None,MeanFunction=None,MFFixed=None,MFFixed_params=None,MF_args=None):
  """
  Generalised log likelihood function for the GP model. Easy to add priors etc, but
  best done by writing your own likelihood function...
  
  """
  #print p,X,f,KernelFunction,HPFixed,HPFixed_params,MeanFunction,MFFixed,MFFixed_params,MF_args 

  if MeanFunction == None:
    r = f[:] #r is just the target values then
    hp = p[:] #get hyperparams
  
  if MeanFunction != None: #split parameters into hyperpars and mf pars
    n_hp = (np.array(HPFixed)==0).sum() #get number of (variable) h-parameters
    mfp = p[n_hp:] #get function parameters
    hp = p[:n_hp] #get hyperparameters
    
    #subtract mean function from target values to calculate likelihood function
#    r = f[:] - MeanFunction(mfp,MF_args,fixed=MFFixed,fixed_par=MFFixed_params)
    r = f[:] - CallFunc_FixedPar(MeanFunction,mfp,MF_args,fixed=MFFixed,fixed_par=MFFixed_params)
  
  #ensure r is an (n x 1) column vector
  r = np.matrix(np.array(r).flatten()).T

  #create the covariance matrix
#  K = CovarianceMatrix(hp,X,fixed_par=HPFixed_params,fixed=HPFixed,KernelFunction=KernelFunction)
  K = CallFunc_FixedPar(CovarianceMatrix,hp,X,fixed=HPFixed,fixed_par=HPFixed_params,KernelFunction=KernelFunction)
  
  #get log det and invert the covariance matrix (these need optimised!)
  sign,logdetK = np.linalg.slogdet( K ) #sign shouldn't matter for positive semi-def matrix! (ie is always positive)
  
  #Get precision matrix
  #PrecMatrix = np.linalg.inv( K ) #old method

  #get prec matrix via lu decomposition
#  PrecMatrix = LA.lu_solve(LA.lu_factor(K),np.eye(K[0].size))

  #get prec matrix via cholesky decomposition
#  PrecMatrix = LA.cho_solve(LA.cho_factor(K),np.eye(K[0].size))
  
#  np.mat(LA.lu_solve(LA.lu_factor(K),r))
  
  #calculate the log likelihood function
#  logP = -0.5 * r.T * PrecMatrix * r - 0.5 * logdetK - (r.size/2.) * np.log(2*np.pi)
  #solve Kx = r (via lu decomposition) to get x = (K^-1)*r - much faster than an inverse
  logP = -0.5 * r.T *  np.mat(LA.lu_solve(LA.lu_factor(K),r)) - 0.5 * logdetK - (r.size/2.) * np.log(2*np.pi)

  return np.array(logP).flatten()


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

