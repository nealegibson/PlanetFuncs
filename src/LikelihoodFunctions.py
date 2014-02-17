"""
Docstring...
"""

import numpy as np
#from GaussianProcesses import CovarianceMatrix
#import MyMatrix
import scipy.linalg as LA
import time

###############################################################################
#'main' log likelihood functions here
def LogLikelihood_iid(r,beta,sig=1.):
  """
  logP function from residuals and noise vector (sig)
  - when sig = 1. it does not contribute to the noise budget,
  hence beta is then the global noise parameter, otherwise a rescale
  """
  
  N = r.size
  
  logP = - (1./2.) * ( r**2 / (beta*sig)**2 ).sum() - np.log(sig).sum() - N*np.log(beta) - N/2.*np.log(2*np.pi)
  
  return logP

def LogLikelihood_iid_mf(p,func,x,y,sig=1.):
  """
  logP function which takes a mean funciton and params to create residuals
  last param must be beta noise (rescale) parameter, rest must be mean_func parameters
  """
  
  r = y - func(p[:-1],x)
  
  return LogLikelihood_iid(r,p[-1],sig=sig)

def LogLikelihood_invK(r,invK,logdetK):
  """
  Generalised log likelihood function from inverse convariance matrix, and log
  determinatnt. Uses fixed noise characteristics for type-II ML.
  """
  
  invK = np.matrix(invK)

#  r = np.matrix(np.array(r)).T

#  print r.shape
#  print invK.shape, type(invK)
#  print logdetK
#  print r.shape, invK.shape, logdetK.shape
  logP = -0.5 * r.T * invK * r - 0.5 * logdetK - (r.size/2.) * np.log(2*np.pi)
#  print logP

  return np.float(logP)

def LogLikelihood_invK_mf(p,func,x,y,invK,logdetK):
  """
  logP function which takes a mean funciton and params to create residuals
  calls the generalised log likelihood function, LogLikelihood_invK
  """
#  print p, func, x.size, y.size, invK.shape, logdetK
#  print x
  r = y - func(p,x)
  r = np.matrix(np.array(r).flatten()).T #ensure r is an (n x 1) matrix

  return LogLikelihood_invK(r,invK,logdetK)

def LogLikelihood_times_invK_mf(p,func,x,y,invK,logdetK):
  """
  logP function which takes a mean funciton and params to create residuals
  calls the generalised log likelihood function, LogLikelihood_invK
  
  Analogous to the invK_mf_times_kf with a fixed covariance matrix
  **not yet tested**
  
  """
  
#  print p, func, x.size, y.size, invK.shape, logdetK
#  print x
  r = y / func(p,x) - 1.
  r = np.matrix(np.array(r).flatten()).T #ensure r is an (n x 1) matrix
  
  return LogLikelihood_invK(r,invK,logdetK)

def LogLikelihood_invK_mf_and_kf(p,t,kf,kf_args,mf,mf_args,n_hp):
#   """
#   p - parameters
#   t - target values
#   kf - kernel function
#   kf_args - arguments to kernel function
#   mf - mean function
#   mf_args - arguments to mean function
#   """
#
  
  hpar = p[:n_hp]
  mf_par = p[n_hp:]
  
#   print "testing:"
#   print "n", n_hp
#   print "t",t.size
#   print "kf_args",kf_args.shape
#   print "mf_args",mf_args.shape
#   print "p", p
#   print "hp:",hpar
#   print "mfp:",mf_par
#   print "kf:",kf
#   print "mf:", mf
  
  
#   if MeanFunction == None:
#     r = f[:] #r is just the target values then
#     hp = p[:] #get hyperparams
#   
#   if MeanFunction != None: #split parameters into hyperpars and mf pars
#     n_hp = (np.array(HPFixed)==0).sum() #get number of (variable) h-parameters
#     mfp = p[n_hp:] #get function parameters
#     hp = p[:n_hp] #get hyperparameters
    
    #subtract mean function from target values to calculate likelihood function
#    r = f[:] - MeanFunction(mfp,MF_args,fixed=MFFixed,fixed_par=MFFixed_params)

  r = t - mf(mf_par,mf_args)
#  r = t - CallFunc_FixedPar(MeanFunction,mfp,mf_args,fixed=MFFixed,fixed_par=MFFixed_params)
  
  #ensure r is an (n x 1) column vector
  r = np.matrix(np.array(r).flatten()).T

  #create the covariance matrix
#  K = CovarianceMatrix(hp,X,fixed_par=HPFixed_params,fixed=HPFixed,KernelFunction=KernelFunction)
#  K = CallFunc_FixedPar(CovarianceMatrix,hp,X,fixed=HPFixed,fixed_par=HPFixed_params,KernelFunction=KernelFunction)
  
  K = CovarianceMatrix(hpar,kf_args,KernelFunction=kf)
  
  #get log det and invert the covariance matrix (these need optimised!)
  sign,logdetK = np.linalg.slogdet( K )

  #'normal' scipy matrix inversion
#   invK = np.linalg.inv( K )
#   return LogLikelihood_invK(r,invK,logdetK)

  #use cholesky factorisation to get inverse?
#  print K
#  invK = LA.cho_solve(LA.cho_factor(K),np.eye(K[0].size))
#  return LogLikelihood_invK(r,invK,logdetK)

  #use lu factorisation to get inverse
#  print K
#  invK = LA.lu_solve(LA.lu_factor(K),np.eye(K[0].size))
#  return LogLikelihood_invK(r,invK,logdetK)
  
  #solve for the inverse inside the logP equation (quicker?)
#  invK_x_R = np.mat(LA.lu_solve(LA.lu_factor(K),r))
#  logP = -0.5 * r.T * np.mat(LA.lu_solve(LA.lu_factor(K,overwrite_a=1),r)) - 0.5 * logdetK - (r.size/2.) * np.log(2*np.pi)
  logP = -0.5 * r.T * np.mat(LA.cho_solve(LA.cho_factor(K,overwrite_a=1),r)) - 0.5 * logdetK - (r.size/2.) * np.log(2*np.pi)
  return np.float(logP)
  
  #invert matrix via my cholesky-gsl code - very slow - probably memory issues...
#  invK = np.copy(K) #first copy K?
#  MyMatrix.SolveChol(K) #invert matrix *in place*
#  return LogLikelihood_invK(r,K,logdetK)  

###############################################################################
def LogLikelihood_invK_mf_plus_GP(p,t,kf,kf_args,mf,mf_args,n_hp):
  """
  Same as LogLikelihood_invK_mf_and_kf (with a more informative name cf times)
  """
  
  #create hyperparameter and mean function vectors
  hpar = p[:n_hp]
  mf_par = p[n_hp:]
    
  #residuals from mean function
  #start = time.clock()
  r = t - mf(mf_par,mf_args)
  #print 'model time:', time.clock() - start,'s'
  
  #ensure r is an (n x 1) column vector
  r = np.matrix(np.array(r).flatten()).T

  #create the covariance matrix
  #start = time.clock()
  K = CovarianceMatrix(hpar,kf_args,KernelFunction=kf)
  #print 'K time:', time.clock() - start,'s'
  
  #get log det and invert the covariance matrix (these need optimised!)
#  sign,logdetK = np.linalg.slogdet( K )
  
#  logP = -0.5 * r.T * np.mat(LA.lu_solve(LA.lu_factor(K),r)) - 0.5 * logdetK - (r.size/2.) * np.log(2*np.pi)
#  logP = -0.5 * r.T * np.mat(LA.lu_solve(LA.lu_factor(K,overwrite_a=1),r)) - 0.5 * logdetK - (r.size/2.) * np.log(2*np.pi)
#  logP = -0.5 * r.T * np.mat(LA.cho_solve(LA.cho_factor(K,overwrite_a=1),r)) - 0.5 * logdetK - (r.size/2.) * np.log(2*np.pi)

  #alternative approach - is this faster?
  #start = time.clock()
  CI = LA.cho_factor(K,overwrite_a=1)
  logdetK = (2*np.log(np.diag(CI[0])).sum()) #compute log determinant
  logP = -0.5 * np.dot(r.T,LA.cho_solve(CI,r)) - 0.5 * logdetK - (r.size/2.) * np.log(2*np.pi)
    #logP = -0.5 * r.T * np.mat(LA.cho_solve(CI,r)) - 0.5 * logdetK - (r.size/2.) * np.log(2*np.pi)
  #ch_time = (time.clock() - start)
#   print 'solve time CHOL:', ch_time,'s', logP

  #start = time.clock()
#   LU = LA.lu_factor(K,overwrite_a=1)
#   logdetK = (np.log(np.diag(LU[0])).sum()) #compute log determinant
#   logP = -0.5 * np.dot(r.T,LA.lu_solve(LU,r)) - 0.5 * logdetK - (r.size/2.) * np.log(2*np.pi)  
#   lu_time = (time.clock() - start)
#   print 'solve time LU:', lu_time,'s', logP
#   print "DIFF:", lu_time - ch_time
#  print
  
#  return np.float64(logP)
  return logP
###############################################################################

def LogLikelihood_invK_mf_times_GP(p,t,kf,kf_args,mf,mf_args,n_hp):
  """
  Similar to LogLikelihood_invK_mf_and_kf, although the data is *divided* by
  the mean function before fitting a GP(1,C)
  """
  
  #create hyperparameter and mean function vectors
  hpar = p[:n_hp]
  mf_par = p[n_hp:]
    
  #residuals from *division* by mean function (and subtract 1)
  r = t / mf(mf_par,mf_args) - 1.
  
  #ensure r is an (n x 1) column vector
  r = np.matrix(np.array(r).flatten()).T

  #create the covariance matrix
  K = CovarianceMatrix(hpar,kf_args,KernelFunction=kf)
  
  #get log det and invert the covariance matrix (these need optimised!)
#   sign,logdetK = np.linalg.slogdet( K )
#   
#   logP = -0.5 * r.T * np.mat(LA.lu_solve(LA.lu_factor(K),r)) - 0.5 * logdetK - (r.size/2.) * np.log(2*np.pi)

  CI = LA.cho_factor(K,overwrite_a=1)
  logdetK = (2*np.log(np.diag(CI[0])).sum()) #compute log determinant
  logP = -0.5 * np.dot(r.T,LA.cho_solve(CI,r)) - 0.5 * logdetK - (r.size/2.) * np.log(2*np.pi)

  return np.float64(logP)
  
###############################################################################
#neg chi squared / 2 functions for input to generic MCMC routine
#take exactly the same arguments as log likelihood but removes the constant terms
#from the ratio
def NegChi2_2_iid(r,beta,sig=1.):
  """
  Neg chi2 *over 2* function from residuals and noise vector (sig).
  Takes exactly the same arguments as log likelihood functions but terms which
  cancel out in the likelihood ratio (assuming fixed noise) are removed.
  
  r - residuals from model fit
  beta - noise when sig = 1., or rescale factor when sig array is provided
  sig - array of noise (stddev) - when = 1. does not contribute to the noise budget,
    hence beta is then the global noise parameter.
  """

  c2_2 = - (1./2.) * ( r**2 / (beta*sig)**2 ).sum()

  return c2_2

def NegChi2_2_iid_mf(p,func,x,y,sig=1.):
  """
  Neg chi2 *over 2* function from residuals and noise vector (sig).
  First subtracts the mean funciton from the data to get the residuals, then
  passes the residuals to NegChi2_2_iid functon
  """

  r = y - func(p[:-1],x)

  return NegChi2_2_iid(r,p[-1],sig=sig)

###############################################################################
#copied over from GaussianProcess module
def CovarianceMatrix(theta,X,KernelFunction=None,dtype=np.float64):
  """
  X - input matrix (n x D)
  theta - hyperparameter array/list
  K - (n x n) covariance matrix
  """
  
  #allow fixed parameters and variable params to be passed
#   if fixed_par==None:
#     theta = par[:]
#   else:
#     fixed = np.array(fixed) #ensure fixed is a np array
#     theta = np.empty(fixed.size) #create empty pars array
#     #assign parameters to normal param vector
#     theta[np.where(fixed==True)] = fixed_par
#     theta[np.where(fixed!=True)] = par
    
  K = KernelFunction(X,X,theta,white_noise=True)
  
  return np.matrix(K)

###############################################################################
