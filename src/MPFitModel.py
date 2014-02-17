
import numpy as np
import mpfit

###############################################################################

def ErrFunc(p, x = None, y = None, err = None, fjac = None, func = None):
  """Light curve error function to minimize for fitting.""" 
  
#   print func
#   print err
#   print x.size, y.size, p
  
  if func == None:
    print "func must be provided!"
    return
  if err == None: return [0, func(p, x) - y, None]
  return [0, (func(p, x) - y) / err, None]

###############################################################################

def MPFitFunc(x,y,func,guess_params,fixed=None,limited=None,limits=None,\
  err=None,return_res=True,plot=False,save=False):
  """
  Wrapper to fit transit lightcurves using mpfit.
  (don't trust errors yet - maybe bug in mpfit?)
  
  time - numpy time array
  flux - numpy flux array
  guess_params - list of guess parameters
  fixed - numpy array the same length as guess_params
    set to 0 to vary and 1 to fix for each corresponding parameter
  return_res - set to true to return the residuals
  """
  
  #test
#   import pylab
#   pylab.figure(13)
#   pylab.clf()
#   pylab.plot(x,y,'.')
#   pylab.plot(x,func(guess_params,x),'-')
  
  #call mpfit to fit light curve
  fa = {'x': x, 'y': y, 'err': err, 'func': func}
  input_params = guess_params[:]
  m = mpfit.mpfit(ErrFunc,input_params,functkw = fa,fixed=fixed,\
    quiet=True,limited=limited,limits=limits)
  print "test: mpfit done!"
  
  #check for errors
  if (m.status <= 0): print 'error message = ', m.errmsg
  fit_params = np.copy(m.params) #get parameters
  err_params = np.copy(m.pcerror) #get errors
  if return_res:
    residuals = y - func(fit_params,x) #make residuals

  #test
#   pylab.plot(x,func(fit_params,x),'-')
#   raw_input()
  
  #return parameters (and residuals if return_res==True)
  if return_res:
    return fit_params[:],err_params,residuals
  else:
    return fit_params[:],err_params

###############################################################################
