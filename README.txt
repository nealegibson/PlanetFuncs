
Neale Gibson
ngibson@eso.org
nealegibby@gmail.com

MyFuncs python module - collection of python functions, most usefully including functions to compute exoplanet transit light curves (quadratic and nonlinear), and some likelihood functions useful for MCMC and GaussianProcesses (some also in the Infer module). This has been designed for easy use alongside the Infer module. It contains install instructions and many examples, please contact me if you have any problems or suggestions.

Changes prior to first github commit:

Version 1.4
17/02/2014 - Neale Gibson
Minor edits and uploaded to github for version control

Version 1.3dev
03/07/2012 - Neale Gibson
Added the optimiser function. This is a wrapper for the fmin scipy functions, using a similar syntax to MCMC and allowing parameters to be fixed during the
fit.

Version 1.3
01/05/2012 - Neale Gibson

Added support for wavelet likelihood functions adapted from WaveletMCMC module. This incorporates the wavelet based likelihood function of Carter & Winn 2009

Version 1.2
01/05/2012 - Neale Gibson

Added simple radial velocity functions (from TransitFit).

Version 1.1
26/04/2012 - Neale Gibson

Added eccentric light curve functions (in EccLightCurveModels.py). These use PlanetOrbit to calculate the 3d coordinates using a C extension. A bug in PlanetOrbit was fixed to stop memory leaks (caused by using PyBuild_Value to return python objects, which increases the ref count and should be decremented first).

Version 1.0
23/02/2012 - Neale Gibson

1.0 - first 'release' version. This is mainly exoplanet transit light curve functions and likelihood functions used for MyMCMC and GPs. The C extensions are codes to compute Mandel & Agol (2002) light curves, both quadratic and non-linear. These use the gsl scientific library. LibHyp is basically a hack into the scipy 0.8 source code, which calculates the Gauss hypergeometric function reasonably quickly compared to most other codes I tried (required for the non-linear limb darkening model). The appell hypergeometric function is slow for some input parameters, so more analytic continuations probably need to be added to it or to the appell hypergeometric function in flux_func.c. To install with just the quadratic law ignore LibHyp and comment out the import of FluxNonlin and any functions that require it in LightCurveModels.py.