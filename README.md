
Neale Gibson
n.gibson@tcd.ie
nealegibby@gmail.com

PlanetFuncs python module - collection of python functions for exoplanet calculations,
written in C for fast implementation where necessary.

This code is a streamlined and updated version of MyFuncs containing only planet funcs,
and removed inference functions which are now in inferno. I have also simplified
installation from depending on the libhyp hack into scipy code - now all the relevant
code is contained in the TransitFlux submodule.

While there is support for computing RV curves, eccentric orbits etc, I don't use these
regularly and they are neither well tested nor optimised.

***

### INSTALLATION

This module first needs the GSL scientific C library installed.

The simple version is to install via conda
`$ conda install gsl`

Then the module can be installed as normal
`$ python setup.py build`
`$ python setup.py install`
or
`$ python setup.py install --prefix=<INSTALL_DIR>`

***

### INSTALLATION of GSL from source

Alternatively you can install GSL from source, and make sure it's available on your CPATH
and LIBRARY_PATH. It's quite easy, in principle could speed up the code if it's better
tuned to your machine.

You can either download the latest version of GSL from http://www.gnu.org/software/gsl/
or use one of the versions included in this distribution.

It is easy to install, just download, unpack, enter directory:
$ tar -xzvf gsl-1.16.tar.gz
$ cd gsl-1.16

Then run configure
$ ./configure
or first configure to a different directory? - see gsl install instructions
$ ./configure --prefix=<INSTALL_DIR>
or for your home directory
$ ./configure --prefix=$HOME

Then need to build and install with make
$ make
$ make install

Finally, to complete the installation your C compilers need to know where to look for the
headers and code. Depending on the path, they may already be automatically found. If in
doubt, you'll need to add to your CPATH and LIBRARY_PATH:
e.g. in bash:
$ export CPATH="${HOME}/include:$CPATH"
$ export LIBRARY_PATH="${HOME}/lib:$LIBRARY_PATH"
or for a non-standard location:
$ export CPATH="<INSTALL_DIR>/include:$CPATH"
$ export LIBRARY_PATH="<INSTALL_DIR>/lib:$LIBRARY_PATH"

Then you can run python setup.py install as above.

***

### Simple example

```
from PlanetFuncs import transit_quad, transit_nonlin

t = np.linspace(-0.1,0.1,1000)

#define parameter vector [T0,P,a/Rs,Rp/Rs,b,limb_darkening...,f_oot,T_grad]
p_quad = [0,3,10,0.1,0.,0.1,0.1,1.,0]
p_nonlin = [0,3,10,0.1,0.,0.1,0.1,0.1,0.1,1.,0]

flux_quad = transit_quad(p_quad,t)
flux_nonlin = transit_nonlin(p_nonlin,t)

ax = plt.subplot(xlabel='time',ylabel='rel flux')
ax.plot(t,flux_quad,'r-')
ax.plot(t,flux_nonlin,'k-')

```

***

Some old README stuff from MyFuncs that I'll leave here...

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
