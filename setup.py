
from setuptools import setup, Extension
import os
import numpy as np

setup(
  #metadata
  name = "PlanetFuncs", version = "1.0",
  description='Module containing useful planet functions taken from PlanetFuncs, plus C python extension to calculate fast Mandel & Agol 2002 transit light curves', 
  author='Neale Gibson',
  author_email='n.gibson@tcd.ie',
  
  #define pure python package(s)
  packages=['PlanetFuncs','PlanetFuncs.TransitFlux','PlanetFuncs.PlanetOrbit','PlanetFuncs.RadialVelocity'],
  package_dir={'PlanetFuncs':'src','PlanetFuncs.TransitFlux':'src/TransitFlux','PlanetFuncs.PlanetOrbit':'src/PlanetOrbit','PlanetFuncs.RadialVelocity':'src/RadialVelocity'},
  
  #define extension module(s)
  ext_modules = [
    Extension("PlanetFuncs.TransitFlux.TransitFlux",sources=["src/TransitFlux/TransitFlux.c","src/TransitFlux/flux_quad.c","src/TransitFlux/flux_nonlin.c","src/TransitFlux/flux_func.c",
      "src/TransitFlux/hyp2f1.c","src/TransitFlux/gamma.c","src/TransitFlux/psi.c","src/TransitFlux/round.c","src/TransitFlux/const.c","src/TransitFlux/fabs.c","src/TransitFlux/polevl.c","src/TransitFlux/mtherr.c"],libraries=['gsl','cblas','m'],
#      library_dirs=['./LibHyp/lib/'], #libhyp.a must be compiled first, and in a directory ./LibHyp/lib (obviously this can be changed if necessary)
#      include_dirs=[np.get_include(),'./LibHyp/include/']
      ),
    Extension("PlanetFuncs.PlanetOrbit.PlanetOrbit",sources=["src/PlanetOrbit/PlanetOrbit.c","src/PlanetOrbit/PlanetOrbit_functions.c"]),
    Extension("PlanetFuncs.RadialVelocity.RadialVelocity",sources=["src/RadialVelocity/RadialVelocity.c"]),
#     Extension("PlanetFuncs.TransitFlux.TransitFlux_ctypes",sources=["src/TransitFlux/TransitFlux_ctypes.c","src/TransitFlux/flux_quad.c","src/TransitFlux/flux_nonlin.c","src/TransitFlux/flux_func.c",
#       "src/TransitFlux/hyp2f1.c","src/TransitFlux/gamma.c","src/TransitFlux/psi.c","src/TransitFlux/round.c","src/TransitFlux/const.c","src/TransitFlux/fabs.c","src/TransitFlux/polevl.c","src/TransitFlux/mtherr.c"],libraries=['gsl','cblas','m'],
#       ),
    ],
  include_dirs=[np.get_include(),],      
  )
