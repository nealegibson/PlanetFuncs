
from setuptools import setup, Extension
import os
import numpy as np
import sys


extensions = [
    Extension("PlanetFuncs.TransitFlux.TransitFlux",sources=["src/TransitFlux/TransitFlux.c","src/TransitFlux/flux_quad.c","src/TransitFlux/flux_nonlin.c","src/TransitFlux/flux_func.c",
      "src/TransitFlux/hyp2f1.c","src/TransitFlux/gamma.c","src/TransitFlux/psi.c","src/TransitFlux/round.c","src/TransitFlux/const.c","src/TransitFlux/fabs.c","src/TransitFlux/polevl.c","src/TransitFlux/mtherr.c"],libraries=['gsl',],
      include_dirs=[np.get_include(),'src/TransitFlux/']
      ),
    Extension("PlanetFuncs.PlanetOrbit.PlanetOrbit",sources=["src/PlanetOrbit/PlanetOrbit.c","src/PlanetOrbit/PlanetOrbit_functions.c"]),
    Extension("PlanetFuncs.RadialVelocity.RadialVelocity",sources=["src/RadialVelocity/RadialVelocity.c"]),
    ]

if sys.platform != 'win32':
  extensions.append(
      Extension("PlanetFuncs.TransitFlux.TransitFlux_ctypes",sources=["src/TransitFlux/TransitFlux_ctypes.c","src/TransitFlux/flux_quad.c","src/TransitFlux/flux_nonlin.c","src/TransitFlux/flux_func.c",
        "src/TransitFlux/hyp2f1.c","src/TransitFlux/gamma.c","src/TransitFlux/psi.c","src/TransitFlux/round.c","src/TransitFlux/const.c","src/TransitFlux/fabs.c","src/TransitFlux/polevl.c","src/TransitFlux/mtherr.c"],libraries=['gsl','cblas','m'],)
      )
    
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
  ext_modules = extensions,
  include_dirs=[np.get_include(),],      
  )

