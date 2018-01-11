
from numpy.distutils.core import setup, Extension
import os

#first need to build a C library 
os.system('make -f build_libhyp')

setup(
  #metadata
  name = "MyFuncs", version = "1.4",
  description='Module containing useful functions, plus C python extension to calculate fast Mandel & Agol 2002 transit light curves', 
  author='Neale Gibson',
  author_email='nealegibby@gmail.com',
  
  #define pure python package(s)
  packages=['MyFuncs','MyFuncs.TransitFlux','MyFuncs.PlanetOrbit','MyFuncs.RadialVelocity','MyFuncs.WaveletLikelihood','MyFuncs.WaveletTransform'],
  package_dir={'MyFuncs':'src','MyFuncs.TransitFlux':'src/TransitFlux','MyFuncs.PlanetOrbit':'src/PlanetOrbit','MyFuncs.RadialVelocity':'src/RadialVelocity','MyFuncs.WaveletLikelihood':'src/WaveletLikelihood','MyFuncs.WaveletTransform':'src/WaveletTransform'},
  
  #define extension module(s)
  ext_modules = [
    Extension("MyFuncs.TransitFlux.TransitFlux",sources=["src/TransitFlux/TransitFlux.c","src/TransitFlux/flux_quad.c","src/TransitFlux/flux_nonlin.c","src/TransitFlux/flux_func.c"],libraries=['gsl','gslcblas','hyp'],
      library_dirs=['./LibHyp/lib/'], #libhyp.a must be compiled first, and in a directory ./LibHyp/lib (obviously this can be changed if necessary)
      include_dirs=['./LibHyp/include/']),
    Extension("MyFuncs.PlanetOrbit.PlanetOrbit",sources=["src/PlanetOrbit/PlanetOrbit.c","src/PlanetOrbit/PlanetOrbit_functions.c"]),
    Extension("MyFuncs.RadialVelocity.RadialVelocity",sources=["src/RadialVelocity/RadialVelocity.c"]),
    Extension("MyFuncs.WaveletLikelihood.WaveletLikelihood",sources=["src/WaveletLikelihood/WaveletLikelihood.c","src/WaveletLikelihood/FWT.c"],libraries=['gsl','gslcblas']),
    Extension("MyFuncs.WaveletTransform.WaveletTransform",sources=["src/WaveletTransform/WaveletTransform.c","src/WaveletTransform/FWT.c"],libraries=['gsl','gslcblas']),
    Extension("MyFuncs.TransitFlux.TransitFlux_ctypes",sources=["src/TransitFlux/TransitFlux_ctypes.c","src/TransitFlux/flux_quad.c","src/TransitFlux/flux_nonlin.c","src/TransitFlux/flux_func.c"],libraries=['gsl','gslcblas','hyp'],
      library_dirs=['./LibHyp/lib/'], #libhyp.a must be compiled first, and in a directory ./LibHyp/lib (obviously this can be changed if necessary)
      include_dirs=['./LibHyp/include/']),
    ]
  )
