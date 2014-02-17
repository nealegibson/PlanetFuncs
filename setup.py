from numpy.distutils.core import setup, Extension

setup(
  #metadata
  name = "MyFunctions", version = "1.3dev",
  description='Module containing useful functions, plus C python extension to calculate fast Mandel & Agol 2002 transit light curves', 
  author='Neale Gibson',
  author_email='Neale.Gibson@astro.ox.ac.uk',

  #define pure python package(s)
  packages=['MyFunctions','MyFunctions.TransitFlux','MyFunctions.PlanetOrbit','MyFunctions.RadialVelocity','MyFunctions.WaveletLikelihood','MyFunctions.WaveletTransform'],
  package_dir={'MyFunctions':'src','MyFunctions.TransitFlux':'src/TransitFlux','MyFunctions.PlanetOrbit':'src/PlanetOrbit','MyFunctions.RadialVelocity':'src/RadialVelocity','MyFunctions.WaveletLikelihood':'src/WaveletLikelihood','MyFunctions.WaveletTransform':'src/WaveletTransform'},
  
  #define extension module(s)
  ext_modules = [
    Extension("MyFunctions.TransitFlux.TransitFlux",sources=["src/TransitFlux/TransitFlux.c","src/TransitFlux/flux_quad.c","src/TransitFlux/flux_nonlin.c","src/TransitFlux/flux_func.c"],libraries=['gsl','gslcblas','hyp'],
      library_dirs=['./LibHyp/lib/'], #libhyp.a must be compiled first, and in a directory ./LibHyp/lib (obviously this can be changed if necessary)
      include_dirs=['./LibHyp/include/']),
    Extension("MyFunctions.PlanetOrbit.PlanetOrbit",sources=["src/PlanetOrbit/PlanetOrbit.c","src/PlanetOrbit/PlanetOrbit_functions.c"]),
    Extension("MyFunctions.RadialVelocity.RadialVelocity",sources=["src/RadialVelocity/RadialVelocity.c"]),
    Extension("MyFunctions.WaveletLikelihood.WaveletLikelihood",sources=["src/WaveletLikelihood/WaveletLikelihood.c","src/WaveletLikelihood/FWT.c"],libraries=['gsl','gslcblas']),
    Extension("MyFunctions.WaveletTransform.WaveletTransform",sources=["src/WaveletTransform/WaveletTransform.c","src/WaveletTransform/FWT.c"],libraries=['gsl','gslcblas'])
    ]
    
  )
