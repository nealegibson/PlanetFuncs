
Neale Gibson (n.gibson@tcd.ie, nealegibby@gmail.com)

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

```
$ conda install gsl
```

Then the module can be installed as normal
```
$ python setup.py build
$ python setup.py install
```
or
```
$ python setup.py install --prefix=<INSTALL_DIR>
```

***

### INSTALLATION on Windows 11

This was a little painful...


### INSTALLATION of GSL from source

Alternatively you can install GSL from source, and make sure it's available on your CPATH
and LIBRARY_PATH. It's quite easy, in principle could speed up the code if it's better
tuned to your machine.

You can either download the latest version of GSL from http://www.gnu.org/software/gsl/
or use one of the versions included in this distribution.

It is easy to install, just download, unpack, enter directory:
```
$ tar -xzvf gsl-1.16.tar.gz
$ cd gsl-1.16
```
Then run configure
```
$ ./configure
#or first configure to a different directory? - see gsl install instructions
$ ./configure --prefix=<INSTALL_DIR>
#or for your home directory
$ ./configure --prefix=$HOME
```

Then need to build and install with make
```
$ make
$ make install
```

Finally, to complete the installation your C compilers need to know where to look for the
headers and code. Depending on the path, they may already be automatically found. If in
doubt, you'll need to add to your CPATH and LIBRARY_PATH:

```
#e.g. in bash:
$ export CPATH="${HOME}/include:$CPATH"
$ export LIBRARY_PATH="${HOME}/lib:$LIBRARY_PATH"
#or for a non-standard location:
$ export CPATH="<INSTALL_DIR>/include:$CPATH"
$ export LIBRARY_PATH="<INSTALL_DIR>/lib:$LIBRARY_PATH"
```

Then you can run python setup.py install as above.

***

### Simple example

```
from PlanetFuncs import transit_quad, transit_nonlin
import matplotlib.pyplot as plt

#define time array (units just need to be consistent with T0 and P)
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
