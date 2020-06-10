
#set the CPATH to include python include dir and numpy include dir
PYCPATH:=$(shell python -c 'from distutils.sysconfig import get_python_inc; print(get_python_inc())')
NUMCPATH:=$(shell python -c 'import numpy as np; print(np.get_include())')
CPATH:=$(PYCPATH):$(NUMCPATH):$(CPATH)

build:
	@echo $(CPATH)
	make -C src

install:
	@mkdir -p lib
	@mkdir -p include
	mv src/libhyp.a lib/.
	cp src/*.h include/.

clean:
	make -C src clean
	@rm -vf libhyp.a
	@rm -rvf lib include
	
