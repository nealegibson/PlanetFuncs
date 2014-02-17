
OBJ = hyp2f1.o gamma.o psi.o round.o const.o fabs.o polevl.o mtherr.o

#INSTALL_DIR = .

#INCLUDE = -I/home/GibsonN/Library/Python/2.6/site-packages/numpy-1.5.1-py2.6-macosx-10.6-universal.egg/numpy/core/include/:/System/Library/Frameworks/Python.framework/Versions/2.6/include/python2.6/

#CFLAGS = -g -O2 -Wall $(INCLUDE) -fno-builtin
CFLAGS = -g -O2 -Wall -fno-builtin -fpic
#CC = gcc -arch i386 -arch x86_64
#CC = gcc #this needs to be removed to properly install on by linux box (where CC = gcc -m32 to force 32 bit)

#create only static library
all: libhyp.a

libhyp.a: $(OBJ)
	ar rv libhyp.a $(OBJ)
#	$(CC) -shared -o libhyp.so $(OBJ)

#install: libhyp.a
#	cp libhyp.a $(INSTALL_DIR)

clean:
	@rm -f *.o
	@rm -f libhyp.a
#	rm -f libhyp.so
