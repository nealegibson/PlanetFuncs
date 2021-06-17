
#define PY_SSIZE_T_CLEAN //this redefines length of input args I think - recommended to call before Python.h
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION //this switches off the ability to use the old API
#include <Python.h>
#include <numpy/arrayobject.h> //for PyArray_Type objects
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//macros
#define SQ(x) ((x)*(x))

//main function
double Rad_vel(double time, double e, double p, double w, double K);
double ecc_anom(double e, double M);
