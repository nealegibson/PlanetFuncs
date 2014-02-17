
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
