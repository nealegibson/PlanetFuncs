
/********************************************************************************/

#define PY_SSIZE_T_CLEAN //this redefines length of input args I think - recommended to call before Python.h
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION //this switches off the ability to use the old API
#include <Python.h>
#include <numpy/arrayobject.h> //for PyArray_Type objects
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

/*********************************************************************************
Constants
*********************************************************************************/
#define AU 1.4959787e11 //m
#define G 6.6742e-11 //m3kg-1s-2
#define MSUN 1.989e30 //kg
#define MJUP 1.9e27 //kg
#define RSUN 6.95e8 //m
#define RJUP 7.1492e7 //m

/*********************************************************************************
Macros
*********************************************************************************/
#define SQ(x) ((x)*(x))

/*********************************************************************************
function prototypes
*********************************************************************************/

static PyObject * get_x_py(PyObject *self, PyObject *args, PyObject *keywds);
static PyObject * get_y_py(PyObject *self, PyObject *args, PyObject *keywds);
static PyObject * get_z_py(PyObject *self, PyObject *args, PyObject *keywds);
static PyObject * get_norm_py(PyObject *self, PyObject *args, PyObject *keywds);

/*********************************************************************************
docstrings
*********************************************************************************/

/*
#define MCMCFUNC_DOCSTR \
"Function docstring...\n\
...\n"

#define MODULE_DOCSTR \
"Module docstring...\n\
...\n"
*/