
/********************************************************************************/

#include <Python.h>
#include <numpy/arrayobject.h> //for PyArray_Type objects
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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