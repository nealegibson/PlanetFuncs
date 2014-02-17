/*********************************************************************************
Neale Gibson 15/3/2007 flux_quad.c
This function uses the analytic equations of Mandel and Agol (2000) to compute the
flux of a star during occultation given a normalized separation z, ratio of planet
to star radii, and the two quadratic limb darkening coefficients gam_1 and gam_2.

Python Wrapper Aug 2009. flux_quad_python.c
*********************************************************************************/

#include <Python.h>
#include <numpy/arrayobject.h> //for PyArray_Type objects
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf.h>

/*********************************************************************************
Initialise the module and define docstrings etc
*********************************************************************************/
//python function prototypes
static PyObject * Py_flux_quad_np(PyObject *self, PyObject *args);
static PyObject * Py_flux_nonlin_np(PyObject *self, PyObject *args);
static PyObject * Py_flux_quad_single(PyObject *self, PyObject *args);
static PyObject * scale(PyObject *self, PyObject *args);

//define methods in the module
PyMethodDef TransitFluxMethods[] = {
    {"flux_quad_single", Py_flux_quad_single, METH_VARARGS, "Returns flux as a function of z"},
    {"FluxQuad", Py_flux_quad_np, METH_VARARGS, "Returns flux as a function of z"},
    {"FluxNonlin", Py_flux_nonlin_np, METH_VARARGS, "Returns flux as a function of z"},
    {"scale", scale, METH_VARARGS, "test"},
    {NULL, NULL, 0, NULL}
};

//define the scipy.special.hyp2f1 function as global
//PyObject *pName, *pModule, *pFunc;
//double qqq = 234.987;
//static PyObject * pName = PyString_FromString("scipy.special");
//Py_DECREF(pName);
//static PyObject * pModule = PyImport_Import(pName);
//static PyObject * pFunc = PyObject_GetAttrString(pModule, "hyp2f1");

//initialise the module
PyMODINIT_FUNC 
initTransitFlux(void)
{
    (void) Py_InitModule3("TransitFlux", TransitFluxMethods, "Flux_docstring...");
    import_array(); //must be present for numpy stuff
}
/*********************************************************************************
*********************************************************************************/

//macros
//#define MODE GSL_PREC_DOUBLE
//#define SQ(x) ((x)*(x))


// /*********************************************************************************
// c function prototypes
// *********************************************************************************/

//main function
double flux_quad(double z, double p, double gam_1, double gam_2);
double flux_nonlin(double z, double p, double c_1, double c_2, double c_3, double c_4);

// 
// //functions to compute uniform source (lambda_e)
// double lambda_e_pi_funct(double z, double p);
// double kappa_0(double z, double p);
// double kappa_1(double z, double p);
// 
// //lambda_d functions
// double lambda_1(double z, double p,double a, double b, double q, double k);
// double lambda_2(double z, double p,double a, double b, double q, double k);
// double lambda_3(double p, double k);
// double lambda_4(double p, double k);
// double lambda_5(double p);
// double lambda_6(double p);
// 
// //eta_d functions
// double eta_1(double z, double p, double a, double b);
// double eta_2(double z, double p);
// 
// //further functions
// double mod(double x);

/*********************************************************************************
Python functions
*********************************************************************************/
//test wrapper function to pass and return a numpy array of arbitrary dimensions
//and do a simple scaling calculation
static PyObject * scale(PyObject *self, PyObject *args)
{
	//define array object
	PyObject * input;
	PyArrayObject * array;
	PyArrayObject * output;
	
	//
	double scale;
	
	//read in object
	if (!PyArg_ParseTuple(args, "Od", &input, &scale))
  	  return NULL;
	
	//convert to contiguous array of doubles - will accept up to 10d arrays
	array = (PyArrayObject*) PyArray_ContiguousFromObject(input, PyArray_DOUBLE, 1, 10);
	if(array == NULL) return NULL;
	
	//calculate no of data points in array
	int n,size=1;
	for(n=0;n<array->nd;n++)
		size *= array->dimensions[n];
		
	printf("size = %d\n", size);

	//create output array
	output = (PyArrayObject*) PyArray_SimpleNew(array->nd, array->dimensions, PyArray_DOUBLE);
	
	//create data for output array from input array
	for(n=0;n<size;n++)
		*(((double*)output->data) + n) = *(((double*)array->data) + n) * scale;
	
	//must destroy reference to array
	Py_DECREF(array);
	
	//return array
	return PyArray_Return(output);
}

//wrapper for flux_quad c function - only accepts single value of z
static PyObject * Py_flux_quad_single(PyObject *self, PyObject *args)
{
	//pass arguments from python to function
	double z, p, gam_1, gam_2;
	if (!PyArg_ParseTuple(args, "dddd", &z, &p, &gam_1, &gam_2))
  	  return NULL;
	
	return Py_BuildValue("d", flux_quad(z,p,gam_1,gam_2));	
}

static PyObject * Py_flux_quad_np(PyObject *self, PyObject *args)
{
	//define array object
	PyObject * z;
	PyArrayObject * array;
	PyArrayObject * output;
	
	double p, gam_1, gam_2;
	
	//read in object
	if (!PyArg_ParseTuple(args, "Oddd", &z, &p, &gam_1, &gam_2))
  	  return NULL;
	
	//convert to contiguous array of doubles - will accept up to 10d arrays
	array = (PyArrayObject*) PyArray_ContiguousFromObject(z, PyArray_DOUBLE, 1, 10);
	if(array == NULL) return NULL;
	
	//calculate no of data size of array
	int n,size=1;
	for(n=0;n<array->nd;n++)
		size *= array->dimensions[n];
	
	//create output array
	output = (PyArrayObject*) PyArray_SimpleNew(array->nd, array->dimensions, PyArray_DOUBLE);
	
	//create data for output array from input array
	for(n=0;n<size;n++)
		*(((double*)output->data) + n) = flux_quad(*(((double*)array->data) + n),p,gam_1,gam_2);
	
	//must destroy reference to array
	Py_DECREF(array);
	
	//return array
	return PyArray_Return(output);
}

static PyObject * Py_flux_nonlin_np(PyObject *self, PyObject *args)
{
	//define array object
	PyObject * z;
	PyArrayObject * array;
	PyArrayObject * output;
	
	double p, c1, c2, c3, c4;
	
	//read in object
	if (!PyArg_ParseTuple(args, "Oddddd", &z, &p, &c1, &c2, &c3, &c4))
  	  return NULL;
	
	//convert to contiguous array of doubles - will accept up to 10d arrays
	array = (PyArrayObject*) PyArray_ContiguousFromObject(z, PyArray_DOUBLE, 1, 10);
	if(array == NULL) return NULL;
	
	//calculate no of data size of array
	int n,size=1;
	for(n=0;n<array->nd;n++)
		size *= array->dimensions[n];
	
	//create output array
	output = (PyArrayObject*) PyArray_SimpleNew(array->nd, array->dimensions, PyArray_DOUBLE);
	
	//create data for output array from input array
	for(n=0;n<size;n++)
		*(((double*)output->data) + n) = flux_nonlin(*(((double*)array->data) + n),p,c1,c2,c3,c4);
	
	//must destroy reference to array
	Py_DECREF(array);
	
	//return array
	return PyArray_Return(output);
}
