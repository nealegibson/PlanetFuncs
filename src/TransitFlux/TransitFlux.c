/*********************************************************************************
Neale Gibson 15/3/2007 flux_quad.c
This function uses the analytic equations of Mandel and Agol (2000) to compute the
flux of a star during occultation given a normalized separation z, ratio of planet
to star radii, and the two quadratic limb darkening coefficients gam_1 and gam_2.

Python Wrapper Aug 2009. flux_quad_python.c
*********************************************************************************/

#define PY_SSIZE_T_CLEAN //this redefines length of input args I think - recommended to call before Python.h
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION //this switches off the ability to use the old API
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
static PyObject * Py_flux_quad_np(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject * Py_flux_nonlin_np(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject * Py_flux_quad_single(PyObject *self, PyObject *args);

//define methods in the module
static PyMethodDef mymethods[] = { //funcs with kwargs need to be cast as a PyCFunction
    {"flux_quad_single", Py_flux_quad_single, METH_VARARGS, "Returns flux as a function of z"},
    {"FluxQuad", (PyCFunction) Py_flux_quad_np, METH_VARARGS|METH_KEYWORDS, "Returns flux as a function of z"},
    {"FluxNonlin", (PyCFunction) Py_flux_nonlin_np, METH_VARARGS|METH_KEYWORDS, "Returns flux as a function of z"},
    {NULL, NULL, 0, NULL}
};

//define the scipy.special.hyp2f1 function as global
//PyObject *pName, *pModule, *pFunc;
//double qqq = 234.987;
//static PyObject * pName = PyString_FromString("scipy.special");
//Py_DECREF(pName);
//static PyObject * pModule = PyImport_Import(pName);
//static PyObject * pFunc = PyObject_GetAttrString(pModule, "hyp2f1");

#if PY_MAJOR_VERSION >= 3
//initialise the module for python3
static struct PyModuleDef Methods = {
    PyModuleDef_HEAD_INIT,
    "TransitFlux",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    mymethods
};

PyMODINIT_FUNC PyInit_TransitFlux(void)
{
    PyObject *module = PyModule_Create(&Methods);
    import_array(); //must be present for numpy stuff
    return module;
}

#else
//initialise the module for python2
PyMODINIT_FUNC 
initTransitFlux(void)
{
    (void) Py_InitModule3("TransitFlux", mymethods, "Flux_docstring...");
    import_array(); //must be present for numpy stuff
}
#endif


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
// This uses old numpy API so will remove

// static PyObject * scale(PyObject *self, PyObject *args)
// {
// 	//define array object
// 	PyObject * input;
// 	PyArrayObject * array;
// 	PyArrayObject * output;
// 	
// 	//
// 	double scale;
// 	
// 	//read in object
// 	if (!PyArg_ParseTuple(args, "Od", &input, &scale))
//   	  return NULL;
// 	
// 	//convert to contiguous array of doubles - will accept up to 10d arrays
// 	array = (PyArrayObject*) PyArray_ContiguousFromObject(input, NPY_DOUBLE, 1, 10);
// 	if(array == NULL) return NULL;
// 	
// 	//calculate no of data points in array
// 	int n,size=1;
// 	for(n=0;n<array->nd;n++)
// 		size *= array->dimensions[n];
// 		
// 	printf("size = %d\n", size);
// 
// 	//create output array
// 	output = (PyArrayObject*) PyArray_SimpleNew(array->nd, array->dimensions, NPY_DOUBLE);
// 	
// 	//create data for output array from input array
// 	for(n=0;n<size;n++)
// 		*(((double*)output->data) + n) = *(((double*)array->data) + n) * scale;
// 	
// 	//must destroy reference to array
// 	Py_DECREF(array);
// 	
// 	//return array
// 	return PyArray_Return(output);
// }

//wrapper for flux_quad c function - only accepts single value of z
static PyObject * Py_flux_quad_single(PyObject *self, PyObject *args)
{
	//pass arguments from python to function
	double z, p, gam_1, gam_2;
	if (!PyArg_ParseTuple(args, "dddd", &z, &p, &gam_1, &gam_2))
  	  return NULL;
	
	return Py_BuildValue("d", flux_quad(z,p,gam_1,gam_2));	
}

// static PyObject * Py_flux_quad_np(PyObject *self, PyObject *args)
// {
// 	//define array object
// 	PyObject * z;
// 	PyArrayObject * array;
// 	PyArrayObject * output;
// 	
// 	double p, gam_1, gam_2;
// 	
// 	//read in object
// 	if (!PyArg_ParseTuple(args, "Oddd", &z, &p, &gam_1, &gam_2))
//   	  return NULL;
// 	
// 	//convert to contiguous array of doubles - will accept up to 10d arrays
// 	array = (PyArrayObject*) PyArray_ContiguousFromObject(z, NPY_DOUBLE, 1, 10);
// 	if(array == NULL) return NULL;
// 	
// 	//calculate no of data size of array
// 	int n,size=1;
// 	for(n=0;n<array->nd;n++)
// 		size *= array->dimensions[n];
// 	
// 	//create output array
// 	output = (PyArrayObject*) PyArray_SimpleNew(array->nd, array->dimensions, NPY_DOUBLE);
// 	
// 	//create data for output array from input array
// 	for(n=0;n<size;n++)
// 		*(((double*)output->data) + n) = flux_quad(*(((double*)array->data) + n),p,gam_1,gam_2);
// 	
// 	//must destroy reference to array
// 	Py_DECREF(array);
// 	
// 	//return array
// 	return PyArray_Return(output);
// }

// static PyObject * Py_flux_nonlin_np(PyObject *self, PyObject *args)
// {
// 	//define array object
// 	PyObject * z;
// 	PyArrayObject * array;
// 	PyArrayObject * output;
// 	
// 	double p, c1, c2, c3, c4;
// 	
// 	//read in object
// 	if (!PyArg_ParseTuple(args, "Oddddd", &z, &p, &c1, &c2, &c3, &c4))
//   	  return NULL;
// 	
// 	//convert to contiguous array of doubles - will accept up to 10d arrays
// 	array = (PyArrayObject*) PyArray_ContiguousFromObject(z, NPY_DOUBLE, 1, 10);
// 	if(array == NULL) return NULL;
// 	
// 	//calculate no of data size of array
// 	int n,size=1;
// 	for(n=0;n<array->nd;n++)
// 		size *= array->dimensions[n];
// 	
// 	//create output array
// 	output = (PyArrayObject*) PyArray_SimpleNew(array->nd, array->dimensions, NPY_DOUBLE);
// 	
// 	//create data for output array from input array
// 	for(n=0;n<size;n++)
// 		*(((double*)output->data) + n) = flux_nonlin(*(((double*)array->data) + n),p,c1,c2,c3,c4);
// 	
// 	//must destroy reference to array
// 	Py_DECREF(array);
// 	
// 	//return array
// 	return PyArray_Return(output);
// }

/**
New versions below with updated numpy API
**/

static PyObject * Py_flux_quad_np(PyObject *self, PyObject *args, PyObject *kwargs)
{
  //for keywords, define keyword list - should always be NULL terminated
  static char *kwlist[] = {"z","p","c1","c2","out",NULL};

	//define variables for input
	PyArrayObject * z=NULL, *out=NULL; //initiate as null pointers as some are optional
 	//define values passed from python
 	double p, c1, c2;

  //read in the inputs - args and keywords
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Oddd|O", kwlist, &z, &p, &c1, &c2, &out))
    return NULL;

  //define C pointers to data arrays - initiate to NULL
  double *z_data = NULL, *out_data = NULL;
	//define addition variables
	int n, out_new=0;    
  
  //check if output array passed
  if ( ((PyObject*)out) == Py_None) out = NULL;
  
  //check inputs are correct array types
  if (!PyArray_ISCARRAY(z) || !(PyArray_TYPE(z) == NPY_DOUBLE) || !(PyArray_NDIM(z) == 1))
    {PyErr_SetString(PyExc_TypeError, "Py_flux_nonlin_np error: z must be contiguous C 1D array of doubles"); return NULL;}
  
  //create output array if not given
  if (out==NULL){ //create output array if not given
    //just create array with same shape as z
    out_new = 1; //store flag so I can pass the correct array type to PyBuildVal
    out = (PyArrayObject*) PyArray_SimpleNew(1, PyArray_SHAPE(z), NPY_DOUBLE);
    }      
  //otherwise check if shape is ok
  else{
    if (!PyArray_ISCARRAY(out) || !(PyArray_TYPE(out) == NPY_DOUBLE) || !(PyArray_NDIM(out) == 1))
      {PyErr_SetString(PyExc_TypeError, "Py_flux_nonlin_np error: out must be contiguous C 1D array of doubles"); return NULL;}
    if (*PyArray_SHAPE(out) != *PyArray_SHAPE(z))
      {PyErr_SetString(PyExc_TypeError, "Py_flux_nonlin_np error: out must be same shape/size as z"); return NULL;}
    }
  
  //assign double pointer to the input z and output out
  z_data = (double*)PyArray_DATA(z); 
  out_data = (double*)PyArray_DATA(out);
  
	//now compute the output array
	for(n=0;n<*PyArray_SHAPE(z);n++)
		out_data[n] = flux_quad(z_data[n],p,c1,c2);
  
  //finally return the object
	if(out_new) return Py_BuildValue("N", out); //must use N for any new objects, otherwise I'll get a memory leak/
	else return Py_BuildValue("O", out); //must use O for any existing objects, which increases the ref count
}

static PyObject * Py_flux_nonlin_np(PyObject *self, PyObject *args, PyObject *kwargs)
{

  //for keywords, define keyword list - should always be NULL terminated
  static char *kwlist[] = {"z","p","c1","c2","c3","c4","out",NULL};

	//define variables for input
	PyArrayObject * z=NULL, *out=NULL; //initiate as null pointers as some are optional
 	//define values passed from python
 	double p, c1, c2, c3, c4;

  //read in the inputs - args and keywords
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Oddddd|O", kwlist, &z, &p, &c1, &c2, &c3, &c4, &out))
    return NULL;

  //define C pointers to data arrays - initiate to NULL
  double *z_data = NULL, *out_data = NULL;
	//define addition variables
	int n, out_new=0;    
  
  //check if output array passed
  if ( ((PyObject*)out) == Py_None) out = NULL;
  
  //check inputs are correct array types
  if (!PyArray_ISCARRAY(z) || !(PyArray_TYPE(z) == NPY_DOUBLE) || !(PyArray_NDIM(z) == 1))
    {PyErr_SetString(PyExc_TypeError, "Py_flux_nonlin_np error: z must be contiguous C 1D array of doubles"); return NULL;}
  
  //create output array if not given
  if (out==NULL){ //create output array if not given
    //just create array with same shape as z
    out_new = 1; //store flag so I can pass the correct array type to PyBuildVal
    out = (PyArrayObject*) PyArray_SimpleNew(1, PyArray_SHAPE(z), NPY_DOUBLE);
    }      
  //otherwise check if shape is ok
  else{
    if (!PyArray_ISCARRAY(out) || !(PyArray_TYPE(out) == NPY_DOUBLE) || !(PyArray_NDIM(out) == 1))
      {PyErr_SetString(PyExc_TypeError, "Py_flux_nonlin_np error: out must be contiguous C 1D array of doubles"); return NULL;}
    if (*PyArray_SHAPE(out) != *PyArray_SHAPE(z))
      {PyErr_SetString(PyExc_TypeError, "Py_flux_nonlin_np error: out must be same shape/size as z"); return NULL;}
    }
  
  //assign double pointer to the input z and output out
  z_data = (double*)PyArray_DATA(z); 
  out_data = (double*)PyArray_DATA(out);
  
	//now compute the output array
	for(n=0;n<*PyArray_SHAPE(z);n++)
		out_data[n] = flux_nonlin(z_data[n],p,c1,c2,c3,c4);
  
  //finally return the object
	if(out_new) return Py_BuildValue("N", out); //must use N for any new objects, otherwise I'll get a memory leak/
	else return Py_BuildValue("O", out); //must use O for any existing objects, which increases the ref count
  
}
