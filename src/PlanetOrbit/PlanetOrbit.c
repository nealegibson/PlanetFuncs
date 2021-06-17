/*********************************************************************************
*********************************************************************************/

#include "PlanetOrbit.h"
#include "PlanetOrbit_functions.h"

/*********************************************************************************
Initialise the module and define docstrings etc
*********************************************************************************/

//define methods in the module
PyMethodDef PlanetOrbitMethods[] = {
    {"get_x", (PyCFunction)get_x_py, METH_VARARGS, NULL},
    {"get_y", (PyCFunction)get_y_py, METH_VARARGS, NULL},
    {"get_z", (PyCFunction)get_z_py, METH_VARARGS, NULL},
    {"get_norm", (PyCFunction)get_norm_py, METH_VARARGS, NULL},
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
//initialise the module for python3

static struct PyModuleDef Methods = {
    PyModuleDef_HEAD_INIT,
    "PlanetOrbit",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    PlanetOrbitMethods
};

PyMODINIT_FUNC PyInit_PlanetOrbit(void)
{
    PyObject *module = PyModule_Create(&Methods);
    import_array(); //must be present for numpy stuff
    return module;
}

#else
//initialise the module for python2
PyMODINIT_FUNC 
initPlanetOrbit(void)
{
    (void) Py_InitModule3("PlanetOrbit", PlanetOrbitMethods, "PlanetOrbit docstring...");
    import_array(); //must be present for numpy stuff
}
#endif

/********************************************************************************/
static PyObject * get_x_py(PyObject *self, PyObject *args, PyObject *keywds)
{
	//for numpy array
	PyObject * M_object;
	PyArrayObject * M, * x_coord;
	
	//planet parameter inputs
	double a_rstar,ecc,peri; //Mean_anom, a/rstar, ecc, arg of periastron
	
	//PyArrayObject * output;
	int size,i;
	
	//read in python arguments
	if (!PyArg_ParseTuple(args, "Oddd", &M_object, &a_rstar, &ecc, &peri))
	    {
		printf("Problem loading args...\n");
		return NULL;
		}

	//convert numpy object to contiguous array of doubles - will accept up to 10d arrays
	M = (PyArrayObject*) PyArray_ContiguousFromObject(M_object, NPY_DOUBLE, 1, 10);
	if(M == NULL) return NULL;
	
  //get no of data points in array
	size=PyArray_SIZE(M);
	
	//create output double array
	x_coord = (PyArrayObject*) PyArray_SimpleNew(PyArray_NDIM(M), PyArray_SHAPE(M), NPY_DOUBLE);
	
	//do calculation on numpy array
    for(i=0;i<size;i++){
    	*(((double*)PyArray_DATA(x_coord)) + i) = get_x(*(((double*)PyArray_DATA(M))+i),a_rstar,ecc,peri);
		//printf("a: %lf %lf %lf %lf %lf\n", get_x(*(((double*)PyArray_DATA(M))+i),a_rstar,ecc,peri),*(((double*)PyArray_DATA(M))+i),a_rstar,ecc,peri);
		}
	
	//destroy reference to array
	Py_DECREF(M);
	
	//return array
	return PyArray_Return(x_coord);
}
/********************************************************************************/
static PyObject * get_y_py(PyObject *self, PyObject *args, PyObject *keywds)
{
	//for numpy array
	PyObject * M_object;
	PyArrayObject * M, * y_coord;
	
	//planet parameter inputs
	double a_rstar,ecc,peri,incl; //Mean_anom, a/rstar, ecc, arg of periastron
	
	//PyArrayObject * output;
	int size,i;
	
	//read in python arguments
	if (!PyArg_ParseTuple(args, "Odddd", &M_object, &a_rstar, &ecc, &peri, &incl))
	    {
		printf("Problem loading args...\n");
		return NULL;
		}

	//convert numpy object to contiguous array of doubles - will accept up to 10d arrays
	M = (PyArrayObject*) PyArray_ContiguousFromObject(M_object, NPY_DOUBLE, 1, 10);
	if(M == NULL) return NULL;
	
  //get no of data points in array
	size=PyArray_SIZE(M);
	
	//create output double array
	y_coord = (PyArrayObject*) PyArray_SimpleNew(PyArray_NDIM(M), PyArray_SHAPE(M), NPY_DOUBLE);

	//do calculation on numpy array
    for(i=0;i<size;i++){
    	*(((double*)PyArray_DATA(y_coord)) + i) = get_y(*(((double*)PyArray_DATA(M))+i),a_rstar,ecc,peri,incl);
		//printf("a: %lf %lf %lf %lf %lf\n", get_x(*(((double*)PyArray_DATA(M))+i),a_rstar,ecc,peri),*(((double*)PyArray_DATA(M))+i),a_rstar,ecc,peri);
		//printf("%lf\n", *(((double*)y_coord->data) + i));
		}
	
	//destroy reference to array
	Py_DECREF(M);
	
	//return array
	return PyArray_Return(y_coord);
}
/********************************************************************************/
static PyObject * get_z_py(PyObject *self, PyObject *args, PyObject *keywds)
{
	//for numpy array
	PyObject * M_object;
	PyArrayObject * M, * z_coord;
	
	//planet parameter inputs
	double a_rstar,ecc,peri,incl; //Mean_anom, a/rstar, ecc, arg of periastron
	
	//PyArrayObject * output;
	int size,i;
	
	//read in python arguments
	if (!PyArg_ParseTuple(args, "Odddd", &M_object, &a_rstar, &ecc, &peri, &incl))
	    {
		printf("Problem loading args...\n");
		return NULL;
		}

	//convert numpy object to contiguous array of doubles - will accept up to 10d arrays
	M = (PyArrayObject*) PyArray_ContiguousFromObject(M_object, NPY_DOUBLE, 1, 10);
	if(M == NULL) return NULL;
	
  //get no of data points in array
	size=PyArray_SIZE(M);
	
	//create output double array
	z_coord = (PyArrayObject*) PyArray_SimpleNew(PyArray_NDIM(M), PyArray_SHAPE(M), NPY_DOUBLE);

	//do calculation on numpy array
    for(i=0;i<size;i++){
    	*(((double*)PyArray_DATA(z_coord)) + i) = get_z(*(((double*)PyArray_DATA(M))+i),a_rstar,ecc,peri,incl);
		//printf("a: %lf %lf %lf %lf %lf\n", get_x(*(((double*)PyArray_DATA(M))+i),a_rstar,ecc,peri),*(((double*)PyArray_DATA(M))+i),a_rstar,ecc,peri);
		//printf("%lf\n", *(((double*)z_coord->data) + i));
		}
	
	//destroy reference to array
	Py_DECREF(M);
	
	//return array
	return PyArray_Return(z_coord);
}
/********************************************************************************/
static PyObject * get_norm_py(PyObject *self, PyObject *args, PyObject *keywds)
{
	//for numpy array
	PyObject * M_object;
	PyArrayObject * M, * norm;
	
	//planet parameter inputs
	double a_rstar,ecc,peri,incl; //Mean_anom, a/rstar, ecc, arg of periastron
	
	//PyArrayObject * output;
	int size,i;
	
	//read in python arguments
	if (!PyArg_ParseTuple(args, "Odddd", &M_object, &a_rstar, &ecc, &peri, &incl))
	    {
		printf("Problem loading args...\n");
		return NULL;
		}

	//convert numpy object to contiguous array of doubles - will accept up to 10d arrays
	M = (PyArrayObject*) PyArray_ContiguousFromObject(M_object, NPY_DOUBLE, 1, 10);
	if(M == NULL) return NULL;
	
  //get no of data points in array
	size=PyArray_SIZE(M);
	
	//create output double array
	norm = (PyArrayObject*) PyArray_SimpleNew(PyArray_NDIM(M), PyArray_SHAPE(M), NPY_DOUBLE);

	//do calculation on numpy array
    for(i=0;i<size;i++){
    	*(((double*)PyArray_DATA(norm)) + i) = get_norm(*(((double*)PyArray_DATA(M))+i),a_rstar,ecc,peri,incl);
		//printf("a: %lf %lf %lf %lf %lf\n", get_x(*(((double*)PyArray_DATA(M))+i),a_rstar,ecc,peri),*(((double*)PyArray_DATA(M))+i),a_rstar,ecc,peri);
		//printf("%lf\n", *(((double*)z_coord->data) + i));
		}
	
	//destroy reference to array
	Py_DECREF(M);
	
	//return array
	return PyArray_Return(norm);
}
/********************************************************************************/
