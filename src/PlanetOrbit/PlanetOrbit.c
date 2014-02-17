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

//initialise the module
PyMODINIT_FUNC initPlanetOrbit(void)
{
    (void) Py_InitModule3("PlanetOrbit", PlanetOrbitMethods, NULL);
    import_array(); //must be present for numpy stuff
}

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
	M = (PyArrayObject*) PyArray_ContiguousFromObject(M_object, PyArray_DOUBLE, 1, 10);
	if(M == NULL) return NULL;
	
    //calculate no of data points in array
	size=1;
	for(i=0;i<M->nd;i++) size *= M->dimensions[i];
	
	//create output double array
	x_coord = (PyArrayObject*) PyArray_SimpleNew(M->nd, M->dimensions, PyArray_DOUBLE);
	
	//do calculation on numpy array
    for(i=0;i<size;i++){
    	*(((double*)x_coord->data) + i) = get_x(*(((double*)M->data)+i),a_rstar,ecc,peri);
		//printf("a: %lf %lf %lf %lf %lf\n", get_x(*(((double*)M->data)+i),a_rstar,ecc,peri),*(((double*)M->data)+i),a_rstar,ecc,peri);
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
	M = (PyArrayObject*) PyArray_ContiguousFromObject(M_object, PyArray_DOUBLE, 1, 10);
	if(M == NULL) return NULL;
	
    //calculate no of data points in array
	size=1;
	for(i=0;i<M->nd;i++) size *= M->dimensions[i];
	
	//create output double array
	y_coord = (PyArrayObject*) PyArray_SimpleNew(M->nd, M->dimensions, PyArray_DOUBLE);

	//do calculation on numpy array
    for(i=0;i<size;i++){
    	*(((double*)y_coord->data) + i) = get_y(*(((double*)M->data)+i),a_rstar,ecc,peri,incl);
		//printf("a: %lf %lf %lf %lf %lf\n", get_x(*(((double*)M->data)+i),a_rstar,ecc,peri),*(((double*)M->data)+i),a_rstar,ecc,peri);
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
	M = (PyArrayObject*) PyArray_ContiguousFromObject(M_object, PyArray_DOUBLE, 1, 10);
	if(M == NULL) return NULL;
	
    //calculate no of data points in array
	size=1;
	for(i=0;i<M->nd;i++) size *= M->dimensions[i];
	
	//create output double array
	z_coord = (PyArrayObject*) PyArray_SimpleNew(M->nd, M->dimensions, PyArray_DOUBLE);

	//do calculation on numpy array
    for(i=0;i<size;i++){
    	*(((double*)z_coord->data) + i) = get_z(*(((double*)M->data)+i),a_rstar,ecc,peri,incl);
		//printf("a: %lf %lf %lf %lf %lf\n", get_x(*(((double*)M->data)+i),a_rstar,ecc,peri),*(((double*)M->data)+i),a_rstar,ecc,peri);
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
	M = (PyArrayObject*) PyArray_ContiguousFromObject(M_object, PyArray_DOUBLE, 1, 10);
	if(M == NULL) return NULL;
	
    //calculate no of data points in array
	size=1;
	for(i=0;i<M->nd;i++) size *= M->dimensions[i];
	
	//create output double array
	norm = (PyArrayObject*) PyArray_SimpleNew(M->nd, M->dimensions, PyArray_DOUBLE);

	//do calculation on numpy array
    for(i=0;i<size;i++){
    	*(((double*)norm->data) + i) = get_norm(*(((double*)M->data)+i),a_rstar,ecc,peri,incl);
		//printf("a: %lf %lf %lf %lf %lf\n", get_x(*(((double*)M->data)+i),a_rstar,ecc,peri),*(((double*)M->data)+i),a_rstar,ecc,peri);
		//printf("%lf\n", *(((double*)z_coord->data) + i));
		}
	
	//destroy reference to array
	Py_DECREF(M);
	
	//return array
	return PyArray_Return(norm);
}
/********************************************************************************/
