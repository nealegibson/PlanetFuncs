/*********************************************************************************
Neale Gibson Jan 2009

Python wrapper to calculate radial velocites from python and return numpy array.
*********************************************************************************/

#include "RadialVelocity.h"

/*********************************************************************************
Initialise the module and define docstrings etc
*********************************************************************************/
//python function prototypes
static PyObject * RadVel_np(PyObject *self, PyObject *args);

//define methods in the module
PyMethodDef RadialVelocityMethods[] = {
    {"RadVel", RadVel_np, METH_VARARGS, "Returns the radial velocity for params..."},
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
//initialise the module for python3

static struct PyModuleDef Methods = {
    PyModuleDef_HEAD_INIT,
    "RadialVelocity",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    RadialVelocityMethods
};

PyMODINIT_FUNC PyInit_RadialVelocity(void)
{
    PyObject *module = PyModule_Create(&Methods);
    import_array(); //must be present for numpy stuff
    return module;
}

#else
//initialise the module for python2
PyMODINIT_FUNC 
initRadialVelocity(void)
{
    (void) Py_InitModule3("RadialVelocity", RadialVelocityMethods, "RadialVelocity...");
    import_array(); //must be present for numpy stuff
}
#endif

/*********************************************************************************
Python functions
*********************************************************************************/
//numpy wrapper for radial velocity function
static PyObject * RadVel_np(PyObject *self, PyObject *args)
{
	//define array object
	PyObject * t;
	PyArrayObject * array;
	PyArrayObject * output;
	
	double e,P,w,K;
	
	//read in object
	if (!PyArg_ParseTuple(args, "Odddd", &t, &e, &P, &w, &K))
  	  return NULL;
	
	//convert to contiguous array of doubles - will accept up to 10d arrays
	array = (PyArrayObject*) PyArray_ContiguousFromObject(t, PyArray_DOUBLE, 1, 10);
	if(array == NULL) return NULL;
	
	//calculate no of data size of array
	int n,size=1;
	for(n=0;n<array->nd;n++)
		size *= array->dimensions[n];
	
	//create output array
	output = (PyArrayObject*) PyArray_SimpleNew(array->nd, array->dimensions, PyArray_DOUBLE);
	
	//create data for output array from input array
	for(n=0;n<size;n++)
		*(((double*)output->data) + n) = Rad_vel(*(((double*)array->data) + n),e,P,w,K);
	
	//must destroy reference to array
	Py_DECREF(array);
	
	//return array
	return PyArray_Return(output);
}

/*********************************************************************************
*********************************************************************************/
//function to calculate radial velocity given time since periapsis, ecc,
//period, argument of periapsis w, and semi-amplitude K
double Rad_vel(double M, double e, double p, double w, double K)
{
	double E, f;

	//correct so that M lies in 0 > 2PI range
	while(M > 2*M_PI){
		M -= 2*M_PI;}
	while(M < 0){
		M += 2*M_PI;}

	//calculate eccentric anomoly from M and e
	E = ecc_anom(e, M);

	//now can calculate true anomoly from E
	f = acos( (cos(E) - e) / (1 - e*cos(E)) );
	
	//need to correct f however as it is only returned in range 0<f<PI, should
	//be 2PI - f if M is in the range PI < M < 2PI
	if(M >= M_PI && M < 2*M_PI)
		f = 2*M_PI - f;
	
	//now calculate and return radial velocity
	return K * (e*cos(w) + cos(f + w));
}
/***************************************************************************
***************************************************************************/
//Solves the Kepler equation  numerically using a Newton-Raphson iteration
//method, given the eccentricity e and the mean anomoly M. Method outlined 
//in Solar System Dynamics book (Murray & Dermott)
double ecc_anom(double e, double M)
{
	double E,f_E,fder_E;
	//initial guess for E0
	if(sin(M)>0)
		E = M + 0.85*e;
	else
		E = M - 0.85*e;
	
	do{
		f_E = E - e*sin(E) - M;
		fder_E = 1 - e*cos(E);
		E = E - f_E/fder_E;
	}while(f_E > 0.000001 || f_E < -0.000001);
	
	return E;
}
/***************************************************************************
***************************************************************************/
