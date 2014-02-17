/*********************************************************************************
*********************************************************************************/

#include "WaveletTransform.h"

/*********************************************************************************
Initialise the module and define docstrings etc
*********************************************************************************/

//define methods in the module
PyMethodDef WaveletTransformMethods[] = {
    {"FWT", (PyCFunction) WaveletTransform, METH_VARARGS|METH_KEYWORDS, NULL},
    {"IFWT", (PyCFunction) InvWaveletTransform, METH_VARARGS|METH_KEYWORDS, NULL},
    {NULL, NULL, 0, NULL}
};

//initialise the module
PyMODINIT_FUNC 
initWaveletTransform(void)
{
    (void) Py_InitModule3("WaveletTransform", WaveletTransformMethods, NULL);
    import_array(); //must be present for numpy stuff
}

/********************************************************************************/
static PyObject * WaveletTransform(PyObject *self, PyObject *args, PyObject *keywds)
{
	//define array object
	PyArrayObject * input, * output;
	int filter = 0, dorder = 4;
	
    static char *kwlist[] = {"input","filter","dorder", NULL};
	
	//read in object
	if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|ii", kwlist, &input, &filter, &dorder))
  	  return NULL;
	
	//copy input to output - first create the array
	output = (PyArrayObject*) PyArray_SimpleNew(input->nd, input->dimensions, PyArray_DOUBLE);
	if(output == NULL) return NULL;
		
	//calculate no of data points in array
	int n,size=1;
	for(n=0;n<output->nd;n++)
		size *= output->dimensions[n];

	//copy input data to output data
	for(n=0;n<size;n++)
		*(((double*)output->data) + n) = *(((double*)input->data) + n);
	
	//check that data is a power of 2!
	if(!IsPow2(size)){
		fprintf(stderr, "DoFWT Error: data for DoIFWT is not a power of 2. Exiting...\n");
		return PyFloat_FromDouble(1.);
		}
	
	//do wavelet transform **in place**
	DoFWT((double*)output->data,size,dorder);
	
	//return the transformed array
	return PyArray_Return(output);
}

/********************************************************************************/
static PyObject * InvWaveletTransform(PyObject *self, PyObject *args, PyObject *keywds)
{
	//define array object
	PyArrayObject * input, * output;
	int filter = 0, dorder = 4;
	
    static char *kwlist[] = {"input","filter","dorder", NULL};
	
	//read in object
	if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|ii", kwlist, &input, &filter, &dorder))
  	  return NULL;
	
	//copy input to output - first create the array
	output = (PyArrayObject*) PyArray_SimpleNew(input->nd, input->dimensions, PyArray_DOUBLE);
	if(output == NULL) return NULL;
		
	//calculate no of data points in array
	int n,size=1;
	for(n=0;n<output->nd;n++)
		size *= output->dimensions[n];

	//copy input data to output data
	for(n=0;n<size;n++)
		*(((double*)output->data) + n) = *(((double*)input->data) + n);
	
	//check that data is a power of 2!
	if(!IsPow2(size)){
		fprintf(stderr, "DoIFWT Error: data for DoIFWT is not a power of 2. Exiting...\n");
		return PyFloat_FromDouble(1.);
		}
	
	//do inv wavelet transform **in place**
	DoIFWT((double*)output->data,size,dorder);
	
	//return the transformed array
	return PyArray_Return(output);
}

/********************************************************************************/
//function to return the next highest power of two after N
size_t MakePowerOf2(size_t N)
{
	size_t x = 1;
	
	while(x<N) x*=2;
	
	return x;
}
/********************************************************************************/
//function to return the next highest power of two after N
int LogTwo(int N)
{
	int x = 0;
	
	while(N != 1){
		N/=2;
		x+=1;
		}
		
	return x;
}
/********************************************************************************/
