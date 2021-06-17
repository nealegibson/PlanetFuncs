
#include <stdio.h>
#include "flux.h"

//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
//#include <Python.h>
//#include <numpy/arrayobject.h> //for PyArray_Type objects

//function to loop over array and assign results from flux quad
int ctypes_flux_quad_np(double*z, double*f, int N, double p, double g1, double g2)
{

  int i;
  
  //simply loop over the array and assign value from flux quad
  for (i=0;i<N;i++)
    *(f+i) = flux_quad(*(z+i),p,g1,g2);

  return 1;

}

//function to loop over array and assign results from flux nonlin
int ctypes_flux_nonlin_np(double*z, double*f, int N, double p, double c1, double c2, double c3, double c4)
{

  int i;
  
  //simply loop over the array and assign value from flux nonlin
  for (i=0;i<N;i++)
    *(f+i) = flux_nonlin(*(z+i),p,c1,c2,c3,c4);

  return 1;

}

//function to loop over array and assign results from flux quad
// double * ctypes_flux_quad_np2(double*z, int size, double p, double g1, double g2)
// {
//   
//   int i;
//   
//   //define output array
//   double *output;
//   output = malloc(size*sizeof(double));
//   
//   //simply loop over the array and assign value from flux quad
//   for (i=0;i<size;i++){
//     printf("%d\n", i);
//     *(output+i) = flux_quad(*(z+i),p,g1,g2);
//     }
//     
//   //finally return numpy array
//   return output;
//   //return 1;
// 
// }
