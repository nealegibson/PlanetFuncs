
#ifndef __FLUX_H__
#define __FLUX_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "mconf.h"
#include "gsl/gsl_math.h" 
#include "gsl/gsl_sf_ellint.h" //for elliptic integrals
#include "gsl/gsl_sf_gamma.h" //for gamma func, beta func, factorials, Pochhammer symbol
#include "gsl/gsl_sf_hyperg.h" //for hypergeometric funcitons

#include <Python.h>

//scipy/cephes headers
#include "mconf.h"

//#define DEBUG
//#define MODE GSL_PREC_DOUBLE
#define MODE GSL_PREC_SINGLE //can safely use single precision
#define SQ(x) ((x)*(x))
#define EPSILON 2e-12
#define DOUBLE_EQ(x,v) (((v - EPSILON) < x) && (x <( v + EPSILON)))
#define DOUBEQ(a,b) (fabs(a-b) < 1E-8)

//scipy/cephes Gauss hypergeometric function
double hyp2f1( double a, double b, double c, double x );

//call scipy hyp2f1 from python
//double scipy_hyp2f1(double a, double b, double c, double x);

//parameters for the hypergeometric function loop
#define ETOL_HYPER 1e-6 // default 1e-8
#define MAX_ITER_HYPER 30 // default 50

//main quad flux functions, include this prototype in any other c prog calling this
double flux_quad(double z, double p, double gam_1, double gam_2);
double flux_nonlin(double z, double p, double gam_1, double gam_2, double c_3, double c_4);

//functions to compute uniform source (lambda_e)
double lambda_e_pi_funct(double z, double p);
double kappa_0(double z, double p);
double kappa_1(double z, double p);

//lambda_d functions
double lambda_1(double z, double p,double a, double b, double q, double k);
double lambda_2(double z, double p,double a, double b, double q, double k);
double lambda_3(double p, double k);
double lambda_4(double p, double k);
double lambda_5(double p);
double lambda_6(double p);

//eta_d functions
double eta_1(double z, double p, double a, double b);
double eta_2(double z, double p);

//hypergeometric functions (for flux_nonlin)
double Appell_Hyper(double a,double b1,double b2,double c,double x,double y);

//N and M funcitons for flux_nonlin
double NF(int n, double a, double b, double z, double p);
double MF(int n, double a, double b, double z, double p);
double Omega(double * C);

//further functions
//double mod(double x); //modulus function
//double RF(double a, double b); //rising factorial function
void swap_double(double *x1, double *x2); //swap doubles

#endif