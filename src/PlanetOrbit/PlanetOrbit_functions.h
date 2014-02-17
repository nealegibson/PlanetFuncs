
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SQ(x) ((x)*(x))

double get_x(double M, double a_rstar, double ecc, double peri);
double get_y(double M, double a_rstar, double ecc, double peri, double i);
double get_z(double M, double a_rstar, double ecc, double peri, double i);
double get_norm(double M, double a_rstar, double ecc, double peri, double i);
double ecc_anom(double e, double M);
double true_anom(double e, double E);
