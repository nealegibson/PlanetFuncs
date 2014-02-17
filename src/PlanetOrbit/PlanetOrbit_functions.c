
#include "PlanetOrbit_functions.h"

/**************************************************************************/
//West-East co-ord
double get_x(double M, double a_rstar, double ecc, double w)
{
	
	double f,E,r;
	
	//ensure M lies in 0-2pi range
	while(M > 2*M_PI){
		M -= 2*M_PI;
		}
	while(M < -2*M_PI){
		M += 2*M_PI;
		}
	
	//get eccentric anomaly
	E = ecc_anom(ecc,M);

	//get true anomaly
	f = true_anom(ecc,E);
	
	//calculate r as a function of a,e,f (dist between star and planet)
	r = a_rstar * (1 - SQ(ecc)) / (1 + ecc*cos(f));
	
	//calculate and return x coordinate
	//return r * ( cos(f)*cos(w) - sin(f)*sin(w) );
	return r * cos(f + w); //simplified version
}

/**************************************************************************/
//into/out of 'page' co-ord
double get_y(double M, double a_rstar, double ecc, double w, double i)
{
	double f,E,r;
	
	//ensure M lies in 0-2pi range
	while(M > 2*M_PI){
		M -= 2*M_PI;
		}
	while(M < -2*M_PI){
		M += 2*M_PI;
		}
	
	//get eccentric anomaly
	E = ecc_anom(ecc,M);

	//get true anomaly
	f = true_anom(ecc,E);
	
	//calculate r as a function of a,e,f (dist between star and planet)
	r = a_rstar * (1 - SQ(ecc)) / (1 + ecc*cos(f));
	
	//calculate and return x coordinate
	//return r * ( sin(f)*cos(w)*sin(i) + cos(f)*sin(w)*sin(i) );
	return r * sin(f + w) * sin(i); //simplified version
}

/**************************************************************************/
//up-down co-ord
double get_z(double M, double a_rstar, double ecc, double w, double i)
{
	double f,E,r;
	
	//ensure M lies in 0-2pi range
	while(M > 2*M_PI){
		M -= 2*M_PI;
		}
	while(M < -2*M_PI){
		M += 2*M_PI;
		}
	
	//get eccentric anomaly
	E = ecc_anom(ecc,M);

	//get true anomaly
	f = true_anom(ecc,E);
	
	//calculate r as a function of a,e,f (dist between star and planet)
	r = a_rstar * (1 - SQ(ecc)) / (1 + ecc*cos(f));
	
	//calculate and return x coordinate
	//return r * ( sin(f)*cos(w)*cos(i) + cos(f)*sin(w)*cos(i) );
	return r * sin(f + w) * cos(i); //simplified version
}

/**************************************************************************/
//normalised separation according to observer
double get_norm(double M, double a_rstar, double ecc, double w, double i)
{
	
	double f,E,r;
	double x,z,norm;
	
	//ensure M lies in 0-2pi range
	while(M > 2*M_PI){
		M -= 2*M_PI;
		}
	while(M < -2*M_PI){
		M += 2*M_PI;
		}
	
	//get eccentric anomaly
	E = ecc_anom(ecc,M);

	//get true anomaly
	f = true_anom(ecc,E);
	
	//calculate r as a function of a,e,f (dist between star and planet)
	r = a_rstar * (1 - SQ(ecc)) / (1 + ecc*cos(f));
	
	//get x and y co-ord
	//x = r * ( cos(f)*cos(w) - sin(f)*sin(w) );
	//z = r * ( sin(f)*cos(w)*cos(i) + cos(f)*sin(w)*cos(i) );
	//simplified versions
	x = r * cos(f + w);
	z = r * sin(f + w) * cos(i);
	
	//and calculate normalised separation
	norm = sqrt( SQ(x) + SQ(z) );
	
	//printf("%lf %lf %lf\n",x,z,norm);
	
	//calculate and return x coordinate
	return norm;
}

/**************************************************************************/

//Solves the Kepler equation  numerically using a Newton-Raphson iteration
//method, given the eccentricity e and the mean anomoly M. Method outlined 
//in Solar System Dynamics book (Murray & Dermott)
double ecc_anom(double e, double M)
{
	double E;
	//initial guess for E0
	if(sin(M)>0)
		E = M + 0.85*e;
	else
		E = M - 0.85*e;
		
	double f_E, fder_E;
	
	do{
		f_E = E - e*sin(E) - M;
		fder_E = 1 - e*cos(E);
		E = E - f_E/fder_E;
	}while(f_E > 0.000001 || f_E < -0.000001);
	
	return E;
}

/**************************************************************************/

double true_anom(double e, double E)
{
	double f = acos( (cos(E) - e) / (1. - e*cos(E)) );
	
	//correct for pi < E < 2pi and -pi < E < 0 (due to arccos calc)
	if (E > M_PI && E <= 2*M_PI) f = 2*M_PI - f;
	if(E < 0 && E >= -M_PI) f = 2*M_PI - f;
	
	return f;
}

/**************************************************************************/
