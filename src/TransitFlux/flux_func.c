
#include "flux.h"

/**************************************************************************************************/
//Functions for flux_quad - elliptic integrals etc
double lambda_e_pi_funct(double z, double p)
{
	double f;
	
	f = sqrt( (4.0*SQ(z) - SQ(1.0 + SQ(z) - SQ(p)) )/(4.0) );
	
	return ( (1.0/M_PI) * (SQ(p)*kappa_0(z,p) + kappa_1(z,p) - f) );
}

double kappa_1(double z, double p)
{
	return acos( (1.0 - SQ(p) + SQ(z)) / (2.0*z) );
}

double kappa_0(double z, double p)
{
	return acos( (SQ(p) + SQ(z) - 1.0) / (2.0*p*z) );
}

double lambda_1(double z, double p,double a, double b, double q, double k)
{
	return ((1.0/(9.0*M_PI*sqrt(p*z)))*( ((1.0-b)*(2.0*b+a-3.0) - 3.0*q*(b-2.0))*gsl_sf_ellint_Kcomp(k, MODE) +
		4.0*p*z*(SQ(z) + 7.0*SQ(p) - 4.0)*gsl_sf_ellint_Ecomp(k, MODE) -
		3.0*(q/a)* (gsl_sf_ellint_Kcomp( k, MODE ) 
		- ((1.0/a-1.0)/3.0)*gsl_sf_ellint_RJ (0.0, 1.0-SQ(k), 1, 1/a, MODE))       ));
}

double lambda_2(double z, double p,double a, double b, double q, double k)
{
	return ((2.0/(9.0*M_PI*sqrt(1.0-a)))*( (1.0-5.0*SQ(z)+SQ(p)+SQ(q))*gsl_sf_ellint_Kcomp((1.0/k), MODE) +
		(1.0-a)*(SQ(z)+7*SQ(p)-4.0)*gsl_sf_ellint_Ecomp((1.0/k), MODE) -
		3.0*(q/a)*( gsl_sf_ellint_Kcomp( 1/k, MODE ) 
	    - ( (b/a-1)/3.0 )*gsl_sf_ellint_RJ (0.0, 1-SQ(1/k), 1, b/a, MODE)) ) );
}

double lambda_3(double p, double k)
{
	return (1.0/3.0 + ((16.0*p)/(9.0*M_PI))*(2.0*SQ(p)-1.0)*gsl_sf_ellint_Ecomp((1.0/(2.0*k)), MODE) - 
		((32.0*SQ(p)-20.0*p+3.0)/(9.0*M_PI*p))*gsl_sf_ellint_Kcomp((1.0/2.0*k), MODE) );
}

double lambda_4(double p, double k)
{
	return (1.0/3.0 + (2.0/(9.0*M_PI))*( 4.0*(2.0*SQ(p)-1.0)*gsl_sf_ellint_Ecomp(1./k, MODE) +
			(1.0-4.0*SQ(p))*gsl_sf_ellint_Kcomp(1./k, MODE) ) );
}

double lambda_5(double p)
{
	return ( (2.0/(3.0*M_PI))*(acos(1.0-2.0*p)) - (4.0/(9.0*M_PI))*sqrt(p*(1-p))*(3.0+2*p-8.0*SQ(p)) );
}

double lambda_6(double p)
{
	return ( -(2.0/3.0)*pow((1-SQ(p)),(3.0/2.0)) );
}

double eta_1(double z, double p, double a, double b)
{
	return ( (1.0/(2.0*M_PI))*( kappa_1(z,p) + 2.0*eta_2(z, p)*kappa_0(z,p) - 
			(1.0/4.0)*(1+5*SQ(p)+SQ(z))*sqrt( (1.0-a)*(b-1.0) )  ) );
}

double eta_2(double z, double p)
{
	return ( (SQ(p)/2.0)*(SQ(p)+2*SQ(z)) );
}

// double mod(double x)
// {
// //	if (x==0.0) return 0.0;
// 	if(x<0) return (-x);
// 	else return x;
// }

/**************************************************************************************************/
//functions for flux non-lin

double Appell_Hyper(double a,double b1,double b2,double c,double x,double y)
//This can be written in terms of Gauss hypergeometric functions, which gsl
//has implementations of.
//but gsl version quite slow, and does not work for x>1, therefore I've used the scipy version and adapted the C source
//see http://www.eafit.edu.co/revistas/ingenieria-ciencia/Documents/revista10/productInvHyp.pdf for the equations
//also looked at mpmath.py module - some useful transformations etc, but can probably still be significantly
//speeded up with more analytic continuations/transformations for certain ranges of inputs
{
	
	//int m = 0;//, m_max = 90;//, n = 0;
	double ind_term,sum = 0.;

	//printf("Appell hyper (%lf,%lf,%lf,%lf,%lf,%lf):\n", a,b1,b2,c,x,y); fflush(stdout);
	
	//swap terms so that x<y - the outer loop should have a smaller value for greater accuracy and stability
	if(fabs(x)>fabs(y)){
		swap_double(&x,&y); //ie x smaller
		swap_double(&b1,&b2);
		}

	if(fabs(x) >= 0.99){ //coord transformation if x >= 0.99 (speeds things up?)
		double u = (x-y)/(x-1);
		if(fabs(u) <= 0.99) //if u > 0.99 ignore this transformation
			return pow(1.-x,-b1) * pow(1.-y,c-a-b2) * Appell_Hyper(c-a,b1,c-b1-b2,c,u,y);
		}

	//calculate appell hypergeometric func...
	//using eqs frp, http://www.eafit.edu.co/revistas/ingenieria-ciencia/Documents/revista10/productInvHyp.pdf
	//eq A.5
	double t = 1.,h; //first term for r=0
	int r=0;
	while(1){
		//sum up the individual terms from r=0->inf
		h = hyp2f1(a+r,b2,c+r,y);
		
		if(isnan(h)){
		    printf("warning in Appell_Hyper: h in sum is nan! breaking from loop early\nCheck precision of light curves!\n");
		    break;
		    }
		  
		ind_term = t*h;
		sum += ind_term;
		
		//break if max iterations or tolerance reached
		if(fabs(ind_term)<ETOL_HYPER && fabs(h)>10.*ETOL_HYPER) break;
		if(r==MAX_ITER_HYPER){
			//printf("warning: appellF1 eval - max iterations reached\n");
			break;}
        
		//increment r and calculate first term in series
		//just need to multiply by these terms to get the new first term (t)
		//(ie don't need to recalculate the pochammer vals)
		r++; t *= x*(a+r-1)*(b1+r-1)/(c+r-1)/r;
		
		}
	
	//return f1
	return sum;
}

void swap_double(double *x1, double *x2)
{
    double t = *x1;
    *x1 = *x2;
    *x2 = t;
}

// double RF(double a, double b)
// //Rising factorial function
// {
// 	return gsl_sf_poch(a,b);	
// }

double NF(int n, double a, double b, double z, double p)
//the 'N' function
{
	double N;
	
// 	printf("input pars = %d %lf %lf %lf %lf\n",n,a,b,z,p);
// 	printf("gsl_sf_hyperg_2F1 pars = 0.5 0.5 %lf %lf\n", (n+10.)/4.,(1.-a)/(b-a));
// 	printf("appell = %.10f\n", Appell_Hyper(0.5,1.,0.5,(n+10.)/4.,(a-1.)/a,(1.-a)/(b-a)));

	N = pow(1.-a,(n+6.)/4.) / pow(b-a,0.5);
	N *= gsl_sf_beta((n+8.)/4.,0.5);
	N *= ( (SQ(z)-SQ(p))/a*Appell_Hyper(0.5,1.,0.5,(n+10.)/4.,(a-1.)/a,(1.-a)/(b-a))
		- hyp2f1(0.5,0.5,(n+10.)/4.,(1.-a)/(b-a)) );
	
	return N;
}

double MF(int n, double a, double b, double z, double p)
//the 'M' function
{
	double M;

	M = pow((1.-a),(n+4.)/4.);
	M *= ( (SQ(z)-SQ(p))/a*Appell_Hyper(0.5,-(n+4.)/4.,1.,1.,(b-a)/(1.-a),(a-b)/a)
		- hyp2f1(-(n+4.)/4.,0.5,1.,(b-a)/(1.-a)) );
	
//	printf("Appell test (%lf,%lf,%lf,%lf,%lf,%lf) = %lf\n",0.5,-(n+4.)/4.,1.,1.,(b-a)/(1.-a),(a-b)/a,Appell_Hyper(0.5,-(n+4.)/4.,1.,1.,(b-a)/(1.-a),(a-b)/a));
//	printf("gsl_sf_hyperg_2F1 test (%lf,%lf,%lf,%lf) = %lf\n",-(n+4.)/4.,0.5,1.,(b-a)/(1.-a),gsl_sf_hyperg_2F1(-(n+4.)/4.,0.5,1.,(b-a)/(1.-a)));	
	
	return M;
}

double Omega(double * C)
{
	double Om = 0.;
	int n;

	for(n=0.;n<5;n++) Om += C[n] / (n+4.);
	
	return Om;
}

// double scipy_hyp2f1(double a, double b, double c, double x)
// //call scipy hyp2f1 directly from python
// //assumes python already initialised
// //works but too slow - ridiculous overheads when repeatedly calling python functions
// {
//     PyObject *pName, *pModule, *pFunc; //*pDict
//     PyObject *pArgs, *pValue;
//     double result;
// 	
// 	//start the python interpreter
// //    Py_Initialize();
// //    printf("python initialised...\n");
//     
//     //set the python module name
//     pName = PyString_FromString("scipy.special");
//     /* Error checking of pName left out */
// 	
// 	//import the module
//     pModule = PyImport_Import(pName);
//     Py_DECREF(pName); //decref
// 	
// 	//check the module is imported ok
//     if (pModule == NULL) {
//     	printf("error: module not imported.\n"); fflush(stdout);
//     	exit(1);}
// //    printf("module imported...\n");
//     
//     //define the function
//     pFunc = PyObject_GetAttrString(pModule, "hyp2f1");
//     /* pFunc is a new reference */
// //    printf("func defined...\n");
// 
// 	//check function is callable
//     if (!pFunc || !PyCallable_Check(pFunc)){
//     	printf("error: function not callable/available.\n"); fflush(stdout);
//     	exit(1);
//     	}
//     		
// 	//create the args as a python tuple
// 	pArgs = Py_BuildValue("dddd", a,b,c,x);
// 	
// 	//call the function
// 	pValue = PyObject_CallObject(pFunc, pArgs);
// 		
// 	//check the function call is ok
// 	if (pValue == NULL){
// 		Py_DECREF(pFunc);
// 		Py_DECREF(pModule);
//            printf("error: call failed.\n"); fflush(stdout);
// 		exit(1);
// 		}
// 		
// 	//convert answer back to c double
// 	result = PyFloat_AsDouble(pValue);
// 		
// 	//print result
// //	printf("result = %lf\n", result);
//             
//     //decref the arguments
// 	Py_DECREF(pArgs);
//                
// 	//decrease the ref for func and module        
//     Py_XDECREF(pFunc);
//     Py_DECREF(pModule);
// 
// 	//finish python interpreter
// //    Py_Finalize();
//     
//     return result;
// }

/**************************************************************************************************/
