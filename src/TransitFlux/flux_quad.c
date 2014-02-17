/*********************************************************************************
Neale Gibson 15/3/2007 flux_quad.c
This function uses the analytic equations of Mandel and Agol (2000) to compute the
flux of a star during occultation given a normalized separation z, ratio of planet
to star radii, and the two quadratic limb darkening coefficients gam_1 and gam_2.
*********************************************************************************/

#include "flux.h"

double flux_quad(double z, double p, double gam_1, double gam_2)
{
	#ifdef DEBUG
	printf("z = %lf p = %lf c1 = %lf c2 = %lf\n", z, p, gam_1, gam_2);
	#endif
	
	if (p <= 0.0) { printf("p <= 0 for flux calculation !!\n"); exit(-1); }

	//define and work out a, b and q
	double a = SQ(z-p), b = SQ(z+p), q = (SQ(p)-SQ(z));
	double k = sqrt((1-a)/(4*p*z));

	#ifdef DEBUG
	printf("a = %lf b = %lf q = %lf k = %lf\n", a, b, q, k);	
	#endif
	
	//define and work out limb darkening coeffs
	double c_2 = (gam_1+2.0*gam_2), c_4 = (-gam_2);
	
	//define the main function values
	double lambda_e=0, lambda_d=0, eta_d=0, flux;
	
	int _case = -1;
	//now need to decide the situation based on p and z, then work out the flux accordingly
	
	//for p less than 0.5
	if (p < 0.5) {
		
		//case 1
		if (z >= (1+p) ){
			#ifdef DEBUG
			printf("p < 0.5, case 1\n");
			#endif
			lambda_e = 0;  
			lambda_d = 0;
			eta_d = 0;			
			_case = 1;
			}
		
		//case 4 - still seems to be a bug were case 2 is selected instead when z = p!
		//added the final 'else' to account for this
		//else if (z == (1.00-p)){
		else if (DOUBLE_EQ(z,(1.-p))){
			#ifdef DEBUG
			printf("p < 0.5, case 4\n");
			#endif
			lambda_e = SQ(p);  
			lambda_d = lambda_5(p);
			eta_d = eta_2(z, p);
			_case = 4;
			}
			
		//case 2
		//had to add the .00..1 to make sure the case is not selected instead if case 4
		//else if (z > (1.0000000001-p) && z < (1+p) ){
		else if (z > (1.-p) && z < (1+p) ){
			#ifdef DEBUG
			printf("p < 0.5, case 2\n");
			#endif
			lambda_e = lambda_e_pi_funct(z, p);  
			lambda_d = lambda_1(z, p, a, b, q, k);
			eta_d = eta_1(z,p,a,b);
			_case = 2;
			}
			
		//case 3
		else if (z > p+EPSILON && z < (1-p)){
			#ifdef DEBUG
			printf("p < 0.5, case 3\n");
			#endif
			lambda_e = SQ(p);  
			lambda_d = lambda_2(z, p, a, b, q, k);
			eta_d = eta_2(z,p);
			_case = 3;			
			}

		//case 5
		//else if (z == (p) ){
		else if (DOUBLE_EQ(z,p)){
			#ifdef DEBUG
			printf("p < 0.5, case 5\n");
			#endif
			lambda_e = SQ(p);  
			//k = 0.5/k;
			lambda_d = lambda_4(p, k);
			eta_d = eta_2(z,p);
			_case = 5;		
			}
			
		//case 9
		else if (z > 0 && z < p ){
			#ifdef DEBUG
			printf("p < 0.5, case 9\n");
			#endif
			lambda_e = SQ(p);  
			lambda_d = lambda_2(z, p, a, b, q, k) ;
			eta_d = eta_2(z,p);
			_case = 9;					
			}
		
		//case 10
		//else if (z == 0){
		else if (DOUBLE_EQ(z,0.)){
			#ifdef DEBUG
			printf("p < 0.5, case 10\n");
			#endif
			lambda_e = SQ(p);  
			lambda_d = lambda_6(p);
			eta_d = eta_2(z,p);
			_case = 10;						
			}
			
		//case 4 again!
		else {
			#ifdef DEBUG
			printf("p < 0.5, case 4 (again!)\n");
			#endif
			lambda_e = SQ(p);  
			lambda_d = lambda_5(p);
			eta_d = eta_2(z, p);
			_case = 4;
			}
	}
	
	//now add in for p >= 0.5 situations
	else{
	//case 1
		if (z >= (1+p) ){
			#ifdef DEBUG
			printf("p >= 0.5, case 1\n");
			#endif
			lambda_e = 0;  
			lambda_d = 0;
			eta_d = 0;			
			_case = 1;
			}
		
		//case 2
		//had to add the .00..1 to make sure the case is not selected instead if case 4
		else if (z > (p) && z < (1+p) ){
			#ifdef DEBUG
			printf("p >= 0.5, case 2\n");
			#endif
			lambda_e = lambda_e_pi_funct(z, p);  
			lambda_d = lambda_1(z, p, a, b, q, k);
			eta_d = eta_1(z,p,a,b);
			_case = 2;
			}
			
		//case 6
		else if (z == 0.5 && p == 0.5){
			#ifdef DEBUG
			printf("p >= 0.5, case 6\n");
			#endif
			lambda_e = SQ(p);  
			lambda_d = 1.0/3.0 - 4.0/(9.0*M_PI);
			eta_d = 3.0/32.0;
			_case = 6;			
			}

		//case 7
		else if (z == p){	
			#ifdef DEBUG
			printf("p >= 0.5, case 7\n");
			#endif
			lambda_e = lambda_e_pi_funct(z, p);  
			lambda_d = lambda_3(p, k);
			eta_d = eta_1(z,p,a,b);
			_case = 7;			
			}
		
		//case 8
		else if (z > fabs(1-p) && z < p){
			#ifdef DEBUG
			printf("p >= 0.5, case 8\n");
			#endif
			lambda_e = lambda_e_pi_funct(z, p);  
			lambda_d = lambda_1(z, p, a, b, q, k);
			eta_d = eta_1(z,p,a,b);
			_case = 8;		
			}
		
		//case 9
		else if (z > 0 && z < (1-p) ){
			#ifdef DEBUG
			printf("p >= 0.5, case 9\n");
			#endif
			lambda_e = SQ(p);  
			lambda_d = lambda_2(z, p, a, b, q, k) ;
			eta_d = eta_2(z,p);
			_case = 9;					
			}
		
		//case 10
		else if (z == 0 && p<1.0){
			#ifdef DEBUG
			printf("p >= 0.5, case 10\n");
			#endif
			lambda_e = SQ(p);  
			lambda_d = lambda_6(p);
			eta_d = eta_2(z,p);
			_case = 10;						
			}
		
		//case 11
		else if (p>1 && z < p-1){
			#ifdef DEBUG
			printf("p >= 0.5, case 11\n");
			#endif
			lambda_e = 1;  
			lambda_d = 1;
			eta_d = 1;
			flux = 0;
			_case = 11;	
			//printf("%lf\t%lf\tcase %d\t%.10lf %.10lf %.10lf\t%lf\n", z, p, _case, lambda_e, lambda_d, eta_d, flux);
			return flux;					
			}

	}
	
	//this accounts for the mysterious theta function
	if(z < p) lambda_d += 2.0/3.0;
	
	//calculate the flux now that the lambda and eta functions are evaluated
	flux = 1.0 - (1/(1 - gam_1/3 - gam_2/6)) * ( (1-c_2)*(lambda_e)
			 + c_2*(lambda_d)// + 2/3)
			 - c_4*(eta_d) );
	
	//uniform source output flux
	//flux = 1.0 - lambda_e;	
	
	//print some output if DEBUG is on
	#ifdef DEBUG			
	printf("%.15lf\t%.15lf\tcase %d\t%.15lf %.15lf %.15lf\t%lf\n", z, p, _case, lambda_e, lambda_d, eta_d, flux);
	#endif
	
	return flux;
		
}

