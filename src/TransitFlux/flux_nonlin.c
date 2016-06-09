/*********************************************************************************
Neale Gibson Aug 2011 flux_nonlin.c
This function uses the analytic equations of Mandel and Agol (2000) to compute the
flux of a star during occultation given a normalized separation z, ratio of planet
to star radii, and the four non-linear limb darkening coefficients.
*********************************************************************************/

#include "flux.h"

double flux_nonlin(double z, double p, double c1, double c2, double c3, double c4)
{
	#ifdef DEBUG
	printf("FLUX NONLIN INPUT: z = %lf p = %lf c1 = %lf c2 = %lf c3 = %lf c4 = %lf\n", z, p, c1, c2, c3, c4); fflush(stdout);
	#endif
	
	if (p <= 0.0) { printf("p <= 0 for flux calculation !!\n"); exit(-1); }

	//define commonly used parameters
	double flux=0., sum=0.;
	int n;
	
	//define ld params and c0 in an array
	double C[] = {1.-c1-c2-c3-c4,c1,c2,c3,c4};
	
	//now need to decide the situation based on p and z, then work out the flux accordingly
	int Case = -1;
	
	//make sure p is positive
	if (p <= 0.0){
		printf("p <= 0 for flux calculation !!\n"); exit(-1);
		}
	
	//Case 1 - planet does not occult star
	if(z >= (1.+p)){
		Case = 1;
		#ifdef DEBUG
		printf("FLUX NONLIN: Case %d\n", Case);
		#endif
		flux = 1.;
		}

	//Cases 4 and 5 upper and lower limits of 3, respectively
	//Case 4 - planet entirely within star, but does not touch/cross the centre, and touches limb
	//same as Case 3 for nonlin law
	else if( p < 0.5 && DOUBEQ(z,(1.-p)) ){
		Case = 4;
		#ifdef DEBUG
		printf("FLUX NONLIN: Case %d\n", Case);
		#endif
		
		//calculate sum term
		for(n=1,sum=0.;n<4;n++){
			sum += MF(n,SQ(z-p),SQ(z+p),z,p) * C[n] / (n+4.);
			}

		double L = SQ(p) * (1.-SQ(p)/2.-SQ(z));

		flux = 1. - pow(4.*Omega(C),-1.) * ( (C[0] * SQ(p)) + 2.*sum + C[4] * L );
		}
	
	//Case 5 - planet entirely within star, touches centre of star
	else if( p < 0.5 && DOUBEQ(z,p)){
		Case = 5;
		#ifdef DEBUG
		printf("FLUX NONLIN: Case %d\n", Case);
		#endif

		for(n=0,sum=0.;n<5;n++){
			sum += C[n] / (n+4.) * hyp2f1(0.5,-(n+4.)/4.,1.,4.*SQ(p));
			}

		flux = 0.5 + pow(2.*Omega(C),-1.) * sum;
		}
	
	//Case 6 - planet touches star centre and disk, therefore p = 0.5
	else if( DOUBEQ(p,0.5) && DOUBEQ(z,0.5) ){
		Case = 6;
		#ifdef DEBUG
		printf("FLUX NONLIN: Case %d\n", Case);
		#endif
		
		for(n=0,sum=0.;n<5;n++){
			sum += C[n]/(n+4.) * gsl_sf_gamma(1.5+n/4.) / gsl_sf_gamma(2.+n/4.);
			}
		
		flux = 0.5 + 1./(2.*M_PI*Omega(C)) * sum;
		
		flux = 0.;
		}
	
	//Case 7 - planet touches stellar centre, but not entirely within disk
	else if( p > 0.5 && DOUBEQ(z,p) ){
		Case = 7;
		#ifdef DEBUG
		printf("FLUX NONLIN: Case %d\n", Case);
		#endif
		
		for(n=0,sum=0.;n<5;n++){
			sum += C[n]/(n+4.) * gsl_sf_beta(0.5,(n+8.)/4.) * hyp2f1(0.5,0.5,2.5+n/4.,1./(4.*SQ(p)));
			}
				
		flux = 0.5 + 1./(4.*p*M_PI*Omega(C)) * sum;
		
		flux = 0.;
		}
	
	//Case 2 - planet on limb of star but doesn't touch centre
	else if( z > (0.5+fabs(p-0.5)) && z < (1.+p)){
		Case = 2;
		#ifdef DEBUG
		printf("FLUX NONLIN: Case %d\n", Case);
		#endif
	
		for(n=0,sum=0.;n<5;n++){
			sum += NF(n,SQ(z-p),SQ(z+p),z,p) * C[n] / (n+4.);
		  #ifdef DEBUG
		  printf("NF = %.15lf\n", NF(n,SQ(z-p),SQ(z+p),z,p));
		  #endif
			}
		#ifdef DEBUG
		printf("sum = %.15lf\n", sum);
		#endif
	
		flux = 1. - pow(2*M_PI*Omega(C),-1.) * sum;

// 		#ifdef DEBUG
// 		printf("FLUX NONLIN: Case %d\n", Case);
// 		#endif
		}

	//Case 3 - planet entirely within star, but does not touch/cross the centre, or touch limb
	else if( p < 0.5 && z > p && z < (1.-p)){
		Case = 3;
		#ifdef DEBUG
		printf("FLUX NONLIN: Case %d\n", Case);
		#endif

		for(n=1,sum=0.;n<4;n++){
			sum += MF(n,SQ(z-p),SQ(z+p),z,p) * C[n] / (n+4.);
			}

		double L = SQ(p) * (1.-SQ(p)/2.-SQ(z));

		flux = 1. - pow(4.*Omega(C),-1.) * ( C[0]*SQ(p) + 2.*sum + C[4] * L );
		}
		
	//Case 8 - planet covers the centre and limb of the star
	else if( p > 0.5 && z >= fabs(1.-p) && z < p){
		Case = 8;
		#ifdef DEBUG
		printf("FLUX NONLIN: Case %d\n", Case);
		#endif
		
		for(n=0,sum=0.;n<5;n++){
			sum += C[n] * NF(n,SQ(z-p),SQ(z+p),z,p) * pow(n+4.,-1.);
			}
		
		flux = -1./(2.*M_PI*Omega(C)) * sum;
		}

	//Case 9 - planet entirely withing star, and covers the centre
	//if( (p < 1.) && (z > 0.) && (z < (0.5 - abs(p-0.5)))){
	//had to change this from MA02 paper - it allows situations with p<0.5 and the centre is not covered?
	else if( (p < 1.) && (z > 0.) && (z < p) ){
		Case = 9;
		#ifdef DEBUG
		printf("FLUX NONLIN: Case %d\n", Case);
		#endif

		//define L
		double L = SQ(p) * (1.-SQ(p)/2.-SQ(z));

		for(n=1,sum=0.;n<4;n++){
			sum += (C[n] / (n+4.) * MF(n,SQ(z-p),SQ(z+p),z,p));
			}
		
		flux = pow((4.*Omega(C)),-1.) * ( C[0]*(1.-SQ(p)) + C[4]*(0.5-L) -2.*sum );
		}
	
	//Case 10 - planet centred on star, p<1.
	else if( p < 1. && z == 0.){
		Case = 10;
		#ifdef DEBUG
		printf("FLUX NONLIN: Case %d\n", Case);
		#endif

		for(n=0,sum=0.;n<5;n++){
			sum += C[n] * pow(1.-SQ(p),(n+4.)/4.) / (n+4.);
			}

		flux = pow(Omega(C),-1.) * sum;
		}
	
	//Case 11 - `planet' completely occults star
	else if( p >= 1. && z <= p-1.){
		Case = 11;
		#ifdef DEBUG
		printf("FLUX NONLIN: Case %d\n", Case);
		#endif
		
		flux = 0.;
		}
	
	else{
		Case = -1;
		#ifdef DEBUG
		printf("FLUX NONLIN: Case %d\n", Case);
		#endif
		printf("flux nonlin warning: no Case selected!\n");
 		printf("  (z,p,Case) = %.10lf %.10lf %d\n", z,p,Case); fflush(stdout);
		flux = 0.;
		}

	//print some output if DEBUG is on
	#ifdef DEBUG
	printf("FLUX NONLIN OUTPUT: %.15lf\t%.15lf\tcase %d\t%lf\n", z, p, Case, flux);
	#endif
	
// 	if(flux>1.){
// 		printf("warning: flux > 1.\n");
// 		printf("  (z,p,Case) = %.10lf %.10lf %d\n", z,p,Case); fflush(stdout);
// 		}
// 	if(flux<0.2){
// 		printf("warning: flux < 0.2\n");
// 		printf("  (z,p,Case) = %.10lf %.10lf %d\n", z,p,Case); fflush(stdout);
// 		}
			
	return flux;
}
