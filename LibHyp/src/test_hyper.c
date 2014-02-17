#include <stdio.h>

//#include "protos.h"
//double hyp2f1 ( double a, double b, double c, double x );

#include "mconf.h"

int main()
{

	double a,b,c,x;
	
	a = 0.1;
	b = 0.1;
	c = 0.1;
	x = 0.1;
	
	for(x=0.;x<1.8;x+=0.01)	printf("Hyper2F1(%lf,%lf,%lf,%lf), res = %lf\n", a,b,c,x,hyp2f1(a,b,c,x));
	
	return 1;
	
}