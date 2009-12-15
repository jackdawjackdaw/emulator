#include "stdio.h"
#include "math.h"
#include "gsl/gsl_rng.h"


/* generate a gaussian or more in 2d and smush them together into 
 * some overall function, then try to emulate them 
 */

double gauss(double x, double y, double x0, double y0, double a, double b, double c){
	double f = a*pow(x-x0,2.0) + 2*b*(x-x0)*(y-y0) + c*pow(y-y0,2.0);
	return(exp(-f));
}

const gsl_rng_type *T;
gsl_rng *r;

int main (void){
	// the density of points to evaluate 
	int nx = 16;
	int ny = 16;
	double world[nx][ny];
	int i,j;
	 
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	
	double x =0;
	double y = 0;
	double m1x = 0.5;
	double m1y = 0.5;
	double a1 = 100;
	double a2 = 0.0;
	double a3 = 100;


	for (i = 0; i < nx; i++){
		for(j= 0; j < ny; j++){
			world[i][j] = 0;
			x = ((double) i)/((double) nx);
			y = ((double) j)/((double) nx);
			world[i][j] += gauss(x,y, m1x, m1y, a1,a2,a3);
			//world[i][j] += gauss(x,y, 0.1, 0.1, 75,0,75);
		}
	}
	
	for (i = 0; i < nx; i++){
		for(j= 0; j < ny; j++){
			x = ((double) i)/((double) nx);
			y = ((double) j)/((double) nx);
			if(world[i][j] > 1E-5){
				fprintf(stdout, "%g\t%g\t%g\n", x,y, world[i][j]);
			} else {
				fprintf(stdout, "%g\t%g\t%g\n", x,y, 0.01*gsl_rng_uniform(r));
			}
		}
	}

	gsl_rng_free(r);
	return(0);
}
	
	
