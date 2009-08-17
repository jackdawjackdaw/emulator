//! A C implementation of the BFGS algorithm as sketched out in mathematica
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "assert.h"

/** 
 * This should provide an additional maximisation method for libEmu
 */

/**
 * function pointers
 * float fn (float a, float b)
 * 
 * void foo( float (*ptr)(float, float)){
 * result = ptr(a, b);
 * }
 *
 * foo(&fn);
 */

/* test functions */
double fRosenbrock( gsl_vector* x, int nparams); 
double fSphere(gsl_vector *x, int nparams);
double fPowers( gsl_vector*x, int nparams);

/* actual bfgs functions */
void getGradientNumeric(double(*fn)(gsl_vector*, int), gsl_vector* xk, gsl_vector* gradient, int nparams);
void obtainStep(gsl_vector* step, double (*fn)(gsl_vector*, int), gsl_vector* xk, int nparams, gsl_matrix* bkInv);
void getYk(gsl_vector* yk, double(*fn)(gsl_vector*, int), gsl_vector* xk, gsl_vector* xk1, int nparams);
void getNewBInverse(gsl_matrix *bprev,  gsl_vector *s, gsl_vector *y, int nparams);
void getNewB(gsl_matrix* b, gsl_vector* s, gsl_vector* y, int nparams);
double lineSearch(double(*fn)(gsl_vector*, int), gsl_vector* direction, gsl_vector* position, int nparams);
int armGold(double(*fn)(gsl_vector*, int), gsl_vector* direction, gsl_vector* position, double a, int nparams);

double fRosenbrock( gsl_vector* x, int nparams){
	double xp, yp;
	assert(nparams == 2);
	xp = gsl_vector_get(x, 0);
	yp = gsl_vector_get(x, 1);
	return(pow((1-xp), 2.0) + 100.0*pow((yp-pow(xp, 2.0)), 2.0));
}

double fSphere(gsl_vector *x, int nparams){
	double result = 0.0;
	int i;
	for(i = 0; i < nparams; i++){
		result += pow(gsl_vector_get(x, i), 2.0);
	}
	return(result);
}

// this one works in any order...
double fPowers( gsl_vector*x, int nparams){
	double result = 0.0;
	int i;
	for(i = 0; i < nparams; i++){
		result += fabs(pow(gsl_vector_get(x, i), i));
	}
	return(result);
}

// first order finite diffs, not very accurate
void getGradientNumeric(double(*fn)(gsl_vector*, int), gsl_vector* xk, gsl_vector* gradient, int nparams){
	double stepsize = 1.0E-8;
	int i, j;
	gsl_vector* temp = gsl_vector_alloc(nparams);
	gsl_vector* localGradient = gsl_vector_alloc(nparams);
	double xtemp = 0.0;
	double grad = 0.0;

	for(i = 0;i < nparams; i++){
		gsl_vector_memcpy(temp, xk);
		xtemp = gsl_vector_get(xk, i) + stepsize;
		gsl_vector_set(temp, i, xtemp);
		grad = (1.0/stepsize)*(fn(temp, nparams) - fn(xk, nparams));
		gsl_vector_set(localGradient, i, grad);
	}
	
	gsl_vector_memcpy(gradient, localGradient);	
	gsl_vector_free(localGradient);
	gsl_vector_free(temp);
}
	

/**
 * work out the step vector
 * fn is the thing we're trying to optimise
 * xk is the current position in the parameter space,
 * step is a vector of the direction we want to move in,
 * nparams is the number of dimensions in step and xk
 * bkInv is the approximate inverse of the Hessian matrix
 */
void obtainStep(gsl_vector* step, double (*fn)(gsl_vector*, int), gsl_vector* xk, int nparams, gsl_matrix* bkInv){
	double norm = 0.0;
	int i;
	gsl_vector* gradient = gsl_vector_alloc(nparams);
	
	// calculate the gradient of the function
	getGradientNumeric(fn, xk, gradient, nparams);
	
	// step = -Binv.grad
	gsl_blas_dgemv(CblasNoTrans, -1.0, bkInv, gradient, 0.0, step);

	// now we should normalise the step
	for(i = 0; i < nparams; i++){
		norm += pow(gsl_vector_get(step, i), 2.0);
	}
	norm = (1.0)/(sqrt(norm));

	for (i = 0; i < nparams; i++)
		gsl_vector_set(step, i, norm*gsl_vector_get(step, i)); 

	gsl_vector_free(gradient);
}

//! get the change in gradient between two locations
void getYk(gsl_vector* yk, double(*fn)(gsl_vector*, int), gsl_vector* xk, gsl_vector* xk1, int nparams){
	int i;
	gsl_vector* gradXk = gsl_vector_alloc(nparams);
	gsl_vector* gradXk1 = gsl_vector_alloc(nparams);
	
	getGradientNumeric(fn, xk, gradXk, nparams);
	getGradientNumeric(fn, xk1, gradXk1, nparams);
	
	for(i = 0; i < nparams; i++)
		gsl_vector_set(yk, i, gsl_vector_get(gradXk1, i) - gsl_vector_get(gradXk, i));
 
	free(gradXk);
	free(gradXk1);
}

//! calculates a new hessian
/**
 * calculates a new hessian using the previous one, b, the stepsize s, and the change in gradient y
 */ 
void getNewB(gsl_matrix* b, gsl_vector* s, gsl_vector* y, int nparams){
	double product1 = 0.0;
	double product2 = 0.0;
	double temp1, temp2;
	gsl_vector * temp_vec = gsl_vector_alloc(nparams);
	// bnew = b + (y.y/y.s) - ((B.s.(B.s))/(s.B.s))
	//                   ^ product 1
	//                                  ^ product 2
	gsl_blas_ddot(y, y, &temp1);
	gsl_blas_ddot(y, s, &temp2);	
	product1= temp1/temp2;
	
	// calculate the second part
	gsl_blas_dgemv(CblasNoTrans, 1.0, b, s, 0.0, temp_vec);
	// topline of product 2
	gsl_blas_ddot(temp_vec, temp_vec, &temp1);
	gsl_blas_ddot(s, temp_vec, &temp2);
	product2= temp1/temp2;

	// now b is set to bnew
	gsl_matrix_add_constant(b, product1+product2);
	gsl_vector_free(temp_vec);
}

//! approximates the Inverse  hessian
/** 
 * approximate a new inverse hessian	
 * @return bprev is set to the new value 
 */
void getNewBInverse(gsl_matrix *bprev,  gsl_vector *s, gsl_vector *y, int nparams){
	double sdoty = 0.0;
	double sdots = 0.0;
	double ydotbpdoty = 0.0;
	double temp = 0.0;
	gsl_vector * temp_vec = gsl_vector_alloc(nparams);
	
	//b = bp + (s.s)*(s.y - y.bp.y)/((s.y)^2) 
	gsl_blas_ddot(s, y, &sdoty);
	gsl_blas_ddot(s, s, &sdots);
	gsl_blas_dgemv(CblasNoTrans, 1.0, bprev, s, 0.0, temp_vec); 
	gsl_blas_ddot(s, temp_vec, &ydotbpdoty);
	if(sdoty != 0.0){
		temp = sdots*(sdoty - ydotbpdoty) / (pow(sdoty, 2.0));
		// set the answer
		gsl_matrix_add_constant(bprev, sdots); 
	} else {
		// this is likely to be annoying
		fprintf(stderr, "getNewBinverse had a singluar sdoty\n");
		// don't do anything
	}
	
}	
		
//! do a backtracking line search
double lineSearch(double(*fn)(gsl_vector*, int), gsl_vector* direction, gsl_vector* position, int nparams){
	double a = 5.0; // the max stepsize
	double tau = 0.85; // the shrink factor
	while(armGold(fn, direction, position, a, nparams) == 0){
		a = tau *a;
	}
	return(a);
}


//! check to see if the stepsize is too long or not
// this is WAY simpler in mathematica
int armGold(double(*fn)(gsl_vector*, int), gsl_vector* direction, gsl_vector* position, double a, int nparams){
	double c1 = 1E-3;
	double c2 = 0.9;
	double foffset = 0.0;
	double fposition = 0.0;
	double dirdotgradient  = 0.0;
	double dirdotgradientOffset  = 0.0;
	int retval = 0;
	gsl_vector *temp = gsl_vector_alloc(nparams);
	gsl_vector *offset = gsl_vector_alloc(nparams);
	gsl_vector *gradient = gsl_vector_alloc(nparams);
	gsl_vector *gradientOffset = gsl_vector_alloc(nparams);
	

	// set offset -> x + stepsize*position
	gsl_vector_memcpy(offset, position);
	gsl_vector_scale(temp, a);
	gsl_vector_add(offset, temp);

	foffset = fn(offset, nparams);
	fposition = fn(position, nparams);
	
	getGradientNumeric(fn, position, gradient, nparams);
	gsl_blas_ddot(direction, gradient, &dirdotgradient);	
	getGradientNumeric(fn, offset, gradientOffset, nparams);
	gsl_blas_ddot(direction, gradientOffset, &dirdotgradientOffset);

	if( (foffset <= fposition + c1*a*dirdotgradient)  && (dirdotgradientOffset >= c2 * dirdotgradient)){
		retval = 1;
	} else {
		retval = 0;
	}
	
	gsl_vector_free(gradientOffset);
	gsl_vector_free(gradient);
	gsl_vector_free(temp);
	gsl_vector_free(offset);
	return(retval);
}
		
void doSimpleBFGS( double(*fn)(gsl_vector*, int), gsl_vector* xkInit, gsl_matrix* Binit, int nparams, int nsteps){
	int count = 0;
	gsl_vector *xk = gsl_vector_



// just to test
int main (void){
	int nparams = 2.0;
	int i;
	gsl_vector* gradient = gsl_vector_alloc(nparams);
	gsl_vector* xTest = gsl_vector_alloc(nparams);

	gsl_vector_set(gradient, 0, 0.0);
	gsl_vector_set(gradient, 1, 0.0);

	gsl_vector_set(xTest, 0, 1.3);
	gsl_vector_set(xTest, 1, 1.5);

	getGradientNumeric(&fRosenbrock, xTest, gradient, nparams);

	for(i = 0; i < nparams; i++)
		fprintf(stderr, "%g\n", gsl_vector_get(gradient, i));


	gsl_vector_free(gradient);
	gsl_vector_free(xTest);
	return(0);
}
