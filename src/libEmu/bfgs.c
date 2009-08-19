#include "bfgs.h"
void print_vec(gsl_vector* x, int n);

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
	double stepsize = 1.0E-10;
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
void obtainStep(gsl_vector* step, double (*fn)(gsl_vector*, int), \
								void (*gradientFn)( double (*fn)(gsl_vector*, int), gsl_vector*, gsl_vector*, int), \
								gsl_vector* xk, int nparams, gsl_matrix* bkInv){
	double norm = 0.0;
	int i;
	gsl_vector* gradient = gsl_vector_alloc(nparams);
	
	// calculate the gradient of the function
	//getGradientNumeric(fn, xk, gradient, nparams);
	gradientFn(fn, xk, gradient, nparams);
	
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
void getYk(gsl_vector* yk, double(*fn)(gsl_vector*, int),\
					 void (*gradientFn)( double (*fn)(gsl_vector*, int), gsl_vector*, gsl_vector*, int), \
					 gsl_vector* xk, gsl_vector* xk1, int nparams){
	int i;
	gsl_vector* gradXk = gsl_vector_alloc(nparams);
	gsl_vector* gradXk1 = gsl_vector_alloc(nparams);
	
	//getGradientNumeric(fn, xk, gradXk, nparams);
	gradientFn(fn, xk, gradXk, nparams);
	//getGradientNumeric(fn, xk1, gradXk1, nparams);
	gradientFn(fn, xk1, gradXk1, nparams);
	
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
double lineSearch(double(*fn)(gsl_vector*, int),\
					 void (*gradientFn)( double (*fn)(gsl_vector*, int), gsl_vector*, gsl_vector*, int), \
									gsl_vector* direction, gsl_vector* position, int nparams){
	double a = 2.0; // the max stepsize
	double tau = 0.35; // the shrink factor
	while(armGold(fn, gradientFn, direction, position, a, nparams) == 0){
		a = tau *a;
		//printf("%g\n", a); 
		if( a < 1E-10){
			printf("a is very very small: %g\n", a);
			// just use it anyway!
			return(a);
		}
	}
	return(a);
}


//! check to see if the stepsize is too long or not
// this is WAY simpler in mathematica
int armGold(double(*fn)(gsl_vector*, int), \
					 void (*gradientFn)( double (*fn)(gsl_vector*, int), gsl_vector*, gsl_vector*, int), \
						gsl_vector* direction, gsl_vector* position, double a, int nparams){
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
	
	gsl_vector_set_zero(temp);
	gsl_vector_set_zero(offset);

	// set offset -> x + stepsize*position
	gsl_vector_memcpy(offset, position);
	gsl_vector_memcpy(temp, direction);
	gsl_vector_scale(temp, a);
	gsl_vector_add(offset, temp);

	foffset = fn(offset, nparams);
	fposition = fn(position, nparams);
	
	//getGradientNumeric(fn, position, gradient, nparams);
	gradientFn(fn, position, gradient, nparams);
	gsl_blas_ddot(direction, gradient, &dirdotgradient);	
	//getGradientNumeric(fn, offset, gradientOffset, nparams);
	gradientFn(fn, offset, gradientOffset, nparams);
	gsl_blas_ddot(direction, gradientOffset, &dirdotgradientOffset);
	
/* 	if( foffset <= fposition + c1*a*dirdotgradient){ */
/* 		printf("cond1 == TRUE\n"); */
/* 	} */

/* 	if( dirdotgradientOffset >= c2*dirdotgradient){ */
/* 		printf("cond2 == TRUE\n"); */
/* 	} */

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
		
// Binit should be the id matrix 
void doSimpleBFGS( double(*fn)(gsl_vector*, int),\
									 void (*gradientFn)( double (*fn)(gsl_vector*, int), gsl_vector*, gsl_vector*, int), \
									 gsl_vector* xkInit, gsl_vector* xkFinal, gsl_matrix* Binit, int nparams, int nsteps){
	int count = 0;
	int i;
	double stepsize = 0.0;
	double absGrad = 0.0;
	double absGradPrev = 0.0;
	double convergedValue = 1e-8;
	gsl_vector *xk = gsl_vector_alloc(nparams);
	gsl_vector *xk1 = gsl_vector_alloc(nparams);
	gsl_vector *temp = gsl_vector_alloc(nparams);
	gsl_vector *yk = gsl_vector_alloc(nparams);
	gsl_vector *step = gsl_vector_alloc(nparams);
	gsl_matrix *bk = gsl_matrix_alloc(nparams, nparams);
	gsl_matrix *bkInv = gsl_matrix_alloc(nparams, nparams);
	
	gsl_vector_set_zero(xk1);
	gsl_vector_set_zero(yk);
	gsl_vector_set_zero(step);
	gsl_vector_set_zero(temp);

	// set xk to the initial value, and bk
	gsl_vector_memcpy(xk, xkInit);

	gsl_matrix_memcpy(bk, Binit);
	// gonna cheat and just copy the init matrix to the inverse too,
	gsl_matrix_memcpy(bkInv, Binit);
	


	while(count < nsteps){

		//void obtainStep(gsl_vector* step, double (*fn)(gsl_vector*, int), gsl_vector* xk, int nparams, gsl_matrix* bkInv){
		obtainStep(step, fn, gradientFn, xk, nparams, bkInv);

		/* for(i = 0; i<nparams;i++) */
/* 			printf("%g ", gsl_vector_get(step, i)); */
/* 		printf("\n"); */
		
		// find the stepsize
		//double lineSearch(double(*fn)(gsl_vector*, int), gsl_vector* direction, gsl_vector* position, int nparams){		
		stepsize = lineSearch(fn, gradientFn, step, xk, nparams);
		

		//printf("%g\n", stepsize);
		
		// find the new position
		gsl_vector_memcpy(xk1, xk);
		gsl_vector_memcpy(temp, step);
		gsl_vector_scale(temp, stepsize);
		gsl_vector_add(xk1, temp);
		//void getYk(gsl_vector* yk, double(*fn)(gsl_vector*, int), gsl_vector* xk, gsl_vector* xk1, int nparams){	
		getYk(yk, fn, gradientFn, xk, xk1, nparams);
		// find the new Bk
		// void getNewB(gsl_matrix* b, gsl_vector* s, gsl_vector* y, int nparams){		
		getNewB(bk, step, yk, nparams);
		// and find the new approximate inverse
		//void getNewBInverse(gsl_matrix *bprev,  gsl_vector *s, gsl_vector *y, int nparams){
		getNewBInverse(bkInv, step, yk, nparams);

		

		gsl_vector_set_zero(temp);
		print_vec(xk1, nparams);
		
		count++;

		// copy the new step back into xk
		gsl_vector_memcpy(xk, xk1);

		// for convergence we monitor the difference between the current gradient 
		// and that of the last step, if they get small enough we're stuck 
		// void getGradientNumeric(double(*fn)(gsl_vector*, int), gsl_vector* xk, gsl_vector* gradient, int nparams){		
		gradientFn(fn, xk1, temp, nparams);
		for(i = 0; i < nparams; i++)
			absGrad += pow(gsl_vector_get(temp, i),2.0);
		absGrad = sqrt(absGrad);
		
		if(fabs(absGrad - absGradPrev) < convergedValue){
			fprintf(stderr, "converged after %d steps", count);
			break;
		}
		
		// save the gradient as the previous one
		absGradPrev = absGrad;

	}
	
	
	gsl_vector_memcpy(xkFinal, xk1);

	
	gsl_vector_free(step);
	gsl_vector_free(xk);
	gsl_vector_free(xk1);
	gsl_vector_free(temp);
	gsl_vector_free(yk);
	gsl_matrix_free(bk);
	gsl_matrix_free(bkInv);
}



// just to test
#ifdef EXECUTE
int main (void){
	int nparams = 2.0;
	int i;

	double stepsize = 0.0;
	gsl_vector* xTest = gsl_vector_alloc(nparams);
	gsl_vector* xFinal = gsl_vector_alloc(nparams);
	gsl_matrix* Binit = gsl_matrix_alloc(nparams, nparams);

	gsl_vector_set(xTest, 0, 0.3);
	gsl_vector_set(xTest, 1, 0.35);

	gsl_matrix_set_identity(Binit);

	//int armGold(double(*fn)(gsl_vector*, int), gsl_vector* direction, gsl_vector* position, double a, int nparams){	
	//printf("%d\n", armGold( &fRosenbrock, step, xTest, 0.02, nparams));

	//printf("calling BFGS\n");

	//void doSimpleBFGS( double(*fn)(gsl_vector*, int), gsl_vector* xkInit, gsl_vector* xkFinal, gsl_matrix* Binit, int nparams, int nsteps){
	doSimpleBFGS(&fRosenbrock, &getGradientNumeric, xTest, xFinal, Binit, nparams, 100000);
	print_vec(xFinal, nparams);
	
	gsl_matrix_free(Binit);
	gsl_vector_free(xTest);
	gsl_vector_free(xFinal);
	return(0);
}

void print_vec(gsl_vector* x, int n){
	int i;
	for(i =0; i < n;i++){
		fprintf(stdout, "%g ", gsl_vector_get(x, i));
	}
	fprintf(stdout, "\n");
}
#endif
