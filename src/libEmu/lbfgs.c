#include "lbfgs.h"


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

double fRosenbrock( double* x, int nparams, void* args){
	double xp, yp;
	assert(nparams == 2);
	xp = x[0];
	yp = x[1];
	return(pow((1-xp), 2.0) + 100.0*pow((yp-pow(xp, 2.0)), 2.0));
}

// don't need the args field here
double fSphere(double *x, int nparams, void *args){
	double result = 0.0;
	int i;
	for(i = 0; i < nparams; i++){
		result += pow(x[i], 2.0);
	}
	return(result);
}

// this one works in any order...
double fPowers( double *x, int nparams, void *args){
	double result = 0.0;
	int i;
	for(i = 0; i < nparams; i++){
		result += fabs(pow(x[i], (double)i+1));
	}

	return(result);
}


// adpated from the bfgs.c version
void getGradientNumericLBFGS(double(*fn)(double*, int, void*), double* xk, double* gradient, int nparams, void* args){
	double stepsize = 1.0E-10;
	int i;
	int j;
	double* temp = malloc(sizeof(double)*nparams);
	double* localGradient = malloc(sizeof(double)*nparams);
	double xtemp = 0.0;
	double grad = 0.0;

	for(i = 0;i < nparams; i++){
		//gsl_vector_memcpy(temp, xk);
		//temp = memcpy(temp, xk, nparams);
		for(j = 0; j < nparams; j++){
			temp[j] = xk[j];
		}
		xtemp = xk[i] + stepsize;
		temp[i] = xtemp;
		grad = (1.0/stepsize)*(fn(temp, nparams, args) - fn(xk, nparams, args));
		localGradient[i] = grad;
	}

	for(j = 0; j < nparams; j++){
		gradient[j] = localGradient[j];
	}

	//gradient = memcpy(gradient, localGradient, nparams);
	free(localGradient);
	free(temp);
}



void doBoundedBFGS( double(*fn)(double*, int, void*),													\
										void(*gradientFn)(double(*fn)(double*, int, void*), double*, double*, int, void*), \
										gsl_matrix* ranges, 
										gsl_vector *xkInit, gsl_vector* xkFinal, int nparams, int nsteps, void* args){
	int count = 0;
	int go_flag = 1;
	int i;
	// these are set by the fortran routine
	int nmax = 1024;
	int mmax = 17;
	int memsize = 5; // how many corrections to store (5 in driver1 but can be up to mmax)
	int wasize = 2*mmax*nmax+4*nmax + 12*mmax*mmax + 12*mmax; // see driver1.f
	int *nbd = malloc(sizeof(int)*nparams);
	int *iwa = malloc(sizeof(int)*3*nmax);
	int isave[44];
	long int lsave[4];
	
	int stringlength = 61;
	char task[61], csave[61];
	double fnval; // the value of the fn, passed to the fortran routine
	double *wa = 	malloc(sizeof(double)*wasize);
	double *grad = malloc(sizeof(double)*nparams); // holds the gradient
	double *lower = malloc(sizeof(double)*nparams); // holds the lower bounds
	double *upper = malloc(sizeof(double)*nparams); // holds the upper
	double *xvalue = malloc(sizeof(double)*nparams); // holds the current x value

	double dsave[29]; // don't know what this does yet
	int iprint = -1; // turn off output
	/**
	 * factor and gradtol control the accuracy level of the algorithm, 
	 * you can read the lbfgs source/paper but it's probably not worth worrying about
	 *
	 */
	double factor = 1E7;
	double gradtol = 1E-12;

	/**
	 * the fortran code needs it's "internal memory" setup very carefully or it'll 
	 * barf, of particular importance is the task string which is used for flow control.
	 * 
	 * don't mess with this unless something is really messed up.
	 */
	
	set_zero(grad, nparams);
	set_zero(lower,nparams);
	set_zero(upper,nparams);
	//set_zero(wa, wasize);
	/* set_zero(dsave, 29); */
	for(i = 0; i < nparams; i++) nbd[i] = 0;
	for(i = 0; i < nparams; i++) xvalue[i] = 0.0;
	for(i = 0; i < 3*nmax; i++) iwa[i] = 0;
	for(i = 0; i < wasize; i++) wa[i] = 0.0;
	for(i = 0; i < 29; i++) dsave[i] = 0.0;
	for(i = 0; i < 44; i++) isave[i] = 0;
	for(i = 0; i < 4; i++) lsave[i] = 0;
	for(i = 0; i < stringlength; i++){
		task[i] = ' ';
		csave[i] = ' ';
	}
		


	set_bounds(lower, ranges, 0, nparams);
	set_bounds(upper, ranges, 1, nparams);
	set_nbd(nbd, nparams, nparams); // we have upper and lower so this just sets nbd[i] == 2

	copy_gslvec_vec(xkInit, xvalue, nparams); 

	/* // now we have to set the task string to "START" */
	sprintf(task, "%s", "START"); // this might not be the right way
	for(i = 5; i < stringlength; i++){
		task[i] = ' ';
	}

	/* this is all debug */
	/* printf("starting...\n"); */
	/* printf("wasize = %d\n", wasize); */
 	/* printf("nmax = %d\n", nmax); */
	/* printf("mmax = %d\n", mmax); */
	/* printf("factor = %g\n", factor); */
	/* printf("gradtol = %g\n", gradtol); */
	
	fnval = 0.0;

	while(go_flag == 1){
		
		// call the fortran driver
		setulb_(&nparams, &memsize, xvalue, lower, upper, nbd, &fnval, grad, &factor, &gradtol, wa, iwa, task, &iprint, csave, lsave, isave ,dsave, stringlength, stringlength);

		//printf("%s\n", task);

		//printf("%f %f %f \n", dsave[1], dsave[4], dsave[12]);

		/* printf("xval = "); */
		/* for(i = 0; i < nparams; i++){ */
		/* 	printf("%f\t", xvalue[i]); */
		/* } */
		/* printf("\n"); */


		
		if(strncmp(task, "FG",2) == 0){
			//fprintf(stderr,"task is go!\n");
			/* evaluate the eval-fn at the point xvalue */
			fnval = fn(xvalue, nparams, args);  // this is just for testing, will be a bit more complicated... 
			/* this is perhaps not the best idea */
			/* if(isnan(fnval)){			  */
			/* 	fprintf(stderr, "nan! in lbfgs, stopping\n"); */
			/* 	go_flag = 0; */
			/* 	break; */
			/* } else if(isinf(fnval)){ */
			/* 	fprintf(stderr, "inf in lbfgs, stopping\n"); */
			/* 	go_flag = 0; */
			/* 	break; */
			/* } */
			//fprintf(stderr,"fnval = %g\n", fnval);
			/* evaluate the gradient here */
			gradientFn(fn, xvalue, grad, nparams, args);
			/* fprintf(stderr,"grad = "); */
			/* for(i = 0; i < nparams; i++){ */
			/* 	fprintf(stderr,"%f\t", grad[i]); */
			/* } */
			/* fprintf(stderr,"\n"); */

			// done
		}  else if(strncmp(task,"NEW_X",5) == 0){
			//fprintf(stderr,"%s - new_x\n", task);
			// we have a new iterate and we're going to continue
		} else {
			//fprintf(stderr, "%s", task);
			// we didn't have a new_x and we don't need the grad so
			// this is the end, beautiful friend
			go_flag = 0;
		}
		count++;
		//fprintf(stderr, "%s\n", task);
		if(count > nsteps){ // also stop if we go too long
			go_flag = 0;
			fprintf(stderr, "stopping early in lbfgs\n");
		}
	}
		
	copy_vec_gslvec(xvalue, xkFinal, nparams);



	free(nbd); 
	free(iwa); 
	free(grad); 
	free(lower); 
	free(upper); 
	//free(xvalue); /* causes some kind of double free ? */
	free(wa);
}


/* these are very confusing, work out something better!*/
void copy_gslvec_vec(gsl_vector* gsl_vec, double* vec, int n){
	int i;
	for(i = 0; i < n; i++){
		vec[i] = gsl_vector_get(gsl_vec, i);
	}
}


void copy_vec_gslvec(double* vec, gsl_vector* gsl_vec, int n){
	int i;
	for(i = 0; i < n; i++){
		gsl_vector_set(gsl_vec, i, vec[i]);
	}
}


	
void print_vec(gsl_vector* x, int n){
	int i;
	fprintf(stderr, "\n\n");
	for(i =0; i < n;i++){
		fprintf(stdout, "%f ", gsl_vector_get(x, i));
	}
	fprintf(stdout, "\n");
}


										
void set_zero( double *vec, int size){
	int i;
	for(i =0; i < size; i++){
		vec[i] = 0.0;
	}
}


void set_bounds(double *vec, gsl_matrix *ranges, int index, int nparams){
	int i;
	assert(index == 1 || index == 0);
	for( i = 0; i < nparams; i++){
		// get the indexth col of the ranges matrix
		vec[i] = gsl_matrix_get(ranges, i, index);
	}
}
		
void set_nbd(int *vec, int nparams, int nmax){
	int i;
	for( i = 0; i < nparams; i++){
		// get the indexth col of the ranges matrix
		vec[i] = 2;
	}
	for(i = nparams; i < nmax; i++){
		vec[i] = 0;
	}
}


// just to test
#ifdef EXECUTE
int main (void){
	int nparams = 2;
	int i;

	gsl_vector* xTest = gsl_vector_alloc(nparams);
	gsl_vector* xFinal = gsl_vector_alloc(nparams);
	gsl_matrix* ranges = gsl_matrix_alloc(nparams, 2);
	
	for(i =0; i < nparams; i++){
		gsl_matrix_set(ranges, i, 0, -3.0);
		gsl_matrix_set(ranges, i, 1, 3.0);
		gsl_vector_set(xTest, i, 0.8);
	}

	gsl_vector_set_zero(xFinal);

	doBoundedBFGS(&fRosenbrock, &getGradientNumericLBFGS, ranges, xTest, xFinal,  nparams, 1000, NULL);

	print_vec(xFinal, nparams);
	
	gsl_matrix_free(ranges);
	gsl_vector_free(xTest);
	gsl_vector_free(xFinal);
	return(0);
}
#endif
