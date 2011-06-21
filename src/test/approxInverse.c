#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_linalg.h"
#include <sys/time.h>

extern void invertMatrixVectorCuda(gsl_matrix* Cn, gsl_vector* u, gsl_vector *answer, double theta3);

double covFn(double x1, double x2, double theta1, double theta2, double theta3);
void invertMatrixVector(gsl_matrix* Cn, gsl_vector* u, gsl_vector *answer, double theta3);
double computeQValue(gsl_matrix* Cn, gsl_vector* u, gsl_vector *y);
double computeQmax(gsl_matrix* Cn, gsl_vector* u, gsl_vector *y, double theta3);
double computeLambda(gsl_matrix* Cn, gsl_vector* h, gsl_vector* g);
double computeGamma(gsl_vector *gnext, gsl_vector *g);
void computeNextG(gsl_matrix* Cn, gsl_vector* h, gsl_vector* g, double lambda, gsl_vector* gnext);
void computeNextH(gsl_vector* h, gsl_vector* gnext, double gamma, gsl_vector* hnext);
void computeNextY(gsl_vector *y, gsl_vector *h, double lambda, gsl_vector *ynext);

int main (void){
	int npoints = NPOINTS;
	int i,j;
	double theta1, theta2, theta3;
	double ctemp;
	double qmax;
	double texact, tapprox, tcuda, mindev,maxdev;
	double mindevCuda, maxdevCuda;
	struct timeval start, stop, diff;
	FILE *fptr;
	gsl_vector *temp = gsl_vector_calloc(npoints);
	gsl_vector *uvector = gsl_vector_calloc(npoints);
	gsl_vector *yexact = gsl_vector_calloc(npoints);
	gsl_vector *yapprox = gsl_vector_calloc(npoints);
	gsl_vector *diffVec = gsl_vector_calloc(npoints);

	gsl_vector *yapproxCuda = gsl_vector_calloc(npoints);
	gsl_vector *diffVecCuda = gsl_vector_calloc(npoints);

	
	gsl_matrix *cmatrix = gsl_matrix_alloc(npoints, npoints);

	gsl_matrix *cinverse = gsl_matrix_alloc(npoints, npoints);
	
	// we'll suppose a simple 1d model to test
	// space the points evenly on 1..10
	
	for(i = 0; i <npoints;i++){
		gsl_vector_set(uvector, i, 0.1*((double)i));
		//printf("%lf ", gsl_vector_get(uvector, i));
	}
	//printf("\n");

	theta1 = 1;
	theta2 = 0.2;
	theta3 = 0.01;

	for(i = 0; i <npoints;i++){
		for(j = 0; j < npoints; j++){
			ctemp = covFn(gsl_vector_get(uvector, i), gsl_vector_get(uvector, j), theta1, theta2, theta3);
			gsl_matrix_set(cmatrix, i, j, ctemp);
			//printf("%lf ", ctemp);
		}
		//printf("\n");
	}

	
	/* gettimeofday(&start, NULL); // get the time  */
	/* // so now we'll do the inversion exactly */
	/* gsl_matrix_memcpy(cinverse, cmatrix); */
	/* gsl_linalg_cholesky_decomp(cinverse); */
	/* gsl_linalg_cholesky_invert(cinverse); */

	/* //gsl_blas_dsymv(CblasNoTrans, 1.0,Cn, h, 0.0, temp); */
	/* // now form the product */
	/* gsl_blas_dsymv(CblasUpper, 1.0, cinverse, uvector, 0.0, yexact); */
	/* gettimeofday(&stop, NULL); */

	/* timersub(&stop, &start, &diff); */
	/* texact = diff.tv_sec + diff.tv_usec/1000000.0; */
	/* printf("exact time: %lf\n", texact); */

	// we can also compute Qmax to check what's going on
	// this is basically cheating
	/* qmax = computeQValue(cmatrix, uvector, yexact); */
	/* printf("qmax: %lf\t\n", qmax); */

	gettimeofday(&start, NULL); // get the time 	
	invertMatrixVector(cmatrix, uvector, yapprox, theta3);
	gettimeofday(&stop, NULL);

	timersub(&stop, &start, &diff);
	tapprox = diff.tv_sec + diff.tv_usec/1000000.0;
	printf("approx time: %lf\n", tapprox);

	gettimeofday(&start, NULL); // get the time 	
	invertMatrixVectorCuda(cmatrix, uvector, yapproxCuda, theta3);
	gettimeofday(&stop, NULL);

	timersub(&stop, &start, &diff);
	tcuda = diff.tv_sec + diff.tv_usec/1000000.0;
	printf("cuda time: %lf\n", tcuda);
	

	
	/* for(i = 0; i < npoints; i++){ */
	/* 	//printf("%d %lf %lf %lf\n", i, gsl_vector_get(yexact, i), gsl_vector_get(yapprox, i), fabs(gsl_vector_get(yexact,i) - gsl_vector_get(yapprox,i))); */
	/* 	gsl_vector_set(diffVec, i, fabs(gsl_vector_get(yexact,i) - gsl_vector_get(yapprox,i)) / gsl_vector_get(yexact, i)); */
	/* 	gsl_vector_set(diffVecCuda, i, fabs(gsl_vector_get(yexact,i) - gsl_vector_get(yapproxCuda,i)) / gsl_vector_get(yexact, i)); */
	/* } */
	
	/* maxdev = gsl_vector_max(diffVec); */
	/* mindev = gsl_vector_min(diffVec); */

	/* maxdevCuda = gsl_vector_max(diffVecCuda); */
	/* mindevCuda = gsl_vector_min(diffVecCuda); */


	/* printf("max deviation %lf cuda %lf: \n", maxdev, maxdevCuda); */
	/* printf("min deviation %lf cuda %lf: \n", mindev, mindevCuda); */
		
	gsl_vector_free(uvector);
	gsl_vector_free(yexact);
	gsl_vector_free(yapprox);
	gsl_vector_free(diffVec);

	gsl_vector_free(yapproxCuda);
	gsl_vector_free(diffVecCuda);

	gsl_matrix_free(cmatrix);
	gsl_matrix_free(cinverse);
	
	fptr = NULL;
	fptr = fopen("test-results.dat", "a");
	if(fptr == NULL){
		fprintf(stderr, "couldn't open results file!\n");
		return EXIT_FAILURE;
	}
	fprintf(fptr, "%d %lf %lf %lf\n", NPOINTS, tapprox, tcuda);
					

	return EXIT_SUCCESS;
}

double covFn(double x1, double x2, double theta1, double theta2, double theta3){
	double cov;
	double rsq = (x1-x2)*(x1-x2);
	double thetasq = theta2*theta2;
	cov = theta1*exp(-0.5*rsq/thetasq);
	if(fabs(x1 - x2) < 0.00001){
		cov += theta3;
	}
	return cov;
}

/**
 * given a matrix Cn and a vector u we 
 * want to compute Cn^{-1}.u 
 */
void invertMatrixVector(gsl_matrix* Cn, gsl_vector* u, gsl_vector *answer, double theta3){
	int nsteps = 100; // for testing...
	int i;
	gsl_vector* y1, *g1, *h1;
	gsl_vector* ynext, *gnext, *hnext;
	double gamma, lambda;
	double upperBound, lowerBound;
	double lowerPrev;
	double controlRatio;

	double stopDiff = 0.0001;
	
	y1 = gsl_vector_calloc(u->size); // note that y1 is set to zero
	g1 = gsl_vector_alloc(u->size);
	h1 = gsl_vector_alloc(u->size);

	ynext = gsl_vector_alloc(u->size);
	gnext = gsl_vector_alloc(u->size);
	hnext = gsl_vector_alloc(u->size);

	// initial conditions
	gsl_vector_memcpy(g1, u);
	gsl_vector_memcpy(h1, u);

	
	// a single step
	lowerPrev = lowerBound = 0;
	for(i = 0; i < nsteps; i++){		
		lambda = computeLambda(Cn, h1, g1);		
		computeNextG(Cn, h1, g1, lambda, gnext);
		gamma = computeGamma(gnext, g1);
		computeNextH(h1, gnext, gamma, hnext);
		computeNextY(y1, h1, lambda, ynext);
		
		//printf("l:%lf\tg:%lf\n", lambda, gamma);

		// now copy the next values back into the next ones and we do it again
		gsl_vector_memcpy(y1, ynext);
		gsl_vector_memcpy(g1, gnext);
		gsl_vector_memcpy(h1, hnext);

		//upperBound = computeQmax(Cn, u, y1, theta3);
		lowerBound = computeQValue(Cn, u, y1);
		if((lowerBound - lowerPrev) < stopDiff)
			break;
		//controlRatio = fabs((upperBound - lowerBound))/(lowerBound);
		//fprintf(stderr, "upper: %lf\tlower: %lf\tcontrolRatio: %lf\n", upperBound, lowerBound, controlRatio);

		lowerPrev = lowerBound;
	}
	
	printf("nsteps: %d\n", i);
	
	gsl_vector_memcpy(answer, y1);



	gsl_vector_free(y1);
	gsl_vector_free(g1);	
	gsl_vector_free(h1);	
}

////////// CONJUGATE GRADIENT BOUNDING //////////////

/** 
 * Q(y) can be used for the upper and lower bounds
 * Q(y) = y^T u - 1/2 * y^T C_N y
 */
double computeQValue(gsl_matrix* Cn, gsl_vector* u, gsl_vector *y){
	double qval, Atemp;
	gsl_vector *temp = gsl_vector_alloc(u->size);
	gsl_blas_ddot(y, u, &qval);
	
	//printf("computeQ: qval1 = %lf\n", qval);

	//-1/2 * C_n * y 
	gsl_blas_dsymv(CblasUpper, -0.5, Cn, y, 0.0, temp);
	gsl_blas_ddot(y, temp, &Atemp);

	qval = qval + Atemp;

	gsl_vector_free(temp);
	return(qval);
}
	
double computeQmax(gsl_matrix* Cn, gsl_vector* u, gsl_vector *y, double theta3){
	int i,j;
	double qstar, qmax, inner;
	gsl_vector* v = gsl_vector_alloc(u->size);
	gsl_vector* x = gsl_vector_alloc(u->size);
	gsl_matrix* A = gsl_matrix_alloc(u->size,u->size);
	gsl_matrix* temp = gsl_matrix_alloc(u->size,u->size);
	

	gsl_matrix_memcpy(A, Cn);
	gsl_vector_memcpy(v, u);
	gsl_vector_memcpy(x, y);

	// Cn = A + \theta3 I
	// so to get A we need to remove Theta3 from the diagonal
	for(i = 0; i < u->size; i++)
		gsl_matrix_set(A, i, i, (gsl_matrix_get(A, i, i) - theta3));

	/* printf("pre decomp \n\n"); */
	/* for(i = 0;i < u->size; i++){ */
	/* 	for(j = 0; j < u->size; j++){ */
	/* 		printf("%lf ", gsl_matrix_get(A, i, j)); */
	/* 	} */
	/* 	printf("\n"); */
	/* } */

	/* printf("post decomp \n\n"); */

	// v = A^{1/2}*u
	// x = A^{1/2}*y
	// A is now L L^T upper / lower triangular
	gsl_linalg_cholesky_decomp(A);

	
	
	/* for(i = 0;i < u->size; i++){ */
	/* 	for(j = 0; j < u->size; j++){ */
	/* 		//printf("%lf ", gsl_matrix_get(A, i, j)); */
	/* 		if(i <= j){ */
	/* 			gsl_matrix_set(temp, i, j, gsl_matrix_get(A,i, j)); */
	/* 		} */
	/* 		printf("%lf ", gsl_matrix_get(temp, i, j)); */
	/* 	} */
	/* 	printf("\n"); */
	/* } */

	gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, A, v);
	gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, A, x);
	qstar = computeQValue(Cn, v, x);

	gsl_blas_ddot(u, u, &inner);	
	/* printf("Qstar: %lf\n", qstar); */
	/* printf("1/2*inner: %lf\n", 0.5*inner); */
	//printf("1/2*inner - qstar: %lf\n", 0.5*inner - qstar);


	qmax = (1/theta3)*(0.5*inner - qstar);
	//qmax = (1/theta3)*qstar;

	gsl_vector_free(v);
	gsl_vector_free(x);
	gsl_matrix_free(A);
	return(qmax);
}


////////// CONJUGATE GRADIENT SCALARS //////////////

/** 
 * the stepsize at step k
 * \lambda_k = \frac{g_k^{T} g_k}{g_k^{T} C_n h_k}
 */
double computeLambda(gsl_matrix* Cn, gsl_vector* h, gsl_vector* g){
	double lambda;
	double num, denom;	
	gsl_vector * temp = gsl_vector_alloc(h->size);
	gsl_blas_ddot(g, g, &num);
	// we can use the symmetric product since C is symm by defn
	// y = alpha*opcode(A).x + beta
	// opcode, alpha, A, x, beta, y
	gsl_blas_dsymv(CblasUpper, 1.0,Cn, h, 0.0, temp);
	gsl_blas_ddot(g, temp, &denom);
	
	lambda = num / denom;

	//printf("number = %f\tdenom = %f\tlambda = %f\n", num, denom, lambda);
	
	gsl_vector_free(temp);
	return(lambda);
}

/**
 * the gradient conjugate scale
 *
 * gamma_k  = \frac{ g_{k+1}. g_{k+1}}{ g_k . g_k }
 */
double computeGamma(gsl_vector *gnext, gsl_vector *g){
	double gamma;
	double num, denom;
	
	gsl_blas_ddot(gnext, gnext, &num);
	gsl_blas_ddot(g, g, &denom);
	
	gamma = num / denom;
	//printf("g: num = %f\tdenom = %f\tgamma = %f\n", num, denom, num/denom);
	return(gamma);
}
	

////////// CONJUGATE GRADIENT VECTORS //////////////

/**
 * compute the next g vector g_k+1 from Cn, hk, gk and lambdak
 *
 * gnext should be supplied as an allocated gsl_vector
 * 
 * g_{k+1} = g_{k} - \lambda_k C_N h_k
 */
void computeNextG(gsl_matrix* Cn, gsl_vector* h, gsl_vector* g, double lambda, gsl_vector* gnext){
	gsl_vector *temp = gsl_vector_alloc(g->size);

	gsl_vector_memcpy(gnext, g); // copy original g
	gsl_blas_dsymv(CblasUpper, -1.0*lambda, Cn, h, 0.0, temp);
	
	// gnext' = gnext + temp
	gsl_vector_add(gnext, temp); 
	gsl_vector_free(temp);
}

/**
 * compute the next h vector h_k+1 from Cn, hk, gk+1, and gammak
 * 
 * h_{k+1} = g_{k+1} + \gamma_k h_k
 */
void computeNextH(gsl_vector* h, gsl_vector* gnext, double gamma, gsl_vector* hnext){
	gsl_vector_memcpy(hnext, h);
	gsl_vector_scale(hnext, gamma);
	gsl_vector_add(hnext, gnext);
}

/**
 * compute the next actual y-vector from lambda and h
 *
 * y_{k+1} = y_k + \lambda_k h_k
 */
void computeNextY(gsl_vector *y, gsl_vector *h, double lambda, gsl_vector *ynext){
	gsl_vector *temp = gsl_vector_alloc(h->size);
	gsl_vector_memcpy(ynext, y);
	gsl_vector_memcpy(temp, h);
	gsl_vector_scale(temp, lambda);
	gsl_vector_add(ynext, temp);

	gsl_vector_free(temp);
}
