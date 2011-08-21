/* Includes, system */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "math.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_linalg.h"
#include <sys/time.h>


/* Includes, cuda */
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <shrQATest.h>

using namespace std;

// a macro to use i,j indexs on a linear array
# define IDX2C(i,j,ld) (((j)*(ld))+(i))

void *mallocChecked(size_t size);
void cudaMallocChecked(float** devptr, size_t size);
void cudaSetVectorChecked(int npoints, int elSize, float* hvec, int stridex, float* devVec, int stridey);

void handleCudaStatus(cublasStatus_t status, char* message);

void invertMatrixVectorCuda(gsl_matrix* Cn, gsl_vector* u, gsl_vector *answer, double theta3);


float computeLambdaCuda(cublasHandle_t handle, float* d_cn, float* d_h, float* d_g, float* d_temp, int npoints);
float computeQValueCuda(cublasHandle_t handle, float* d_cn, float* d_u, float* d_y, float* d_temp, int npoints);
float computeGammaCuda(cublasHandle_t handle, float* d_gnext, float* d_g, int npoints);
void computeNextGCuda(cublasHandle_t handle, float* d_cn, float* d_h, float* d_g, float* d_gnext, float* d_temp, float lambda, int npoints);
void computeNextHCuda(cublasHandle_t handle, float* d_h, float* d_gnext, float* d_hnext, float* d_temp, float gamma, int npoints);
void computeNextYCuda(cublasHandle_t handle, float *d_y, float* d_ynext, float* d_h, float* d_temp, float lambda, int npoints);


void invertMatrixVectorCuda(gsl_matrix* Cn, gsl_vector* u, gsl_vector *answer, double theta3){
	cublasStatus_t status;
	cublasHandle_t handle;
	cudaError_t err;
	int nsteps = 100;
	int i, j;
	double gamma, lambda;
	double lowerBound, lowerPrev;
	int npoints;
	size_t vectorSize;
	
	lowerBound = lowerPrev = 0;

	const double stopDiff = 0.0001;
	// h_thing is a host array (or vector)
	// d_thing is a device array (or vector)
	float *h_cn, *h_u, *h_answer;
	float *h_g1, *h_h1, *h_y1; // we need host copies of g1, h1 and y1 since they get init
	float *d_cn; // this is a matrix
	float *d_u; // we need this for computing Q
	float *d_y1, *d_g1, *d_h1, *d_temp; 
	float *d_ynext, *d_gnext, *d_hnext; 
	float *tempPtr;
	
	npoints = u->size; 
	vectorSize = npoints*sizeof(float);

	// allocate host memory...
	// which we have to cast because c++ is annoying
	h_cn = (float*)mallocChecked(npoints*vectorSize);
	h_u = (float*)mallocChecked(vectorSize);
	h_answer = (float*)mallocChecked(vectorSize);
	h_g1 = (float*)mallocChecked(vectorSize);
	h_h1 = (float*)mallocChecked(vectorSize);
	h_y1 = (float*)mallocChecked(vectorSize);
	
	// init host memory (stupid)
	for(i = 0; i < npoints; i++){
		for(j =0; j < npoints; j++){
			h_cn[IDX2C(i,j,npoints)] = gsl_matrix_get(Cn, i, j);
		}
	}

	for(i = 0; i < npoints;i++){
		h_y1[i] = 0.0;
		h_u[i] = gsl_vector_get(u, i);
		h_answer[i] = 0.0;
		h_g1[i] = gsl_vector_get(u, i);
		h_h1[i] = gsl_vector_get(u, i);
	}
	
	status = cublasCreate(&handle);
	if(status != CUBLAS_STATUS_SUCCESS){
		fprintf(stderr, "cublas init failed :(\n");
		exit(EXIT_FAILURE);
	} 
	printf("cublas init ok :D \n");

	// alloc device memory
	cudaMallocChecked(&d_cn, vectorSize*npoints);
	cudaMallocChecked(&d_temp, vectorSize);
	cudaMallocChecked(&d_y1, vectorSize);	
	cudaMallocChecked(&d_g1, vectorSize);	
	cudaMallocChecked(&d_h1, vectorSize);	
	cudaMallocChecked(&d_u, vectorSize);	
	cudaMallocChecked(&d_ynext, vectorSize);	
	cudaMallocChecked(&d_gnext, vectorSize);	
	cudaMallocChecked(&d_hnext, vectorSize);	

	
	// init dev memory from host
	printf("alloc vectors\n");
	cudaSetVectorChecked(npoints, sizeof(float), h_y1, 1, d_y1,1);
	cudaSetVectorChecked(npoints, sizeof(float), h_g1, 1, d_g1,1);
	cudaSetVectorChecked(npoints, sizeof(float), h_h1, 1, d_h1,1);
	cudaSetVectorChecked(npoints, sizeof(float), h_u, 1, d_u,1);

	printf("alloc matrix\n");
	status = cublasSetMatrix(npoints, npoints, sizeof(float), h_cn, npoints, d_cn, npoints);
	if(status != CUBLAS_STATUS_SUCCESS){
		fprintf(stderr, "cublas matrix init failed\n");
		if(status == CUBLAS_STATUS_INVALID_VALUE){
			fprintf(stderr, "CUBLAS_STATUS_INVALID_VALUE\n");
		} else if(status == CUBLAS_STATUS_MAPPING_ERROR){
			fprintf(stderr, "CUBLAS_STATUS_MAPPING_ERROR\n");
		}
		exit(EXIT_FAILURE);
	}

	
	for(i = 0; i < nsteps; i++){
		lambda = computeLambdaCuda(handle, d_cn, d_h1, d_g1, d_temp, npoints);

		computeNextGCuda(handle, d_cn, d_h1, d_g1, d_gnext, d_temp, lambda, npoints);
		gamma = computeGammaCuda(handle, d_gnext, d_g1, npoints);

		computeNextHCuda(handle, d_h1, d_gnext, d_hnext, d_temp, gamma, npoints);
		computeNextYCuda(handle, d_y1, d_ynext, d_h1, d_temp, lambda, npoints);

		// swap the pointers
		tempPtr = d_y1;	d_y1 = d_ynext;	d_ynext = tempPtr;
		tempPtr = d_g1; d_g1 = d_gnext; d_gnext = tempPtr;
		tempPtr = d_h1; d_h1 = d_hnext; d_hnext = tempPtr;
		// status = cublasScopy(handle, npoints, d_ynext, 1, d_y1, 1);
		// handleCudaStatus(status, "scopy failed d_ynext->d_y1\n");
		// status = cublasScopy(handle, npoints, d_gnext, 1, d_g1, 1);
		// handleCudaStatus(status, "scopy failed d_gnext->d_g1\n");
		// status = cublasScopy(handle, npoints, d_hnext, 1, d_h1, 1);
		// handleCudaStatus(status, "scopy failed d_hnext->d_h1\n");

		lowerBound = computeQValueCuda(handle, d_cn, d_u, d_y1, d_temp, npoints);

		//printf("l:%lf\tg:%lf\tlower:%f\n", lambda, gamma, lowerBound);

		if((lowerBound - lowerPrev) < stopDiff){
			break;
		}


		lowerPrev = lowerBound;
	}

	printf("cuda: nstep %d, tol %f\n", i, stopDiff);
	status = cublasGetVector(npoints, sizeof(float),  d_y1, 1, h_answer, 1);
	for(i = 0; i < npoints; i++)
		gsl_vector_set(answer, i, h_answer[i]);
	

	cudaFree(d_cn);
	cudaFree(d_y1);	
	cudaFree(d_g1);	
	cudaFree(d_h1);	
	cudaFree(d_u);	
	cudaFree(d_temp);
	cudaFree(d_ynext);	
	cudaFree(d_gnext);	
	cudaFree(d_hnext);	

	status = cublasDestroy(handle);
	if (status != CUBLAS_STATUS_SUCCESS) {
		fprintf (stderr, "!!!! shutdown error (A)\n");
		exit( EXIT_FAILURE);
	}


}

/** 
 * Q(y) can be used for the upper and lower bounds
 * Q(y) = y^T u - 1/2 * y^T C_N y
 */
float computeQValueCuda(cublasHandle_t handle, float* d_cn, float* d_u, float* d_y, float* d_temp, int npoints){
	cublasStatus_t status;
	float qval, temp;
	float alpha = -0.5;
	float beta = 0.0;
	
	status = cublasSdot(handle, npoints, d_y, 1, d_u, 1, &qval);
	handleCudaStatus(status, "computeQCuda: d_y d_u ddot failed\n");
	
	//printf("computeQ: qval1 = %lf\n", qval);

	// -1/2 * C_n * y
	status = cublasSsymv(handle, CUBLAS_FILL_MODE_UPPER, 
											 npoints, &alpha, 
											 d_cn, npoints, 
											 d_y, 1,
											 &beta, d_temp, 1);
	handleCudaStatus(status, "computeQ: dsymv failed\n");
	
	status = cublasSdot(handle, npoints, d_y, 1, d_temp, 1, &temp);
	handleCudaStatus(status, "computeQCuda: d_y d_temp ddot failed\n");

	return(qval + temp);
}

	

/** 
 * the stepsize at step k
 * \lambda_k = \frac{g_k^{T} g_k}{g_k^{T} C_n h_k}
 */
float computeLambdaCuda(cublasHandle_t handle, float* d_cn, float* d_h, float* d_g, float* d_temp, int npoints){
	cublasStatus_t status;
	float number, denom;
	float alpha = 1.0;
	float beta = 0.0;
	// compute g.g (easy one call)
	status = cublasSdot(handle, npoints, d_g, 1, d_g, 1, &number);
	handleCudaStatus(status, "computeLambda: gg ddot failed\n");

	// now we want to compute the product g^{t} C_n h
	// d_temp = C_n %*% h
	status = cublasSsymv(handle, CUBLAS_FILL_MODE_UPPER, 
											 npoints, &alpha, 
											 d_cn, npoints, 
											 d_h, 1,
											 &beta, d_temp, 1);
	handleCudaStatus(status, "computeLambda: dsymv failed\n");

	// compute d_g^t * d_temp
	status = cublasSdot(handle, npoints, d_g, 1, d_temp, 1, &denom);
	handleCudaStatus(status, "computeLambda: g*temp ddot failed\n");

	//printf("l:number = %f\tdenom = %f\n", number, denom);
	
	return (number/denom);
}

/**
 * the gradient conjugate scale
 *
 * gamma_k  = \frac{ g_{k+1}. g_{k+1}}{ g_k . g_k }
 */
float computeGammaCuda(cublasHandle_t handle, float* d_gnext, float* d_g, int npoints){
	cublasStatus_t status;
	float num, denom;
	status = cublasSdot(handle, npoints, d_gnext, 1, d_gnext,1, &num);
	handleCudaStatus(status, "computeGamma: d_gnext*d_gnext ddot failed\n");

	status = cublasSdot(handle, npoints, d_g, 1, d_g, 1, &denom);
	handleCudaStatus(status, "computeGamma: d_g*d_g ddot failed\n");
	
	//printf("g: num = %f\tdenom = %f\n", num, denom);
	
	return(num/denom);
}

/**
 * compute the next g vector g_k+1 from Cn, hk, gk and lambdak
 *
 * gnext should be supplied as an allocated gsl_vector
 * 
 * g_{k+1} = g_{k} - \lambda_k C_N h_k
 */
void computeNextGCuda(cublasHandle_t handle, float* d_cn, float* d_h, float* d_g, float* d_gnext, float* d_temp, float lambda, int npoints){
	cublasStatus_t status;
	float alpha = -1*lambda;
	float beta = 0.0;
	// copy g in to gnext
	status = cublasScopy(handle, npoints, d_g, 1, d_gnext, 1);
	handleCudaStatus(status, "computeNextGCuda: cublaSsCOPY failed\n");

	status = cublasSsymv(handle, CUBLAS_FILL_MODE_UPPER, 
											 npoints, &alpha, 
											 d_cn, npoints, 
											 d_h, 1,
											 &beta, d_temp, 1);
	handleCudaStatus(status, "computeLambda: dsymv failed\n");

	alpha = 1.0;
	// Saxpy(alpha, x, y) :-> y = alpha*x + y
	// so we will add temp and get our answer back into gnext
	status = cublasSaxpy(handle, npoints, &alpha, d_temp, 1, d_gnext, 1);
	handleCudaStatus(status, "computeLambda: Saxpy failed\n");
	
}


/**
 * compute the next h vector h_k+1 from Cn, hk, gk+1, and gammak
 * 
 * h_{k+1} = g_{k+1} + \gamma_k h_k
 * 
 * we need to do this without nuking g_next and h etc
 */
void computeNextHCuda(cublasHandle_t handle, float* d_h, float* d_gnext, float* d_hnext, float* d_temp, float gamma, int npoints){
	cublasStatus_t status;

	// copy d_gnext->d_temp
	status = cublasScopy(handle, npoints, d_gnext, 1, d_temp, 1);
	handleCudaStatus(status, "computeNextHCuda; Scopy d_gnext->d_temp failed\n");
	// d_temp = gamma * d_h + d_temp
	status = cublasSaxpy(handle, npoints, &gamma, d_h, 1, d_temp, 1);
	handleCudaStatus(status, "computeNextHCuda; Saxpy d_temp = gamma*d_h + d_temp failed\n");
	// copy d_temp ->d_hnext
	status = cublasScopy(handle, npoints, d_temp, 1, d_hnext, 1);
	handleCudaStatus(status, "computeNextHCuda; Scopy d_temp->d_hnext failed\n");
}


/**
 * compute the next actual y-vector from lambda and h
 *
 * y_{k+1} = y_k + \lambda_k h_k
 */
void computeNextYCuda(cublasHandle_t handle, float *d_y, float* d_ynext, float* d_h, float* d_temp, float lambda, int npoints){
	cublasStatus_t status;
	// copy d_y->d_ynext
	status = cublasScopy(handle, npoints, d_y, 1, d_ynext, 1);
	handleCudaStatus(status, "computeNextYCuda; Scopy d_y->d_ynext failed\n");

	// d_ynext = lambda * d_h + d_ynext(which is y)
	status = cublasSaxpy(handle, npoints, &lambda, d_h, 1, d_ynext, 1);
	handleCudaStatus(status, "computeNextYCuda; Saxpy d_ynext = lambda*d_h + d_ynext failed\n");
	
}


///////////////////////////////////////////// WRAPPED Routines ////////////////////////////////////////

//! checks for null
/**
 * but what is null in cpp? 0?
 */
void *mallocChecked(size_t size){
	void *r = malloc(size);
	if( r == NULL){
		fprintf(stderr
						, "memory wasn't allocated");
		exit(EXIT_FAILURE);
	}
	return(r);
}

void handleCudaStatus(cublasStatus_t status, char* message){
	if(status != CUBLAS_STATUS_SUCCESS){
		fprintf(stderr,"%s\n", message);
		exit(EXIT_FAILURE);
	}
}

void cudaMallocChecked(float** devptr, size_t size){
	if(cudaMalloc((void**)devptr, size) != cudaSuccess){
		fprintf (stderr, "!!!! device memory allocation error\n");
		exit(EXIT_FAILURE);
	}
}


void cudaSetVectorChecked(int npoints, int elSize, float* hvec, int stridex, float* devVec, int stridey){
	cublasStatus_t stat;
	stat = cublasSetVector(npoints, elSize, hvec, stridex, devVec, stridey);
	if(stat != CUBLAS_STATUS_SUCCESS){
		fprintf(stderr, "! cudaSetVectorChecked failed\n");

		if(stat == CUBLAS_STATUS_INVALID_VALUE){
			fprintf(stderr, "CUBLAS_STATUS_INVALID_VALUE\n");
		} else if(stat == CUBLAS_STATUS_MAPPING_ERROR){
			fprintf(stderr, "CUBLAS_STATUS_MAPPING_ERROR\n");
		}

		exit(EXIT_FAILURE);
	}
}

