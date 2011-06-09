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

// a macro to use i,j indexs on a linear array
# define IDX2C(i,j,ld) (((j)*(ld))+(i))

void *mallocChecked(size_t size);
void cudaMallocChecked(float** devptr, size_t size);
void cudaSetVectorChecked(int npoints, int elSize, float* hvec, int stridex, float* devVec, int stridey);

void handleCudaStatus(cublasStatus_t status, char* message);

void invertMatrixVectorCuda(gsl_matrix* Cn, gsl_vector* u, gsl_vector *answer, double theta3);
float computeLambdaCuda(cublasHandle_t handle, float* d_cn, float* d_h, float* d_g, float* d_temp, int npoints);

void invertMatrixVectorCuda(gsl_matrix* Cn, gsl_vector* u, gsl_vector *answer, double theta3){
	cublasStatus_t status;
	cublasHandle_t handle;
	cudaError_t err;
	int nsteps;
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
	float *d_y1, *d_g1, *d_h1, *d_temp; 
	float *d_ynext, *d_gnext, *d_hnext; 
	
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


	// err = cudaMalloc((void**)&d_cn, npoints*npoints*sizeof(d_cn[0]));
	// if(err != cudaSuccess){
	// 	fprintf (stderr, "!!!! device memory allocation error (allocate d_cn)\n");
	// 	exit(EXIT_FAILURE);
	// } 

	// // init the matrix

		
	// exit(1);
	
	// alloc device memory
	cudaMallocChecked(&d_cn, vectorSize*npoints);
	cudaMallocChecked(&d_temp, vectorSize);
	cudaMallocChecked(&d_y1, vectorSize);	
	cudaMallocChecked(&d_g1, vectorSize);	
	cudaMallocChecked(&d_h1, vectorSize);	
	cudaMallocChecked(&d_ynext, vectorSize);	
	cudaMallocChecked(&d_gnext, vectorSize);	
	cudaMallocChecked(&d_hnext, vectorSize);	

	
	// init dev memory from host
	printf("alloc vectors\n");
	cudaSetVectorChecked(npoints, sizeof(float), h_y1, 1, d_y1,1);
	cudaSetVectorChecked(npoints, sizeof(float), h_g1, 1, d_g1,1);
	cudaSetVectorChecked(npoints, sizeof(float), h_h1, 1, d_h1,1);

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


	lambda = computeLambdaCuda(handle, d_cn, d_h1, d_g1, d_temp, npoints);
	// for(i = 0; i < nsteps; i++){
		
	// }
	printf("lambda = %lf\n", lambda);


	cudaFree(d_cn);
	cudaFree(d_y1);	
	cudaFree(d_g1);	
	cudaFree(d_h1);	
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

	printf("number = %f\tdenom = %f\n", number, denom);
	
	return (number/denom);
}

/**
 * the gradient conjugate scale
 *
 * gamma_k  = \frac{ g_{k+1}. g_{k+1}}{ g_k . g_k }
 */
float computeGammaCuda(cublasHandle_t handle, float* d_gnext, float* d_g, int npoints){
	cublasStatus_t status;
	double num, denom;
	status = cublasSdot(handle, npoints, d_gnext, d_gnext, &num);
	handleCudaStatus(status, "computeGamma: d_gnext*d_gnext ddot failed\n");

	status = cublasSdot(handle, npoints, d_g, d_g, &denom);
	handleCudaStatus(status, "computeGamma: d_g*d_g ddot failed\n");
	
	return(num/denom);
}


void computeNextGCuda(cublasHandle_t handle, float* d_h, float* d_g, float* d_gnext, float lambda, int npoints){
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
		fprintf(stderr, message);
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
