#ifndef __INC_MULTIMODELSTRUCT__
#define __INC_MULTIMODELSTRUCT__


#include "modelstruct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

struct modelstruct;

/**
 * \struct multimodelstruct
 * handle a model with t-variate output i.e y_m = {y_1, ..., y_t}^{T}
 * 
 * this contains some r <= t modelstructs representing the emulators in the PCA projected
 * space
 */
typedef struct multi_modelstruct{
	int nt; // number of observables in OBSERVED space
	int nr; // number of pca-observables kept in the ROTATED space
	int nparams; // number of dimensions in the PARAMETER SPACE
	int nmodel_points;  // number of points in the design
	int cov_fn_index;
	int regression_order;
	double min_length_scale; // (<=0) don't use a min scale (>0) min scale for cov fn
	
	/** 
	 * the design (rows:nmodel_points) (cols:nparams)
	 */
	gsl_matrix *xmodel;
	/**
	 * the full training set in the original space
	 * (rows:nmodel_points) (cols:nt)
	 */
	gsl_matrix *training_matrix;
	/**
	 * a t length vector of the the mean values of the cols of training_matrix 
	 */
	gsl_vector *training_mean; 
	
	/** array of (r) pca decomposed models
	 * 
	 * pointers from this back to xmodel, is this bad?
	 */
	modelstruct** pca_model_array;
	
	/**
	 * the eigenvalues and vectors from the pca decomp (t xt) , not saved for now...
	 */
	//gsl_vector *pca_eigenvalues;
	//gsl_matrix *pca_eigenvectors;

	gsl_vector *pca_evals_r;	/** just the first r evals */
	gsl_matrix *pca_evecs_r;  /** the first r evecs (t x nr)*/
	

	/** 
	 * an (r x nmodel_points) matrix of the pca-transformed training points, these are 
	 * used to init each of the r entries of pca_model_array */
	gsl_matrix *pca_zmatrix; /** z <- (1/sqrt(pca_evals_r)) %*% t(pca_evecs_r) %*% (training_matrix - training_mean) */
	
} multi_modelstruct;

multi_modelstruct* alloc_multimodelstruct(gsl_matrix *xmodel_in, gsl_matrix *training_matrix_in, 
																					int cov_fn_index, int regression_order,
																					double varfrac, double min_length_scale);

void gen_pca_decomp(multi_modelstruct *m, double vfrac);
void gen_pca_model_array(multi_modelstruct *m);

void dump_multi_modelstruct(FILE* fptr, multi_modelstruct *m );
multi_modelstruct *load_multi_modelstruct(FILE* fptr);


double vector_elt_sum(gsl_vector* vec, int nstop);

void free_multimodelstruct(multi_modelstruct *m);

#endif
