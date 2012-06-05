#ifndef __INC_MULTIMODELSTRUCT__
#define __INC_MULTIMODELSTRUCT__



#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

struct modelstruct

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
	int nmodel_points;  // number of points in the design
	int nparams; // number of dimensions in the PARAMETER SPACE
	
	/** 
	 * the design (rows:nmodel_points) (cols:nparams)
	 */
	gsl_matrix *xmodel;

	/**
	 * the full training set in the original space
	 * (rows:nmodel_points) (cols:nt)
	 */
	gsl_matrix *training_matrix;
	
	/** the array of (r) pca decomposed models
	 */
	modelstruct* pca_model_array;

	/**
	 * the eigenvals and vectors from the pca decomp
	 */
	gsl_vector *eigenvalues;
	gsl_matrix *eigenvectors;
	
} multi_modelstruct;

multi_modelstruct* alloc_multmodelstruct(gsl_matrix *model_in, gsl_matrix *training_matrix_in, int cov_fn_index, int regression_order);

void free_multimodelstruct(multi_modelstruct *m);

#endif
