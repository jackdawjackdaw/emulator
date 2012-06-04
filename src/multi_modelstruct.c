#include "multi_modelstruct.h"
#include "multivar_support.h"
#include "assert.h"

multi_modelstruct* alloc_multmodelstruct(gsl_matrix *model_in, 
																				 gsl_matrix *training_matrix_in, 
																				 int cov_fn_index, 
																				 int regression_order){

	assert(training_matrix->size1 == xmodel->size1);
	assert(training_matrix->size1 > 0);
	assert(training_matrix->size2 > 0);
	assert(xmodel->size2 > 0);
	
	int nt = training_matrix->size2;
	int nmodel_points = xmodel->size1;
	int nparams = xmodel->size2;

	/* use default if out of range */
	if (regression_order < 0 || regression_order > 3)
		regression_order = 0;

	/* ntheta is a function of cov_fn_index and nparams */
	int nthetas;
	if ((cov_fn_index == MATERN32) || (cov_fn_index == MATERN52)) {
		nthetas = 3;
	} else if (cov_fn_index == POWEREXPCOVFN) {
		nthetas = nparams + 2;
	} else {
		cov_fn_index = POWEREXPCOVFN;
		nthetas = nparams + 2;
	}
	
	multi_modelstruct *model = (multi_modelstruct*)malloc(sizeof(multi_modelstruct));

	// fill in
	model->nt = nt;
	mdoel->nr = 0; // init at zero
	model->nmodel_points = nmodel_points;
	model->nparams = nparams;
	model->xmodel = xmodel_in;
	model->training_matrix = training_matrix_in;
	
	/* carry out the pca decomp on this model, this is defined in multivar_support for now
	 * this will fill in nr, eigenvalues, eigenvectors, and pca_model_array
	 */
	pca_decomp(model);

	return model;
}

/**
 * this free's everything in m
 */
void free_multimodelstruct(multi_modelstruct *m){
	int i;
	gsl_matrix_free(m->xmodel);
	gsl_matrix_free(m->training_matrix);
	for(i = 0; i < m->nr; i++){
		free_modelstruct_2(m->pca_model_array[i]);
		free_modelstruct(m->pca_model_array);
	}
	free(m->pca_model_array);
	gsl_matrix_free(m->eigenvalues);
	gsl_matrix_free(m->eigenvectors);
}
