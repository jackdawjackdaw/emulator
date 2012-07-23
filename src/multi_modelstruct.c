#include "multi_modelstruct.h"
#include "multivar_support.h"

#include <math.h>
#include <assert.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>

/** 
 * ccs, 05.06.2012
 * \todo this all needs to be checked carefully against an example that works in R
 * \todo valgrind for mem-leaks
 * \bug no error checking in linalg
 * \bug no error checking in alloc
 * \bug is there error checking in file handles?
 *
 * ccs, 19.06.2012
 * if you allocate a multi_modelstruct m
 * then m->pca_model_array[i]->xmodel == m->xmodel
 * the value in the pca_model_array is just a pointer to m, do we really want this?
 */

/**
 * allocates a multi_modelstruct, like alloc_modelstruct_2, 
 * but for multivariate models with t output values at each location
 *
 * @param model_in: (n x d) matrix of the design
 * @param training_matrix: (n x t) matrix of the values of the training values at each of the n locations
 * @param cov_fn_index:  POWEREXPCOVFN, MATERN32, or MATERN52
 * @param regression_order:  0, 1, 2, or 3
 * @param varfrac: the minimum fractional variance that should be retained during the PCA decomp
 * 
 * applies a pca decomp to training_matrix to reduce the dimensionality
 * 
 */
multi_modelstruct* alloc_multimodelstruct(gsl_matrix *xmodel_in, 
																				 gsl_matrix *training_matrix_in, 
																				 int cov_fn_index, 
																				 int regression_order, double varfrac)
{
	assert(training_matrix_in->size1 == xmodel_in->size1);
	assert(training_matrix_in->size1 > 0);
	assert(training_matrix_in->size2 > 0);
	assert(xmodel_in->size2 > 0);
	

	int i;
	double mean_temp = 0.0;
	int nt = training_matrix_in->size2;
	int nmodel_points = xmodel_in->size1;
	int nparams = xmodel_in->size2;
	gsl_vector_view col_view;
	

	/* use default if out of range */
	if (regression_order < 0 || regression_order > 3)
		regression_order = 0;

	/* use a sensible default for the variance fraction */
	if(varfrac < 0 || varfrac > 1)
		varfrac = 0.95;

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
	
	// this doesn't seem to be allocating correctly, seems to be a broken defn of MallocChecked
	// strangely, code will build in this case...
	//multi_modelstruct * model = (multi_modelstruct*)MallocChecked(sizeof(multi_modelstruct));
	multi_modelstruct * model = (multi_modelstruct*)malloc(sizeof(multi_modelstruct));

	// fill in
	model->nt = nt;
	model->nr = 0; // init at zero
	model->nmodel_points = nmodel_points;
	model->nparams = nparams;
	model->xmodel = xmodel_in;
	model->training_matrix = training_matrix_in;
	model->training_mean = gsl_vector_alloc(nt);
	model->regression_order = regression_order;
	model->cov_fn_index = cov_fn_index;
	
	/* fill in the mean vector, should probably sum this more carefully... */
	for(i = 0; i < nt; i++){
		col_view = gsl_matrix_column(model->training_matrix, i);
		mean_temp = vector_elt_sum(&col_view.vector, nmodel_points);
		//printf("%lf\n", (mean_temp/((double)nmodel_points)));
		gsl_vector_set(model->training_mean, i, (mean_temp/((double)nmodel_points)) );
	}

	/* carry out the pca decomp on this model, this is defined in multivar_support for now
	 * this will fill in nr, pca_eigenvalues, pca_eigenvectors, pca_evals_r, pca_evecs_r
	 * 
	 * this is making a mess if nt = 1
	 */
	gen_pca_decomp(model, varfrac);

	/* fill in pca_model_array */
	gen_pca_model_array(model);
	
	return model;
}

/**
 * fill in pca_model_array
 * requires: 
 * - the pca decomp to have been calculated, so nr and all the pca_... fields are allocated and filled in
 * - m to have been allocated up to: nt, nmodel_points, and xmodel
 *
 * is it possible we could start running out of memory doing all these allocs?
 */
void gen_pca_model_array(multi_modelstruct *m)
{
	int nr = m->nr;
	int i;
	gsl_vector_view col_view; 
	gsl_vector* temp_train_vector = gsl_vector_alloc(m->nmodel_points);
	//gsl_matrix* temp_xmodel = gsl_matrix_alloc(m->nmodel_points, m->nparams);

	//gsl_matrix_mempcy(temp_xmodel, m->xmodel);
	
	// alloc the array of nr model structs
	//m->pca_model_array = (modelstruct**)MallocChecked(sizeof(modelstruct*)*nr);
	m->pca_model_array = (modelstruct**)malloc(sizeof(modelstruct*)*nr);
	// fill in the modelstructs correctly
	for(i = 0; i < nr; i++){
		col_view = gsl_matrix_column(m->pca_zmatrix, i);
		gsl_vector_memcpy(temp_train_vector, &(col_view.vector));

		// this isn't copying in the training vector correctly for somereason
		m->pca_model_array[i] = alloc_modelstruct_2(m->xmodel, temp_train_vector,
																								m->cov_fn_index, m->regression_order);

		// see if brute forcing it will work
		m->pca_model_array[i]->training_vector = gsl_vector_alloc(m->nmodel_points);
		gsl_vector_memcpy(m->pca_model_array[i]->training_vector, temp_train_vector);
	}
	gsl_vector_free(temp_train_vector);
}


/**
 * carries out a pca decomp on m->training_matrix;
 * setting m->nr, m->pca_eigenvalues, m->pca_eigenvectors and
 * initializing m->pca_model_array
 * 
 * the pca decomp is pretty simple: 
 * let Y_ij = m->y_training
 * let mu_i = (1/nmodel_points) sum_{j=1}^{nmodel_points} Y_{ij}  // the sample mean
 * let ysub_i = Y_i - rep(mu_i,nmodel_points)  //subtract the sample means from each column 
 * let sigma_ij = (1/nmodel_points) ( ysub * t(ysub)) // this is the covariance matrix
 * 
 * then all we need to do is compute the eigendecomposition
 * 
 * sigma_ij = U^{-1} Lambda U
 * 
 * where Lambda is a diagonal matrix of t eigenvalues and U is a t x t matrix with the eigenvectors as columns
 *
 * requires: 
 * - 0 < vfrac < 1 (0.95 is a good value)
 * - m to have allocated and filled out, nt, nmodel_points, trianing_matrix, (doesn't need xmodel explicitly)
 */
void gen_pca_decomp(multi_modelstruct *m, double vfrac)
{
	FILE *fptr; // for debug output
	int i,j;
	int nt = m->nt;
	int retval;
	double total_variance = 0.0, frac  = 0.0;
	gsl_matrix *pca_zmatrix_temp;
	gsl_matrix *y_sub_mat = gsl_matrix_alloc(m->nmodel_points, nt);
	gsl_matrix *y_temp_mat = gsl_matrix_alloc(m->nmodel_points, nt);

	gsl_matrix *y_cov_mat = gsl_matrix_alloc(nt, nt);
	
	gsl_vector *evals_temp = gsl_vector_alloc(nt);
	gsl_matrix *evecs_temp = gsl_matrix_alloc(nt,nt);
	
	gsl_eigen_symmv_workspace *ework = gsl_eigen_symmv_alloc(nt);
	
	gsl_vector_view col_view;

	// why is this here?
	gsl_matrix_memcpy(y_sub_mat, m->training_matrix);
	

	// subtract out the mean 
	for(i = 0; i < nt; i++){
		//col_view = gsl_matrix_column(y_sub_mat, i);
		printf("# y(%d) mean: %lf\n", i, gsl_vector_get(m->training_mean, i));
		for(j = 0; j < m->nmodel_points; j++){
			gsl_matrix_set(y_sub_mat, j, i, gsl_matrix_get(y_sub_mat, j, i) - gsl_vector_get(m->training_mean, i));
		}
	}
	
	// compute the sample-variance, by multiplying y_sub_mat, with itself transposed
	gsl_matrix_memcpy(y_temp_mat, y_sub_mat);
	gsl_matrix_set_zero(y_cov_mat);

	/** — Function: int gsl_blas_dgemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha, const gsl_matrix * A, const gsl_matrix * B, double beta, gsl_matrix * C)  (always forget this one)*/
	/* want C (nt x nt ) so we need to do: (nt x nmodel_points) * (nmodel_points x nt) */
	/** 
	 * this is strange, the upper triangle is not obviously correct
	 */
	//retval = gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, y_temp_mat, y_sub_mat, 0.0, y_cov_mat);
	retval = gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, y_sub_mat, y_temp_mat, 0.0, y_cov_mat);
	if(retval){
		printf("# gen_pca_decomp:gsl_blas_dgemm %s\n", gsl_strerror(retval));
		exit(EXIT_FAILURE);
	}

	
	gsl_matrix_scale(y_cov_mat, (1.0/((double)m->nmodel_points)));

	/**
	 * debug output
	 */
	#ifdef DEBUGPCA
	fptr = fopen("pca-debug.dat","w");
	fprintf(fptr, "# ycov:\n");
	for(j = 0; j < nt; j++){
		for(i = 0; i < nt; i++)
			fprintf(fptr, "%lf ", gsl_matrix_get(y_cov_mat, i, j));
		fprintf(fptr, "\n");
	}
	#endif

	/** now the eigendecomp 
	 * y_cov_mat is symmetric and better be real so we can use gsl_eigen_symmv,
	 * note that the ?upper? triangle of y_cov_mat is borked during this process
	 * also: the evals are not sorted by order, but the evectors are in the order of the evalues.
	 * so we need to sort the evalues and the evectors correctly before we can use them
	 * using
	 * — Function: int gsl_eigen_symmv_sort (gsl_vector * eval, gsl_matrix * evec, gsl_eigen_sort_t sort_type)
	 * sort them into descending order (biggest first)
	 */
	gsl_eigen_symmv(y_cov_mat, evals_temp, evecs_temp, ework); 
	gsl_eigen_symmv_sort(evals_temp, evecs_temp, GSL_EIGEN_SORT_VAL_DESC);
	/**
	 * eigenvectors are stored in columns of pca_evecs*/
	total_variance = vector_elt_sum(evals_temp, nt);

	
#ifdef DEBUGPCA
	fprintf(fptr, "# evals:\n");
	for(i = 0; i < nt; i++)
		fprintf(fptr, "%lf %lf\n", gsl_vector_get(evals_temp, i), gsl_vector_get(evals_temp, i) / total_variance);

	fprintf(fptr, "# evecs:\n");
	for(j = 0; j < nt; j++){
		for(i = 0; i < nt; i++)
			fprintf(fptr, "%lf ", gsl_matrix_get(evecs_temp, i, j));
		fprintf(fptr, "\n");
	}

	#endif

	i=0;
	while( frac < vfrac && (i+1) < nt){
		frac = (1.0/total_variance) * vector_elt_sum(evals_temp, i);
		i++;
	}
	m->nr = i;
	
	if(nt == 1){
		//printf("# 1d case, nr=1\n");
		m->nr = 1;
	}

	m->pca_evals_r = gsl_vector_alloc(m->nr);
	m->pca_evecs_r = gsl_matrix_alloc(m->nt, m->nr);
	// debug...
	fprintf(stderr, "# nr: %d frac: %lf\n", m->nr, frac);
	
	for(i = 0; i < m->nr; i++){
		gsl_vector_set(m->pca_evals_r, i, gsl_vector_get(evals_temp, i));
		col_view = gsl_matrix_column(evecs_temp, i);
		gsl_matrix_set_col(m->pca_evecs_r, i, &col_view.vector);
	}

	// fill in pca_zmatrix
	m->pca_zmatrix = gsl_matrix_alloc(m->nmodel_points, m->nr);
	pca_zmatrix_temp = gsl_matrix_alloc(m->nmodel_points, m->nr);
	// zmat: (nmodel_points x nr) = (nmodel_points x nt) * ( nt x nr ) 
	/** — Function: int gsl_blas_dgemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha, const gsl_matrix * A, const gsl_matrix * B, double beta, gsl_matrix * C)  (always forget this one)*/
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0,  y_sub_mat, m->pca_evecs_r, 0.0, m->pca_zmatrix);

	#ifdef DEBUGPCA
	fprintf(fptr, "# evecs R:\n");
	for(i = 0; i < m->nt; i++){
		for(j = 0; j < m->nr; j++)
			fprintf(fptr, "%lf ", gsl_matrix_get(m->pca_evecs_r, i, j));
		fprintf(fptr, "\n");
	}
	#endif
	
	gsl_matrix_free(y_temp_mat);
	y_temp_mat = gsl_matrix_alloc(m->nr, m->nr);
	for(i = 0; i < m->nr; i++) // scale the diagonal by the evalue */
		gsl_matrix_set(y_temp_mat, i, i, 1.0/(sqrt(gsl_vector_get(m->pca_evals_r, i))));

	//print_matrix(y_temp_mat, m->nr, m->nr);
	
	// if nr != nt this won't work!
	//gsl_matrix_memcpy(y_sub_mat, m->pca_zmatrix);
	gsl_matrix_memcpy(pca_zmatrix_temp, m->pca_zmatrix);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, pca_zmatrix_temp, y_temp_mat, 0.0, m->pca_zmatrix);

	#ifdef DEBUGPCA
	fprintf(fptr, "# zmat:\n");
	for(i = 0; i < 5; i++){
		for(j = 0; j < m->nr; j++)
			fprintf(fptr, "%lf ", gsl_matrix_get(m->pca_zmatrix, i, j));
		fprintf(fptr, "\n");
	}
	fclose(fptr);
	#endif

	

	gsl_vector_free(evals_temp);
	gsl_matrix_free(evecs_temp);
	gsl_matrix_free(pca_zmatrix_temp);

	gsl_eigen_symmv_free(ework);
	gsl_matrix_free(y_sub_mat);
	gsl_matrix_free(y_temp_mat);
	gsl_matrix_free(y_cov_mat);
}

/**
 * dump the multimodelstruct to fptr, follows from dump_modelstruct_2
 * we first dump the new fields and then iterate through the nr additional models which are dumped as before
 * so we dump a lot of the same info, but this is probably ok, the advantage is that each section defining
 * a model can be pulled out and worked on separately...
 */
void dump_multi_modelstruct(FILE* fptr, multi_modelstruct *m){
	assert(fptr);

	int i,j;
	int nt = m->nt;
	int nr = m->nr;
	int nparams = m->nparams;
	int nmodel_points = m->nmodel_points;
	int cov_fn_index = m->cov_fn_index;
	int regression_order = m->regression_order;
	
	fprintf(fptr, "%d\n", nt);
	fprintf(fptr, "%d\n", nr);
	fprintf(fptr, "%d\n", nparams);
	fprintf(fptr, "%d\n", nmodel_points);
	fprintf(fptr, "%d\n", cov_fn_index);
	fprintf(fptr, "%d\n", regression_order);

	// multimodel thetas are inside pca_model_array...

	for(i = 0; i < nmodel_points; i++){
		for(j = 0; j < nparams; j++)
			fprintf(fptr, "%.17lf ", gsl_matrix_get(m->xmodel, i, j));
		fprintf(fptr, "\n");
	}
	
	for(i = 0; i < nmodel_points; i++){
		for(j = 0; j < nt; j++)
			fprintf(fptr, "%.17lf ", gsl_matrix_get(m->training_matrix, i, j));
		fprintf(fptr, "\n");
	}
	
	// now the rest of the pca information
	for(i = 0; i < nr; i++)
			fprintf(fptr, "%.17lf ", gsl_vector_get(m->pca_evals_r, i));
		fprintf(fptr, "\n");

	for(i = 0; i < nt; i++){
		for(j = 0; j < nr; j++)
			fprintf(fptr, "%.17lf ", gsl_matrix_get(m->pca_evecs_r, i, j));
		fprintf(fptr, "\n");
	}

	for(i = 0; i < nmodel_points; i++){
		for(j = 0; j < nr; j++)
			fprintf(fptr, "%.17lf ", gsl_matrix_get(m->pca_zmatrix, i, j));
		fprintf(fptr, "\n");
	}

	// now the pca_model_array
	// these could go in separate files...
	for(i = 0; i < nr; i++){
		dump_modelstruct_2(fptr, m->pca_model_array[i]);
	}
	
}

/**
 * loads a multivariate modelstructure from fptr
 */
multi_modelstruct *load_multi_modelstruct(FILE* fptr){
	multi_modelstruct *m = (multi_modelstruct*)malloc(sizeof(multi_modelstruct));
	int i,j;
	int nt, nr;
	int nparams, nmodel_points;
	int cov_fn_index;
	int regression_order;

	double mean_temp;
	gsl_vector_view col_view;

	fscanf(fptr, "%d%*c", & nt);
	fscanf(fptr, "%d%*c", & nr);
	fscanf(fptr, "%d%*c", & nparams);
	fscanf(fptr, "%d%*c", & nmodel_points);
	fscanf(fptr, "%d%*c", & cov_fn_index);
	fscanf(fptr, "%d%*c", & regression_order);
	
	m->nt = nt;
	m->nr = nr;
	m->nparams = nparams;
	m->nmodel_points = nmodel_points;
	m->cov_fn_index = cov_fn_index;
	m->regression_order = regression_order;

	// now we can allocate everything in m
	m->xmodel = gsl_matrix_alloc(nmodel_points, nparams);
	m->training_matrix = gsl_matrix_alloc(nmodel_points, nt);
	m->training_mean = gsl_vector_alloc(nt); // do we need this? (yes!)
	m->pca_model_array = (modelstruct**)malloc(sizeof(modelstruct*)*nr);
	m->pca_evals_r = gsl_vector_alloc(nr);
	m->pca_evecs_r = gsl_matrix_alloc(nt, nr);
	m->pca_zmatrix = gsl_matrix_alloc(nmodel_points, nr);


	for(i = 0; i < nmodel_points; i++)
		for(j = 0; j < nparams; j++)
			fscanf(fptr, "%lf%*c", gsl_matrix_ptr(m->xmodel, i, j));
	
	for(i = 0; i < nmodel_points; i++)
		for(j = 0; j < nt; j++)
			fscanf(fptr, "%lf%*c", gsl_matrix_ptr(m->training_matrix, i, j));
	
	// now the rest of the pca information
	for(i = 0; i < nr; i++)
			fscanf(fptr, "%lf%*c", gsl_vector_ptr(m->pca_evals_r, i));

	for(i = 0; i < nt; i++)
		for(j = 0; j < nr; j++)
			fscanf(fptr, "%lf%*c", gsl_matrix_ptr(m->pca_evecs_r, i, j));

	for(i = 0; i < nmodel_points; i++)
		for(j = 0; j < nr; j++)
			fscanf(fptr, "%lf%*c", gsl_matrix_ptr(m->pca_zmatrix, i, j));
	
	for(i = 0; i < nr; i++)
		m->pca_model_array[i] = load_modelstruct_2(fptr);

	/* fill in the mean vector */
	for(i = 0; i < nt; i++){
		col_view = gsl_matrix_column(m->training_matrix, i);
		mean_temp = vector_elt_sum(&col_view.vector, nmodel_points);
		gsl_vector_set(m->training_mean, i, (mean_temp/((double)nmodel_points)) );
	}

	return m;
}

/**
 * return the sum of the elements of vec from 0:nstop
 */
double vector_elt_sum(gsl_vector* vec, int nstop)
{
	assert(nstop >= 0); 
	assert((unsigned)nstop <= vec->size);
	int i;
	double sum = 0.0;
	for(i = 0; i < nstop; i++){
		sum += gsl_vector_get(vec, i);
	}
	return(sum);
}


/**
 * this free's everything in m
 */
void free_multimodelstruct(multi_modelstruct *m)
{
	int i;
	gsl_vector_free(m->training_mean);
	for(i = 0; i < m->nr; i++){
		free_modelstruct_2(m->pca_model_array[i]);
		//gsl_matrix_free(m->pca_model_array[i]->xmodel);
	}
	free(m->pca_model_array);
	gsl_matrix_free(m->xmodel);
	gsl_matrix_free(m->training_matrix);
	gsl_vector_free(m->pca_evals_r);
	gsl_matrix_free(m->pca_evecs_r);
	gsl_matrix_free(m->pca_zmatrix);
}
