#include "multivar_support.h"
#include "multi_modelstruct.h"
#include "emulator_struct.h"
#include "modelstruct.h"
#include "libEmu/estimate_threaded.h"
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>



/**
 * estimates the thetas for an allocated multvariate model and outputs the 
 * info to outfp
 * 
 * requries: 
 * m must have been succesfully allocated and pca'd
 * outfp is an open fileptr
 */
void estimate_multi(multi_modelstruct *m, FILE* outfp)
{
	int i;
	for(i = 0; i < m->nr; i++){
		estimate_thetas_threaded(m->pca_model_array[i], m->pca_model_array[i]->options);
	}
	// dump the trained modelstruct, this maybe should go somewhere else
	dump_multi_modelstruct(outfp, m);
}

multi_emulator *alloc_multi_emulator(multi_modelstruct *m)
{
	int i;
	multi_emulator *e = (multi_emulator*)malloc(sizeof(multi_emulator));
	e->nt = m->nt;
	e->nr = m->nr;
	e->nparams = m->nparams;
	e->nmodel_points = m->nmodel_points;
	// have to read these out of the pca array
	e->nregression_fns = m->pca_model_array[0]->options->nregression_fns;
	e->nthetas = m->pca_model_array[0]->options->nthetas;
	e->model = m;
	
	// allocate the nr emulators we need
	e->emu_struct_array = (emulator_struct**)malloc(sizeof(emulator_struct*)*e->nr);

	for(i = 0; i < e->nr; i++)
		e->emu_struct_array[i] = alloc_emulator_struct(m->pca_model_array[i]);
	
	return e;
}


void free_multi_emulator(multi_emulator *e)
{
	int i;
	free_multimodelstruct(e->model);
	for(i = 0; i < e->nr; i++)
		free_emulator_struct(e->emu_struct_array[i]);
	free(e->emu_struct_array);
}


/**
 * emulate the values of a point in a multivariate model
 * 
 * requries:
 * emu has been allocated from a multi_modelstruct
 * 
 * @argument emu: a multi_emulator structure created from an estimated multi_modelstruct with alloc_multi_emulator
 * @argument the_point: the location in parameter space (nparams) length vector
 * @argument the_mean: the emulated mean values (t length)
 * @argument the_variance: the emulated variance values (t length)
 */
void emulate_point_multi(multi_emulator *emu, gsl_vector *the_point,
												 gsl_vector *the_mean, gsl_vector *the_variance)
{
	int i, j;
	int nr = emu->nr;
	int nt = emu->nt;
	double vec_mat_sum;
	// the mean and variance in the PCA space
	gsl_vector *mean_pca = gsl_vector_alloc(nr);
	gsl_vector *var_pca = gsl_vector_alloc(nr); 

	// first we sample the nr emulators
	for(i = 0; i < emu->nr; i++)
		emulate_point(emu->emu_struct_array[i], the_point, gsl_vector_ptr(mean_pca, i), gsl_vector_ptr(var_pca, i));

	/**
	 *  now we have to project back into the REAL space from the PCA space
	 * mean_real (nt) = training_mean (nt) + pca_evec_r %*% diag(sqrt(pca_evals_R)) %*% mean_pca
	 */
	gsl_vector *mean_real = gsl_vector_alloc(nt);
	gsl_vector *temp = gsl_vector_alloc(nt);
	gsl_vector_memcpy(mean_real, emu->model->training_mean);

	vec_mat_sum = 0.0;
	for(i = 0; i < nt; i++){
		for(j = 0; j < nr; j++)
			vec_mat_sum += gsl_matrix_get(emu->model->pca_evecs_r, i, j) * gsl_vector_get(mean_pca, j);
		// save the sum scaled by the sqrt of the eval
		gsl_vector_set(temp, i, vec_mat_sum*sqrt(gsl_vector_get(emu->model->pca_evals_r, i)));
		vec_mat_sum = 0.0;
	}
	
	gsl_vector_add(mean_real, temp);
	gsl_vector_memcpy(the_mean, mean_real); // save the final mean

	/**
	 * project back the variance 
	 * yVar_i = Sum_k=0^{nr} ( U_i_k * U_i_k * lambda_k * V_k)
	 */
	gsl_vector_set_zero(temp);
	vec_mat_sum = 0.0;
	for(i = 0; i < nt; i++){
		for(j = 0; j < nr; j++)
			vec_mat_sum += pow(gsl_matrix_get(emu->model->pca_evecs_r, i, j), 2.0) *
				gsl_vector_get(emu->model->pca_evals_r, j) * gsl_vector_get(var_pca, j);
		gsl_vector_set(the_variance, i, vec_mat_sum);
		vec_mat_sum = 0.0;
	}
	
	gsl_vector_free(mean_real);
	gsl_vector_free(temp);
	gsl_vector_free(mean_pca);
	gsl_vector_free(var_pca);
}


