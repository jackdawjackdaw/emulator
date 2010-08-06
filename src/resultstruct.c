#include "resultstruct.h"

/**
 * free a resultstruct 
 */
void free_resultstruct(resultstruct *res){
	gsl_matrix_free(res->new_x);
	gsl_vector_free(res->emulated_mean);
	gsl_vector_free(res->emulated_var);
}

/** 
 * allocate the resultstruct from the options
 */
void alloc_resultstruct(resultstruct *res, optstruct* opts){
	res->new_x = gsl_matrix_alloc(opts->nemulate_points, opts->nparams);
	res->emulated_mean = gsl_vector_alloc(opts->nemulate_points);
	res->emulated_var = gsl_vector_alloc(opts->nemulate_points);
	res->options = opts;
}

/**
 * src -> dst
 */
void copy_resultstruct(resultstruct *dst, resultstruct *src){
	gsl_matrix_memcpy(dst->new_x, src->new_x);
	gsl_vector_memcpy(dst->emulated_mean, src->emulated_mean);
	gsl_vector_memcpy(dst->emulated_var, src->emulated_var);
	dst->options = src->options;
}



