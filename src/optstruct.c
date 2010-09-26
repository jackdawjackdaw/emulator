
#include "optstruct.h"

void free_optstruct(optstruct *opts){
	gsl_matrix_free(opts->grad_ranges);
}

void copy_optstruct(optstruct *dst, optstruct* src){
	dst->nthetas = src->nthetas;
	dst->nparams = src->nparams;
	dst->nmodel_points = src->nmodel_points;
	dst->nemulate_points = src->nemulate_points;
	dst->nregression_fns = src->nregression_fns;
	strcpy(dst->filename, src->filename);
	strcpy(dst->outputfile, src->outputfile);
	dst->cov_fn_alpha = src->cov_fn_alpha;
	dst->covariance_fn = src->covariance_fn;
	
	// this needs to be sized by the number of hyperparams
	dst->grad_ranges = gsl_matrix_alloc(src->nthetas, 2);
	gsl_matrix_memcpy(dst->grad_ranges, src->grad_ranges);
}
