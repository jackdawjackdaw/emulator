#include "string.h"
#include "assert.h"
#include "math.h"
#include "modelstruct.h"
#include "ioread.h"
#include "libEmu/regression.h"
#include "libEmu/emulator.h"
#include "libEmu/maxmultimin.h"

/**
 * allocate a modelstruct from the params in optstruct
 */
void alloc_modelstruct(modelstruct* the_model, optstruct* options){
	// row x column
	the_model->xmodel = gsl_matrix_alloc(options->nmodel_points, options->nparams);
	the_model->training_vector = gsl_vector_alloc(options->nmodel_points);
	the_model->thetas = gsl_vector_alloc(options->nthetas);
	the_model->sample_scales = gsl_vector_alloc(options->nparams);
}

/**
 * free a modelstruct
 */
void free_modelstruct(modelstruct* the_model){
	gsl_matrix_free(the_model->xmodel);
	gsl_vector_free(the_model->training_vector);
	gsl_vector_free(the_model->thetas);
	gsl_vector_free(the_model->sample_scales);
}

/**
 * copy a modelstruct from src->dst
 */
void copy_modelstruct(modelstruct* dst, modelstruct* src){
	gsl_matrix_memcpy(dst->xmodel, src->xmodel);
	gsl_vector_memcpy(dst->training_vector, src->training_vector);
	gsl_vector_memcpy(dst->thetas, src->thetas);
	gsl_vector_memcpy(dst->sample_scales, src->sample_scales);
}
	

/**
 * dump a model struct to fptr, in ascii each field is dumped in order they 
 * are defined.
 * vectors take a single line and matricies take nrows lines
 * we use the optstruct to get the sizes of everything
 */
void dump_modelstruct(FILE *fptr, modelstruct* the_model, optstruct *opts){
	int i,j;
	int nparams = opts->nparams;
	int nmp = opts->nmodel_points;
	int nthetas = opts->nthetas;

	for(i = 0; i < nmp; i++){
		for(j = 0; j < nparams; j++){
			fprintf(fptr, "%lf ", gsl_matrix_get(the_model->xmodel, i, j));
		}
		fprintf(fptr, "\n");
	}

	for(i = 0; i < nmp; i++)
		fprintf(fptr, "%lf ", gsl_vector_get(the_model->training_vector, i));
	
	for(i = 0; i < nthetas; i++)
		fprintf(fptr, "%lf ", gsl_vector_get(the_model->thetas, i));
	
	for(i = 0; i < nparams; i++)
		fprintf(fptr, "%lf ", gsl_vector_get(the_model->sample_scales, i));
	
}

/**
 * load a model struct from fptr, we use the optstruct to allocate 
 * the fields in the supplied modelstruct before filling them
 */
void load_modelstruct(FILE* fptr, modelstruct* the_model, optstruct* opts){
	int nparams = opts->nparams;
	int nmp = opts->nmodel_points;
	int nthetas = opts->nthetas;
	int i,j;
	
	// allocate everything first
	the_model->xmodel = gsl_matrix_alloc(nmp, nparams);
	the_model->training_vector = gsl_vector_alloc(nmp);
	the_model->thetas = gsl_vector_alloc(nthetas);
	the_model->sample_scales = gsl_vector_alloc(nparams);


	for(i = 0; i < nmp; i++){
		for(j = 0; j < nparams; j++){
			scan_double(fptr,gsl_matrix_ptr(the_model->xmodel, i, j));
		}
	}

	for(i = 0; i < nmp; i++){
		scan_double(fptr,gsl_vector_ptr(the_model->training_vector, i));
	}
	
	for(i = 0; i < nthetas; i++){
		scan_double(fptr,gsl_vector_ptr(the_model->thetas, i));
	}
	
	for(i = 0; i < nparams; i++){
		scan_double(fptr,gsl_vector_ptr(the_model->sample_scales, i));
	}

}



/**
 * fill a modelstruct, from a big happy array of chars from the stdin
 */
void fill_modelstruct(modelstruct* the_model, optstruct* options, char** input_data, int number_lines){
	int i,j;
	double temp_value;
	char* split_string;
	char buffer[128];
	gsl_vector *differences = gsl_vector_alloc(options->nmodel_points - 1);
	double min_value, average_value;

	// there's a bug, this can't handle empty lines at the end of the input!
	for(i = 0; i < options->nmodel_points; i++){
		split_string = strtok(input_data[i], "\t ");		
		for(j=0; j < options->nparams; j++){
			//printf("%s\n", split_string);
			// split string into tab or space tokens
			// each time you do it split_string is pointed to the next block
			// it will come up null when you're done
			assert(split_string != NULL);
			sscanf(split_string, "%lg", &temp_value);
			//fprintf(stderr,"param: %s\n", split_string);
			gsl_matrix_set(the_model->xmodel, i, j, temp_value);
			split_string = strtok(NULL, "\t ");
		}
		assert(split_string != NULL);
		sscanf(split_string,"%lg", &temp_value);
		//fprintf(stderr,"train: %s\n", split_string);
		gsl_vector_set(the_model->training_vector, i, temp_value);
	}

	// compute the average separations
	for(i = 0; i < options->nparams; i++){
		average_value = 0;
		for(j = 0; j < (options->nmodel_points-1); j++){
			gsl_vector_set(differences, j, fabs(gsl_matrix_get(the_model->xmodel, j+1, i) - 
																					gsl_matrix_get(the_model->xmodel, j, i))) ;
			average_value += gsl_vector_get(differences, j);
		}
		// compute the min difference
		min_value = gsl_vector_min(differences);
		// compute the average difference
		average_value /= (options->nmodel_points-1);
		gsl_vector_set(the_model->sample_scales, i, min_value);
		fprintf(stderr, "# param %d min-value %lf average %lf\n", i, min_value, average_value);
	}

	/* turn this off */
	/* fprintf(stderr, "read the following input matrix: %d x %d\n", options->nmodel_points, options->nparams); */
	/* message(buffer, 2); */
	/* print_matrix(the_model->xmodel, options->nmodel_points, options->nparams); */
	/* fprintf(stderr, "the training data is:\n"); */
	/* print_vector_quiet(the_model->training_vector, options->nmodel_points); */

	gsl_vector_free(differences);
}


/*********************************************************************
copied from src/libRbind/rbind.c
Fills in model->sample_scales vector, based on model->xmodel.
*********************************************************************/
gsl_vector * fill_sample_scales(gsl_matrix* xmodel) {
	if (xmodel == NULL)
		return NULL;
	int i, j;
	int nmodel_points = xmodel->size1;
	int nparams = xmodel->size2;
	gsl_vector * sample_scales = gsl_vector_alloc(nparams);
	double min_value, value;
	for(i = 0; i < nparams; i++) {
		min_value = fabs(
			gsl_matrix_get(xmodel, 1, i) -
			gsl_matrix_get(xmodel, 0, i));
		for(j = 1; j < (nmodel_points - 1); j++) {
			value = fabs(
				gsl_matrix_get(xmodel, j + 1, i) -
				gsl_matrix_get(xmodel, j,     i));
			if (value < min_value)
				min_value = value;
		}
		if(min_value < 1.0e-5)
			min_value = 1.0e-5;
		gsl_vector_set(sample_scales, i, min_value);
	}
	return sample_scales;
}

/*********************************************************************
Set some global variables: makeHVector, covariance_fn, and
makeGradMatLength.  Copied from optstruct.c setup_cov_fn() and
setup_regression().
*********************************************************************/
void set_global_ptrs(modelstruct * model) {
	switch (model->options->regression_order) {
	case 1:
		makeHVector = &(makeHVector_linear);
		model->makeHVector = &(makeHVector_linear);
		break;
	case 2:
		makeHVector = &(makeHVector_quadratic);
		model->makeHVector = &(makeHVector_quadratic);
		break;
	case 3:
		makeHVector = &(makeHVector_cubic);
		model->makeHVector = &(makeHVector_cubic);
		break;
	default:
		makeHVector = &(makeHVector_trivial);
		model->makeHVector = &(makeHVector_trivial);
	}
	switch(model->options->cov_fn_index){
	case MATERN32:
		model->covariance_fn = &(covariance_fn_matern_three);
		model->makeGradMatLength = &(derivative_l_matern_three);
		covariance_fn = &(covariance_fn_matern_three);
		makeGradMatLength = &(derivative_l_matern_three);
		break;
	case MATERN52:
		model->covariance_fn = &(covariance_fn_matern_five);
		model->makeGradMatLength = &(derivative_l_matern_five);
		covariance_fn = &(covariance_fn_matern_five);
		makeGradMatLength = &(derivative_l_matern_five);
		break;
	default:
		model->covariance_fn = &(covariance_fn_gaussian);
		model->makeGradMatLength = &(derivative_l_gauss);
		covariance_fn = &(covariance_fn_gaussian);
		makeGradMatLength = &(derivative_l_gauss);
	}
}

/*********************************************************************
Inspired by the functions in src/libRbind/rbind.c, but simplified.

Allocates and populates both modelstruct and optstruct.

@param xmodel: (n x d) matrix containing the training points.
@param training_vector: n-size vector containing the training values.
@param cov_fn_index:  POWEREXPCOVFN, MATERN32, or MATERN52
@param regression_order:  0, 1, 2, or 3

Does not estimate the thetas, since that is a labor-intensive.

Sets global variables.  I'd like to eliminate those globals and move
that information into the options structure.  Global variables means
that we can't have two models in use at once.

ccs, the fnptrs are now also in the modelstruct, the rub is that changing the 
estimation process to use the fnptrs will break the Rlibrary which is not ideal
so estimation remains non-thread-safe but sampling the mean/variance at different locations 
is safe if you use the emulator_struct form
*********************************************************************/
modelstruct * alloc_modelstruct_2(
		gsl_matrix * xmodel,
		gsl_vector * training_vector,
		int cov_fn_index,
		int regression_order) {
	assert(training_vector->size == xmodel->size1);
	assert(training_vector->size > 0);
	assert(xmodel->size2 > 0);

	/* Read parameters from dimensions of xmodel */
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

	modelstruct * model = (modelstruct*) malloc(sizeof(modelstruct));
	model->options = (optstruct*) malloc(sizeof(optstruct));

	model->options->nparams = nparams;
	model->options->nmodel_points =  nmodel_points;
	model->options->nthetas = nthetas;
	model->options->cov_fn_index = cov_fn_index;
	model->options->regression_order = regression_order;
	model->options->grad_ranges = gsl_matrix_alloc(nthetas, 2);
	model->options->nregression_fns = 1 + (regression_order * nparams);


	/* Set some global variables: makeHVector, covariance_fn, and makeGradMatLength */
	/* still need to set at least makeGradMatLength at the global scale
	 * or estimation will crash
	 */
	set_global_ptrs(model);

	
	/* alloc_modelstruct replacement code */
	model->xmodel = xmodel;
	model->training_vector = training_vector;
	model->thetas = gsl_vector_alloc(nthetas);

	model->sample_scales = fill_sample_scales(model->xmodel);

	setup_optimization_ranges(model->options, model);
	return model;
}


/*********************************************************************
@param model: pointer to the modelstruct to be freed.

Does not free model->xmodel or model->training_vector
since alloc_modelstruct_2() doesn't take "ownership" of those
data structures.
*********************************************************************/
void free_modelstruct_2(modelstruct * model) {
	/* gsl_matrix_free(model->xmodel); */
	/* gsl_vector_free(model->training_vector); */
	gsl_vector_free(model->thetas);
	gsl_vector_free(model->sample_scales);
	gsl_matrix_free(model->options->grad_ranges);
	free((void *)(model->options));
	free((void *)model);
}


/*********************************************************************
Dump a modelstruct+optstruct to fptr in ASCII.  Inverse of
load_modelstruct_2.
*********************************************************************/
void dump_modelstruct_2(FILE *fptr, modelstruct* the_model){
	int i,j;
	int nparams = the_model->options->nparams;
	int nmodel_points = the_model->options->nmodel_points;
	int nthetas = the_model->options->nthetas;

	fprintf(fptr, "%d\n", nthetas);
	fprintf(fptr, "%d\n", nparams);
	fprintf(fptr, "%d\n", nmodel_points);
	fprintf(fptr, "%d\n", the_model->options->nemulate_points);
	fprintf(fptr, "%d\n", the_model->options->regression_order);
	fprintf(fptr, "%d\n", the_model->options->nregression_fns);
	fprintf(fptr, "%d\n", the_model->options->fixed_nugget_mode);
	fprintf(fptr, "%.17lf\n", the_model->options->fixed_nugget);
	fprintf(fptr, "%d\n", the_model->options->cov_fn_index);
	fprintf(fptr, "%d\n", the_model->options->use_data_scales);
	for(i = 0; i < nthetas; i++)
		fprintf(fptr, "%.17lf %.17lf\n",
			gsl_matrix_get(the_model->options->grad_ranges, i, 0),
			gsl_matrix_get(the_model->options->grad_ranges, i, 1));
	for(i = 0; i < nmodel_points; i++){
		for(j = 0; j < nparams; j++)
			fprintf(fptr, "%.17lf ", gsl_matrix_get(the_model->xmodel, i, j));
		fprintf(fptr, "\n");
	}
	for(i = 0; i < nmodel_points; i++)
		fprintf(fptr, "%.17lf ", gsl_vector_get(the_model->training_vector, i));
	fprintf(fptr, "\n");
	for(i = 0; i < nthetas; i++)
		fprintf(fptr, "%.17lf ", gsl_vector_get(the_model->thetas, i));
	fprintf(fptr, "\n");
	for(i = 0; i < nparams; i++)
		fprintf(fptr, "%.17lf ", gsl_vector_get(the_model->sample_scales, i));
	fprintf(fptr, "\n");
}


/*********************************************************************
Load a modelstruct+optstruct from fptr. Inverse of dump_modelstruct_2.

Sets global variables.  I'd like to eliminate those globals and move
that information into the options structure.  Global variables means
that we can't have two models in use at once.
*********************************************************************/
modelstruct* load_modelstruct_2(FILE *fptr) {
	modelstruct* model = (modelstruct*)malloc(sizeof(modelstruct));
	model->options = (optstruct*)malloc(sizeof(optstruct));

	int i,j;
	int nparams, nmodel_points, nthetas;

	
	scan_int(fptr, & nthetas);
	scan_int(fptr, & nparams);
	scan_int(fptr, & nmodel_points);

	model->options->nparams = nparams;
	model->options->nmodel_points = nmodel_points;
	model->options->nthetas = nthetas;

	scan_int(fptr, &(model->options->nemulate_points));
	scan_int(fptr, &(model->options->regression_order));
	scan_int(fptr, &(model->options->nregression_fns));
	scan_int(fptr, &(model->options->fixed_nugget_mode));
	scan_double(fptr, &(model->options->fixed_nugget));
	scan_int(fptr, &(model->options->cov_fn_index));
	scan_int(fptr, &(model->options->use_data_scales));

	model->options->grad_ranges = gsl_matrix_alloc(nthetas, 2);
	for(i = 0; i < nthetas; i++) {
		scan_double(fptr, gsl_matrix_ptr(model->options->grad_ranges,i,0));
		scan_double(fptr, gsl_matrix_ptr(model->options->grad_ranges,i,1));
	}

	model->xmodel = gsl_matrix_alloc(nmodel_points, nparams);
	for(i = 0; i < nmodel_points; i++)
		for(j = 0; j < nparams; j++)
			scan_double(fptr, gsl_matrix_ptr(model->xmodel, i, j));

	model->training_vector = gsl_vector_alloc(nmodel_points);
	for(i = 0; i < nmodel_points; i++)
		scan_double(fptr, gsl_vector_ptr(model->training_vector, i));

	model->thetas = gsl_vector_alloc(nthetas);
	for(i = 0; i < nthetas; i++)
		scan_double(fptr, gsl_vector_ptr(model->thetas, i));

	model->sample_scales = gsl_vector_alloc(nparams);
	for(i = 0; i < nparams; i++)
		scan_double(fptr, gsl_vector_ptr(model->sample_scales, i));

	set_global_ptrs(model);
	return model;
}
