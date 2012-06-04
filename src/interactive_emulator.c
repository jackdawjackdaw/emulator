/*********************************************************************
INTERACTIVE N-D GAUSSIAN PROCESS INTERPOLATOR (MODEL EMULATOR)
Copyright 2012, The University of North Carolina at Chapel Hill.

DESCRIPTION:
  Meant to be included in the GP-Emulator program
  <https://github.com/jackdawjackdaw/emulator> by C.Coleman-Smith
  <cec24@phy.duke.edu>, Copyright 2009-2011 Duke University.

ACKNOWLEDGMENTS:
  This software was written in 2012 by Hal Canary <cs.unc.edu/~hal>,
  based off of code written by Christopher Coleman-Smith
  <cec24@phy.duke.edu> in 2010-2012 while working for the MADAI project
  <http://madai.us/>.

LICENSE:
  (Since we link to the GNU Scientific Library, we MUST use GPL.)
  This program is free software: you can redistribute it and/or modify
  it under the terms of version 3 of the GNU General Public License as
  published by the Free Software Foundation.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  A copy of version 3 of the GNU General Public License can be found at
  <http://www.gnu.org/licenses/gpl-3.0-standalone.html>.

USE:
  useage:
    interactive_emulator estimate_thetas INPUT_MODEL_FILE MODEL_SNAPSHOT_FILE [OPTIONS]
  or
    interactive_emulator interactive_mode MODEL_SNAPSHOT_FILE

  INPUT_MODEL_FILE can be "-" to read from standard input.

  The input MODEL_SNAPSHOT_FILE for interactive_mode must match the format of
  the output MODEL_SNAPSHOT_FILE from estimate_thetas.

  Options for estimate_thetas can include:
    --regression_order=0
    --regression_order=1
    --regression_order=2
    --regression_order=3
    --covariance_fn=POWER_EXPONENTIAL
    --covariance_fn=MATERN32
    --covariance_fn=MATERN52
  The defaults are regression_order=0 and covariance_fn=POWER_EXPONENTIAL.

  These options will be saved in MODEL_SNAPSHOT_FILE.

INPUT_MODEL_FILE FORMAT:
  BEGIN EXAMPLE
    nparams
    nmodel_points
    X[0,0]
    ...
    X[0,nparams-1]
    X[1,0]
    ...
    X[nmodel_points,nparams-1]
    y[0]
    ...
    y[nmodel_points-1]
  END EXAMPLE

  nparams and nmodel_points should be positive integers.  X[i,j] and
  y[i] will be read as double-precision floats.

  You don't have to use newlines to separate values.  If fscanf(fp,
  "%lf%*c", ptr) can read the value, it will work.

BUGS:
  Writing MODEL_SNAPSHOT_FILE to stdout is not allowed because
  GPEmulatorLib is excessively verbose.  This should be fixed when we
  refactor.

  This program probably should be two seperate executables, but for
  now it is self-contained in one source file.

  Multiple models can't be run simultaniously due to the use of global
  variables in GPEmulatorLib.

  "interactive_emulator estimate_thetas ..." is mostly redundant
  against "estimator", except that the input file format and
  command-line arguments are more sane (by Hal's standards).

  The functions fill_sample_scales, alloc_modelstruct_2,
  free_modelstruct_2, dump_modelstruct_2, and load_modelstruct_2
  should all be moved into modelstruct.h/.c.

BIGGEST BUG:
  What do I do when I have multiple training vectors (Y's) i.e. I'm
  trying to emulate several functions defined on the same space and
  sampled on the same training points?  We need to modify this to
  efficiently do that.

Possible TODO:
  allow the design x[0,0]...x[nmodel_points,nparams-1] to be read separately from the 
  training points
   
*********************************************************************/

/* #define BINARY_INTERACTIVE_MODE */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#undef NDEBUG
#include <assert.h>

#include "main.h"
#include "resultstruct.h"
#include "emulator_struct.h"
#include "libEmu/maxmultimin.h"

/********************************************************************/
static const char useage [] =
	"useage:\n"
	"  interactive_emulator estimate_thetas INPUT_MODEL_FILE MODEL_SNAPSHOT_FILE [OPTIONS]\n"
	"or\n"
	"  interactive_emulator interactive_mode MODEL_SNAPSHOT_FILE\n"
	"\n"
	"INPUT_MODEL_FILE can be \"-\" to read from standard input.\n"
	"\n"
	"The input MODEL_SNAPSHOT_FILE for interactive_mode must match the format of\n"
	"the output MODEL_SNAPSHOT_FILE from estimate_thetas.\n"
	"\n"
	"Options for estimate_thetas can include:\n"
	"  --regression_order=0\n"
	"  --regression_order=1\n"
	"  --regression_order=2\n"
	"  --regression_order=3\n"
	"  --covariance_fn=POWER_EXPONENTIAL\n"
	"  --covariance_fn=MATERN32\n"
	"  --covariance_fn=MATERN52\n"
	"The defaults are regression_order=0 and covariance_fn=POWER_EXPONENTIAL.\n";


/********************************************************************/
int perr(const char * s) {
	fprintf(stderr,"%s\n",s);
	return EXIT_FAILURE;
}


/*********************************************************************
copied from src/libRbind/rbind.c
Fills in model->sample_scales vector, based on model->xmodel.
*********************************************************************/
void fill_sample_scales(modelstruct* model) {
	int i, j;
	int nparams = model->options->nparams;
	int nmodel_points = model->options->nmodel_points;
	double min_value, value;
	for(i = 0; i < nparams; i++) {
		min_value = fabs(
			gsl_matrix_get(model->xmodel, 1, i) -
			gsl_matrix_get(model->xmodel, 0, i));
		for(j = 1; j < (nmodel_points - 1); j++) {
			value = fabs(
				gsl_matrix_get(model->xmodel, j + 1, i) -
				gsl_matrix_get(model->xmodel, j,     i));
			if (value < min_value)
				min_value = value;
		}
		if(min_value < 1.0e-5)
			min_value = 1.0e-5;
		gsl_vector_set(model->sample_scales, i, min_value);
	}
}

/*********************************************************************
Set some global variables: makeHVector, covariance_fn, and
makeGradMatLength.  Copied from optstruct.c setup_cov_fn() and
setup_regression().
*********************************************************************/
void set_global_ptrs(int regression_order, int cov_fn_index) {
	switch (regression_order) {
	case 1:
		makeHVector = &(makeHVector_linear);
		break;
	case 2:
		makeHVector = &(makeHVector_quadratic);
		break;
	case 3:
		makeHVector = &(makeHVector_cubic);
		break;
	default:
		makeHVector = &(makeHVector_trivial);
	}
	switch (cov_fn_index) {
	case MATERN32:
		covariance_fn = &(covariance_fn_matern_three);
		makeGradMatLength = &(derivative_l_matern_three);
		break;
	case MATERN52:
		covariance_fn = &(covariance_fn_matern_five);
		makeGradMatLength = &(derivative_l_matern_five);
		break;
	default:
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


	/** 
	 *	ccs: why do we need to set makeHVector here and then also in set_global_ptrs?
	 *  these fn pts are now in modelstruct 
	 */
	if (model->options->regression_order == 0){
		makeHVector = &(makeHVector_trivial);
		model->makeHVector = &(makeHVector_trivial);
	}
	else if  (model->options->regression_order == 1){
		makeHVector = &(makeHVector_linear);
		model->makeHVector = &(makeHVector_linear);
	}
	else if (model->options->regression_order == 2){
		makeHVector = &(makeHVector_quadratic);
		model->makeHVector = &(makeHVector_quadratic);
	}
	else if (model->options->regression_order == 3){
		makeHVector = &(makeHVector_cubic);
		model->makeHVector = &(makeHVector_cubic);
	}
	
	/**
	 * ccs: set the covfn ptrs
	 */
	switch(cov_fn_index){
	case MATERN32:
		model->covariance_fn = &(covariance_fn_matern_three);
		model->makeGradMatLength = &(derivative_l_matern_three);
		break;
	case MATERN52:
		model->covariance_fn = &(covariance_fn_matern_five);
		model->makeGradMatLength = &(derivative_l_matern_five);
		break;
	default:
		model->covariance_fn = &(covariance_fn_gaussian);
		model->makeGradMatLength = &(derivative_l_gauss);
	}


	
	/* Set some global variables: makeHVector, covariance_fn, and makeGradMatLength */
	/* still need to set at least makeGradMatLength at the global scale
	 * or estimation will crash
	 */
	set_global_ptrs(regression_order, cov_fn_index);
	

	
	/* alloc_modelstruct replacement code */
	model->xmodel = xmodel;
	model->training_vector = training_vector;
	model->thetas = gsl_vector_alloc(nthetas);
	model->sample_scales = gsl_vector_alloc(nparams);

	fill_sample_scales(model);
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

	fscanf(fptr, "%d%*c", & nthetas);
	fscanf(fptr, "%d%*c", & nparams);
	fscanf(fptr, "%d%*c", & nmodel_points);

	model->options->nparams = nparams;
	model->options->nmodel_points = nmodel_points;
	model->options->nthetas = nthetas;

	fscanf(fptr, "%d%*c", &(model->options->nemulate_points));
	fscanf(fptr, "%d%*c", &(model->options->regression_order));
	fscanf(fptr, "%d%*c", &(model->options->nregression_fns));
	fscanf(fptr, "%d%*c", &(model->options->fixed_nugget_mode));
	fscanf(fptr, "%lf%*c", &(model->options->fixed_nugget));
	fscanf(fptr, "%d%*c", &(model->options->cov_fn_index));
	fscanf(fptr, "%lf%*c", &(model->options->use_data_scales));

	model->options->grad_ranges = gsl_matrix_alloc(nthetas, 2);
	for(i = 0; i < nthetas; i++) {
		fscanf(fptr, "%lf%*c", gsl_matrix_ptr(model->options->grad_ranges,i,0));
		fscanf(fptr, "%lf%*c", gsl_matrix_ptr(model->options->grad_ranges,i,1));
	}

	model->xmodel = gsl_matrix_alloc(nmodel_points, nparams);
	for(i = 0; i < nmodel_points; i++)
		for(j = 0; j < nparams; j++)
			fscanf(fptr, "%lf%*c", gsl_matrix_ptr(model->xmodel, i, j));

	model->training_vector = gsl_vector_alloc(nmodel_points);
	for(i = 0; i < nmodel_points; i++)
		fscanf(fptr, "%lf%*c", gsl_vector_ptr(model->training_vector, i));

	model->thetas = gsl_vector_alloc(nthetas);
	for(i = 0; i < nthetas; i++)
		fscanf(fptr, "%lf%*c", gsl_vector_ptr(model->thetas, i));

	model->sample_scales = gsl_vector_alloc(nparams);
	for(i = 0; i < nparams; i++)
		fscanf(fptr, "%lf%*c", gsl_vector_ptr(model->sample_scales, i));

	set_global_ptrs(model->options->regression_order, model->options->cov_fn_index);
	return model;
}


/*********************************************************************
Reads a file containing xmodel and training_vector info and create
data structures to hold that information.
@param input_filename: name of file to open.
  can be "-" or "stdin" to represent standard input
@param xmodel_ptr: returns a newly allocated gsl_matrix
@param training_vector_ptr: returns a newly allocated gsl_vector
*********************************************************************/
int open_model_file(char * input_filename,
		gsl_matrix ** xmodel_ptr, gsl_vector ** training_vector_ptr) {
	FILE * input_file;
	if ((0 == strcmp(input_filename, "-")) || (0 == strcmp(input_filename, "stdin")))
		input_file = stdin;
	else
		input_file = fopen(input_filename, "r");
	if (input_file == NULL)
		return 0; /* failure */
	int i, j, nparams, nmodel_points;
	fscanf(input_file,"%d%*c", & nparams);
	fscanf(input_file,"%d%*c", & nmodel_points);
	gsl_matrix * xmodel = gsl_matrix_alloc(nmodel_points, nparams);
	gsl_vector * training_vector = gsl_vector_alloc(nmodel_points);
	for (i = 0; i < nmodel_points; i++)
		for (j = 0; j < nparams; j++)
			fscanf(input_file, "%lf%*c", gsl_matrix_ptr(xmodel,i,j));
	for (i = 0; i < nmodel_points; i++)
		fscanf(input_file, "%lf%*c", gsl_vector_ptr(training_vector, i));
	if (input_file != stdin)
		fclose(input_file);
	*xmodel_ptr = xmodel;
	*training_vector_ptr = training_vector;
	return 1;
}


/*********************************************************************
return mean and variance at a point.
thread-safe, uses fnptrs defined in e->modelstruct
*********************************************************************/
void emulate_point(emulator_struct* e, gsl_vector * point,  double * mean, double * variance){
	gsl_vector * kplus = gsl_vector_alloc(e->nmodel_points); //allocate variable storage
	gsl_vector * h_vector = gsl_vector_alloc(e->nregression_fns); //allocate variable storage

	makeKVector_es(kplus, point, e);
	/* use the fnptr in modelstruct */
	e->model->makeHVector(h_vector, point, e->nparams);
	
	(*mean) = makeEmulatedMean(e->cinverse, e->model->training_vector,
		kplus, h_vector, e->h_matrix, e->beta_vector, e->nmodel_points);
	double kappa = covariance_fn(point, point, e->model->thetas, e->nthetas, e->nparams);
	(*variance) = makeEmulatedVariance(e->cinverse, kplus, h_vector,
		e->h_matrix, kappa, e->nmodel_points, e->nregression_fns);
	gsl_vector_free(kplus);
	gsl_vector_free(h_vector);
	return;
}



/*********************************************************************
Return true only if the beginning of s1 matches s2.
*********************************************************************/
#define starts_with(s1,s2) (strncmp ((s1), (s2), strlen(s2)) == 0)


/*********************************************************************
Return true only if s1 equals s2.
*********************************************************************/
#define str_equal(s1,s2) (strcmp ((s1), (s2)) == 0)


/*********************************************************************
If main called with "estimate_thetas" argument, this happens.
*********************************************************************/
int estimate_thetas (int argc, char ** argv) {
	/* read a file into a GSL_matrix and a GSL_vector */
	gsl_matrix * xmodel;
	gsl_vector * training_vector;

	if (argc < 2)
		return perr("Not enough arguments\n");

	if (! open_model_file(argv[0], &xmodel, &training_vector))
		return perr("File read failed.");

	/*
	FILE * outfp;
	if ((argc < 2) || (0 == strcmp(argv[1], "-")) || (0 == strcmp(argv[1], "stdout")))
		outfp = stdout;
	else
		outfp = fopen(argv[1],"w"); */
	FILE * outfp = fopen(argv[1], "w");
	if (outfp == NULL)
		return perr("Opening output file failed.");

	int cov_fn_index = POWEREXPCOVFN; /* POWEREXPCOVFN: 1, MATERN32: 2, MATERN52: 3 */
	int regression_order = 0; /* 0, 1, 2, or 3 */
	argc -= 2;
	argv += 2;
	while (argc > 0) {
		if (starts_with(argv[0], "--covariance_fn=")) {
			if (starts_with(&(argv[0][16]), "POWER_EXPONENTIAL")) {
				cov_fn_index = POWEREXPCOVFN; /* 1 */
			} else if (starts_with(&(argv[0][16]), "MATERN32")) {
				cov_fn_index = MATERN32; /* 2 */
			} else if (starts_with(&(argv[0][16]), "MATERN52")) {
				cov_fn_index = MATERN52;  /* 3 */
			} else {
				fprintf(stderr, "Unknown covariance_function \"%s\".\n", &(argv[0][16]));
				return perr(useage);
			}
		} else if (starts_with(argv[0], "--regression_order=")) {
			if (str_equal(&(argv[0][19]), "0")) {
				regression_order = 0;
			} else if (str_equal(&(argv[0][19]), "1")) {
				regression_order = 1;
			} else if (str_equal(&(argv[0][19]), "2")) {
				regression_order = 2;
			} else if (str_equal(&(argv[0][19]), "3")) {
				regression_order = 3;
			} else {
				fprintf(stderr, "Invalid regression_order: \"%s\"\n", &(argv[0][19]));
				return perr(useage);
			}
		}
		else {
			fprintf(stderr, "invalid option: \"%s\"\n", argv[0]);
			return perr(useage);
		}
		argc--;
		argv++;
	}

	modelstruct * model = alloc_modelstruct_2(xmodel, training_vector,
		cov_fn_index, regression_order);

	/* actually do the estimation using libEmu */
	estimate_thetas_threaded(model, model->options);

	/* write to file */
	dump_modelstruct_2(outfp, model);
	fclose(outfp);

	gsl_matrix_free(xmodel);
	gsl_vector_free(training_vector);
	free_modelstruct_2(model);
	return 0;
}

/*********************************************************************
If main called with "interactive_mode" argument, this happens.
*********************************************************************/
int interactive_mode (int argc, char** argv) {
	FILE * interactive_input = stdin;
	FILE * interactive_output = stdout;
	int i;
	double the_mean, the_variance;

	FILE * fp = fopen(argv[0],"r");
	if (fp == NULL)
		return perr("Error opening file");
	modelstruct* model = load_modelstruct_2(fp);
	fclose(fp);

	emulator_struct * the_emulator = alloc_emulator_struct(model);

	int nparams = model->options->nparams;
	gsl_vector * the_point = gsl_vector_alloc(nparams);

	#ifdef BINARY_INTERACTIVE_MODE
		fwrite(&nparams, sizeof(int), 1, interactive_output);
	#else
		fprintf(interactive_output,"%d\n",nparams);
	#endif
	fflush(interactive_output);
	while (! feof(interactive_input)) {
		for(i =0; i < nparams; i++) {
			#ifdef BINARY_INTERACTIVE_MODE
				fread(gsl_vector_ptr(the_point, i), sizeof(double), 1, interactive_input);
			#else
				fscanf(interactive_input, "%lf%*c", gsl_vector_ptr(the_point, i));
			#endif
		}
		emulate_point(the_emulator, the_point, &the_mean, &the_variance);
		#ifdef BINARY_INTERACTIVE_MODE
			fwrite(&the_mean, sizeof(double), 1, interactive_output);
			fwrite(&the_variance, sizeof(double), 1, interactive_output);
		#else
			fprintf(interactive_output, "%.17f\n%.17f\n", the_mean, the_variance);
		#endif
		fflush(interactive_output);
	}
	free_emulator_struct(the_emulator);
	gsl_matrix_free(model->xmodel);
	model->xmodel = NULL;
	gsl_vector_free(model->training_vector);
	model->training_vector = NULL;
	free_modelstruct_2(model);
	gsl_vector_free(the_point);
	return 0;
}


/*********************************************************************
main
*********************************************************************/
int main (int argc, char ** argv) {
	if (argc < 3)
		return perr(useage);
	if (0 == strcmp(argv[1], "estimate_thetas"))
		return estimate_thetas(argc - 2, argv + 2);
	if (0 == strcmp(argv[1], "interactive_mode"))
		return interactive_mode(argc - 2, argv + 2);
	return perr(useage);
}
