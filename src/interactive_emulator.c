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
    --covariance_fn=0 (POWER_EXPONENTIAL)
    --covariance_fn=1 (MATERN32)
    --covariance_fn=2 (MATERN52)
  The defaults are regression_order=0 and covariance_fn=POWER_EXPONENTIAL.

  These options will be saved in MODEL_SNAPSHOT_FILE.

INPUT_MODEL_FILE (With multi-output) FORMAT:
  Models with multivariate output values y = y_1...y_t which we will think of as rows of a matrix, in the following spec
  number_outputs = t
  BEGIN EXAMPLE
    number_outputs
    number_params 
    number_model_points
    X[0,0]
    ...
    X[0,number_params-1]
    X[1,0]
    ...
    X[number_model_points-1,number_params-1]
    Y[0, 0]
    Y[0, 1]
    ...
    Y[0, number_outputs-1]
    ...
    Y[1, 1]
    ... 
    Y[number_model_points-1, number_outputs-1]
   END EXAMPLE

  number_outputs, number_params and number_model_points should be
  positive integers.  X[i,j] and Y[i,j] will be read as
  double-precision floats.

  You don't have to use newlines to separate values.  If fscanf(fp,
  "%lf%*c", ptr) can read the value, it will work.

BUGS:
  Writing MODEL_SNAPSHOT_FILE to stdout is not allowed because
  GPEmulatorLib is excessively verbose.  This should be fixed when we
  refactor.

  This program probably should be two seperate executables, but for
  now it is self-contained in one source file.

  Multiple models can't be run simultaniously due to the use of global
  variables in GPEmulatorLib.  This is being fixed.

  "interactive_emulator estimate_thetas ..." is mostly redundant
  against "estimator", except that the input file format and
  command-line arguments are different.

BIGGEST BUG:
  What do I do when I have multiple training vectors (Y's) i.e. I'm
  trying to emulate several functions defined on the same space and
  sampled on the same training points?  We need to modify this to
  efficiently do that.
  
  this is somewhat supported now
  

 TODO:
  - allow the design x[0,0]...x[nmodel_points,nparams-1] to be read separately from the 
  training points
  - check for scaling of design and training values, warn if they different columns have very different scales
  - make estimation process respect the quiet flag?
  

   
*********************************************************************/

/* #define BINARY_INTERACTIVE_MODE */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>

#include "main.h"
#include "modelstruct.h"
#include "emulator_struct.h"
#include "multi_modelstruct.h"
#include "multivar_support.h"
#include "libEmu/maxmultimin.h"



/** 
 * contents will be set by global_opt_parse, can be passed into estimate_thetas and interactive_mode 
 * to be used appropriately
 */
struct cmdLineOpts{
	int regOrder;  /* -r --regression_order = */
	int covFn;     /* -c --covariance_fn = */
	int quietFlag; /* -q --quiet */
	int pcaOutputFlag; /* -z --pca_output  turns on pca only output */
	double pca_variance; /* -v --pca_variance =  fractional value for the pca decomp */
	
	// add additional flags as needed here

	// what mode to run in
	char* run_mode; /* arg 1 */
	char *inputfile; /* arg 2 (maybe) */
	// the state filename
	char* statefile; /* arg 2 or 3 */
} cmdLineOpts;


struct cmdLineOpts* global_opt_parse(int argc, char** argv);

/********************************************************************/
static const char useage [] =
	"useage:\n"
	"  interactive_emulator estimate_thetas INPUT_MODEL_FILE MODEL_SNAPSHOT_FILE [OPTIONS]\n"
	"or\n"
	"  interactive_emulator interactive_mode MODEL_SNAPSHOT_FILE [OPTIONS]\n"
	"\n"
	"INPUT_MODEL_FILE can be \"-\" to read from standard input.\n"
	"\n"
	"The input MODEL_SNAPSHOT_FILE for interactive_mode must match the format of\n"
	"the output MODEL_SNAPSHOT_FILE from estimate_thetas.\n"
	"\n"
	"Options which only influence estimate_thetas:\n"
	"  --regression_order=0 (const)\n"
	"  --regression_order=1 (linear)\n"
	"  --regression_order=2 (quadratic)\n"
	"  --regression_order=3 (cubic)\n"
	"  --covariance_fn=0 (POWER_EXPONENTIAL)\n"
	"  --covariance_fn=1 (MATERN32)\n"
	"  --covariance_fn=2 (MATERN52)\n"
	"  (-v=FRAC) --pca_variance=FRAC : sets the pca decomp to keep cpts up to variance fraction frac\n"
	"options which incluence interactive_mode:\n"
	"  (-q) --quiet: run without any extraneous output\n"
	"  (-z) --pca_output: emulator output is left in the pca space\n"
	"general options:\n"
	"  -h -? print this dialogue\n"
	"The defaults are regression_order=0 and covariance_fn=POWER_EXPONENTIAL.\n";

/**
 * Convienence function for exiting a function.
 */
int perr(const char * s) {
	fprintf(stderr,"%s\n",s);
	return EXIT_FAILURE;
}


/**
 * Reads a file containing xmodel and training_matrix info and create
 * data structures to hold that information.
 * @param input_filename: name of file to open.
 *   can be "-" or "stdin" to represent standard input
 * @param xmodel_ptr: returns a newly allocated gsl_matrix (nmodel_points x nparams)
 * @param training_matrix_ptr: returns a newly allocated gsl_matrix (nmodel_points x nt)
 * @returns 1 on sucess and 0 on error.
 * 
 * Example of use:
 * 	gsl_matrix * xmodel;
 * 	gsl_matrix * training_matrix;
 * 	open_model_file(filename, &xmodel, &training_matrix)
 */
int open_model_file(char * input_filename,
		gsl_matrix ** xmodel_ptr, gsl_matrix ** training_matrix_ptr) {
	FILE * input_file;
	if ((0 == strcmp(input_filename, "-")) || (0 == strcmp(input_filename, "stdin")))
		input_file = stdin;
	else
		input_file = fopen(input_filename, "r");
	if (input_file == NULL)
		return 0; /* failure */
	int i, j, number_outputs, number_params, number_model_points;
	fscanf(input_file,"%d%*c", & number_outputs);
	fscanf(input_file,"%d%*c", & number_params);
	fscanf(input_file,"%d%*c", & number_model_points);

	assert(number_outputs > 0);
	assert(number_params > 0);
	assert(number_model_points > 0);

	gsl_matrix * xmodel = gsl_matrix_alloc(number_model_points, number_params);
	gsl_matrix * training_matrix = gsl_matrix_alloc(number_model_points, number_outputs);
	
	for (i = 0; i < number_model_points; i++)
		for (j = 0; j < number_params; j++)
			fscanf(input_file, "%lf%*c", gsl_matrix_ptr(xmodel,i,j));
	for (i = 0; i < number_model_points; i++)
		for(j = 0; j < number_outputs; j++)
			fscanf(input_file, "%lf%*c", gsl_matrix_ptr(training_matrix, i, j));
	
	if (input_file != stdin)
		fclose(input_file);
	*xmodel_ptr = xmodel;
	*training_matrix_ptr = training_matrix;
	return 1;
}




/**
 * Return true only if the beginning of s1 matches s2.
 */
#define starts_with(s1,s2) (strncmp ((s1), (s2), strlen(s2)) == 0)

/**
 * Return true only if s1 equals s2.
 */
#define str_equal(s1,s2) (strcmp ((s1), (s2)) == 0)

/**
 * read out the cov-fn and regression order from argc and argv
 * 
 * ccs: can we use getopt_long here?
 */
void parse_arguments_interactive(int* cov_fn_index, int* regression_order, int argc, char ** argv){
	argc -= 2;
	argv += 2;
	while (argc > 0) {
		if (starts_with(argv[0], "--covariance_fn=")) {
			if (starts_with(&(argv[0][16]), "POWER_EXPONENTIAL")) {
				*cov_fn_index = POWEREXPCOVFN; /* 1 */
			} else if (starts_with(&(argv[0][16]), "MATERN32")) {
				*cov_fn_index = MATERN32; /* 2 */
			} else if (starts_with(&(argv[0][16]), "MATERN52")) {
				*cov_fn_index = MATERN52;  /* 3 */
			} else {
				fprintf(stderr, "Unknown covariance_function \"%s\".\n", &(argv[0][16]));
				exit(perr(useage));
			}
		} else if (starts_with(argv[0], "--regression_order=")) {
			if (str_equal(&(argv[0][19]), "0")) {
				*regression_order = 0;
			} else if (str_equal(&(argv[0][19]), "1")) {
				*regression_order = 1;
			} else if (str_equal(&(argv[0][19]), "2")) {
				*regression_order = 2;
			} else if (str_equal(&(argv[0][19]), "3")) {
				*regression_order = 3;
			} else {
				fprintf(stderr, "Invalid regression_order: \"%s\"\n", &(argv[0][19]));
				exit(perr(useage));
			}
		}
		else {
			fprintf(stderr, "invalid option: \"%s\"\n", argv[0]);
			exit(perr(useage));
		}
		argc--;
		argv++;
	}
}


/**
 * if main called with "estimate_thetas_multi"
 */
int estimate_thetas(struct cmdLineOpts* cmdOpts) {
	gsl_matrix *xmodel = NULL;
	gsl_matrix *training_matrix = NULL;
	double varfrac = 0.95; // this could be set by an arg

	/* if (argc < 2) */
	/* 	return perr("Not enough arguments\n"); */

	if (! open_model_file(cmdOpts->inputfile, &xmodel, &training_matrix))
		return perr("Input File read failed.");

	FILE * outfp = fopen(cmdOpts->statefile, "w");
	if (outfp == NULL)
		return perr("Opening statefile failed.");

	int cov_fn_index = POWEREXPCOVFN; /* POWEREXPCOVFN: 1, MATERN32: 2, MATERN52: 3 */
	int regression_order = 0; /* 0, 1, 2, or 3 */
	
	cov_fn_index = cmdOpts->covFn;
	regression_order = cmdOpts->regOrder;

	if(cmdOpts->pca_variance <= 1.0 && cmdOpts->pca_variance > 0){
		varfrac = cmdOpts->pca_variance; // use the option from the cmd line
	}
	
	if(cov_fn_index < 0 || cov_fn_index > 3){
		fprintf(stderr, "#ERROR cov_fn_index %d not supported\n", cov_fn_index);
		exit(EXIT_FAILURE);
	}

	if(regression_order < 0 || regression_order > 3){
		fprintf(stderr, "#ERROR regression_order %d not supported\n", cov_fn_index);
		exit(EXIT_FAILURE);
	}

	/* allocate the multi-model, do the pca decomp
	 * this is a little chatty on stderr
	 */

	multi_modelstruct * model = NULL;
	model = alloc_multimodelstruct(
		xmodel, training_matrix,
		cov_fn_index, regression_order, varfrac);

	if(model == NULL)
		return perr("Failed to allocated multi_modelstruct.\n");
	
	estimate_multi(model, outfp);

	fclose(outfp);

	//gsl_matrix_free(xmodel);
	//gsl_matrix_free(training_matrix);
	free_multimodelstruct(model);
	return EXIT_SUCCESS;
}


/**
 * If main called with "interactive_mode" argument, this happens.
 */
int interactive_mode (struct cmdLineOpts* cmdOpts) {
	FILE * interactive_input = stdin;
	FILE * interactive_output = stdout;
	int i, r, expected_r;

	FILE * fp = fopen(cmdOpts->statefile,"r");
	if (fp == NULL)
		return perr("Error opening file");
	multi_modelstruct *model = load_multi_modelstruct(fp);

	fclose(fp);

	multi_emulator *the_multi_emulator = alloc_multi_emulator(model);

	int number_params = model->nparams;
	gsl_vector * the_point = gsl_vector_alloc(number_params);

	int number_outputs = model->nt;
	gsl_vector *the_mean = gsl_vector_alloc(number_outputs);
	gsl_vector *the_variance = gsl_vector_alloc(number_outputs);



#ifdef BINARY_INTERACTIVE_MODE
	r = expected_r = sizeof(double);
#else
	r = expected_r = 1;
#endif

	if(!cmdOpts->quietFlag){
		fprintf(interactive_output,"%d\n",number_params);
		for(i = 0; i < number_params; i++) {
			/* FIXME we may want parameter identifiers in the future. */
			fprintf(interactive_output,"%s%d\n","param_",i);
		}

	
		int nreturns = 2 * number_outputs ; /* FIXME - also return the joint implausibility. */
		fprintf(interactive_output,"%d\n",nreturns);
		for(i = 0; i < number_outputs; i++) {
			/* FIXME we may want output identifiers in the future. */
			fprintf(interactive_output,"%s_%d\n%s_%d\n","mean",i,"variance",i);
		}

		fflush(interactive_output);
	}
	
	while (! feof(interactive_input)) {
		for(i =0; (i < number_params) && (r == expected_r); i++) {
			#ifdef BINARY_INTERACTIVE_MODE
				r = fread(gsl_vector_ptr(the_point, i), sizeof(double), 1, interactive_input);
			#else
				r = fscanf(interactive_input, "%lf%*c", gsl_vector_ptr(the_point, i));
			#endif
		}
		if (r < expected_r) /* probably eof, otherwise error */
			break;
		if(cmdOpts->pcaOutputFlag == 0){ 
			emulate_point_multi(the_multi_emulator, the_point, the_mean, the_variance);
		} else { /* support output in the pca space */
			emulate_point_multi_pca(the_multi_emulator, the_point, the_mean, the_variance);
		}
		for(i = 0; i < number_outputs; i++) {
			#ifdef BINARY_INTERACTIVE_MODE
				fwrite(gsl_vector_ptr(the_mean,i), sizeof(double), 1, interactive_output);
				fwrite(gsl_vector_ptr(the_variance,i), sizeof(double), 1, interactive_output);
			#else
				fprintf(interactive_output, "%.17f\n", gsl_vector_get(the_mean,i));
				fprintf(interactive_output, "%.17f\n", gsl_vector_get(the_variance,i));
			#endif
		}
		fflush(interactive_output);
	}

	free_multi_emulator(the_multi_emulator);
	// this is causing segfaults, i guess i dont understand the allocation pattern here properly
	//free_multimodelstruct(model);
	gsl_vector_free(the_point);
	gsl_vector_free(the_mean);
	gsl_vector_free(the_variance);
	return 0;
}





/*********************************************************************
main
*********************************************************************/
int main (int argc, char ** argv) {

	if(argc < 3) { // the minimum possible number of args is 3?
		return perr(useage);
	}
		
	struct cmdLineOpts *opts = global_opt_parse(argc, argv);

	/* just to check that the options work */
	/* printf("regOrder: %d\n", opts->regOrder); */
	/* printf("covfn: %d\n", opts->covFn); */
	/* printf("quiet: %d\n", opts->quietFlag); */
	/* printf("run_mode: %s\n", opts->run_mode); */
	/* printf("pcaOutputFlag: %d\n", opts->pcaOutputFlag); */
	/* printf("pca_variance: %lf\n", opts->pca_variance); */
	/* if(opts->inputfile != NULL){ */
	/* 	printf("inputfile: %s\n", opts->inputfile); */
	/* } */
	/* printf("statefile: %s\n", opts->statefile); */

	if (str_equal(opts->run_mode, "estimate_thetas")){
		estimate_thetas(opts);
	} else if (str_equal(opts->run_mode, "interactive_mode")) {
		interactive_mode(opts);
	} else{
		free(opts);
		return perr(useage);
	}
	
	free(opts);
	return(EXIT_SUCCESS);
}


/**
 * option parsing using getoptlong
 */
struct cmdLineOpts* global_opt_parse(int argc, char** argv)
{

	/* note: flags followed with a colon come with an argument */
	static const char *optString = "r:c:qh?"; 

	// should add a variance option for the pca_decomp
	// and a flag to do return output in pca space
	static const struct option longOpts[] = {
		{ "regression_order", required_argument , NULL, 'r'},
		{ "covariance_fn", required_argument , NULL, 'c'},
		{ "pca_variance", required_argument, NULL , 'v'}, // set the var fraction for the pca decomp
		{ "pca_output",  no_argument , NULL , 'z'},  // output from interactive emulator is left in pca space 
		{ "quiet", no_argument , NULL, 'q'},
		{ "help", no_argument , NULL, 'h'},
		{ NULL, no_argument, NULL, 0} 
	};

	struct cmdLineOpts *opts = (struct cmdLineOpts*) malloc(sizeof(struct cmdLineOpts));
	// init with default values
	opts->regOrder = 0;
	opts->covFn = 0;
	opts->quietFlag = 0;
	opts->pca_variance = 0.99; 
	opts->pcaOutputFlag = 0; // if this is one, output multivar results without rotation
	opts->run_mode = NULL;
	opts->inputfile = NULL;
	opts->statefile = NULL;

	int longIndex;
	int opt;

	// add more options here 
	opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    while( opt != -1 ) {
        switch( opt ) {
				case 'r':
					opts->regOrder = atoi(optarg); 
					break;
				case 'c':
					opts->covFn = atoi(optarg); /* this expects the cov fn to be 0,1,2? */
					break;
				case 'v':
					opts->pca_variance = atof(optarg); /* expect the var to be a float < 1 > 0 */
					if(opts->pca_variance < 0.0 || opts->pca_variance > 1.0){
						fprintf(stderr, "# err pca_variance argument given incorrect value: %lf\n", opts->pca_variance);
						opts->pca_variance = 0.95;
						fprintf(stderr, "# using default value: %lf\n", opts->pca_variance);
					}
				case 'z':
					opts->pcaOutputFlag = 1;
				case 'q':
					opts->quietFlag = 1;
					break;
				case 'h':   /* fall-through is intentional */
				case '?':
					//display_usage();
					exit(perr(useage));
					break;
				default:
					/* You won't actually get here. */
					break;
				}
				opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    }

		// set the remaining fields
		opts->run_mode = argv[optind];
		if(str_equal(opts->run_mode, "estimate_thetas")){
			opts->inputfile = argv[optind+1];
			opts->statefile = argv[optind+2];
		} else {
			opts->statefile = argv[optind+1];
		}

		return opts;
}
