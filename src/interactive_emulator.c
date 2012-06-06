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
  variables in GPEmulatorLib.  This is being fixed.

  "interactive_emulator estimate_thetas ..." is mostly redundant
  against "estimator", except that the input file format and
  command-line arguments are different.

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
	scan_int(input_file, & nparams);
	scan_int(input_file, & nmodel_points);
	gsl_matrix * xmodel = gsl_matrix_alloc(nmodel_points, nparams);
	gsl_vector * training_vector = gsl_vector_alloc(nmodel_points);
	for (i = 0; i < nmodel_points; i++)
		for (j = 0; j < nparams; j++)
			scan_double(input_file, gsl_matrix_ptr(xmodel,i,j));
	for (i = 0; i < nmodel_points; i++)
		scan_double(input_file, gsl_vector_ptr(training_vector, i));
	if (input_file != stdin)
		fclose(input_file);
	*xmodel_ptr = xmodel;
	*training_vector_ptr = training_vector;
	return 1;
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
	int i, r, expected_r;
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
		r = expected_r = sizeof(double);
	#else
		r = expected_r = 1;
	#endif
	fprintf(interactive_output,"%d\n",nparams);
	for(i = 0; i < nparams; i++) {
		/* FIXME we may want parameter identifiers in the future. */
		fprintf(interactive_output,"%s%d\n","param_",i);
	}
	int nreturns = 2; /* More models would mean more returns. */
	fprintf(interactive_output,"%d\n",nreturns);
	fprintf(interactive_output,"%s\n%s\n","mean","variance");
	fflush(interactive_output);
	while (! feof(interactive_input)) {
		for(i =0; (i < nparams) && (r == expected_r); i++) {
			#ifdef BINARY_INTERACTIVE_MODE
				r = fread(gsl_vector_ptr(the_point, i), sizeof(double), 1, interactive_input);
			#else
				r = fscanf(interactive_input, "%lf%*c", gsl_vector_ptr(the_point, i));
			#endif
		}
		if (r < expected_r) /* probably eof, otherwise error */
			break;
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
