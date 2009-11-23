#ifndef __INC_MAIN__
#define __INC_MAIN__

#include "stdio.h"
#include "stdlib.h"
#include "unistd.h"
#include "libEmu/estimator.h"
#include "libEmu/emulator.h"
#include "libEmu/maximise.h"
#include "libEmu/maxbfgs.h"
#include "ioread.h"
#include "sys/time.h"
#include "useful.h"
#include "estimate_threaded.h"


#define NELDER

/**
 * @file
 * @author Chris Coleman-Smith cec24@phy.duke.edu
 * @version 0.2
 * @section Description
 * 
 * The main apparatus to run the emulator, reads in options from the command line and
 * reads the actual data points from the given filename
 * 
 * threaded the estimator call function
 */
 

/* 
 * 1 -> read command line parameters
 * 2 -> read and process stream data
 * 3 -> estimate thetas
 * 4 -> calculate new means and variances
 * 5 -> output
 */


/* 
 * in the 1d case
 *  now there are only 3 hyperparams by default 
 *  -> vertical-scale theta0
 *  -> nugget theta1
 *  -> length-scale theta2...theta(Nparams-2⎈)
 */
#define NTHETASDEFAULT 4
#define NPARAMSDEFAULT 1
#define NEMULATEDEFAULT 4096
#define EMULATEMINDEFAULT 0.0
#define EMULATEMAXDEFAULT 1.0

//! holds command line options
/** 
 * designed to hold basic command line 
 * options
 */
typedef struct optstruct{
	int nthetas;
	int nparams;
	int nmodel_points;
	int nemulate_points;
	double emulate_min;
	double emulate_max;
	char  filename[128];
	char outputfile[128];
} optstruct;



void print_usage(void);
void print_vector_quiet(gsl_vector* vec, int npts);
void parse_arguments(int argc, char** argv, optstruct* options);
void emulate_model(gsl_matrix* xmodel, gsl_vector* training, gsl_vector*thetas, optstruct* options);
void estimate_thetas(gsl_matrix* xmodel_input, gsl_vector* training_vector, gsl_vector* thetas, optstruct* options);
void read_input_bounded(gsl_matrix* model, gsl_vector* training, optstruct * options);
void read_input_fromfile(gsl_matrix *xmodel, gsl_vector *training, optstruct *options);


#endif
