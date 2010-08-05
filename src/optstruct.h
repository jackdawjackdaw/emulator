#ifndef __INC_OPTSTRUCT__
#define __INC_OPTSTRUCT__

#include <gsl/gsl_matrix.h>
#include <string.h>

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
	int nregression_fns;
	double emulate_min;
	double emulate_max;
	char  filename[128];
	char outputfile[128];
	// this holds the ranges for the optimisation routine
	gsl_matrix* grad_ranges;
	// the covariance function
	double cov_fn_alpha;
	double (*covariance_fn)(gsl_vector*, gsl_vector*, gsl_vector*, int, int, double);
} optstruct;

void free_optstruct(optstruct *opts);
void copy_optstruct(optstruct *dst, optstruct* src);


#endif
