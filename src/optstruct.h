#ifndef __INC_OPTSTRUCT__
#define __INC_OPTSTRUCT__

#include <gsl/gsl_matrix.h>

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
} optstruct;

#endif
