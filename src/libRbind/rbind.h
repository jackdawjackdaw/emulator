#ifndef _INC_RBIND_
#define _INC_RBIND_

#include "libEmu/emulator.h"
#include "libEmu/estimator.h"
#include "libEmu/emulate-fns.h"
#include "libEmu/regression.h"
#include "libEmu/estimate_threaded.h"

#include "useful.h"
#include "stdio.h"

#include "optstruct.h"
#include "modelstruct.h"
#include "resultstruct.h"

#include "../defaults.h"

/**
 * \file rbind.h
 * \brief defines the interface to R (messy)
 */



void callEstimate(double* xmodel_in, int* nparams_in, double* training_in, int *nmodelpts, int *nthetas_in, double* final_thetas, 
									int* use_fixed_nugget, double* fixed_nugget_in,
									int* cov_fn_index_in, int* regression_order_in);

void callEmulateAtList(double *xmodel_in, int *nparams_in, double* points_in, int *nemupoints, double* training_in,
											 int *nmodelpts, double* thetas_in, int *nthetas_in, double* final_emulated_y, 
											 double* final_emulated_variance, int* cov_fn_index_in, int* regression_order_in);

void callEmulateAtPt(double* xmodel_in, int* nparams_in, double* point_in, double* training_in,
										 int* nmodelpts, double* thetas_in, int* nthetas_in, double* final_emulated_y,
										 double* final_emulated_variance, int* cov_fn_index_in, int*regression_order_in);

void callEvalLhoodList(double *xmodel_in, int *nparams_in, double *pointList_in,
											 int *nevalPoints_in, double *training_in, int *nmodelPoints_in,
											 int *nthetas_in, double *answer, int*cov_fn_index_in, int* regression_order_in);

void fill_sample_scales(modelstruct* the_model, optstruct* options);

// conversion from R arrays to C

void convertDoubleToMatrix(gsl_matrix* the_matrix, double* input, int ny, int nx);
void convertDoubleToVector(gsl_vector* the_vec, double* input, int nx);
void testConvert(double* matrix, int *nx_in, int *ny_in);

// deprecated
/* void callEmulator(double* xmodel_in, int *nparams_in, double* training_in, int *nmodelpts, int *nthetas, double* final_x, int* nemupts, \ */
/* 									double* finaly, double* finalvar, double* rangemin, double* rangemax); */


/* void callEmulate(double* xmodel_in, int* nparams_in, double* training_in, int* nmodelpts, double* thetas_in, int* nthetas_in, double* final_emulated_x, int *nemupts_in, double* final_emulated_y, double* final_emulated_variance, double* range_min_in, double*range_max_in); */

/* void callEvalLikelyhood(double * xmodel_in, int* nparams_in, double* training_in, \ */
/* 													int *nmodelpts_in, int* nthetas_in, double* thetas_in, \ */
/* 												double* answer); */





#endif
