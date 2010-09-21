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

/**
 * \file rbind.h
 * \brief defines the interface to R (messy)
 */


void testConvert(double* matrix, int *nx_in, int *ny_in);
void convertDoubleToMatrix(gsl_matrix* the_matrix, double* input, int ny, int nx);

void convertDoubleToVector(gsl_vector* the_vec, double* input, int nx);

void callEmulateAtPt(double* xmodel_in, int* nparams_in, double* point_in, double* training_in, int* nmodelpts, double* thetas_in, int* nthetas_in, double* final_emulated_y, double* final_emulated_variance);


void callEmulator(double* xmodel_in, int *nparams_in, double* training_in, int *nmodelpts, int *nthetas, double* final_x, int* nemupts, \
									double* finaly, double* finalvar, double* rangemin, double* rangemax);

void callEstimate(double* xmodel_in, int* nparams_in, double* training_in, int *nmodelpts, int *nthetas_in, double* final_thetas);

void callEmulate(double* xmodel_in, int* nparams_in, double* training_in, int* nmodelpts, double* thetas_in, int* nthetas_in, double* final_emulated_x, int *nemupts_in, double* final_emulated_y, double* final_emulated_variance, double* range_min_in, double*range_max_in);

void callEvalLikelyhood(double * xmodel_in, int* nparams_in, double* training_in, \
													int *nmodelpts_in, int* nthetas_in, double* thetas_in, \
												double* answer);


#endif
