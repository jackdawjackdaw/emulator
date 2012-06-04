#ifndef __INC_EMULATORSTRUCT__
#define __INC_EMULATORSTRUCT__


#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/**
 * \@file emulator_struct.h
 * \brief defines the emulator_struct which contains the details of a working emulator at runtime
 */

struct modelstruct;

/*********************************************************************
This data structure holds pointers to the information used to execute
the emulator at a new point.  This is especially helpful since I want
to stream in queries and get results streamed back out.
*********************************************************************/
typedef struct emulator_struct {
	int nparams;
	int nmodel_points;
	int nregression_fns;
	int nthetas;
	struct modelstruct * model;
	gsl_matrix * cinverse;
	gsl_vector * beta_vector;
	gsl_matrix * h_matrix;
} emulator_struct;

#include "modelstruct.h"

emulator_struct * alloc_emulator_struct(modelstruct * model);
void free_emulator_struct(emulator_struct * e);

void makeHMatrix_es(gsl_matrix *h_matrix, emulator_struct *e);
void makeCovMatrix_es(gsl_matrix *cov_matrix, emulator_struct *e);
void makeKVector_es(gsl_vector *kvector, gsl_vector *point, emulator_struct *e);
void estimateBeta_es(gsl_vector *beta_vector, emulator_struct *e);
void emulate_point(emulator_struct* e, gsl_vector * point,  double * mean, double * variance);
#endif
