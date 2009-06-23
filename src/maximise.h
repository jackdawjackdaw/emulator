#ifndef __INC__MAXIMISE__
#define __INC__MAXIMISE__
#include "emulator.h"
#include "estimator.h"
#include "gsl/gsl_rng.h"

void gradDesc(gsl_rng* rand, int max_tries, int nsteps, int gamma, gsl_matrix* ranges,  gsl_matrix *xmodel, gsl_vector *trainingvector, gsl_vector *thetas, int nmodel_points, int nthetas, int nparams);
void set_random_initial_value(gsl_rng* rand, gsl_vector* x, gsl_matrix* ranges,int  nthetas);
int range_check(gsl_vector* x, gsl_matrix* ranges, int nthetas);
int vector_components_equal(gsl_vector *x, double test_value, int nthetas);

#endif
