#ifndef _INC_MAXLBFGS_
#define _INC_MAXLBFGS_

#include "estimate_threaded.h"
#include "pthread.h" // for debug info
#include "lbfgs.h"
#include "../useful.h"


double evalFnLBFGS(double *xinput, int nthetas, void* args);

void getGradientExactGauss(double *xinput, double* gradient, int nparams, void* args);
double getGradientCn(gsl_matrix * dCdtheta, gsl_matrix *cinvese,  gsl_vector* training_vector,int nmodel_points, int nthetas);
void setupdCdThetaLength(gsl_matrix *dCdTheta, gsl_matrix *covsub, gsl_matrix* xmodel, double thetaLength, int index, int nmodel_points);

//! wrapper for evalFnLBFGS
double evalLikelyhoodLBFGS_struct(struct estimate_thetas_params *p);

//! sets the vector x to some set of random initial values within the set of ranges given
void set_random_initial_value(gsl_rng* rand, gsl_vector* x, gsl_matrix* ranges,int  nthetas);

void maxWithLBFGS(struct estimate_thetas_params *params);


#endif
