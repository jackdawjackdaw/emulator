#ifndef _INC_MAXLBFGS_
#define _INC_MAXLBFGS_

#include "estimate_threaded.h"
#include "pthread.h" // for debug info
#include "lbfgs.h"
#include "regression.h"
#include "emulator.h"
#include "estimator.h"
#include "../useful.h"

double evalFnLBFGS(double *xinput, int nthetas, void* args);

/**
 * compute the gradient matrix for the length setting theta values
 * dC/dTheta = (C-nugget) * (1/2)*(x_i - x_j)^(alpha) * alpha / (thetaLength) 
 * 
 * this is a fn-ptr which will be set when the cov function is setup
 * the different target fns are in emulator.c called derivative_l_<covfnname>
 * 
 * 
 * @param covsub the covariance matrix with the nugget term subtracted
 * @param thetaLength the current value of the length scale we're differentiating wrt
 * @param index, the direction we're looking in
 * @return dCdTheta the matrix of the derivative of C wrt index
 */
void (*setupdCdThetaLength)(gsl_matrix *dCdTheta,  gsl_matrix* xmodel, 
														double thetaLength, int index, int nmodel_points, int nparams);

void getGradientExactGauss(double *xinput, double* gradient, int nparams, void* args);

double getGradientCn(gsl_matrix * dCdtheta, gsl_matrix *cinvese,  
										 gsl_vector* training_vector,int nmodel_points, int nthetas);

//! wrapper for evalFnLBFGS
double evalLikelyhoodLBFGS_struct(struct estimate_thetas_params *p);

//! sets the vector x to some set of random initial values within the set of ranges given
void set_random_initial_value(gsl_rng* rand, gsl_vector* x, gsl_matrix* ranges,int  nthetas);

void maxWithLBFGS(struct estimate_thetas_params *params);


#endif
