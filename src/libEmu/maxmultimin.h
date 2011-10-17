#ifndef _INC_MAXMULTIMIN_
#define _INC_MAXMULTIMIN_

#include "estimate_threaded.h"
#include "pthread.h" // for debug info
#include "regression.h"
#include "emulator.h"
#include "estimator-fns.h"
#include "../useful.h"


#include "gsl/gsl_vector.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_errno.h"

/**
 * compute the gradient matrix for the length setting theta values
 * dC/dTheta = (C-nugget) * (1/2)*(x_i - x_j)^(alpha) * alpha / (thetaLength) 
 * 
 * this is a fn-ptr which will be set when the cov function is setup
 * the different target fns are in emulator.c called derivative_l_<covfnname>
 * 
 * @param covsub the covariance matrix with the nugget term subtracted
 * @param thetaLength the current value of the length scale we're differentiating wrt
 * @param index, the direction we're looking in
 * @return dCdTheta the matrix of the derivative of C wrt index
 *
 */
void (*makeGradMatLength)(gsl_matrix *dCdTheta,  gsl_matrix* xmodel, 
														double thetaLength, int index, int nmodel_points, int nparams);



void maxWithMultiMin(struct estimate_thetas_params *params);

double evalFnMulti(const gsl_vector* theta_vec, void* params_in);

void gradFnMulti(const gsl_vector* theta_vec, void* params_in, gsl_vector * grad_vec);


double getGradientCn(gsl_matrix * dCdtheta, gsl_matrix *cinverse,  
										 gsl_vector* training_vector,int nmodel_points, int nthetas);


void evalFnGradMulti(const gsl_vector* theta_vec, void* params, 
										 double* fnval, gsl_vector * grad_vec);

double estimateSigma(gsl_matrix* cinverse, void* params_in);

double estimateSigmaFull(gsl_vector *thetas, void* params_in);

void doOptimizeMultiMin( double(*fn)(const gsl_vector*, void*),													
												void(*gradientFn)(const gsl_vector*,void*, gsl_vector*),
												void(*fnGradFn)(const gsl_vector*, void*, double*, gsl_vector*),
												gsl_vector *thetaInit, gsl_vector* thetaFinal, void* args);

void set_random_init_value(gsl_rng* rand, gsl_vector* x, gsl_matrix* ranges,int  nthetas);

#endif

