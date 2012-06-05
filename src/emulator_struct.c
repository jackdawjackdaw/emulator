#include "emulator_struct.h"
#include "libEmu/emulator.h"
#include "libEmu/regression.h"
#include "libEmu/emulate-fns.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

/*********************************************************************
Allocates emulator_struct.
*
* ccs, should be thread safe now
*********************************************************************/
emulator_struct * alloc_emulator_struct(modelstruct * model)
{
	emulator_struct * e = (emulator_struct *)malloc(sizeof(emulator_struct));
	double determinant_c	= 0.0;
	e->nparams = model->options->nparams;
	e->nmodel_points = model->options->nmodel_points;
	e->nregression_fns = model->options->nregression_fns;
	e->nthetas = model->options->nthetas;
	e->model = model;
	gsl_matrix * c_matrix = gsl_matrix_alloc(e->nmodel_points, e->nmodel_points);
	e->cinverse = gsl_matrix_alloc(e->nmodel_points, e->nmodel_points);
	e->beta_vector = gsl_vector_alloc(e->nregression_fns);
	e->h_matrix = gsl_matrix_alloc(e->nmodel_points, e->nregression_fns);
	gsl_matrix * temp_matrix = gsl_matrix_alloc(e->nmodel_points, e->nmodel_points);

	makeCovMatrix_es(c_matrix, e);
	gsl_matrix_memcpy(temp_matrix, c_matrix);
	chol_inverse_cov_matrix(model->options, temp_matrix, e->cinverse, &determinant_c);
	makeHMatrix_es(e->h_matrix, e);
	estimateBeta_es(e->beta_vector, e);
	
	gsl_matrix_free(c_matrix);
	gsl_matrix_free(temp_matrix);
	return e;
}


/*********************************************************************
Free emulator_struct.
*********************************************************************/
void free_emulator_struct(emulator_struct * e)
{
	gsl_matrix_free(e->cinverse);
	gsl_matrix_free(e->h_matrix);
	gsl_vector_free(e->beta_vector);
	free((void *)e);
}


/**
 * the following fns are wrapped versions of commonly used fns in emulator.h etc, this 
 * should make manipulations a little cleaner without breaking libRbind
 * they also use the emulator_struct->model->fnptrs instead
 */

/**
 * make the matrix of regression fns HMatrix, (thread safe?)
 * uses modelstruct fnptr: (void)(makeHvector)(gsl_vector*, gsl_vector*, int);
 *
 * although emulator_struct contains a h_matrix, we might want this for somethign else, so it is given as the 1st arg
 * the h_matrix in e is *not set* by this function
 */
void makeHMatrix_es(gsl_matrix *h_matrix, emulator_struct *e)
{
	makeHMatrix_fnptr(h_matrix, 
							e->model->xmodel,
							e->nmodel_points,
							e->nparams,
							e->nregression_fns,
							e->model->makeHVector);
}


/**
 * make a covariance matrix: makeCovMatrix_es, (thread safe?)
 * uses the modelstruct fnptr (double)(*covariance_fn)(gsl_vector*, gsl_vector*, gsl_vector*, int, int)
 */
void makeCovMatrix_es(gsl_matrix *cov_matrix, emulator_struct *e)
{
	makeCovMatrix_fnptr(cov_matrix, 
											e->model->xmodel, 
											e->model->thetas,
											e->nmodel_points,
											e->nthetas, 
											e->nparams,
											e->model->covariance_fn);
}


/**
 * compute a kvector at point k(point, emulator_struct) := { c(point, x_1), ... , c(point, x_n) }^T
 *
 * threadsafe
 *
 * requires: e->model to be init
 * e->model->thetas needs to contain sensible values
 */
void makeKVector_es(gsl_vector *kvector, gsl_vector *point, emulator_struct *e)
{
	makeKVector_fnptr(kvector, e->model->xmodel, point,
										e->model->thetas, e->nmodel_points, e->nthetas, e->nparams, 
										e->model->covariance_fn);
}

/** 
 * estimate betas from an emulator_struct, no need to modify anything to make this threadsafe
 *
 * betas are set into the first vector, not in the emulator_struct
 *
 * @require e->cinverse and e->h_matrix to be init
 */
void estimateBeta_es(gsl_vector *beta_vector, emulator_struct *e)
{
	estimateBeta(beta_vector, e->h_matrix, e->cinverse, e->model->training_vector, 
							 e->nmodel_points, e->nregression_fns);
}

/*********************************************************************
return mean and variance at a point.
thread-safe, uses fnptrs defined in e->modelstruct
*********************************************************************/
void emulate_point(emulator_struct* e, gsl_vector * point,  double * mean, double * variance)
{
	gsl_vector * kplus = gsl_vector_alloc(e->nmodel_points); //allocate variable storage
	gsl_vector * h_vector = gsl_vector_alloc(e->nregression_fns); //allocate variable storage

	makeKVector_es(kplus, point, e);
	/* use the fnptr in modelstruct */
	e->model->makeHVector(h_vector, point, e->nparams);
	
	(*mean) = makeEmulatedMean(e->cinverse, e->model->training_vector,
		kplus, h_vector, e->h_matrix, e->beta_vector, e->nmodel_points);
	double kappa = covariance_fn(point, point, e->model->thetas, e->nthetas, e->nparams);
	(*variance) = makeEmulatedVariance(e->cinverse, kplus, h_vector,
		e->h_matrix, kappa, e->nmodel_points, e->nregression_fns);
	gsl_vector_free(kplus);
	gsl_vector_free(h_vector);

	// can we really return on void?
	return;
}

