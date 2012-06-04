#ifndef __INC_MODELSTRUCT__
#define __INC_MODELSTRUCT__



#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/** 
 * @file modelstruct.h
 * \brief defines modelstruct which contains the raw data to be worked on by emulator and estimator
 */
struct optstruct;
//! holds the main data structures for a single model emulator run

/**
 * \struct modelstruct
 * this structure is everything you need to do run the estimation and the 
 * emulation
 *
 * it contains:
 * training vector
 * model points
 * theta vector (hyperparameters)
 * the options used to define this model 
 * \note: the options used to be kept separately, they are now always contained in the modelstruct
 */
typedef struct modelstruct{
	/** 
	 * \brief x_values (nmodel_points * nparams)
	 * 
	 * the points at which the model (that which is to be emulated) was evaluated, 
	 * each row is a new model point, the columns are the dimensions 
	 */
	gsl_matrix* xmodel;
	/** 
	 * \brief y_values (nmodel_points)
	 * 
	 * the values of the model evaluated at each model point
	 * this is what we aim to interpolate
	 */
	gsl_vector* training_vector;
	/** 
	 * \brief hyperparameters (nthetas)
	 * 
	 * working and eventually best values of the covariance fns (hyper)parameters which 
	 * should force the emulated model to pass through all of the training points 
	 */
	gsl_vector* thetas;
	/**
   * the average distances between samples in the design in each dimension
	 * this sets soft minima on the correlation function length
	 */
	gsl_vector* sample_scales; 
	/**
	 * the options struct used to init this structure
	 */
	struct optstruct* options;
	
	/**
	 * fnptrs which define various aspects of the regression etc
	 * should these be in optstruct?
	 */

	/**
	 * setup the regression H vector, depends upon the regression order 
	 * 
	 * (needed in regression.h)
	 */
	void (*makeHVector)(gsl_vector *h_vector, gsl_vector *x_location, int nparams); 

	/**
	 * points to the covariance function, 
	 * 
	 * (needed in emulator.h)
	 */
	double (*covariance_fn)(gsl_vector*, gsl_vector*, gsl_vector*, int, int);

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
	 * (needed in maxmultimin.h)
	 */
	void (*makeGradMatLength)(gsl_matrix *dCdTheta,  gsl_matrix* xmodel, 
														double thetaLength, int index, int nmodel_points, int nparams);



} modelstruct;


#include "optstruct.h"

void alloc_modelstruct(modelstruct* the_model, optstruct* options);
void free_modelstruct(modelstruct* the_model);
void fill_modelstruct(modelstruct* the_model, optstruct* options, char** input_data, int number_lines);	

void load_modelstruct(FILE* fptr, modelstruct* the_model, optstruct* opts);
void dump_modelstruct(FILE* fptr, modelstruct* the_model, optstruct* opts);

void copy_modelstruct(modelstruct* dst, modelstruct* src);

modelstruct * alloc_modelstruct_2(gsl_matrix* xmodel, gsl_vector* training_vector, int cov_fn_index, int regression_order);
void free_modelstruct_2(modelstruct * model);
void dump_modelstruct_2(FILE *fptr, modelstruct* the_model);
modelstruct* load_modelstruct_2(FILE *fptr);


#endif
