#ifndef __INC_MODELSTRUCT__
#define __INC_MODELSTRUCT__

#include "optstruct.h"
#include "assert.h"
#include "stdlib.h"
#include "useful.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/** 
 * @file modelstruct.h
 * \brief defines modelstruct which contains the raw data to be worked on by emulator and estimator
 */

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
 * a pointer to the options used to init this model (should *not* be free'd)
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
	 * a pointer back to the options struct used to init this structure, just for fun really
	 */
	optstruct*  options;
} modelstruct;


void alloc_modelstruct(modelstruct* the_model, optstruct* options);
void free_modelstruct(modelstruct* the_model);
void fill_modelstruct(modelstruct* the_model, optstruct* options, char** input_data, int number_lines);	
void copy_modelstruct(modelstruct* dst, modelstruct* src);

#endif
