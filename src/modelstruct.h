#ifndef __INC_MODELSTRUCT__
#define __INC_MODELSTRUCT__

#include "optstruct.h"
#include "assert.h"
#include "stdlib.h"
#include "useful.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

//! holds the main data structures for a single model emulator run

/**
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
	gsl_matrix* xmodel;
	gsl_vector* training_vector;
	gsl_vector* thetas;
	optstruct*  options;
} modelstruct;


void alloc_modelstruct(modelstruct* the_model, optstruct* options);
void free_modelstruct(modelstruct* the_model);
void fill_modelstruct(modelstruct* the_model, optstruct* options, char** input_data, int number_lines);	
void copy_modelstruct(modelstruct* dst, modelstruct* src);

#endif
