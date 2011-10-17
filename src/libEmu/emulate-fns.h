#ifndef _INC_EMULATE_FNS_
#define _INC_EMULATE_FNS_

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"


#include "../optstruct.h"
#include "../modelstruct.h"
#include "../resultstruct.h"


void emulate_model_results(modelstruct *the_model, optstruct* options, resultstruct* results);
void emulateAtPoint(modelstruct *the_model, gsl_vector* the_point, optstruct* options, double* the_mean, double* the_variance);
void emulateAtPointList(modelstruct *the_model, gsl_matrix* point_list, optstruct* options,
												double* the_mean, double* the_variance);


/** 
 * to be used with callEmulateMC from rbind.c
 */
void emulateQuick(modelstruct *the_model, gsl_vector* the_point, optstruct* options, 
									 double* mean_out, double* var_out, gsl_matrix *h_matrix, 
									gsl_matrix *cinverse, gsl_vector* beta_vector );



void emulate_ith_location(modelstruct *the_model, optstruct *options, resultstruct *results,int i, gsl_matrix* h_matrix, gsl_matrix* cinverse, gsl_vector *beta_vector);
void chol_inverse_cov_matrix(optstruct* options, gsl_matrix* temp_matrix, gsl_matrix* result_matrix, double* final_determinant_c);
#endif
