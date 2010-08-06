#ifndef _INC_RESULTSTRUCT_
#define _INC_RESULTSTRUCT_

#include "optstruct.h"
#include "modelstruct.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"

typedef struct resultstruct{
	/* the points that the emluator has been evaluated at*/
	gsl_matrix* new_x; 
	/* the emulated mean at these points */
	gsl_vector* emulated_mean;
	/* the emulated variance at these points */
	gsl_vector* emulated_var;
	/* the optstruct which "owns" these results */ 
	optstruct* options;
	/* the modelstruct which led to these results */
	modelstruct* model;
} resultstruct;

void free_resultstruct(resultstruct *res);
void alloc_resultstruct(resultstruct *res, optstruct *opts);
void copy_resultstruct(resultstruct *dst, resultstruct *src);



#endif
