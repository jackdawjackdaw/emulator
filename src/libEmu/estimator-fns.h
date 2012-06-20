#ifndef __ESTIMATE_INC_
#define __ESTIMATE__INC_

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

double getLogLikelyhood(gsl_matrix *cinverse, double det_cmatrix,  gsl_matrix *xmodel, 
												gsl_vector *trainingvector, gsl_vector *thetas, gsl_matrix *h_matrix, int nmodel_points, 
												int nthetas, int nparams, int nregression_fns, 
												void (*makeHVector)(gsl_vector *, gsl_vector*, int));


#endif
