#ifndef __ESTIMATE_INC_
#define __ESTIMATE__INC_
#include "emulator.h"
#include "regression.h"

double getLogLikelyhood(gsl_matrix *cinverse, double det_cmatrix,  gsl_matrix *xmodel, gsl_vector *trainingvector, gsl_vector *thetas, gsl_matrix *h_matrix, int nmodel_points, int nthetas, int nparams, int nregression_fns);


#endif

/* double getGradient(gsl_matrix *cinverse, gsl_matrix *xmodel, gsl_vector *trainingvector, gsl_vector *thetas, int index, int nmodel_points, int nthetas, int nparams); */
/* void  getdCdt(gsl_matrix* dcdt, gsl_matrix* xmodel, gsl_vector* thetas, int index, int nmodel_points, int nthetas, int nparams); */
