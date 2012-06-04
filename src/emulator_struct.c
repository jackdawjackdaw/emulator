#include "emulator_struct.h"
#include "libEmu/emulator.h"
#include "libEmu/regression.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

/*********************************************************************
Allocates emulator_struct.
*********************************************************************/
emulator_struct * alloc_emulator_struct(modelstruct * model) {
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

	makeCovMatrix(c_matrix, model->xmodel, model->thetas, e->nmodel_points,
		e->nthetas, e->nparams);
	gsl_matrix_memcpy(temp_matrix, c_matrix);
	chol_inverse_cov_matrix(model->options, temp_matrix, e->cinverse, &determinant_c);
	makeHMatrix(e->h_matrix, model->xmodel, e->nmodel_points, e->nparams, e->nregression_fns);
	estimateBeta(e->beta_vector, e->h_matrix, e->cinverse, model->training_vector,
							 e->nmodel_points, e->nregression_fns);
	gsl_matrix_free(c_matrix);
	gsl_matrix_free(temp_matrix);
	return e;
}


/*********************************************************************
Free emulator_struct.
*********************************************************************/
void free_emulator_struct(emulator_struct * e) {
	gsl_matrix_free(e->cinverse);
	gsl_matrix_free(e->h_matrix);
	gsl_vector_free(e->beta_vector);
	free((void *)e);
}



