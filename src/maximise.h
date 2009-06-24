#ifndef __INC__MAXIMISE__
#define __INC__MAXIMISE__
#include "emulator.h"
#include "estimator.h"
#include "gsl/gsl_rng.h"

void gradDesc(gsl_rng* rand, int max_tries, int nsteps, double gamma, gsl_matrix* ranges,  gsl_matrix *xmodel, gsl_vector *trainingvector, gsl_vector *thetas, int nmodel_points, int nthetas, int nparams);
void set_random_initial_value(gsl_rng* rand, gsl_vector* x, gsl_matrix* ranges,int  nthetas);
int range_check(gsl_vector* x, gsl_matrix* ranges, int nthetas);
int vector_components_equal(gsl_vector *x, double test_value, int nthetas);

typedef struct evalunit {
	int index;
	double value;
}evalunit;


int compare_evals(const evalunit *a, const evalunit *b);
void nelderMead(gsl_rng *rand, int max_tries, int nsteps, gsl_vector* the_answer,  gsl_matrix* ranges, gsl_matrix* xmodel, gsl_vector *trainingvector, gsl_vector* thetas, int nmodel_points, int nthetas, int nparams);
void make_new_vlist(gsl_matrix* new_verticies, gsl_matrix* verticies, double sigma, int nverticies, int nthetas);
double evalLikelyhood(gsl_vector *vertex, gsl_matrix *xmodel, gsl_vector *trainingvector, int nmodel_points, int nthetas, int nparams);
void sort_vertex_list(gsl_matrix *verticies, evalunit* evalList, int nverticies, int nthetas);
void calc_com(gsl_matrix *verticies, gsl_vector *com, int nverticies, int nthetas);
#endif
