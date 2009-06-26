
#ifndef __INC_MULTIFIT__
#define __INC_MULTIFIT__
#include "emulator.h"
#include "estimator.h"
#include "maximise.h"
#include "gsl/gsl_rng.h"


typedef struct eopts{
	int nmodel_points;
	int nemu_points;
	int nparams;
	int nthetas;
	double range_min;
	double range_max;
	gsl_matrix *xmodel;
	gsl_vector *training;
	gsl_vector *thetas;
	char filename[128];
} eopts;


typedef struct emuResult{
	int nemu_points;
	int nparams;
	gsl_matrix* new_x;
	gsl_vector* new_mean;
	gsl_vector* new_var;
} emuResult;


void emulate_region(gsl_matrix *new_x, gsl_vector* emulated_mean, gsl_vector* emulated_variance , eopts* options);
void estimate_region(eopts* options, gsl_rng *random);
void evaluate_region(emuResult *results, eopts* options, gsl_rng* random);
int is_smooth(double smooth_val, gsl_vector* xemu, gsl_vector* mean_emu, gsl_vector* var_emu, eopts* options);
double get_mse( double mean, double variance);

#endif
