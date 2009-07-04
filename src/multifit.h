
#ifndef __INC_MULTIFIT__
#define __INC_MULTIFIT__
#include "emulator.h"
#include "estimator.h"
#include "maximise.h"
#include "gsl/gsl_rng.h"


//! options for a single run of the emulator
/**
 * The options and data for a single run of th eemulator over 
 * some range. Contains pointers to arrays and vectors which 
 * need to be allocated, and of course free'd at some point.
 * 
 * This should be enough to fully specifity a complete emulation,
 * pass this to the multifit.c methods and create an emuresult 
 */
typedef struct eopts{
	//! The number of points the model is trained at
	int nmodel_points;
	//! How many points the emulator will be evaluated at, these are by default evenly spaced and stupid in >1d
	int nemu_points;
	//! How many parameters the model takes, i.e how wide the xmodel matrix will be 
	int nparams;
	//! How mnay theta's are in the covariance function, by default this is nparams + 3
	int nthetas;
	//! Where to start running the emulator (in the first parameter)
	double range_min;
	//! Where to stop running the emulator (again in the first parameter)
	double range_max;
	//! An array of nmodel_points x nparams input values that the training vec was created at
	gsl_matrix *xmodel;
	//! A vector of nmodel_points created by the model 
	gsl_vector *training;
	//! A vec of nthetas to  be set to the most likely thetas
	gsl_vector *thetas;
	//! A filename, to store the final data. This is poorly checked!
	char filename[128];
} eopts;



//! Store the results from an emulator run.
typedef struct emuResult{
	//! the number of points the emulator was evaluated at, sets the size of the included things
	int nemu_points;
	//! The number of parameters for which the emulator was evaluated
	int nparams;
	//! A matrix of nemu_points x nparams of the points at which the emulator was evaluated
	gsl_matrix* new_x;
	//! The emulated mean at these points, vector nemu_points long
	gsl_vector* new_mean;
	//! The emulated variance at these points, vector new_var long
	gsl_vector* new_var;
} emuResult;


void emulate_region(gsl_matrix *new_x, gsl_vector* emulated_mean, gsl_vector* emulated_variance , eopts* options);
void estimate_region(eopts* options, gsl_rng *random);
void evaluate_region(emuResult *results, eopts* options, gsl_rng* random);
int is_smooth(double smooth_val, gsl_vector* xemu, gsl_vector* mean_emu, gsl_vector* var_emu, eopts* options);
double get_mse( double mean, double variance);

#endif
