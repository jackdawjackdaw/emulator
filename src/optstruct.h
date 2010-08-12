#ifndef __INC_OPTSTRUCT__
#define __INC_OPTSTRUCT__

#include <gsl/gsl_matrix.h>
#include <string.h>

/**
 * @file optstruct.h
 * \brief defins the optstruct which holds all the mundane options and dimensions etc
 */



/**
 * \struct optstruct
 * \brief holds the main parameters etc needed by all the emulator fns
 * 
 * this is the most important structure after the modelstruct, this contains
 * all the dimensions of xmodel, training_vec etc
 * 
 */
typedef struct optstruct{
	/**
	 * the number of parameters in the covariance fn
	 */
	int nthetas;
	/**
	 * the number of parameters in the model to be emulated, this sets the dimensionality
	 * of the estimation space and the "hardness" of the estimation process
	 */
	int nparams;
	/**
	 * the number of locations in the model space that we have training values at
	 * another contributing factor to how complex the problem becomes, the covariance matrix
	 * scales as (nmodel_points**2) and we need to invert this every step of the estimation 
	 * process
	 */
	int nmodel_points;
	/** 
	 * how many points to evaluate the emulated model at
	 */
	int nemulate_points;
	/**
	 * indirectly controls the shape of the regression model, be very careful with this one
	 */
	int nregression_fns;
	/**
	 * sets the min position of the "hyper-cube" over which we evaluted the estimated mode;
	 */
	double emulate_min;
	/**
	 * sets the max position of the "hyper-cube" over which we evaluted the estimated mode;
	 */
	double emulate_max;
	char  filename[128];
	char outputfile[128];
	/** this holds the ranges for the optimisation routine*/
	gsl_matrix* grad_ranges;
	/** sets the smoothness of the gaussian covaraince fn  */
	double cov_fn_alpha;
	/**
	 * the fn ptr to the covariance function, this is the most called function in libEmu
	 * you can change this when you setup the optstruct
	 * WARNING: changing this will probably break the code as there is perhaps a final 
	 * argument to the gaussian cov fn and not the others
	 */
	double (*covariance_fn)(gsl_vector*, gsl_vector*, gsl_vector*, int, int, double);
} optstruct;

void free_optstruct(optstruct *opts);
void copy_optstruct(optstruct *dst, optstruct* src);


#endif
