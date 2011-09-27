#ifndef _INC_RBIND_
#define _INC_RBIND_

#include "stdio.h"

#include "libEmu/emulator.h"
#include "libEmu/estimator.h"
#include "libEmu/emulate-fns.h"
#include "libEmu/regression.h"
#include "libEmu/estimate_threaded.h"

#include "useful.h"

#include "optstruct.h"
#include "modelstruct.h"
#include "resultstruct.h"

#include "../defaults.h"

/** 
 * this file contains #defines set by cmake
 * currently we only export VERSION_NUMBER but there
 * could be other nice things we can set here, like if we want to 
 * include the R headers or not...
 */
#include "../buildConfig.h"


/**
 * \author C.Coleman-Smith cec24@phy.duke.edu
 * \file rbind.h
 * \brief defines an interface between libEmu and R
 * 
 * 
 */

/**
 * we'll alloc an instance of this to hold the various data structures
 * needed by setupEmulateMC, callEmulateMC 
 * allowing many rapid samples to be generated (hopefully)
 */
struct emulateMCData{
	optstruct* options;
	modelstruct* model;
	gsl_matrix* cov_matrix;
	gsl_matrix* cov_matrix_inverse;
	gsl_vector* beta_vector;
	gsl_matrix* h_matrix;
};
	
struct emulateMCData emuMCData;

struct emulateMCData *emuMCDataMulti;


/**
 * functions with names: call<whatever> are to be called from the external process (R etc)
 * 
 * these are: callEstimate, callEmulateAtList, callEmulateAtPt, callEvalLhoodList
 */

void callEstimate(double* xmodel_in, int* nparams_in, double* training_in, int *nmodelpts, int *nthetas_in, double* final_thetas, 
									int* use_fixed_nugget, double* fixed_nugget_in,
									int* cov_fn_index_in, int* regression_order_in);

void callEmulateAtList(double *xmodel_in, int *nparams_in, double* points_in, int *nemupoints, double* training_in,
											 int *nmodelpts, double* thetas_in, int *nthetas_in, double* final_emulated_y, 
											 double* final_emulated_variance, int* cov_fn_index_in, int* regression_order_in);

void callEmulateAtPt(double* xmodel_in, int* nparams_in, double* point_in, double* training_in,
										 int* nmodelpts, double* thetas_in, int* nthetas_in, double* final_emulated_y,
										 double* final_emulated_variance, int* cov_fn_index_in, int*regression_order_in);

void callEvalLhoodList(double *xmodel_in, int *nparams_in, double *pointList_in,
											 int *nevalPoints_in, double *training_in, int *nmodelPoints_in,
											 int *nthetas_in, double *answer, int*cov_fn_index_in, int* regression_order_in);


void callInfo(void);


/** 
 * functions for quick calls for getting samples of the posterior density
 * you need to call setupEmulateMC first
 */
void setupEmulateMC(double* xmodel_in, int* nparams_in, double* training_in, 
										 int* nmodelpts, double* thetas_in, int* nthetas_in, 
										int *cov_fn_index_in, int*regression_order_in);



void callEmulateMC(double* point_in, double* mean_out, double* var_out);
void freeEmulateMC(void);

void setupEmulateMCHelper(struct emulateMCData* emuMCData, double* xmodel_in, 
													int* nparams_in,  double* training_in, 
										 int* nmodelpts, double* thetas_in, int* nthetas_in, 
													int *cov_fn_index_in, int*regression_order_in);

void setupEmulateMCMulti(double* xmodel_in, int* nparams_in,  
												 double* training_in, int* nydims_in,
												 int* nmodelpts_in, double* thetas_in, int* nthetas_in, 
												 int *cov_fn_index_in, int*regression_order_in);

void callEmulateMCMulti(double* point_in, int* nydims_in, double* final_mean, double* final_var);

void freeEmulateMCMulti(int *nydims_in);

/**
 * these are INTERNAL routines, not to be called from outside 
 */



#ifdef APPLE

/* some r header information
 * currently things are working in linux so we'll only 
 * do this linking in apple builds
 */
#include "Rdefines.h" 
#include "R_ext/Rdynload.h"
/** 
 * some boilerlate for registering things with R
 * 
 * i'm not sure if adding this or switching to an "is.loaded" then load
 * setup in testMac made things work, stupid heisenbugs
 *
 * currently this is just turned on in the CMAKELISTS file
 * 
 */

R_NativePrimitiveArgType callEstArgs[10] = {REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP,
																						INTSXP, REALSXP, INTSXP, INTSXP};

R_NativePrimitiveArgType callEmuAtListArgs[12] = {REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, 
																									INTSXP, REALSXP, INTSXP, REALSXP, 
																									REALSXP, INTSXP, INTSXP};

R_NativePrimitiveArgType callEmuAtPtArgs[11] = {REALSXP, INTSXP, REALSXP, REALSXP, INTSXP,
																								REALSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP};

R_NativePrimitiveArgType callEvalListArgs[10] = {REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, INTSXP,
																						 INTSXP, REALSXP, INTSXP, INTSXP};

R_NativePrimitiveArgType setupEmulateMCArgs[8] = {REALSXP, INTSXP, REALSXP,
																									INTSXP, REALSXP, INTSXP,
																									INTSXP, INTSXP};

R_NativePrimitiveArgType callEmulateMCArgs[3] = {REALSXP, REALSXP, REALSXP};


// multidim
R_NativePrimitiveArgType setupEmulateMCMultiArgs[9] = {REALSXP, INTSXP, REALSXP,
																											 INTSXP, INTSXP, REALSXP, INTSXP,
																											 INTSXP, INTSXP};

R_NativePrimitiveArgType callEmulateMCMultiArgs[4] = {REALSXP, INTSXP, REALSXP, REALSXP};

R_NativePrimitiveArgType freeEmulateMCMultiArgs[1] = {INTSXP};



R_CMethodDef cMethods[] = {
	{"callEstimate", (DL_FUNC)&callEstimate, 10, callEstArgs},
	{"callEmulateAtList", (DL_FUNC)&callEmulateAtList, 12, callEmuAtListArgs},
	{"callEmulateAtPt", (DL_FUNC)&callEmulateAtPt, 11, callEmuAtPtArgs},
	{"callEvalLhoodList", (DL_FUNC)&callEvalLhoodList, 10, callEvalListArgs},
	{"setupEmulateMC", (DL_FUNC)&setupEmulateMC, 8, setupEmulateMCArgs},
	{"callEmulateMC", (DL_FUNC)&callEmulateMC, 3, callEmulateMCArgs},
	{"freeEmulateMC", (DL_FUNC)&freeEmulateMC, 0, NULL},
	{"setupEmulateMCMulti", (DL_FUNC)&setupEmulateMCMulti, 9, setupEmulateMCMultiArgs},
	{"callEmulateMC", (DL_FUNC)&callEmulateMCMulti, 4, callEmulateMCMultiArgs},
	{"freeEmulateMC", (DL_FUNC)&freeEmulateMCMulti, 1, freeEmulateMCMultiArgs},

	{NULL, NULL, 0}
};


void R_init_libRBIND(DllInfo *dll){
	R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
}
		
#endif


void fill_sample_scales(modelstruct* the_model, optstruct* options);

/** 
 * conversion from R arrays to C
 * we're switching from row major (R) to column major (C) 
 * with the added wrinkles that R indexs from 1 AND the 2d structures are 
 * interleaved (matrix)
 */

void convertDoubleToMatrix(gsl_matrix* the_matrix, double* input, int ny, int nx);
void convertDoubleToVector(gsl_vector* the_vec, double* input, int nx);
void testConvert(double* matrix, int *nx_in, int *ny_in);





#endif
