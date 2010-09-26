#ifndef _INC_BIN_SUPPORT_
#define _INC_BIN_SUPPORT_

#include "stdio.h"
#include "stdlib.h"
#include "getopt.h"
#include "defaults.h"
#include "optstruct.h"
#include "modelstruct.h"
#include "useful.h"

#include "libEmu/emulator.h"
#include "libEmu/regression.h"


/**
 * @file bin_support.h
 * \brief contains the fns needed for the estimator and emulator binaries, 
 * shouldn't contain anything vital for the actual statistical code 
 */ 

void print_usage_estimator(void);
void print_usage_emulator(void);

/**
 * read the args from the cmdline input, setup sensible defaults from defaults.h
 * and then try and overwrite all the ones we can. This needs some serious improvement
 */
void parse_arguments(int argc, char** argv, optstruct* options);

/** 
 * sets the covariance_fn pointer in the optstruct, this is used deep in libEmu
 */
void setup_cov_fn(optstruct *options);
/**
 * fills in the grad_ranges struct in the optstruct, vital for the estimator
 */
void setup_optimization_ranges(optstruct* options);


#endif
