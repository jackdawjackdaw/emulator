#ifndef __INC_MAIN__
#define __INC_MAIN__

#include "stdio.h"
#include "stdlib.h"
#include "unistd.h"
#include "getopt.h"
#include "libEmu/estimator.h"
#include "libEmu/emulator.h"
#include "libEmu/estimate_threaded.h"
#include "libEmu/regression.h"
#include "ioread.h"
#include "sys/time.h"
#include "useful.h"
#include "optstruct.h"
#include "modelstruct.h"
#include "defaults.h"
#include "bin_support.h"



/**
 * @file
 * @author Chris Coleman-Smith cec24@phy.duke.edu
 * @version 0.2
 * @section Description
 */
 

void write_thetas(char* theta_file, gsl_vector* thetas, optstruct *options);

#endif
