#ifndef __INC_MAIN__
#define __INC_MAIN__

#include "stdio.h"
#include "stdlib.h"
#include "getopt.h"
#include "libEmu/estimator.h"
#include "libEmu/emulator.h"
#include "libEmu/emulate-fns.h"
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
 * @file main.h
 * @author Chris Coleman-Smith cec24@phy.duke.edu
 * @version 0.2
 * @section Description
 *
 * \brief The emulator project is intended to run a gaussian process interpolation 
 * model on a data set. 
 * 
 * This is the main header for the two binaries, emulator and estimator. They share a lot of
 * common boiler plate code. 
 * 
 * Estimator should be run first, this reads data from stdin, and 
 * estimates an optimial set of hyperparameters which then specify entirely the gaussian-process 
 * model (which may be of variable fidelity) of the supplied data. Estimator produces a file
 * thetas.txt containing the optimum values of the hyperparams
 *
 * Emulator reads the thetas.txt file and again reads the model data from stdin, the data
 * is then interpolated/emulated and values of the gaussian-process model of the supplied
 * data are produced in a regular grid over the parameter space. These are written to 
 * emulator-out.
 * 
 * 
 */
 

void write_thetas(char* theta_file, gsl_vector* thetas, optstruct *options);

#endif
