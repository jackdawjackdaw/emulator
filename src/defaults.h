#ifndef _INC_DEFAULTS_
#define _INC_DEFAULTS_

/**
 * @file defaults.h
 * \brief sets some defaults for command line parsing
 *
 * NTHETASDEFAULT number of adjustable values within the covariance fn which estimator is seeking to find the "best" set of
 * NPARAMSDEFAULT dimensions of the input data
 * NEMULATEDEFAULT how many points to run emulator at
 * EMULATEMINDEFAULT & EMULATEMAXDEFAULT set the size of the grid to run the emulator over
 * 
 */

#define NTHETASDEFAULT 4
#define NPARAMSDEFAULT 1
#define NEMULATEDEFAULT 1024
#define EMULATEMINDEFAULT 0.0
#define EMULATEMAXDEFAULT 4.0

#endif
