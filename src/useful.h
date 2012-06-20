/**
 * \date 11-1-09
 * \file useful.h,
 * \brief some basic functions one always uses in scientific code
 */

#include "gsl/gsl_rng.h"
#include "gsl/gsl_vector.h"

void print_vector_quiet(gsl_vector *x, int n);
//******************
// tries to get a seed from /dev/random or otherwise uses the system time
//*****************
unsigned long int get_seed(void);
unsigned long int get_seed_noblock(void);

// errors
void *MallocChecked(size_t size);

