#ifndef __INC_MULTIVARSUPPORT__
#define __INC_MULTIVARSUPPORT__

#include "multi_modelstruct.h"
#include "emulator_struct.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

struct multi_modelstruct;
struct emulator_struct;

typedef struct multi_emulator {
	int nt;
	int nr;
	int nparams;
	int nmodel_points;
	int nregression_fns;
	int nthetas;
	multi_modelstruct *model;
	emulator_struct **emu_struct_array;
	
} multi_emulator;


multi_emulator *alloc_multi_emulator(multi_modelstruct *model);

void free_multi_emulator(multi_emulator *e);

void estimate_multi(multi_modelstruct *m, FILE* outfp);

void emulate_point_multi(multi_emulator *emu, gsl_vector *the_point,
												 gsl_vector *the_mean, 
												 gsl_vector *the_variance);

void emulate_point_multi_covar(multi_emulator *emu, gsl_vector *the_point,
															 gsl_vector *the_mean, gsl_vector *the_variance,
															 gsl_matrix *the_covar);

void emulate_point_multi_pca(multi_emulator *emu, gsl_vector *the_point,
												 gsl_vector *the_mean, 
												 gsl_vector *the_variance);



#endif
