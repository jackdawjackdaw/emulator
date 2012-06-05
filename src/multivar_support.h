#ifndef __INC_MULTIVARSUPPORT__
#define __INC_MULTIVARSUPPORT__

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

struct multi_modelstruct;

typedef struct multi_emulator {
	int nt;
	int nr;
	int nparams;
	int nmodel_points;
	int nregression_fns;
	int nthetas;
	struct multi_modelstruct *model;
	struct emulator_struct **emu_struct_arrary;
	
} multi_emulator;


multi_emulator *alloc_multi_emulator(multi_modelstruct *model);

void free_multi_emulator(multi_emulator *e);

void estimate_multi(multi_modelstruct *m, FILE* outfp);

void emulate_point_multi(multi_emulator *emu, gsl_vector *the_point,
												 gsl_vector *the_mean, gsl_vector *the_variance);


#endif
