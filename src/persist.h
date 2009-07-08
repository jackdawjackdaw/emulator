#ifndef __INC_PERSIST__
#define __INC_PERSIST__

#include "stdio.h"
#include "stdlib.h"
#include "emulator.h"
#include "estimator.h"
#include "maximise.h"
#include "multifit.h"
#include "useful.h"

void dump_emuresult(emuResult* eres, FILE* fptr);
void load_emuresult(emuResult* eres, FILE* fptr);
void dump_eopts(eopts* the_struct, FILE* fptr);
void load_eopts(eopts* the_struct, FILE *fptr);

#endif
