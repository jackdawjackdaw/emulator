#ifndef __INC_PERSIST__
#define __INC_PERSIST__

#include "stdio.h"
#include "stdlib.h"
#include "libEmu/emulator.h"
#include "libEmu/estimator.h"
#include "libEmu/maximise.h"
#include "multi/multifit.h"
#include "useful.h"

void dump_emuresult(emuResult* eres, FILE* fptr);
void load_emuresult(emuResult* eres, FILE* fptr);
void dump_eopts(eopts* the_struct, FILE* fptr);
void load_eopts(eopts* the_struct, FILE *fptr);

#endif
