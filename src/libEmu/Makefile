## makefile for libEmu
LTOOL = libtool
COMPILE = --mode=compile gcc 
CFLAGS = -O -c -O2 -Wextra -pedantic -std=c99  
#CFLAGSDEBUG = -O -c -O0 -g 
#CFLAGS = $(CFLAGSDEBUG)
LFLAGS = -lgsl -lgslcblas -lm
LINK = --mode=link gcc  -no-install -shared

all:lbfgs-test

lbfgsFort: lbfgs/routines.f
	gfortran -c lbfgs/routines.f -g -O0  -o lbfgs/routines.o

lbfgs: lbfgs.c
	gcc -c lbfgs.c -g -O0 $(LFLAGS)  -Wextra -std=c99 -DEXECUTE

lbfgs-test: lbfgs lbfgsFort
	gcc -o lbfgs-test lbfgs.o lbfgs/routines.o  -lgsl -lgslcblas -g  -lgfortran



estimator:
	$(LTOOL) $(COMPILE) $(CFLAGS) estimator.c 

emulator:
	$(LTOOL) $(COMPILE) $(CFLAGS) emulator.c 

maximise:
	$(LTOOL) $(COMPILE) $(CFLAGS) maximise.c 

maxbfgs:
	$(LTOOL) $(COMPILE) $(CFLAGS) maxbfgs.c

bfgs:
	$(LTOOL) $(COMPILE) $(CFLAGS) bfgs.c

libEmu: estimator emulator maximise maxbfgs bfgs
	$(LTOOL) $(LINK) -o libEmu.la estimator.lo emulator.lo maximise.lo maxbfgs.lo bfgs.lo $(CFLAGS) $(LFLAGS)

clean:
	rm lbfgs.o
	rm lbfgs/routines.o
	rm lbfgs-test
	# rm *.la
	# rm *.lo
