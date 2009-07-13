## makefile for libEmu
LTOOL = libtool
COMPILE = --mode=compile gcc 
CFLAGS = -O -c -O2 -Wextra -pedantic -std=c99 
LFLAGS = -lgsl -lgslcblas -lm
LINK = --mode=link gcc 

all: libEmu

estimator:
	$(LTOOL) $(COMPILE) $(CFLAGS) estimator.c 

emulator:
	$(LTOOL) $(COMPILE) $(CFLAGS) emulator.c 

maximise:
	$(LTOOL) $(COMPILE) $(CFLAGS) maximise.c 

libEmu: estimator emulator maximise
	$(LTOOL) $(LINK) -o libEmu.la estimator.lo emulator.lo maximise.lo $(CFLAGS) $(LFLAGS)

clean:
	rm *.la
	rm *.lo
	