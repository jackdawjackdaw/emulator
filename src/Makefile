
CC = gcc
LDFLAGS = -lgsl -lgslcblas -lm
CFLAGS = -O0 -Wall -g
TAGS = etags #for emacs

all: estimatortest

tags:
	$(TAGS) emulator.c 
	$(TAGS) --append estimator.c

estimator: estimator.c estimator.h
	$(CC) -c estimator.c $(LDFLAGS) $(CFLAGS)

emulator: emulator.c emulator.h
	$(CC) -c emulator.c $(LDFLAGS) $(CFLAGS)

testing: emulator 
	$(CC) testing.c -o testing.out emulator.o $(LDFLAGS) $(CFLAGS)

estimatortest:  estimator emulator
	$(CC) estimator-test.c -o estimator-test.out estimator.o emulator.o $(LDFLAGS) $(CFLAGS)

clean:
	rm *.o
	rm *.out