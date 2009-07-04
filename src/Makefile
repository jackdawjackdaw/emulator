
CC = gcc
LDFLAGS = -lgsl -lgslcblas -lm
CFLAGS = -O0 -Wall -g
#CFLAGS = -O3 -Wall
TAGS = etags #for emacs

all: multidrive


multifit:
	$(CC) -c multifit.c $(CFLAGS) $(LDFLAGS)

multidrive: multifit emulator estimator maximise persist
	$(CC) multi-drive.c -o multi-drive.out estimator.o persist.o emulator.o maximise.o multifit.o $(LDFLAGS) $(CFLAGS)

persist:
	$(CC) -c persist.c $(CFLAGS) $(LDFLAGS)


tags:
	$(TAGS) emulator.c 
	$(TAGS) --append estimator.c

maximise: maximise.c maximise.h
	$(CC) -c maximise.c $(LDFLAGS) $(CFLAGS)

gradtest: maximise estimator emulator
	$(CC) grad-testing.c maximise.o estimator.o emulator.o -o grad-test.out $(LDFLAGS) $(CFLAGS)

estimator: estimator.c estimator.h
	$(CC) -c estimator.c $(LDFLAGS) $(CFLAGS)

emulator: emulator.c emulator.h
	$(CC) -c emulator.c $(LDFLAGS) $(CFLAGS)

testing: emulator 
	$(CC) testing.c -o testing.out emulator.o $(LDFLAGS) $(CFLAGS)

estimatortest:  estimator emulator
	$(CC) estimator-test.c -o estimator-test.out estimator.o emulator.o $(LDFLAGS) $(CFLAGS)

main: emulator estimator maximise
	$(CC) main.c -o ../emulator.out estimator.o emulator.o maximise.o $(LDFLAGS) $(CFLAGS)


clean:
	rm *.o
	rm *.out