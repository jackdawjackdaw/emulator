CC = gcc
LDFLAGS = -lgsl -lgslcblas -lm
CFLAGS = -O0 -Wall -g -pedantic -Wextra -std=gnu99 -I../libEmu -I..
#CFLAGS = -O3 -Wall
TAGS = etags #for emacs
LTOOL = libtool
LINK = --mode=link gcc
LIBEMU = ../libEmu/libEmu.a


all: bfgstest

libEmulator:
	$(MAKE) -C ../libEmu

gradtest: 
	$(LTOOL) $(LINK) grad-testing.c   $(LIBEMU) -o grad-test.out $(LDFLAGS) $(CFLAGS)

bfgstest:
	$(LTOOL) $(LINK) test-bfgs.c  ../useful.o ../persist.o ../ioread.o $(LIBEMU) -o test-bfgs.out $(LDFLAGS) $(CFLAGS)

## this is a gsl test mainly
testing:  
	$(LTOOL) $(LINK) testing.c $(LIBEMU) -o testingGSL.out  $(LDFLAGS) $(CFLAGS)

estimatortest:  
	$(LTOOL) $(LINK) estimator-test.c -o estimator-test.out $(LIBEMU) $(LDFLAGS) $(CFLAGS)

test-resize:
	$(CC) test-resize.c -o test-resize.out $(CFLAGS)

array-test:
	$(CC) array-test.c -o array-test.out $(CFLAGS) $(LDFLAGS)

list-test:
	$(CC) list-test.c -o list-test.out $(CFLAGS)

testRbind:
	$(LTOOL) $(LINK) testRbind.c -o testRbind.out ../libRbind.so ../ioread.o ../useful.o $(LIBEMU) $(LDFLAGS) $(CFLAGS)