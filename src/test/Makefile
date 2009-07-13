CC = gcc
LDFLAGS = -lgsl -lgslcblas -lm
CFLAGS = -O0 -Wall -g -pedantic -Wextra -std=gnu99 -I../libEmu
#CFLAGS = -O3 -Wall
TAGS = etags #for emacs
LTOOL = libtool
LINK = --mode=link gcc
LIBEMU = ../libEmu/libEmu.la


all: gradtest testing estimatortest test-resize list-test array-test

libEmulator:
	$(MAKE) -C ../libEmu

gradtest: libEmulator
	$(LTOOL) $(LINK) grad-testing.c $(LIBEMU) -o grad-test.out $(LDFLAGS) $(CFLAGS)

## this is a gsl test mainly
testing:  libEmulator
	$(LTOOL) $(LINK) testing.c $(LIBEMU) -o testingGSL.out  $(LDFLAGS) $(CFLAGS)

estimatortest:  libEmulator
	$(LTOOL) $(LINK) estimator-test.c -o estimator-test.out $(LIBEMU) $(LDFLAGS) $(CFLAGS)

test-resize:
	$(CC) test-resize.c -o test-resize.out $(CFLAGS)

array-test:
	$(CC) array-test.c -o array-test.out $(CFLAGS) $(LDFLAGS)

list-test:
	$(CC) list-test.c -o list-test.out $(CFLAGS)