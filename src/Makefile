CC = gcc
LDFLAGS = -lgsl -lgslcblas -lm
CFLAGS = -O0 -Wall -g -pedantic -Wextra -std=gnu99 -I$(PWD)/libEmu
#CFLAGS = -O3 -Wall
TAGS = etags #for emacs
LTOOL = libtool
LINK = --mode=link gcc



all: main multidrive

multifit:
	$(CC) -c multifit.c $(CFLAGS) $(LDFLAGS)

multidrive: multifit  persist ioread useful libEmulator
	$(LTOOL) $(LINK) multi-drive.c -o multi-drive.out ./libEmu/libEmu.la  persist.o ioread.o useful.o  multifit.o $(LDFLAGS) $(CFLAGS)

libEmulator:
	$(MAKE) -C libEmu

test:
	$(MAKE) -C test

persist:
	$(CC) -c persist.c $(CFLAGS) $(LDFLAGS)

useful:
	$(CC) -c useful.c $(CFLAGS) $(LDFLAGS)

ioread:
	$(CC) -c ioread.c $(CFLAGS) 

tags:
	$(TAGS) emulator.c 
	$(TAGS) --append estimator.c

main:  ioread persist useful libEmulator
	$(LTOOL) $(LINK) main.c -o emulator.out  ./libEmu/libEmu.la ioread.o persist.o useful.o $(LDFLAGS) $(CFLAGS)

clean:
	rm *.o
	rm *.out