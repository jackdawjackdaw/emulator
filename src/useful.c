// some useful things
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "sys/time.h"
#include "useful.h"


//! checks for null
void *MallocChecked(size_t size){
	void *r = malloc(size);
	if( r == NULL)
		unix_error("memory wasn't allocated");
	return(r);
}

//! print a vector to stderr
void print_vector_quiet(gsl_vector *x, int n){
int i;
	for(i =0; i < n; i++){
		fprintf(stderr, "%g\t", gsl_vector_get(x, i));
	}
	fprintf(stderr,"\n");
}

// RNG 
// tries to read from /dev/random, or otherwise uses the system time
unsigned long int get_seed(void){
	unsigned int seed;
	struct timeval tv;
	FILE *devrandom;

	if((devrandom = fopen("/dev/random", "r")) == NULL){
		gettimeofday(&tv, 0);
		seed = tv.tv_sec + tv.tv_usec;
		//fprintf(stderr,"Got seed %u from gettimeofday()\n", seed);
	}
	else {
		fread(&seed, sizeof(seed), 1, devrandom);
		//fprintf(stderr, "Got seed %u from /dev/random\n", seed);
		fclose(devrandom);
	}
	return(seed);
}

// RNG 
// tries to read from /dev/random, or otherwise uses the system time
unsigned long int get_seed_noblock(void){
	unsigned long int seed;
	struct timeval tv;
	FILE *devrandom;

	if((devrandom = fopen("/dev/urandom", "r")) == NULL){
		gettimeofday(&tv, 0);
		seed = tv.tv_sec + tv.tv_usec;
		//fprintf(stderr,"Got seed %u from gettimeofday()\n", seed);
	}
	else {
		fread(&seed, sizeof(seed), 1, devrandom);
		//fprintf(stderr, "Got seed %u from /dev/random\n", seed);
		fclose(devrandom);
	}
	return(seed);
}

