// some useful things
#include "useful.h"


//! takes a print statement and either prints it out, or logs it 
/**
 * \note should finish implementing this, or find a nice logging lib to do it
 */
void message(char* the_message, int level){
	// for now it just prints
	char buffer[128];
	sprintf(buffer, "LOG: %s\n", the_message);
	fprintf(stderr, "%s\n", the_message);
}



//! checks for null
void *MallocChecked(size_t size){
	void *r = malloc(size);

	if( r == NULL)
		unix_error("memory wasn't allocated");
	
	return(r);
}

// unix error fn
void unix_error(char *msg){
	fprintf(stderr, "%s: %s\n", msg, strerror(errno));
	exit(0);
}


void posix_error(int code, char* msg){
	fprintf(stderr, "%s: %s\n", msg, strerror(code));
}

void error(char *msg){
	fprintf(stderr, "%s\n", msg);
	exit(0);
}

//! print a vector to stderr
void print_vector_quiet(gsl_vector *x, int n){
int i;
	for(i =0; i < n; i++){
		fprintf(stderr, "%g\t", gsl_vector_get(x, i));
	}
	fprintf(stderr,"\n");
}


// memory management

/** allocates a 2d array of ints
 * @arg size_x the x dimension of the array
 * @arg size_y the y dimension of the array
 * @arg init_val an initial value to set all entries in the array to
 * @return a pointer to the allocated array, don't forget to free this properly later!
 */
int* alloc_2d_int_array(int size_x, int size_y, int init_val){
	int i;
	int *the_array;
	int the_total_size = size_x*size_y;
	the_array = MallocChecked(sizeof(int)*the_total_size);
	
	for (i = 0; i < the_total_size; i++){
		the_array[i] = init_val;
	}
	
	return(the_array);
}

/** allocates a 2d array of floats
 * @arg size_x the x dimension of the array
 * @arg size_y the y dimension of the array
 * @arg init_val an initial value to set all entries in the array to
 * @return a pointer to the allocated array, don't forget to free this properly later!
 */
double* alloc_2d_double_array( int size_x, int size_y, double init_val){
	double *local_array;
	int i;
	int the_total_size = size_x*size_y;
	local_array = MallocChecked(sizeof(double)*the_total_size);
	
	for(i = 0; i < the_total_size; i++){
		local_array[i] = init_val;
	}
	return(local_array);
}


void print_int_array(int* the_array, int n_x, int n_y){
	int i, j;
	for(i = 0; i < n_x; i++){
		for(j = 0; j < n_y; j++){
			printf("%d ", the_array[i+j*n_x]);
		}
		printf("\n");
	}
}

void print_double_array(double* the_array, int n_x, int n_y){
	int i, j;
	for(i = 0; i < n_x; i++){
		for(j = 0; j < n_y; j++){
			printf("%g ", the_array[i+j*n_x]);
		}
		printf("\n");
	}
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
		fprintf(stderr,"Got seed %u from gettimeofday()\n", seed);
	}
	else {
		fread(&seed, sizeof(seed), 1, devrandom);
		fprintf(stderr, "Got seed %u from /dev/random\n", seed);
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
		fprintf(stderr,"Got seed %u from gettimeofday()\n", seed);
	}
	else {
		fread(&seed, sizeof(seed), 1, devrandom);
		fprintf(stderr, "Got seed %u from /dev/random\n", seed);
		fclose(devrandom);
	}
	return(seed);
}


/* binary io routines */
// input
void in_int(FILE *fptr, int *iptr){
	if (fread((void*) iptr, sizeof(int), 1, fptr) != 1) {		
		error( "in_int: fread failed!"); 
	}
}


void in_double(FILE *fptr, double *dptr){
	if (fread((void*) dptr, sizeof(double), 1, fptr) != 1){
		error("in_double: fread failed!\n");
	}
}


void in_vector(FILE *fptr, double *vec){ // reads in a 3d vector
	if (fread((void*) vec, sizeof(double), NDIM, fptr) != NDIM){
		error( "in_vec: fread failed!\n"); 
	}
}

void in_blob(FILE *fptr, int n_things, double *vec){ // reads in n_things doubles
	if (fread((void*) vec, sizeof(double), n_things, fptr) != (unsigned)n_things){
		error( "in_blob: fread_failed!\n");
	}
}


// output
void out_int(FILE* fptr, int ival){
	if(fwrite((void*) &ival, sizeof(int), 1, fptr) != 1){
		error( "out_int: fwrite failed!\n"); 
	}
}

void out_double(FILE* fptr, double dval){
	if(fwrite((void*) &dval, sizeof(double), 1, fptr) != 1){
		error( "out_double: fwrite failed!\n"); 
	}
}
	
void out_vector(FILE* fptr, double* vec){
	if(fwrite((void*) vec, sizeof(double), NDIM, fptr) != NDIM){
		error( "out_vec: fwrite failed!\n"); 
	}
}

