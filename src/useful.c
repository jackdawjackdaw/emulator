// some useful things
#include "useful.h"

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
	the_array = malloc(sizeof(int)*the_total_size);
	if(the_array == NULL){
		fprintf(stderr, "cannot allocate array\n");
		exit(1);
	}
	
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
	if((local_array = malloc(sizeof(double)*the_total_size)) == NULL){
		fprintf(stderr, "cannot allocate array\n");
		exit(1);
	}
	
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
		printf("Got seed %u from gettimeofday()\n", seed);
	}
	else {
		fread(&seed, sizeof(seed), 1, devrandom);
		printf("Got seed %u from /dev/random\n", seed);
		fclose(devrandom);
	}
	return(seed);
}


/* binary io routines */
// input
void in_int(FILE *fptr, int *iptr){
	if (fread((void*) iptr, sizeof(int), 1, fptr) != 1) {
		fprintf(stderr, "in_int: fread failed!\n"); 
		exit(1);
	}
}


void in_double(FILE *fptr, double *dptr){
	if (fread((void*) dptr, sizeof(double), 1, fptr) != 1){
		fprintf(stderr, "in_double: fread failed!\n"); 
		exit(1);
	}
}


void in_vector(FILE *fptr, double *vec){ // reads in a 3d vector
	if (fread((void*) vec, sizeof(double), NDIM, fptr) != NDIM){
		fprintf(stderr, "in_vec: fread failed!\n"); 
		exit(1);
	}
}

void in_blob(FILE *fptr, int n_things, double *vec){ // reads in n_things doubles
	if (fread((void*) vec, sizeof(double), n_things, fptr) != n_things){
		fprintf(stderr, "in_blob: fread_failed!\n");
		exit(1);
	}
}


// output
void out_int(FILE* fptr, int ival){
	if(fwrite((void*) &ival, sizeof(int), 1, fptr) != 1){
		fprintf(stderr, "out_int: fwrite failed!\n"); 
		exit(1);
	}
}

void out_double(FILE* fptr, double dval){
	if(fwrite((void*) &dval, sizeof(double), 1, fptr) != 1){
		fprintf(stderr, "out_double: fwrite failed!\n"); 
		exit(1);
	}
}
	
void out_vector(FILE* fptr, double* vec){
	if(fwrite((void*) vec, sizeof(double), NDIM, fptr) != NDIM){
		fprintf(stderr, "out_vec: fwrite failed!\n"); 
		exit(1);
	}
}

