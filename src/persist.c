#include "persist.h"

/**
 *  some stuff to persist the emulator to file, and load it again
 * at it should save:
 * - a complete eopts struct including the attached							
 *  + model data
 *  + training data
 *  + thetas 
 * - a copy of the covariance matrix (as is) (saves time doing it again)
 *  
 * - a complete emuResult struct including the attached
 *  + new_x
 *  + new_mean
 *  + new_var
 * 
 * These binary writes are perhaps not very portable.
 */


//! Drop the given emuresult in 
/**
 * dumps the emuresult to a file
 * @param fptr -> stream to write to (open alredy!)
 * @param eres -> the emu result we want to write
 * 
 * Dump takes place in the order the struct is specified, <nemu_points>, <nparams>, 
 * <new_x>, <new_mean>, <new_var>
 * 
 */ 
void dump_emuresult(emuResult* eres, FILE* fptr){
	int i, j;
	assert(fptr != NULL);
	// write out the basic info
	out_int(fptr, eres->nemu_points);
	out_int(fptr, eres->nparams);

	
	// now drop the new_x matrix
	for(i = 0; i < eres->nemu_points; i++){
		for(j = 0; j < eres->nparams; j++){
			out_double(fptr, gsl_matrix_get(eres->new_x, i, j));
		}
	}
	
	// drop the new_mean vector
	for(i = 0; i < eres->nemu_points; i++){
		out_double(fptr, gsl_vector_get(eres->new_mean, i));
	}
	
	// drop the new_variance
	for(i = 0; i < eres->nemu_points; i++){
		out_double(fptr, gsl_vector_get(eres->new_var, i));
	}
}

//! read an emuresult from file, allocating data as you go
/** 
 * read the fptr as an emuresult. 
 * 
 * read takes place in the order the struct is declared. 
 * @return eres is allocated to the nemu_points and nparams given in the 
 * file and takes the assigned values. 
 * 
 * This allocates NEW memory for the fields in the emuresult, withoiut 
 * trying to see they already point to something, don't pass in an already
 * allocated struct (am stupid and don't really know how to check for this 
 * without requiring all unallocated pointers to be null or something. 
 */
void load_emuresult(emuResult* eres, FILE* fptr){
	int i, j;
	int temp_int_val;
	double temp_double_val;
	assert(fptr != NULL);

	in_int(fptr, &temp_int_val);
	assert(temp_int_val > 0);
	eres->nemu_points = temp_int_val;

	in_int(fptr, &temp_int_val); 
	assert(temp_int_val > 0);
	eres->nparams = temp_int_val;

	// allocate the new_x matrix
	eres->new_x = gsl_matrix_alloc(eres->nemu_points, eres->nparams);
	// allocate the new_mean
	eres->new_mean = gsl_vector_alloc(eres->nemu_points);
	eres->new_var = gsl_vector_alloc(eres->nemu_points);

	// assign vals to the new_x matrix
	for(i = 0; i < eres->nemu_points; i++){
		for(j = 0; j < eres->nparams; j++){
			in_double(fptr, &temp_double_val);
			gsl_matrix_set(eres->new_x, i, j, temp_double_val);
		}
	}
	
	// read the new mean
	for(i = 0; i < eres->nemu_points; i++){
		in_double(fptr, &temp_double_val);
		gsl_vector_set(eres->new_mean, i, temp_double_val);
	}

	// read the new variance
	for(i = 0; i < eres->nemu_points; i++){
		in_double(fptr, &temp_double_val);
		gsl_vector_set(eres->new_var, i, temp_double_val);
	}
	
}


//! dump an eopt struct, without the filename field
/** 
 * dump the given struct to the fptr given, as a binary file. The 
 * struct is written in the order it is defined, 
 */
void dump_eopts(eopts* the_struct, FILE* fptr){
	int i, j;
	assert(fptr != NULL); // check it's at least vaguely right

	out_int(fptr, the_struct->nmodel_points);
	out_int(fptr, the_struct->nemu_points);
	out_int(fptr, the_struct->nparams);
	out_int(fptr, the_struct->nthetas);
	out_double(fptr, the_struct->range_min);
	out_double(fptr, the_struct->range_max);
	
	for(i = 0; i < the_struct->nmodel_points; i++){
		for(j = 0; j < the_struct->nparams; j++){
			out_double(fptr, gsl_matrix_get(the_struct->xmodel, i, j));
		}
	}

	for(i = 0; i < the_struct->nmodel_points; i++){
		out_double(fptr, gsl_vector_get(the_struct->training, i));
	}

	for(i = 0; i < the_struct->nthetas;i++){
		out_double(fptr, gsl_vector_get(the_struct->thetas, i));
	}
}

//! load an eopts struct, without the filename
/**
 * load an eopts struct from the given fptr, note that the filename field is not filled
 * but you could use the name of the FILE you've loaded etc. 
 * 
 * This method allocates the memory for xmodel, training and thetas and so 
 * one must free them when you're done.
 */
void load_eopts(eopts* the_struct, FILE* fptr){
	int i,j;
	int temp_int;
	double temp_double;
	assert(fptr != NULL);
	
	// read the nmodel_points
	in_int(fptr, &temp_int);
	assert(temp_int > 0);
	the_struct->nmodel_points = temp_int;
	
	// read the nemu_points
	in_int(fptr, &temp_int);
	assert(temp_int > 0);
	the_struct->nemu_points = temp_int;

	// read the nparams
	in_int(fptr, &temp_int);
	assert(temp_int > 0);
	the_struct->nparams = temp_int;

	// read the nthetas
	in_int(fptr, &temp_int);
	assert(temp_int > 0);
	the_struct->nthetas = temp_int;
	
	// read in range_min
	in_double(fptr, &temp_double);
	the_struct->range_min = temp_double;
	
	// read in range_max
	in_double(fptr, &temp_double);
	the_struct->range_max = temp_double;

	// allocate the data in the_struct
	the_struct->xmodel = gsl_matrix_alloc(the_struct->nmodel_points, the_struct->nparams);
	the_struct->training = gsl_vector_alloc(the_struct->nmodel_points);
	the_struct->thetas = gsl_vector_alloc(the_struct->nthetas);

	// not fast, could read blocks buuuut...
	for(i = 0; i < the_struct->nmodel_points; i++){
		for(j = 0; j < the_struct->nparams; j++){
			in_double(fptr, &temp_double);
			gsl_matrix_set(the_struct->xmodel, i,j, temp_double);
		}
	}
	
	for(i = 0; i < the_struct->nmodel_points; i++){
		in_double(fptr, &temp_double);
		gsl_vector_set(the_struct->training, i,temp_double);
	}

	for(i = 0; i < the_struct->nthetas; i++){
		in_double(fptr, &temp_double);
		gsl_vector_set(the_struct->thetas, i, temp_double);
	}
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
