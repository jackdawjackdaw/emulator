
#include "stdio.h"
#include "stdlib.h"
#include "emulator.h"
#include "estimator.h"
#include "maximise.h"

#include "multifit.h"

int main (void){
	char inputfile[128];
	eopts the_options;
	FILE *fptr;
	sprintf(inputfile, "model-cut.dat");

	fptr = fopen(inputfile, "r");
	
	// we kind of have to hope that the inputfile is sorted, could sort it...
	read_input_fromfile(inputfile, the_options);


	return(0);
}

void read_input_from_file(char* filename, eopts* options){
	int i = 0; 
	int j = 0;
	double temp_value = 0.0;
	FILE *fptr;
	fptr = fopen(options->filename, "r");

	for(i =0; i < options->nmodel_points; i++){
		for(j = 0; j < options->nparams; j++){
			fscanf(fptr, "%lg", &temp_value);
			gsl_matrix_set(xmodel, i, j, temp_value);
		}
		fscanf(fptr, "%lg", &temp_value);
		gsl_vector_set(training, i, temp_value);
	}
	print_matrix(xmodel, options->nmodel_points, options->nparams);
	vector_print(training, options->nmodel_points);
}
