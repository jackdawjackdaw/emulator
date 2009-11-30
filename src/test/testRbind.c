#include "../libEmu/emulator.h"
#include "../libEmu/estimator.h"
#include "../libEmu/maximise.h"
#include "../multifit.h"
#include "../useful.h"
#include "../ioread.h"

void process_input_data(char** input_data, eopts* the_options);

extern void callEmulator(double* xmodel_in, int* nparams_in,  double* training_in, int *nmodelpts, int* nthetas_in,\
 double* final_emulated_x, int *nemupts_in, double* final_emulated_y, double* final_emulated_variance , double* range_min_in,\
 double* range_max_in );

int main (void){
	char filename[128];
	char **input_data;
	int number_lines;
	FILE *fptr;
	eopts the_options;
	int i, j;

	double* xmodel_flat;
	double* training_flat;
	double* final_emulated_x;
	double* final_emulated_y;
	double* final_emulated_variance;
	double range_min = 0.0;
	double range_max = 1.0;	
	int nparams = 2;	
	int nmodelpts;
	int nthetas = 4;
	int nemupts= 100;

	the_options.nthetas = 4;
	the_options.nparams= 1;


	sprintf(filename, "../../gauss.dat");
	
	if((fptr = fopen(filename, "r")) == NULL){
		fprintf(stderr, "couldn't open file: %s\n", filename);
		exit(1);
	}
	
	// read the data
	input_data = unconstrained_read(filename, &number_lines);
	fprintf(stderr, "read in %d lines \n", number_lines);
	the_options.nmodel_points = number_lines;
	
	process_input_data(input_data, &the_options);
	
	// alloc the data we need to use callEMulator
	// what a rigmarole
	nmodelpts = the_options.nmodel_points;
	xmodel_flat = malloc(sizeof(double)*nparams*nmodelpts);
	training_flat = malloc(sizeof(double)*nmodelpts);
	final_emulated_x = malloc(sizeof(double)*nparams*nemupts);
	final_emulated_y = malloc(sizeof(double)*nemupts);
	final_emulated_variance = malloc(sizeof(double)*nemupts);

	
	for(i = 0; i < nemupts*nparams; i++)
		final_emulated_x[i] = 0.0;
	
	for(i = 0; i < nemupts; i++){
		final_emulated_y[i] = 0.0;
		final_emulated_variance[i] = 0.0;
	}

	
	for(j = 0; j < nparams; j++){
		for(i = 0; i < nmodelpts; i++){
			xmodel_flat[i+j*nmodelpts] = gsl_matrix_get(the_options.xmodel, i, j);
		}
	}

	for(i = 0; i < nmodelpts; i++)
		training_flat[i] = gsl_vector_get(the_options.training, i);

	callEmulator(xmodel_flat, &nparams, training_flat, &nmodelpts, &nthetas, final_emulated_x, &nemupts, final_emulated_y,\
							 final_emulated_variance, &range_min, &range_max);
	
	for(i = 0; i < nemupts; i++)
		printf("%g\t%g\t%g\n", final_emulated_x[i], final_emulated_y[i], final_emulated_variance[i]);

	// now call it again, why? because ooh it doesn't work! 
	
	callEmulator(xmodel_flat, &nparams, training_flat, &nmodelpts, &nthetas, final_emulated_x, &nemupts, final_emulated_y,\
							 final_emulated_variance, &range_min, &range_max);
	

	for(i = 0; i < nemupts; i++)
		printf("%g\t%g\t%g\n", final_emulated_x[i], final_emulated_y[i], final_emulated_variance[i]);


	free(xmodel_flat);
	free(training_flat);
	free(final_emulated_x);
	free(final_emulated_variance);
	free(final_emulated_y);
	gsl_matrix_free(the_options.xmodel);
	gsl_vector_free(the_options.training);
	gsl_vector_free(the_options.thetas);

	for(i = 0; i < number_lines; i++)
		free(input_data[i]);
	
	free(input_data);

	fclose(fptr);
	return(0);
}
	

//! splits up the input data 
/**
 * fixed to work with multiple params using strtok
 * input data is split on tabs and spaces
 * 
 * figures out the min and max of the input data and uses this to set the 
 * range_min and range_max fields in the_options
 * 
 * splits up the raw char input data into real things (xmodel etc)
 * this requires that the following fields in the_options be set correctly
 * nparams, nthetas, nmodel_points
 */
void process_input_data(char** input_data, eopts* the_options){
	int i,j;
	// set these off the first index
	double range_min = 0.0;
	double range_max = 0.0;
	char* split_string;
	double temp_value;
	
	assert(the_options->nmodel_points > 0); 
	// first allocate the buffers in options
	the_options->xmodel = gsl_matrix_alloc(the_options->nmodel_points, the_options->nparams);
	the_options->training = gsl_vector_alloc(the_options->nmodel_points);
	the_options->thetas = gsl_vector_alloc(the_options->nthetas);

	for(i = 0; i < the_options->nmodel_points; i++){
		split_string = strtok(input_data[i], "\t ");
		for( j = 0; j < the_options->nparams; j++){
			assert(split_string != NULL);
			sscanf(split_string, "%lg", &temp_value); 		
			gsl_matrix_set(the_options->xmodel, i, j, temp_value);
			split_string = strtok(NULL, "\t ");
			if(i == 0 && j == 0){
				range_max = temp_value;
				range_min = temp_value;
			}
			if(j ==0){ // only use the zero index param to set the max-min ranges
				if(temp_value > range_max){
					range_max = temp_value;
				}
				if(temp_value < range_min){
					range_min = temp_value;
				}
			}
	 }
		assert(split_string != NULL);		
		sscanf(split_string, "%lg", &temp_value);
		gsl_vector_set(the_options->training, i, temp_value);
	}
	
	fprintf(stderr, "read data in range: %g..%g\n", range_min, range_max);
	the_options->range_min = range_min;
	the_options->range_max = range_max;
	
}
