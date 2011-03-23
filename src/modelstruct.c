#include "modelstruct.h"
#include "libEmu/emulator.h"


/**
 * allocate a modelstruct from the params in optstruct
 */
void alloc_modelstruct(modelstruct* the_model, optstruct* options){
	the_model->xmodel = gsl_matrix_alloc(options->nmodel_points, options->nparams);
	the_model->training_vector = gsl_vector_alloc(options->nmodel_points);
	the_model->thetas = gsl_vector_alloc(options->nthetas);
	the_model->options = options;
}

/**
 * free a modelstruct
 */
void free_modelstruct(modelstruct* the_model){
	gsl_matrix_free(the_model->xmodel);
	gsl_vector_free(the_model->training_vector);
	gsl_vector_free(the_model->thetas);
}

/**
 * copy a modelstruct from src->dst
 */
void copy_modelstruct(modelstruct* dst, modelstruct* src){
	gsl_matrix_memcpy(dst->xmodel, src->xmodel);
	gsl_vector_memcpy(dst->training_vector, src->training_vector);
	gsl_vector_memcpy(dst->thetas, src->thetas);
}
	



/**
 * fill a modelstruct, from a big happy array of chars from the stdin
 */
void fill_modelstruct(modelstruct* the_model, optstruct* options, char** input_data, int number_lines){
	int i,j;
	double temp_value;
	char* split_string;
	char buffer[128];

	// there's a bug, this can't handle empty lines at the end of the input!
	for(i = 0; i < options->nmodel_points; i++){
		split_string = strtok(input_data[i], "\t ");		
		for(j=0; j < options->nparams; j++){
			printf("%s\n", split_string);
			// split string into tab or space tokens
			// each time you do it split_string is pointed to the next block
			// it will come up null when you're done
			assert(split_string != NULL);
			sscanf(split_string, "%lg", &temp_value);
			//fprintf(stderr,"param: %s\n", split_string);
			gsl_matrix_set(the_model->xmodel, i, j, temp_value);
			split_string = strtok(NULL, "\t ");
		}
		assert(split_string != NULL);
		sscanf(split_string,"%lg", &temp_value);
		//fprintf(stderr,"train: %s\n", split_string);
		gsl_vector_set(the_model->training_vector, i, temp_value);
	}

	sprintf(buffer, "read the following input matrix: %d x %d\n", options->nmodel_points, options->nparams);
	message(buffer, 2);
	print_matrix(the_model->xmodel, options->nmodel_points, options->nparams);
	fprintf(stderr, "the training data is:\n");
	print_vector_quiet(the_model->training_vector, options->nmodel_points);
	
}
	
