#include "string.h"
#include "assert.h"
#include "math.h"
#include "modelstruct.h"


/**
 * allocate a modelstruct from the params in optstruct
 */
void alloc_modelstruct(modelstruct* the_model, optstruct* options){
	// row x column
	the_model->xmodel = gsl_matrix_alloc(options->nmodel_points, options->nparams);
	the_model->training_vector = gsl_vector_alloc(options->nmodel_points);
	the_model->thetas = gsl_vector_alloc(options->nthetas);
	the_model->sample_scales = gsl_vector_alloc(options->nparams);
}

/**
 * free a modelstruct
 */
void free_modelstruct(modelstruct* the_model){
	gsl_matrix_free(the_model->xmodel);
	gsl_vector_free(the_model->training_vector);
	gsl_vector_free(the_model->thetas);
	gsl_vector_free(the_model->sample_scales);
}

/**
 * copy a modelstruct from src->dst
 */
void copy_modelstruct(modelstruct* dst, modelstruct* src){
	gsl_matrix_memcpy(dst->xmodel, src->xmodel);
	gsl_vector_memcpy(dst->training_vector, src->training_vector);
	gsl_vector_memcpy(dst->thetas, src->thetas);
	gsl_vector_memcpy(dst->sample_scales, src->sample_scales);
}
	

/**
 * dump a model struct to fptr, in ascii each field is dumped in order they 
 * are defined.
 * vectors take a single line and matricies take nrows lines
 * we use the optstruct to get the sizes of everything
 */
void dump_modelstruct(FILE *fptr, modelstruct* the_model, optstruct *opts){
	int i,j;
	int nparams = opts->nparams;
	int nmp = opts->nmodel_points;
	int nthetas = opts->nthetas;

	for(i = 0; i < nmp; i++){
		for(j = 0; j < nparams; j++){
			fprintf(fptr, "%lf ", gsl_matrix_get(the_model->xmodel, i, j));
		}
		fprintf(fptr, "\n");
	}

	for(i = 0; i < nmp; i++)
		fprintf(fptr, "%lf ", gsl_vector_get(the_model->training_vector, i));
	
	for(i = 0; i < nthetas; i++)
		fprintf(fptr, "%lf ", gsl_vector_get(the_model->thetas, i));
	
	for(i = 0; i < nparams; i++)
		fprintf(fptr, "%lf ", gsl_vector_get(the_model->sample_scales, i));
	
}

/**
 * load a model struct from fptr, we use the optstruct to allocate 
 * the fields in the supplied modelstruct before filling them
 */
void load_modelstruct(FILE* fptr, modelstruct* the_model, optstruct* opts){
	int nparams = opts->nparams;
	int nmp = opts->nmodel_points;
	int nthetas = opts->nthetas;
	int i,j;
	double temp;
	
	// allocate everything first
	the_model->xmodel = gsl_matrix_alloc(nmp, nparams);
	the_model->training_vector = gsl_vector_alloc(nmp);
	the_model->thetas = gsl_vector_alloc(nthetas);
	the_model->sample_scales = gsl_vector_alloc(nparams);


	for(i = 0; i < nmp; i++){
		for(j = 0; j < nparams; j++){
			fscanf(fptr, "%lf ", &temp);
			gsl_matrix_set(the_model->xmodel, i, j, temp);
		}
	}

	for(i = 0; i < nmp; i++){
		fscanf(fptr, "%lf ", &temp);
		gsl_vector_set(the_model->training_vector, i, temp);
	}
	
	for(i = 0; i < nthetas; i++){
		fscanf(fptr, "%lf ", &temp);
		gsl_vector_set(the_model->thetas, i, temp);
	}
	
	for(i = 0; i < nparams; i++){
		fscanf(fptr, "%lf ", &temp);
		gsl_vector_set(the_model->sample_scales, i, temp);
	}

}



/**
 * fill a modelstruct, from a big happy array of chars from the stdin
 */
void fill_modelstruct(modelstruct* the_model, optstruct* options, char** input_data, int number_lines){
	int i,j;
	double temp_value;
	char* split_string;
	char buffer[128];
	gsl_vector *differences = gsl_vector_alloc(options->nmodel_points - 1);
	double min_value, average_value;

	// there's a bug, this can't handle empty lines at the end of the input!
	for(i = 0; i < options->nmodel_points; i++){
		split_string = strtok(input_data[i], "\t ");		
		for(j=0; j < options->nparams; j++){
			//printf("%s\n", split_string);
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

	// compute the average separations
	for(i = 0; i < options->nparams; i++){
		average_value = 0;
		for(j = 0; j < (options->nmodel_points-1); j++){
			gsl_vector_set(differences, j, fabs(gsl_matrix_get(the_model->xmodel, j+1, i) - 
																					gsl_matrix_get(the_model->xmodel, j, i))) ;
			average_value += gsl_vector_get(differences, j);
		}
		// compute the min difference
		min_value = gsl_vector_min(differences);
		// compute the average difference
		average_value /= (options->nmodel_points-1);
		gsl_vector_set(the_model->sample_scales, i, min_value);
		fprintf(stderr, "# param %d min-value %lf average %lf\n", i, min_value, average_value);
	}

	/* turn this off */
	/* fprintf(stderr, "read the following input matrix: %d x %d\n", options->nmodel_points, options->nparams); */
	/* message(buffer, 2); */
	/* print_matrix(the_model->xmodel, options->nmodel_points, options->nparams); */
	/* fprintf(stderr, "the training data is:\n"); */
	/* print_vector_quiet(the_model->training_vector, options->nmodel_points); */


	gsl_vector_free(differences);
	
}
	
