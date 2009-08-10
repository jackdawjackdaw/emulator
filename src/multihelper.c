
void dump_result(emuResult *res, FILE *fptr);
void alloc_emuRes(emuResult *thing, eopts *options);
void free_eopts(eopts* options);
void free_emuRes(emuResult *thing);
void print_splits(gsl_matrix* splits, int n);

	
//! setup an emuResult struct from the given options
/**
 * copy the right parts of the options into the result and alloc the data structures inside.
 * after this the emuResult has to be freed with free_emuRes or you'll leak the arrays inside. 
 * 
 * Also note that the gsl allocs are not checked.
 */ 
void alloc_emuRes(emuResult *thing, eopts *options){
	int n = options->nemu_points;
	int np = options->nparams;
	thing->nemu_points = n;
	thing->new_x = gsl_matrix_alloc(n, np);
	thing->new_mean = gsl_vector_alloc(n);
	thing->new_var = gsl_vector_alloc(n);
	thing->nparams = np;
}

//! frees an emures struct
void free_emuRes(emuResult *thing){
	gsl_matrix_free(thing->new_x);
	gsl_vector_free(thing->new_mean);
	gsl_vector_free(thing->new_var);
}

//! frees an eopts struct
void free_eopts(eopts* options){
	gsl_matrix_free(options->xmodel);
	gsl_vector_free(options->training);
	gsl_vector_free(options->thetas);
}

//! debugging
void print_splits(gsl_matrix* splits, int n){
	int i =0;
	for(i = 0; i < n; i++){
		printf("%d (%g..%g)\n", i, gsl_matrix_get(splits, i, 0), gsl_matrix_get(splits, i, 1));
	}
}



void dump_result(emuResult *res, FILE *fptr){
	int i;
	double goodness = 0.0;
	for(i = 0; i < res->nemu_points; i++){
		fprintf(fptr, "%g\t%g\t%g\t", gsl_matrix_get(res->new_x, i, 0), gsl_vector_get(res->new_mean, i), gsl_vector_get( res->new_var, i));
		goodness = 1 / pow(gsl_vector_get(res->new_var, i), 2.0);
		fprintf(fptr, "%g\n", goodness);
	}
}

