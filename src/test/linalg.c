#include "stdio.h"
#include "stdlib.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

const gsl_rng_type *T;
gsl_rng *r;

int main (void){
	int n = 256;
	int i,j;
	gsl_matrix * test_matrix = gsl_matrix_alloc(n,n);
	gsl_matrix * inverse = gsl_matrix_alloc(n,n);
	gsl_permutation *lu_permutation = gsl_permutation_alloc(n);
	
	int lu_signum = 0;
	
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			gsl_matrix_set(test_matrix, i, j, gsl_ran_gaussian(r, 0.1));
		}
	}
	
	gsl_linalg_LU_decomp(test_matrix, lu_permutation, &lu_signum);
	gsl_linalg_LU_invert(test_matrix, lu_permutation, inverse);

	gsl_matrix_free(test_matrix);
	gsl_matrix_free(inverse);

	gsl_rng_free(r);
	return(0);
}
