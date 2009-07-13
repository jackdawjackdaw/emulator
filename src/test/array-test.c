#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h> // high level interface to blas
#include <gsl/gsl_linalg.h> // lapack

void print_matrix(gsl_matrix* m, int nx, int ny);

int main (void) {
	int i, j;
	int signum = 0;
	double determinant_m = 0.0;
	gsl_permutation  *m_LU_permutation = gsl_permutation_alloc(8);
	gsl_matrix *m_test = gsl_matrix_alloc(8, 8);
	gsl_matrix *m_original = gsl_matrix_alloc(8,8);
	gsl_matrix *m_inverted = gsl_matrix_alloc(8,8);
	
	/*for (i = 0; i < 8; i++){
		for(j = 0; j < 8; j++){
			gsl_matrix_set(m_test, i, j, 1);
		}
	}*/
	
	gsl_matrix_set(m_test, 0, 7, 1);
	gsl_matrix_set(m_test, 7,0, 1);
	gsl_matrix_set(m_test, 3,2, 5);
	gsl_matrix_set(m_test, 5,2, 5);


	for( i = 0; i < 8; i++){
		gsl_matrix_set(m_test, i, i, 3);
	}


	// set the oritinal one too
	gsl_matrix_memcpy(m_original, m_test); 	
	// now print the matrix out
	//print_matrix(m_test, 8, 8);
	
	printf("\n");
	print_matrix(m_original, 8,8);

	// simple matrix operations
	//gsl_matrix_add(m_test, m_test);
	//gsl_matrix_add_constant(m_test, 1.05);

	//print_matrix(m_test, 8,8);

	// creates the lu decomp of the matrix
	// need to allocate the permutation
	gsl_linalg_LU_decomp(m_test, m_LU_permutation, &signum);

	// calculate the determinant of m_test from the lu
	determinant_m = gsl_linalg_LU_det(m_test, signum);
	printf("det_m = %g\n", determinant_m);
	// invert the matrix
	gsl_linalg_LU_invert(m_test, m_LU_permutation, m_inverted);

	// multiply the inverse and the matrix together  (can't use the LU one bc it's fucked up)
	// m_test should be I
	gsl_blas_dgemm( CblasNoTrans,CblasNoTrans , 1.0, m_original, m_inverted, 0.0, m_test);

	print_matrix(m_test, 8,8);


	gsl_permutation_free(m_LU_permutation);
	gsl_matrix_free(m_inverted);
	gsl_matrix_free(m_test);
	return(0);
}

void print_matrix(gsl_matrix* m, int nx, int ny){
	int i,j;
	for(i = 0; i < nx; i++){
		for(j= 0; j < ny; j++){
			printf("%g ", gsl_matrix_get(m, i, j));
		}
		printf("\n");
 	}
}
