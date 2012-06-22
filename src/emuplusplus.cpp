#include "emuplusplus.h"
#include <string>
#include <vector>
#include <fstream>

#include <stdio.h> // for the c io to load the statefile
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/* the c-interface to the emulator */
#include "emulatorstruct.h" 
#include "modelstruct.h"
#include "multi_modelstruct.h"

/**
 * 22.06.2012, ccs
 * cec24@phy.duke.edu
 * 
 * a c++ interface to an emulator, for use in MCMC or similar, 
 * requries a trained emulator statefile as produced by interactive_emulator.c
 * 
 * for reference to the interface see parameterset.cc and mcmc.cc in madai-stat
 * 
 * interface should expose the following required fns:
 *
 * emulator(string StateFilePath, bool PcaOnly)
 * @argument filepath -> location of the emulator state file
 * @argument PcaOnly -> if true query emulator will return the results in the pca space (z)
 * without rotating them back to the true output space (y)
 * 
 * creates an emulator object
 * filepath should point to the emulator state file
 * 
 * DestroyEmulator(...)
 * kills everything nicely, the destructor
 * 
 * QueryEmulator(ParameterSet xpoint, vector<double> &Means, vector<double> &Errors);
 * sample the emulator at the point in the parameter space xpoint, 
 * @return: Means is a vector of the emulator means at xpoint
 * @return: Errors is a vector of the emulator error (sqrt(var)) at xpoint
 * 
 * Optional (but useful) functions:
 *
 * GetEmulatorPCADecomp(...)
 * returns the evalues and evectors of the pca decomp for a multivariate emulator, need this to use the emulator 
 * in pca-only mode.
 * 
 */

// default constructor
emulator::emulator(string StateFilePath){
	self->outputPCAValues = false; // default to outputing values in the Y space
	self->StateFilePath = StateFilePath;
	
	self->the_model = 0;
	FILE *fptr = fopen(self->StateFilePath.c_str(), "r");
	// this will load the modelstruct
	self->the_model = load_multi_modelstruct(fptr);
	fclose(fptr);
	
	if(self->the_model == 0){
		cerr << "error opening statefile: " << StateFilePath << endl;
		exit(EXIT_FAILURE); // could throw here?
	}
	
	self->the_emulator = 0;
	self->the_emulator = alloc_multi_emulator(self->the_model);
	if(self->the_emulator == 0){
		cerr << "error initializing multi_emulator" << endl;
		exit(EXIT_FAILURE);
	}

	self->number_outputs = self->the_model->nt;
	self->number_params = self->the_model->nparams;
	
	self->the_emulate_point = gsl_vector_alloc(self->number_params);
	self->the_emulate_mean = gsl_vector_alloc(self->number_outputs);
	self->the_emulate_var = gsl_vector_alloc(self->number_outputs);
}

// same as the default but sets the output flag and changes the number of outputs we expect
emulator::emulator(string StateFilePath, bool PcaOnly):emulator(StateFilePath){
	self->outputPCAValues = PcaOnly;
	self->number_outputs = self->the_model->nr;
	
	// free and re-alloc the mean and var since (nt != nr)
	gsl_vector_free(self->the_emulate_mean);
	gsl_vector_free(self->the_emulate_var);
	self->the_emulate_mean = gsl_vector_alloc(self->number_outputs);
	self->the_emulate_var = gsl_vector_alloc(self->number_outputs);
	
}


/**
 * query the emulator at a given point
 * 
 * if outputPCAValues is true, output is the Z values in PCA space (length is self->the_model->nr)
 * if outputPCAValues is false, output is the Y values in the real space (length is self->the_model->nt)
 * 
 * @argument xpoint: the location in the parameter space where to sample the emulator
 * @return means: filled with the emulator mean for the outputs at the point
 * @return errors: filled with the sqrt of the emulator variance for the outputs at a point
 *
 */
void emulator::QueryEmulator(const void<double> xpoint, vector<double> &Means, vector<double> &Errors){
	if(xpoint.size != self->number_params){
		cerr << "Error::QueryEmulator called with incorrect number of dimensions in xpoint" << endl;
		cerr << "xpoint.length: " << xpoint.size << " emulator->number_params: " << self->number_params << endl;
		exit(EXIT_FAILURE);
	}

	// setup the gsl vectors
	for(int i = 0; i < self->number_params; i++){ // copy in the xpoint cpts into the internal gsl vector
		gsl_vector_set(self->the_emulate_point, i, xpoint[i]);
	}
	
	gsl_vector_set_zero(self->the_emulate_mean);
	gsl_vector_set_zero(self->the_emulate_var);

	// call the emulator
	if(!self->outputPCAValues){
		// do the output in the real space
		emulate_point_multi(self->the_emulator, self->the_emulate_point, self->the_emulate_mean, self->the_emulate_var);
	} else{
		// do the output in the pca space
		emulate_point_multi_pca(self->the_emulator, self->the_emulate_point, self->the_emulate_mean, self->the_emulate_var);		
	}
	
	// now copy out the results
	for(int i = 0; i < self->number_outputs; i++){
		Means.push_back(gsl_vector_get(self->the_emulate_mean, i));
		// note that we're filling the errrors with the square-root of the variance...
		Errors.push_back(sqrt(gsl_vector_get(self->the_emulate_var, i)));
	}
	
}

/**
 * samples the emulator at the location specified by xpoint
 *
 * if self->outputPCAValues is true, output is the Z values in PCA space (length is self->the_model->nr)
 * if self->outputPCAValues is false, output is the Y values in the real space (length is self->the_model->nt)
 *
 *
 * @return Means is the vector of emulated means
 * @return Errors is the vector of emulated s.d's
 */
void emulator::QueryEmulator(const ParameterSet xpoint, vector<double> &Means, vector<double> &Errors){
// just pass through to the simpler method
	QueryEmulator(xpoint.Values, Means, Errors);
}

/**
 * get the pca decomp from the emulator, you'll need this to rotate your 'real' output into 
 * the pca space if you want to compare
 *
 * the PCA decomp on a set Y of training data, each observation being a vector of length nt and 
 * with the set Y being n long:
 * 1) compute the column mean:  pca_mean_j = (1/n) Sum_{i=1}^{n} Y_j
 * where pca_mean is then a vector of length nt.
 * 2) compute the observed covariance matrix (nt x nt)  pca_cov_matrix_ij = (1/n) (Y-pca_mean)_ik (Y-pca_mean)_kj
 * 3) eigendecompose this covarance matrix: pca_cov_matrix = U^{-1} LAMBDA U
 * where U is an (nt x nt) (column) matrix of the eigenvectors and LAMDBA is a diagonal (nt x nt) matrix of the 
 * eigenvalues.
 *
 * @requires: the vectors have been allocated
 * 
 * @return pca_evals is set to the eigenvalues of the pca-decomp, sorted in descending absolute value (is length nr)
 * @return pca_evecs is set to the eigenvectors of the pca-decomp, the eigenvectors are stored in the columns and sorted to match the eigenvalues (is (rows)nt by (columns)nr)
 * @return pca_mean is set to the mean of the training data (and is length nt)
 */

void emulator::getEmulatorPCA(vector<double> *pca_evals, vector< vector<double> > *pca_evecs, 
															vector<double> *pca_mean){
	
	for(int i = 0; i < self->the_model->nr; i++){
		pca_evals.push_back(gsl_vector_get(self->the_model->pca_evals_r, i));
	} 
	
	for(int i = 0; i < self->the_model->nt; i++){
		pca_mean.push_back(gsl_vector_get(self->the_model->pca_evals_r, i));
	}

	// pca_evecs must have been allocated...
	for(int i = 0; i< self->the_model->nt; i++){
		for(int j = 0; j < self->the_model->nr; j++){
			pca_evecs[i][j] = gsl_matrix_get(self->pca_evecs_r, i, j);
		}
	}
		
}


/**
 * default destructor for the emulator
 */
void emulator::~emulator(){
	// free the c-structures
	free_multi_emulator(self->the_emulator);
	free_multimodelstruct(self->the_model); // this is a memory leaker  :(
	gsl_vector_free(self->the_emulate_point);
	gsl_vector_free(self->the_emulate_mean);
	gsl_vector_free(self->the_emulate_var);
}


