#ifndef __INCLUDE_EMUPLUSPLUS__
#define __INCLUDE_EMUPLUSPLUS__



#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <string>
#include <vector>

extern "C"{
 #include "multi_modelstruct.h"
 #include "multivar_support.h"
}


using namespace std;

/**
 * @file emuplusplus.h
 * \brief a simple c++ interface to a trained emulator
 * 
 * Instance of the emulator class are built from a trained emulator, 
 * the emulator can then be sampled at a point in the parameter space by the function QueryEmulator
 * 
 * They must be initialized from an interactive_emulator statefile, this includes all the data for 
 * a trained emulator
 * 
 */

class emulator{
 public:
	emulator(string StateFilePath); // default constructor will return Y values
	emulator(string StateFilePath, bool PcaOnly);
	~emulator();
	
	/** query the emulator with a vector xpoint in the parameter space
	 * could add a method to output the means, errors and a covaraince matrix */
		void QueryEmulator(const vector<double> &xpoint, vector<double> &Means, vector<double> &Errors);

	/**
	 * get the emulator pca decomp
	 */
	void getEmulatorPCA(vector<double> *pca_evals, vector< vector<double> > *pca_evecs, vector<double> *pca_mean);

	int getRegressionOrder(void){return(the_model->regression_order);};
	int getCovFnIndex(void){return(the_model->cov_fn_index);};

	int number_params;
	int number_outputs;
	
 private:

	
	// if true the values are output in the pca space, otherwise they're in the real space
	bool outputPCAValues; 
	string StateFilePath;
	multi_modelstruct *the_model;  // the c structure which defines the model
	multi_emulator *the_emulator; // the c structure which defines the emulator

	gsl_vector *the_emulate_point;
	gsl_vector *the_emulate_mean;
	gsl_vector *the_emulate_var;
	
};


#endif
