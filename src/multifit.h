
typedef struct eopts{
	int nmodel_points;
	int nemu_points;
	int nparams;
	int nthetas;
	double range_min;
	double range_max;
	gsl_matrix *xmodel;
	gsl_vector *training;
	gsl_vector *thetas;
} eopts;


void emulate_region(gsl_vector *new_x, gsl_vector* emulated_mean, gsl_vector* emulated_variance , eopts* options);
void estimate_region(eopts* options, gsl_rng *random);
void evaluate_region(gsl_vector* new_x, gsl_vector* new_mean, gsl_vector* new_variance, eopts* options, gsl_rng * random);
int is_smooth(double smooth_val, gsl_vector* xemu, gsl_vector* mean_emu, gsl_vector* var_emu, eopts* options);
double get_mse( double mean, double variance);
