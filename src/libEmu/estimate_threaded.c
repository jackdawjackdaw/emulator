#include "estimate_threaded.h"


// this lives in libEmu/emulator.c
extern emulator_opts the_emulator_options;


//#define NUMBERTHREADS 2

#ifdef USEMUTEX
// globals for the threads to use
pthread_mutex_t job_counter_mutex;// = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t results_mutex;// = PTHREAD_MUTEX_INITIALIZER;
#else
pthread_spinlock_t job_counter_spin;
pthread_spinlock_t results_spin;
#endif

/* how many lots of thread_level_tries to do */
int ntries = 2; 
/* mutex protected counter to keep track of completed jobs */
int jobnumber = 0; 
/* global spot for the best thetas to be kept in */
gsl_vector *best_thetas;
/* the best likelyhood we find */
double best_likelyhood_val = -1000;

int get_number_cpus(void){
	int ncpus = 0;
	/* now try and find the number of threads by seeing how many cpus we have */
	ncpus = sysconf(_SC_NPROCESSORS_ONLN); // man 3 sysconf
	fprintf(stderr, "NCPUS: %d\n", ncpus);
	return(ncpus);
}

// hack
#define DEBUGMODE

//! threaded estimate thetas 
/** 
 * uses the nelder mead estimator (or the probably broken) bfgs method to 
 * esimate the most likley hyperparams, the number of threads is set to the number of cpus
 * you can switch between mutexes and spinlocks for threadsynch which 
 * by defining USEMUTEX (or not and then using spins). 
 * Spinlocks are slightly faster but they are probably not universally supported...
 */
void estimate_thetas_threaded(gsl_matrix* xmodel_input, gsl_vector* training_vector, gsl_vector* thetas, optstruct* options){
	
	/* thread data */
	int nthreads = get_number_cpus();

	#ifdef DEBUGMODE
	nthreads = 1;
  #endif


	fprintf(stderr, "nthreads = %d\n", nthreads);
	
	/* how many attempts to maximise should we make */
	/* each thread will make this number of tries and then compare its best values 
	 * to the ones in best_thetas, if it wins it will save them
	 * we only care about the *best* so it doesn't matter if we just throw 
	 * the rest out the window... 
	 */
	int thread_level_tries = 5; 
	if(nthreads > 2) {
		thread_level_tries = thread_level_tries / nthreads;		
	}
	fprintf(stderr, "thread_level_tries %d\n", thread_level_tries);

	#ifdef DEBUGMODE
	thread_level_tries = 1;
	#endif

	pthread_t *threads;
	struct estimate_thetas_params *params;

	threads = MallocChecked(sizeof(pthread_t)*nthreads);
	params = MallocChecked(sizeof(struct estimate_thetas_params)*nthreads);

	best_thetas = gsl_vector_alloc(options->nthetas);



	
	// set the jobnumber back to zero otherwise running twice will kill ya
	jobnumber = 0;
	best_likelyhood_val = -1000;

	/* regular stuff */
	const gsl_rng_type *T;
	
	int i; 
	int number_steps = 20;
	gsl_matrix *grad_ranges = gsl_matrix_alloc(options->nthetas, 2);
	T = gsl_rng_default;

	/* set the ranges for the initial values of the NM lookup, 
	 * might want to adjust these as required etc, but whatever */
	/* \TODO replace this this set_likelyhood_ranges ? */
	for(i = 0; i < options->nthetas; i++){
		gsl_matrix_set(grad_ranges, i, 0, 0.13);
		gsl_matrix_set(grad_ranges, i, 1, 2.0);
		gsl_vector_set(best_thetas, i, 0.0);
	}

	if(the_emulator_options.usematern ==0){
		// hackity hack, force the nugget to be small
		gsl_matrix_set(grad_ranges, 1, 0, 0.00001);
		gsl_matrix_set(grad_ranges, 1, 1, 0.1);
	} else {
		// also force the nugget to be small for the matern
		gsl_matrix_set(grad_ranges, 3, 0, 0.000001);
		gsl_matrix_set(grad_ranges, 3, 1, 0.001);
	}
		
	for(i = 0; i < options->nthetas;i++)
		fprintf(stderr, "%d %g %g\n", i, gsl_matrix_get(grad_ranges, i, 0), gsl_matrix_get(grad_ranges, i, 1));
	


	/* setup the thread params */
	for(i = 0; i < nthreads; i++){
		// alloc a rng for each thread
		params[i].random_number = gsl_rng_alloc(T);
		// this is blocking right now (slooow)
		gsl_rng_set(params[i].random_number, get_seed_noblock());
		// not sure about this, perhaps each thread should make 10 tries
		params[i].max_tries = thread_level_tries;
		params[i].thetas = gsl_vector_alloc(options->nthetas);
		params[i].grad_ranges = gsl_matrix_alloc(options->nthetas, 2);		
		params[i].model_input = gsl_matrix_alloc(options->nmodel_points, options->nparams);
		params[i].training_vector = gsl_vector_alloc(options->nmodel_points);
		params[i].nmodel_points = options->nmodel_points;
		params[i].nthetas = options->nthetas;
		params[i].nparams = options->nparams;		
		params[i].number_steps = number_steps;
		params[i].nregression_fns = options->nregression_fns;
		// now actually copy the stuff into the vectors / matrices
		gsl_vector_memcpy(params[i].thetas, thetas);
		gsl_matrix_memcpy(params[i].grad_ranges, grad_ranges);
		gsl_matrix_memcpy(params[i].model_input, xmodel_input);
		gsl_vector_memcpy(params[i].training_vector, training_vector);
		
	}
	
	#ifdef USEMUTEX
	// didn't know you needed to do this?
	pthread_mutex_init(&job_counter_mutex, NULL);
	pthread_mutex_init(&results_mutex, NULL);
	#else 
	pthread_spin_init(&job_counter_spin, 0);
	pthread_spin_init(&results_spin, 0);
	#endif

	// create the threads
	for(i = 0; i < nthreads; i++)
		pthread_create(&threads[i], NULL, &estimate_thread_function, &params[i]);
	
	// wait to rejoin
	for(i = 0; i < nthreads; i++)
		pthread_join(threads[i], NULL);

	#ifdef USEMUTEX
	// now kill the mutexs
	pthread_mutex_destroy(&job_counter_mutex);
	pthread_mutex_destroy(&results_mutex);
	#else 
	pthread_spin_destroy(&job_counter_spin);
	pthread_spin_destroy(&results_spin);
	#endif

	fprintf(stderr, "final best L: %g\n", best_likelyhood_val);
	fprintf(stderr, "THETAS WE WILL USE: \t");
	print_vector_quiet(best_thetas, options->nthetas);

	// tear down the thread params
	for(i = 0; i < nthreads; i++){
		gsl_rng_free(params[i].random_number);
		gsl_matrix_free(params[i].grad_ranges);
		gsl_matrix_free(params[i].model_input);
		gsl_vector_free(params[i].thetas);
		gsl_vector_free(params[i].training_vector);
	}

	// copy the global best_theta into the one provided 
	gsl_vector_memcpy(thetas, best_thetas);
	// now free best_thetas
	gsl_vector_free(best_thetas);
	gsl_matrix_free(grad_ranges);

	free(threads);
	free(params);
	
}	

// THIS USES GLOBAL VARIABES DEFINED ABOVE, WATCH OUT!
// what the threads actually call, put this in here so that the function
// will share the same scope as the rest of the crap here
void* estimate_thread_function(void* args){
	// cast the args back
	struct estimate_thetas_params *p = (struct estimate_thetas_params*) args;
	int next_job;
	unsigned long my_id = pthread_self();
	double my_theta_val = 0.0;
	while(1){
		/* see if we've done enough */
		#ifdef USEMUTEX
		pthread_mutex_lock(&job_counter_mutex);
		#else 
		pthread_spin_lock(&job_counter_spin);
		#endif
		if(jobnumber == ntries){
			next_job = -1;
		} else {
			next_job = jobnumber;
			jobnumber++;
			printf("job: %d by %lu\n", next_job, my_id); 
		}
		/* now we can unlock the job counter */
		#ifdef USEMUTEX		
		pthread_mutex_unlock(&job_counter_mutex);
		#else 
		pthread_spin_unlock(&job_counter_spin);
		#endif
		
		/* we're done so stop */
		if(next_job == -1)
			break;
		
		#ifdef NELDER
		/* else we do the nelder mead stuff */
		nelderMead(p->random_number, p->max_tries, p->number_steps, p->thetas, p->grad_ranges, p->model_input, p->training_vector, p->nmodel_points, p->nthetas, p->nparams, p->nregression_fns);
		#elif BFGS
		maxWithBFGS(p->random_number, p->max_tries, p->number_steps, p->grad_ranges, p->model_input, p->training_vector, p->thetas,	\
								p->nmodel_points, p->nthetas, p->nparams, p->nregression_fns);
		#else 
		maxWithLBFGS(p->random_number, p->max_tries, p->number_steps, p->grad_ranges, p->model_input, p->training_vector, p->thetas, p->nmodel_points, p->nthetas, p->nparams, p->nregression_fns);
		#endif


		// kind of sneakily calling into the maximise.c api (aah well...)
		// won't work without some more fiddling
		my_theta_val = evalLikelyhood(p->thetas, p->model_input, p->training_vector, p->nmodel_points, p->nthetas, p->nparams, p->nregression_fns);

		#ifdef USEMUTEX
		pthread_mutex_lock(&results_mutex);
		#else 
		pthread_spin_lock(&results_spin);
		#endif
		printf("results locked by %lu\n", my_id);
		if(my_theta_val > best_likelyhood_val){
			// this thread has produced better thetas than previously there
			gsl_vector_memcpy(best_thetas, p->thetas); // save them
			// save the new best too
			best_likelyhood_val = my_theta_val;
			printf("thread %lu, won with %g\n", my_id, my_theta_val);
		}
		#ifdef USEMUTEX
		pthread_mutex_unlock(&results_mutex);
		#else 
		pthread_spin_unlock(&results_spin);
		#endif
		printf("results unlocked by: %lu\n", my_id);
	}
	// and relax...
	return NULL;
}
