#include "estimate_threaded.h"
#include "maxmultimin.h"


#define SCREWUPVALUE -20000

#define NUMBERTHREADS 1

#ifdef USEMUTEX
// globals for the threads to use
pthread_mutex_t job_counter_mutex;// = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t results_mutex;// = PTHREAD_MUTEX_INITIALIZER;
#else
pthread_spinlock_t job_counter_spin;
pthread_spinlock_t results_spin;
#endif

/**
 * GLOBALS FOR THIS FILE
 *
 * these are used in the estimate_thread_function, because i'm a dumbass and can't figure out a 
 * better way 
 */

/* how many lots of thread_level_tries to do */
int ntries = 1; 
/* mutex protected counter to keep track of completed jobs */
int jobnumber = 0; 
/* global spot for the best thetas to be kept in */
gsl_vector *best_thetas;
/* the best likelyhood we find 
 * have to init this to a large negative value, since the L may be negative
 */
double best_likelyhood_val = SCREWUPVALUE;

/**
 * end of global section
 */

int get_number_cpus(void){
	int ncpus = 0;
	char buffer[256];
	/* now try and find the number of threads by seeing how many cpus we have */
	ncpus = sysconf(_SC_NPROCESSORS_ONLN); // man 3 sysconf (should be posix?)
	if(ncpus == -1){
		ncpus = 1;
	}
	fprintf(stderr, "NCPUS: %d\n", ncpus);
	return(ncpus);
}


//#define DEBUGMODE

/**
 * setup the params_array structure, this is an array of estimate_thetas_params structs which 
 * each contain a copy of the  modelstruct and an opstruct, to be passed to each thread.
 */
void setup_params(struct estimate_thetas_params *params_array, modelstruct* the_model, optstruct* options, int nthreads, int max_tries){
	int i;
	for(i = 0; i < nthreads; i++){
		copy_optstruct(params_array[i].options, options);
		// have to allocate the memory inside each model for each parameter
		alloc_modelstruct(params_array[i].the_model, options);
		copy_modelstruct(params_array[i].the_model, the_model);
		params_array[i].max_tries = max_tries;		
	}
	
}

//! threaded estimate thetas 
/** 
 * uses the lbfgs estimator, see maxlbfgs.c, routines.f 
 * esimate the most likley hyperparams, the number of threads is set to the number of cpus
 * you can switch between mutexes and spinlocks for thread synch which 
 * by defining USEMUTEX (or not and then using spins). 
 * Spinlocks are slightly faster
 */
void estimate_thetas_threaded(modelstruct* the_model, optstruct* options){
	int i;
	int rc = 0; 
	/* thread data */
	/* \bug
	 * at least on os-x there is some kind of bug where the code gets stuck hanging with
	 * all but one of the threads finished, this only happens with some orderings of execution
	 * and only if there are >= 4 threads (on my MBP). don't want to spend too much time trying 
	 * to fix this right now
	 * 
	 * \bug
	 * further issues with R hanging after calling estimateThetas a second time
	 * 
	 */
	int nthreads = 1; // get_number_cpus();


	
	/* force each thread to do at least one of the tries */
	if(ntries < nthreads){
		ntries = nthreads;
	}
	ntries = 1;

	//fprintf(stderr, "nthreads = %d\tntries = %d\n", nthreads, ntries);

	/* how many attempts to maximise should we make */
	/* each thread will make this number of tries and then compare its best values 
	 * to the ones in best_thetas, if it wins it will save them
	 * we only care about the *best* so it doesn't matter if we just throw 
	 * the rest out the window... 
	 */
	int thread_level_tries = 50; 

	pthread_t *threads;
	struct estimate_thetas_params *params;


	threads = MallocChecked(sizeof(pthread_t)*nthreads);
	params = MallocChecked(sizeof(struct estimate_thetas_params)*nthreads);
	
	for(i=0; i < nthreads; i++){
		params[i].my_best = SCREWUPVALUE;
		params[i].the_model = MallocChecked(sizeof(modelstruct));
		params[i].options = MallocChecked(sizeof(optstruct));
	}

	/*
	 * setup the bulk of the parameter structures, but we need to 
	 * setup the random_number afterwards
	 */
	setup_params(params, the_model, options, nthreads, thread_level_tries);

	best_thetas = gsl_vector_alloc(options->nthetas);
	gsl_vector_set_zero(best_thetas);

	// set the jobnumber back to zero otherwise running twice will kill ya
	jobnumber = 0;
	best_likelyhood_val = SCREWUPVALUE;

	/* regular stuff */
	const gsl_rng_type *T;
	int number_steps = 20;

	T = gsl_rng_default;

	/* 
	 * grad_ranges now live in the optstruct
	 */

	/* setup the thread params */
	for(i = 0; i < nthreads; i++){
		// alloc a rng for each thread
		params[i].random_number = gsl_rng_alloc(T);
		// this is blocking right now (slooow)
		gsl_rng_set(params[i].random_number, get_seed_noblock());
	}
	
	#ifdef USEMUTEX
	pthread_mutex_init(&job_counter_mutex, NULL);
	pthread_mutex_init(&results_mutex, NULL);
	#else 
	pthread_spin_init(&job_counter_spin, 0);
	pthread_spin_init(&results_spin, 0);
	#endif

	// create the threads
	for(i = 0; i < nthreads; i++){
		rc = pthread_create(&threads[i], NULL, &estimate_thread_function, &params[i]);
		if(rc){
			fprintf(stderr, "pthread_create (%d)  err: %d\n", i, rc);
			perror(NULL);
			exit(EXIT_FAILURE);
		}
	}
	
	// wait to rejoin
	for(i = 0; i < nthreads; i++){
		if( (rc = pthread_join(threads[i], NULL)) != 0){
			fprintf(stderr, "pthread_join (%d)  err: %d\n", i, rc);
			perror(NULL);
			exit(EXIT_FAILURE);
		}
	}

	#ifdef USEMUTEX
	// now kill the mutexs
	pthread_mutex_destroy(&job_counter_mutex);
	pthread_mutex_destroy(&results_mutex);
	#else 
	pthread_spin_destroy(&job_counter_spin);
	pthread_spin_destroy(&results_spin);
	#endif

	// check the local best values from the threads
	printf("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
	for(i = 0; i < nthreads; i++)
		printf("%d\t%lf\n", i, params[i].my_best);
	printf("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");

	/* fprintf(stderr, "final best L: %g\n", best_likelyhood_val); */
	/* fprintf(stderr, "THETAS WE WILL USE: \t"); */
	/* print_vector_quiet(best_thetas, options->nthetas); */

	
	// tear down the thread params
	for(i = 0; i < nthreads; i++){
		gsl_rng_free(params[i].random_number);
		free_modelstruct(params[i].the_model);
		gsl_matrix_free(params[i].options->grad_ranges);
		gsl_matrix_free(params[i].h_matrix); // hehe was free'ing this way too early before, gives strange behaviour
		free(params[i].the_model);
		free(params[i].options);
	}

	// copy the global best_theta into the one provided 
	gsl_vector_memcpy(the_model->thetas, best_thetas);
	// now free best_thetas
	gsl_vector_free(best_thetas);

	free(threads);
	free(params);
	
}	

// THIS USES GLOBAL VARIABES DEFINED ABOVE, WATCH OUT!
// what the threads actually call, put this in here so that the function
// will share the same scope as the rest of the crap here
void* estimate_thread_function(void* args){
	// cast the args back
	struct estimate_thetas_params *params = (struct estimate_thetas_params*) args;
	int next_job;
	int rc;
	pthread_t my_id = pthread_self();
	double my_theta_val = 0.0; /* this is the goodness of your current evaluation */
	while(1){
		/* see if we've done enough */
		#ifdef USEMUTEX
		rc = pthread_mutex_lock(&job_counter_mutex);
		if(rc != 0){
			fprintf(stderr, "ERR: pthread_mutex_lock(&job_counter_mutex) rc is %d\n", rc); 
			perror(NULL); 
			exit(EXIT_FAILURE);
		}
		#else 
		pthread_spin_lock(&job_counter_spin);
		#endif
		if(jobnumber == ntries){
			next_job = -1;
		} else {
			next_job = jobnumber;
			jobnumber++;
			/* printf("job: %d by ", next_job);  */
			/* fprintPt(stdout, my_id); */
			/* printf("\n"); */
		}
		//printf("jobnumber = %d\tnext_job = %d\tntries = %d\n", jobnumber, next_job, ntries);
		/* now we can unlock the job counter */
		#ifdef USEMUTEX		
		pthread_mutex_unlock(&job_counter_mutex);
		#else 
		pthread_spin_unlock(&job_counter_spin);
		#endif
		
		/* we're done so stop */
		if(next_job == -1)
			break;

		/* just support LBFGS maximisation now */
		//maxWithLBFGS(params);
		/* switch to multimin for testing */
		maxWithMultiMin(params);
		my_theta_val = params->lhood_current;

		/* printf("my_theta_val = %g\n", my_theta_val); */
		/* printf("params->best = %g\n", params->my_best); */
		/* printf("params->lhood = %g\n", params->lhood_current); */

		/* store you local best value too */
		if(my_theta_val > params->my_best)
			params->my_best = my_theta_val;
		
		
		#ifdef USEMUTEX
		rc =pthread_mutex_lock(&results_mutex);
		if(rc != 0){
			fprintf(stderr, "ERR: pthread_mutex_lock(&results_mutex) rc is %d\n", rc); 
			perror(NULL); 
			exit(EXIT_FAILURE);
		}
		#else 
		pthread_spin_lock(&results_spin);
		#endif
		/* printf("results locked: "); */
		/* fprintPt(stdout, my_id); */
		/* printf("\n"); */

		if(my_theta_val > best_likelyhood_val){
			// this thread has produced better thetas than previously there
			gsl_vector_memcpy(best_thetas, params->the_model->thetas); // save them
			// save the new best too
			best_likelyhood_val = my_theta_val;
			fprintPt(stdout, my_id);
			printf(" won with %g\n",  my_theta_val);
		}
		#ifdef USEMUTEX
		pthread_mutex_unlock(&results_mutex);
		#else 
		pthread_spin_unlock(&results_spin);
		#endif
		/* printf("results unlocked by:"); */
		/* fprintPt(stdout, my_id); */
		/* printf("\n"); */

	}
	// and relax...
	printf("thread: ");
	fprintPt(stdout, my_id);
	printf(" is done\n");
	return NULL;
}

/**
 * from stackoverflow.
 * print pthread id in a nice way?
 * well if you like hex anyway...
 */
void fprintPt(FILE *f, pthread_t pt) {
  unsigned char *ptc = (unsigned char*)(void*)(&pt);
  fprintf(f, "0x");
  for (size_t i=0; i<sizeof(pt); i++) {
    fprintf(f, "%02x", (unsigned)(ptc[i]));
  }
}
