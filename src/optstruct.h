#ifndef __INC_OPTSTRUCT__
#define __INC_OPTSTRUCT__

//! holds command line options
/** 
 * designed to hold basic command line 
 * options
 */
typedef struct optstruct{
	int nthetas;
	int nparams;
	int nmodel_points;
	int nemulate_points;
	double emulate_min;
	double emulate_max;
	char  filename[128];
	char outputfile[128];
} optstruct;

#endif
