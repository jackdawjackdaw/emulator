#include "stdio.h"
#include "stdlib.h"
#include "string.h"


typedef struct region{
	int region_start;
	int region_stop;
	// snip....
} region;

int resize_region_array(region** the_array, int current_length, int grow_length);
void *MallocChecked(size_t size);
void unix_error(char *msg);


//! copy source->target
/**
 * copies the contents of the source region into the target region, useful because i can't make memcpy just take a 
 * block of these guys and do it for me.
 */
void copy_region_array(region* target, region* source, int length){
	int i;
	for (i = 0; i < length; i++){
		target[i].region_start = source[i].region_start;
		target[i].region_stop = source[i].region_stop;
	}	
}

int resize_region_array(region** the_array, int current_length, int grow_length){
	size_t current_size = current_length*sizeof(region);
	size_t new_size = current_size + sizeof(region)*grow_length; // the grown array
	region* buffer = MallocChecked(current_size);
	copy_region_array(buffer, *the_array, current_length);
	free(*the_array);
	*the_array = MallocChecked(new_size);
	copy_region_array(*the_array, buffer, current_length);
	free(buffer);
	// and return the new length
	return(current_length+grow_length);
}

//! checks for null
void *MallocChecked(size_t size){
	void *r = malloc(size);
	
	if( r == NULL){
		fprintf(stderr, "memory wasn't allocated");
		exit(1);
	}
	
	return(r);
}


int main (void){
	int i;
	int initial_array_length = 1;
	int grow_array = 10;
	int array_length = initial_array_length;
	region* region_array = MallocChecked(array_length*sizeof(region));
	

	array_length = resize_region_array(&region_array, array_length, grow_array);
	
	


	free(region_array);
	return(0);
}
