
#include "ioread.h"

/**
 * some useful io stuff
 */



//! read a file of unknown length
/**
 * the file can have comments in it, as long as they are preceeded by #s
 *
 * @return the read in data
 * @param line_count thenumebr of lines we read
 */
char** unconstrained_read(char* filename, int* line_count_final){
	void copy_char_arrays(char** dest, char** source, int lx, int ly);
	int i;
	FILE *fptr;
	char** input_data;
	char** temp_buffer; 
	int init_number_lines = 20;
	int actual_number_lines = init_number_lines;
	int previous_number_lines;
	int line_width = 256; // assume that lines are not wider than this... (right?)
	char temp_line[line_width];
	int line_count = 0;
	char* is_end = 0;
	int buffer_size;
	fptr = fopen(filename, "r");
	if(fptr == NULL){
		fprintf(stderr, "could not open inputfile\n");
		exit(1);
	}



	input_data = malloc(sizeof(char*)*actual_number_lines);
	for(i = 0; i < actual_number_lines; i++){
		input_data[i] = malloc(sizeof(char)*line_width);
	}
	
	buffer_size = (sizeof(char*)*actual_number_lines)*sizeof(char)*line_width;
	fprintf(stderr, "buffer_size is %d\n", buffer_size);

	temp_buffer = malloc(sizeof(char*)*actual_number_lines);
	for(i = 0; i < actual_number_lines; i++){
		temp_buffer[i] = malloc(sizeof(char)*line_width);
	}
		

	do{
		// read line_width chars or up to EOF or EOL
		// if read EOF then is_end == NULL
		is_end = fgets(temp_line, line_width, fptr);
		
		if(strncmp(temp_line, "#", 1) != 0){
			// not a comment so set the input_data part
			memcpy(input_data[line_count], temp_line, line_width);
			line_count++;
		} else {
			fprintf(stderr, "comment!\n");
		}
			

		if(line_count > actual_number_lines-1){
			// i.e next read will drop us off the world
			fprintf(stderr, "allocating more space!\n");
			copy_char_arrays(temp_buffer, input_data, line_width, actual_number_lines);			
			// free the old space, this is a bit tricky since we are trying to use structured data
			free_char_array(input_data, line_width, actual_number_lines);		
			
			previous_number_lines = actual_number_lines;

			actual_number_lines = actual_number_lines + init_number_lines; // grow the size
			buffer_size = (sizeof(char*)*actual_number_lines)*sizeof(char)*line_width;
			// reallocate the buffer
			input_data = malloc(sizeof(char*)*actual_number_lines);
			for(i = 0; i < actual_number_lines; i++){
				input_data[i] = malloc(sizeof(char)*line_width);
			}
			
			// copy the data back in 
			copy_char_arrays(input_data, temp_buffer, line_width, previous_number_lines);
			// finally we have to free and realloc the temp buffer
			free_char_array(temp_buffer, line_width, previous_number_lines);
			// and allocate it again
			temp_buffer = malloc(sizeof(char*)*actual_number_lines);
			for(i = 0; i < actual_number_lines; i++){
				temp_buffer[i] = malloc(sizeof(char)*line_width);
			}
						
		}
	  
	} while (is_end != NULL);
	line_count--; // (reading EOF overcounts by one)

	fprintf(stderr, "read %d\n", line_count);
	free_char_array(temp_buffer, line_width, actual_number_lines);

	// realloc temp to be just big enough
	temp_buffer = malloc(sizeof(char*)*line_count);
	for(i = 0; i < line_count; i++){
		temp_buffer[i]  = malloc(sizeof(char)*line_width);
	}
	// copy in the final data
	copy_char_arrays(temp_buffer, input_data, line_width, line_count); 

	free_char_array(input_data, line_width, actual_number_lines);
	fclose(fptr);
	*line_count_final = line_count;
	return(temp_buffer);
}


//! copy char_arrays of length lx,ly
void copy_char_arrays(char** dest, char** src, int lx, int ly){
	int i;
	for(i = 0; i < ly; i++){
		memcpy(dest[i], src[i], lx);
	}
}

//! free a 2d array
void free_char_array(char** array, int lx, int ly){
	int i;
	for(i = 0; i < ly;i++){
		free(array[i]);
	}
	free(array);
}
