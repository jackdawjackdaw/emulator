#include "stdio.h"
#include "stdlib.h"

/**
 * have another look at search.h
 * this implements a binary search tree, 
 * and includes fns for building, stripping down and various traversal algs
 * the breakdown of regions can of course be rep'd as a bin search tree 
 * sorted on start of region
 *
 * root -> region1 (0,0)
 *      |
 *      -> node1(0.1, 1) -> region2 (0.1, 0.5)
 *                       |
 *                       -> node2(0.5, 1) -> region3(0.5, 0.8)
 *                       |
 *                       -> region4(0.8, 1.0)
 *                       
 * 
 * This is not quite as easy to deal with an unsorted linked list
 * head-> region1 -> region3 -> region2 -> region4 -> end
 * 
 * But it's implicitly sorted by left-rightness so it makes postprocessing easier
 * otherwise we have to sort the linked list , or more simply transfer the final
 * list into a fixed size array and then qsort on that. 
 *
 */

//! for a linked list of regions
typedef struct region{
	int region_start;
	int region_stop;
	double emu_x_start;
	double emu_x_stop;
} region;


typedef struct list{
	region the_region;
	struct list* next;
	struct list* prev;
} list;


list *create_list(void);
int empty_list(list* the_list);
void push_region( region x, list* the_list);
region pop_region(list* the_list);
region top_region(list* the_list);
void remove_list(list* the_list);
list* get_node();
void return_node(list* the_node);
void append_region(region x, list* the_list);

list *create_list(void){
	list* l;
	l = get_node();
	l->next = NULL;
	l->prev = NULL;
	return(l);
}

int empty_list(list* the_list){
	return(the_list->next == NULL);
}

// add a region to the head
void push_region( region x, list* the_list){
	list* temp;
	temp = get_node();
	temp->the_region = x;
	temp->next = the_list->next;
	the_list->next = temp;
}

// add to the tail
void append_region(region x, list* the_list){
	list* temp;
	list* tailval;
	temp = the_list->next;
	while (temp->next != NULL){
		temp = temp->next;
	} 
 	// now it's the last one
	tailval = get_node();
	tailval->next = NULL;
	tailval->the_region = x;
	temp->next = tailval;
}
	
// remove from tail
region tail_pull(list* the_list){
	list* temp;
	region temp_region;
	temp = the_list->next;
	// while you're not the one before the end
	while (temp->next->next != NULL){
		temp = temp->next;
	}
	// get the last region
	temp_region = temp->next->the_region;
	// tell this penulatimate region it's now the end

	// ditch the end
	return_node(temp->next->next);
	temp->next = NULL;
	return(temp_region);
}

int get_length(list* the_list){
	list *temp;
	int count = 0;
	temp = the_list;
	while (temp->next != NULL){
		temp = temp->next;
		count++;
	}
	return(count);
}


region pop_region(list* the_list){
	list* temp;
	region temp_region;
	temp = the_list->next;
	the_list->next = temp->next;
	temp_region = temp->the_region;
	return_node(temp);
	return(temp_region);
}

// what's on top?
region top_region(list* the_list){
	return(the_list->next->the_region);
}

void remove_list(list* the_list){
	list* temp;
	do{
		temp = the_list->next;
		return_node(the_list);
		the_list = temp;
	} while (temp != NULL);
}


list* get_node(){
	list* temp;
	temp = malloc(sizeof(list));
	return(temp);
}

void return_node(list* the_node){
	free(the_node);
}


void walknprint(list* the_list, list* important){
	list* temp = the_list->next;
	while(temp != NULL){
		printf("%d %d\n", temp->the_region.region_start, temp->the_region.region_stop);
		if(temp->the_region.region_start == 2 && temp->the_region.region_stop == 3){
			important = temp;
			printf("0x%x\n", important);
		}
		temp = temp->next;
	}
}

int main (void){
	int i;
	region temp_region;
	list* test  = create_list();
	list* lptr;

	for(i = 0; i < 10; i++){
		temp_region.region_start = i;
		temp_region.region_stop = i+1;	
		push_region(temp_region, test);
	}
	
	temp_region.region_start = 100;
	temp_region.region_stop = 200;
	append_region(temp_region, test);

	//printf("top->%d %d\n", top_region(test).region_start, top_region(test).region_stop);

	walknprint(test, lptr);

 	printf("list is %d long\n", get_length(test));
	
	//printf("found important: 0x%x %d %d\n", lptr, lptr->the_region.region_start, lptr->the_region.region_stop);

	temp_region = tail_pull(test);
	printf("pulled tail: %d %d\n", temp_region.region_start, temp_region.region_stop);
	

	for(i = 0; i < 10; i++){
		temp_region = pop_region(test);
		printf("pop->%d %d\n", temp_region.region_start, temp_region.region_stop);
	}

	remove_list(test);

	return(0);
}
	
	
	
