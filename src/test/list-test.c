#include "stdio.h"
#include "stdlib.h"

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


list *create_list(void){
	list* l;
	l = get_node();
	l->next = NULL;
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

	printf("top->%d %d\n", top_region(test).region_start, top_region(test).region_stop);

	walknprint(test, lptr);
	
	printf("found important: 0x%x %d %d\n", lptr, lptr->the_region.region_start, lptr->the_region.region_stop);



	for(i = 0; i < 10; i++){
		temp_region = pop_region(test);
		printf("pop->%d %d\n", temp_region.region_start, temp_region.region_stop);
	}

	remove_list(test);

	return(0);
}
	
	
	
