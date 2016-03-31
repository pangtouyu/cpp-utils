#ifndef __HEAP_H
#define __HEAP_H

#include "qgeneral.h"

#define DEFAULT_HEAP_INIT_SIZE (128)

typedef double heap_key_t;	// MUST be char, int, float or double, 
				// or else some other type with overloaded
				// >, <, and == operators provided

typedef customer_t heap_ele_t;	// the data stored in the heap

struct heap_t	// this data structure is used by the heap routines,
		// and is returned by Heap::all_elements(int *)
  {
  heap_ele_t val;
  heap_key_t key;
  };

class Heap
  {
  private:
	int size, curr_max_size;
	heap_t *heap;

	int left_of(int node);
	int right_of(int node);
	int parent_of(int node);
	void heap_error(char *msg);
	void down_heap(int node);
	void expand_heap(void);

  public:
	const heap_t *all_elements(int *n_ele=NULL);	
			// a bit dangerous, and breaks ADT framework, but useful
	heap_key_t query_min_key(void);
	heap_ele_t query_min(void);
	int is_empty(void);
	int size_val(void);
	void empty_heap(void);
	void insert(heap_ele_t new_ele, heap_key_t key);
	heap_ele_t delete_min(void);
	Heap(int init_max_size=DEFAULT_HEAP_INIT_SIZE);
	~Heap();

  };

#endif	// __HEAP_H

