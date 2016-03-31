#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "heap.h"
#include "qgeneral.h"

#define GROWTH_FACTOR (1.5)

// a bit dangerous, and breaks ADT framework, but useful
const heap_t *Heap::all_elements(int *n_ele)
  {
  if (n_ele != NULL)
	*n_ele = size;

  return heap;
  }

heap_ele_t Heap::query_min(void)
  {
  if (is_empty())
	heap_error("Somebody tried to query_min while empty.");
	
  return heap[0].val;
  }

heap_key_t Heap::query_min_key(void)
  {
  if (is_empty())
	heap_error("Somebody tried to query_min while empty.");
	
  return heap[0].key;
  }

int Heap::is_empty(void)
  {
  return (size==0);
  }

int Heap::size_val(void)
  {
  return size;
  }

void Heap::empty_heap(void) 
  {
  size=0;

  return;
  }

int Heap::left_of(int node) 
  {
  return 2*node+1;
  }

int Heap::right_of(int node) 
  {
  return 2*node+2;
  }

int Heap::parent_of(int node)
  {
  return (node-1)/2;
  }

void Heap::heap_error(char *msg)
  {
  fflush(stdout);

  fprintf(stderr, "Heap ERROR: %s\n", msg);
  fflush(stderr);

//  abort();
  exit(-1);
  }

void Heap::down_heap(int node)
  {
  int l, r, smallest; 
  heap_t temp;

  l=left_of(node);
  r=right_of(node);

  smallest = node;
  if (l < size && heap[l].key < heap[node].key)
	smallest = l;
  if (r < size && heap[r].key < heap[smallest].key)
	smallest = r;

  if (smallest != node)
	{
	temp = heap[smallest];
	heap[smallest] = heap[node];
	heap[node] = temp;

	down_heap(smallest);
	}

  return;
  }

void Heap::expand_heap(void)
  {
  int new_heap_size;
  int ctr;
  heap_t *h;

  new_heap_size = (int)ceil(curr_max_size * GROWTH_FACTOR);
  if (new_heap_size == curr_max_size)
	new_heap_size++;

//  printf("Heap: expanding from %d to %d \n", curr_max_size, new_heap_size);

  h = (heap_t*)malloc(sizeof(heap_t)*new_heap_size);
  if (h == NULL)
	heap_error("Out of memory while trying to expand heap.");

  for (ctr=0; ctr < curr_max_size; ctr++)
	h[ctr] = heap[ctr];

  free(heap);

  heap = h;
  curr_max_size = new_heap_size;

  return;
  }

void Heap::insert(heap_ele_t new_ele, heap_key_t key)
  {
  int node;
  int ctr;

  if (size >= curr_max_size)
	expand_heap();
  assert(size < curr_max_size);

  node = size;
  size++; 
  while (node > 0 && heap[parent_of(node)].key > key)
	{
	heap[node] = heap[parent_of(node)];
	node = parent_of(node);
	}

  heap[node].key = key;
  heap[node].val = new_ele;

  /* For debugging, and doesn't work due to some earliere changes. */
/*
  printf("Heap: ");
  for (ctr=0; ctr < size; ctr++)
	printf("%f ", heap[ctr]);
  printf("\n");
*/

  return;
  }

heap_ele_t Heap::delete_min(void)
  {
  heap_ele_t min;

  if (is_empty())
	heap_error("Somebody tried to delete while empty.");

  min = query_min();

  size--;
  heap[0] = heap[size];

  down_heap(0);

  return min;
  }

Heap::Heap(int init_max_size)
  {
  curr_max_size = init_max_size;
  size=0;

  heap = (heap_t*) malloc(sizeof(heap_t) * curr_max_size);

  if (heap == NULL)
	heap_error("Out of memory. Could not create heap.");

  return;
  }

Heap::~Heap()
  {
  if (heap != NULL)
	free(heap);

  heap = NULL;

  return;
  }

