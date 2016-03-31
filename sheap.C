#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sheap.h"
#include "splay.h"

int SHeap::is_empty(void)
  {
  if (tree->size_val() == 0)
	{
	assert(tree->is_empty());
	return 1;
	}
  else
	{
	assert(!tree->is_empty());
	return 0;
	}
  }

sheap_ele_t SHeap::query_min(void)
  {
  if (is_empty())
	heap_error("Somebody tried to query_min while empty.");
	
  return tree->query_min();
  }

sheap_ele_t SHeap::delete_min(void)
  {
  if (is_empty())
	heap_error("Somebody tried to delete_min while empty.");
	
  return tree->delete_min();
  }

void SHeap::empty_heap(void)
  {
  while (!is_empty()) 
	tree->delete_min();

  return;
  }

void SHeap::insert(sheap_ele_t new_ele, sheap_key_t key)
  {
  int err;
  
  err=tree->insert(key, new_ele);

  if (err == 1)
	heap_error("insert: Failed to insert.");

  return;
  }

void SHeap::delete_ele(sheap_ele_t ele, sheap_key_t key)
  {
  int err;

  err = tree->delete_node(key, ele);

  if (err == 1)
	heap_error("delete_ele() failed. (Not found).");

  return;
  }

void SHeap::heap_error(char *msg)
  {
  fflush(stdout);

  fprintf(stderr, "Heap ERROR: %s\n", msg);
  fflush(stderr);

  exit(-1);
  }

SHeap::SHeap(void)
  {
  tree = new SplayTree;

  if (tree == NULL)
	heap_error("Out of memory. Could not create heap.");

  return;
  }

SHeap::~SHeap()
  {
  assert(tree != NULL);
  delete tree;

  return;
  }

