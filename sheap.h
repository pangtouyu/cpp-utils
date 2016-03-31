#ifndef _SHEAP_H
#define _SHEAP_H

// sheap.[Ch]: Splay-tree implementation of heaps.

#include "event.h"
#include "qgeneral.h"

typedef double sheap_key_t;	// MUST be char, int, float or double, 
				// or else some other type with overloaded
				// >, <, and == operators provided
typedef event_t sheap_ele_t;	// the data stored in the heap

// for splay.[Ch], which use this header
typedef double splaytree_key_t;
typedef sheap_ele_t splaytree_data_t;

class SplayTree;		// defined in splay.[Ch]

class SHeap
  {
  private:
	SplayTree *tree;

	void heap_error(char *msg);

  public:
	sheap_ele_t query_min(void);
	int is_empty(void);
	void empty_heap(void); 
	void insert(sheap_ele_t new_ele, sheap_key_t key);
	void delete_ele(sheap_ele_t ele, sheap_key_t key);
	sheap_ele_t delete_min(void);
	SHeap(void);
	~SHeap();
  };

#endif	// _SHEAP_H

