#ifndef _SPLAY_H
#define _SPLAY_H

#include "sheap.h"

typedef sheap_key_t splaytree_key_t;
typedef sheap_ele_t splaytree_data_t;

#define SPLAY_ALLOW_SPECIFIC_DELETES	// if allowed, then we need an
					// == operator for splaytree_data_t

struct splaytree_node {
    splaytree_node * left, * right;
    splaytree_key_t key;
    splaytree_data_t val;
};

class SplayTree
  {
  private:
        splaytree_node *root;
        int size;
	int allow_identical_keys;

	splaytree_node *splay (splaytree_key_t key, splaytree_node *t) const;
	void splaytree_error(const char *msg);
	int delete_node_helper(splaytree_key_t key, splaytree_data_t val,
			splaytree_node **t);

  public:
        int insert(splaytree_key_t key, splaytree_data_t val);
        int delete_node(splaytree_key_t key);
	int size_val(void) {return size;}
	int is_empty(void) {return (root == NULL);}
	int delete_node(splaytree_key_t key, splaytree_data_t val);
	splaytree_key_t query_min_key(void);
	splaytree_data_t query_min(void);
	splaytree_data_t delete_min(void);
	SplayTree(int allow_identical_keys_);
	SplayTree(void);
	~SplayTree(void);
  };

#endif 		// _SPLAY_H

