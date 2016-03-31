#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "splay.h"
#include "misc.h"

void SplayTree::splaytree_error(const char *msg)
  {
  fflush(stdin);
  fprintf(stderr,"\n\nSplayTree ERROR: %s\n\n", msg);
  fflush(stderr);

  exit(-1);
  }

// Simple top down splay, not requiring i to be in the tree t.  
splaytree_node *SplayTree::splay(splaytree_key_t key, splaytree_node *t) const {
    splaytree_node N, *l, *r, *y;
    if (t == NULL) return t;
    N.left = N.right = NULL;
    l = r = &N;

    for (;;) {
        if (key < t->key) {
            if (t->left == NULL) break;
            if (key < t->left->key) {
                y = t->left;                           // rotate right 
                t->left = y->right;
                y->right = t;
                t = y;
                if (t->left == NULL) break;
            }
            r->left = t;                               // link right 
            r = t;
            t = t->left;
        } else if (key > t->key) {
            if (t->right == NULL) break;
            if (key > t->right->key) {
                y = t->right;                          // rotate left 
                t->right = y->left;
                y->left = t;
                t = y;
                if (t->right == NULL) break;
            }
            l->right = t;                              // link left 
            l = t;
            t = t->right;
        } else {
            break;
        }
    }
    l->right = t->left;                                // assemble 
    r->left = t->right;
    t->left = N.right;
    t->right = N.left;
    return t;
}

// Insert node the tree, unless it's already there.  
// returns 1 if FAILED to insert. (i.e. if key already there,
// and not allowing duplicate keys)
int SplayTree::insert(splaytree_key_t key, splaytree_data_t val) {
    splaytree_node *new_node;

    new_node = (splaytree_node *)safe_malloc(sizeof (splaytree_node));
    new_node->key = key;
    new_node->val = val;

    if (root == NULL) {
        new_node->left = new_node->right = NULL;
	assert(size == 0);
        size = 1;
        root = new_node;
	return 0;
    }
    root = splay(key,root);
    if (key < root->key) {
        new_node->left = root->left;
        new_node->right = root;
        root->left = NULL;
        size ++;
        root = new_node;
	return 0;
    } else if (key > root->key || allow_identical_keys) {
        new_node->right = root->right;
        new_node->left = root;
        root->right = NULL;
        size++;
        root = new_node;
	return 0;
   } else { // already in the tree, and not allowing identical keys.
            // So don't add it again.
        warn_once("Splay tree: didn't insert since identical key.");
        free(new_node);
        return 1;
   }
}

#ifdef SPLAY_ALLOW_SPECIFIC_DELETES
// Deletes node with specified key from the tree if it's there.
// returns 1 if FAILED to delete
int SplayTree::delete_node(splaytree_key_t key, splaytree_data_t val) 
  {
  return delete_node_helper(key, val, &root);
  }

// t is a pointer to a pointer to a tree. This will try to
// delete the specified element, and modify t if necessary.
int SplayTree::delete_node_helper(splaytree_key_t key, splaytree_data_t val,
			splaytree_node **t) 
  {
  splaytree_node *x;

  if (*t==NULL) 
	return 1;

  *t = splay(key,*t);
  if (key == (*t)->key) 
	{
	// found the key. Check if data's correct.
	// (note we do !(..==..) rather than ..!=.. so that we
	//  only require == to be defined.)
	if (!(val == (*t)->val))
		{
		// data's NOT the same.
		// try delete-left
		if (!delete_node_helper(key, val, &(*t)->left))
			return 0;
		// try delete-right
		if (!delete_node_helper(key, val, &(*t)->right))
			return 0;
		return 1;	// couldn't find it
		}

	// if we get here, that means the data was the same.
	// So, go ahead and remove it!
        if ((*t)->left == NULL) 
	    {
            x = (*t)->right;
            } 
	else 
	    {
	    x = splay(key, (*t)->left);
	    x->right = (*t)->right;
            }
	size--;
	free(*t);
	*t = x;
	return 0;
	}

  return 1;                          // It wasn't there 
  }
#endif 		// SPLAY_ALLOW_SPECIFIC_DELETES

// Deletes node with specified key from the tree if it's there.
// returns 1 if FAILED to delete
int SplayTree::delete_node(splaytree_key_t key) {
    splaytree_node *x;

    if (root==NULL) 
	return 1;

    root = splay(key,root);
    if (key == root->key) {               // found it 
        if (root->left == NULL) {
            x = root->right;
        } else {
            x = splay(key, root->left);
            x->right = root->right;
        }
        size--;
        free(root);
        root = x;
	return 0;
    }
    return 1;                          // It wasn't there 
}

splaytree_key_t SplayTree::query_min_key(void)
  {
  splaytree_node *t;

  if (root == NULL)
	splaytree_error("Tried to query_min_key() on empty tree.");

  for (t=root; t->left != NULL; t = t->left);

  root = splay(t->key, root);		// for logarithmic performance.
					// will splay it to the root.

  return root->key;
  }

splaytree_data_t SplayTree::query_min(void)
  {
  splaytree_node *t;

  if (root == NULL)
	splaytree_error("Tried to query_min() on empty tree.");

  for (t=root; t->left != NULL; t = t->left);

  root = splay(t->key, root);		// for logarithmic performance.
					// will splay it to the root.

  return root->val;
  }

splaytree_data_t SplayTree::delete_min(void)
  {
  splaytree_key_t minkey;
  splaytree_data_t retval;

  minkey = query_min_key();		// should splay it to the root
  assert(root->key == minkey);

  retval = root->val;

  delete_node(minkey);

  return retval;
  }

SplayTree::SplayTree(int allow_identical_keys_)
  {
  allow_identical_keys = allow_identical_keys_;
  root = NULL; 
  size=0;

  return;
  }

SplayTree::SplayTree(void)
  {
  allow_identical_keys = 1;
  root = NULL; 
  size=0;

  return;
  }

SplayTree::~SplayTree()
  {
  assert(size >= 0);

  while (size > 0)
	delete_node(root->key);

  return;
  }

/* Old code used to exercise the splay tree
// A sample use of these functions.  Start with the empty tree,         
// insert some stuff into it, and then delete it                       
int main(void) {
    SplayTree root;
    int i;
    foo_t foo;

    for (i = 0; i < 1024; i++) {
	foo.num = (541*i) & (1023);
        root.insert((541*i) & (1023), foo);
    }
    for (i = 0; i < 1000; i++) {
        root.delete_node((541*i) & (1023));
    }
    printf("size = %d\n", root.size_val());
}
*/
