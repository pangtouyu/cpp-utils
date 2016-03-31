/*************************************************

  assocarr.h - Associative arrays

  Andrew Y. Ng, 1996

**************************************************/

#ifndef _ASSOCARR_H
#define _ASSOCARR_H

#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "str.h"
#include "misc.h"


// If aa is an AssocArray<T,U>, then space is allocated each
// time aa[u] is referenced, for a new u. To "see" the value
// of aa[u] but without allocating new memory, use aa.query_val(u)
// instead.
// 
// An AssocArray<T,U> Maps type U to type T. For each 
// U, we need a int hashit(U,int) that hashes U and returns 
// an integer from 0 to hashsize-1, and a int data_eq(U,U) 
// that returns true iff the 2 arguments it is passed are to
// be considered equal. Also, we assume that "=" works
// on T's for assignment.

int hashit(int x, int hashsize);
int data_eq(int x, int y);
int hashit(double x, int hashsize);
int data_eq(double x, double y);
int hashit(const char *str, int hashsize);
int data_eq(const char *str1, const char *str2);
int hashit(const String &str, int hashsize);
int data_eq(const String &str1, const String &str2);

// data structure used by AssocArray internally to keep
// track of what's mapped to what.
template <class T, class U>
struct aa_ref_t 
  {
  U index;
  T dat;
  };

template <class T, class U>
class AssocArray
  {
  private:
        const double GROWTH_FACTOR;	// initialized in constructor
        const int INIT_SIZE;		// initialized in constructor
        aa_ref_t<T,U> **arr;    // hashtable of size arrsize.
        T default_value;
        int arrsize, nele;

        inline void expand_array(void);
        inline T &subscript_helper(U &index);           // like [], but never expands array

  public:
        inline int nele_val(void) const { return nele; }
	inline T defaultValue_val(void) const {return default_value; }
        inline T &operator [](const U &index);
        inline AssocArray(T default_value_);
        inline ~AssocArray();
        inline int has_key(const U &index) const;
        inline U *new_all_keys(void) const;
        inline T *new_all_values(void) const;
        inline T query_val(const U &index) const;
	inline void forgetAll(int shrinkArray = 0);
  };

template <class T>
class HashedSet
  {
  private:
        const double GROWTH_FACTOR;	// initialized in constructor
        const int INIT_SIZE;		// initialized in constructor
        T **arr;
        int arrsize, nele;

        inline void expand_array(void);

  public:
        inline int nele_val(void) const { return nele; }
        inline void insertElement(const T &ele);	// error if ele is already in the set
        inline void replaceElement(const T &ele);	// okay whether ele is currently there
							// or not
        inline int containsElement(const T &ele);
        inline HashedSet(void);
        inline ~HashedSet();
        inline T *new_allElements(void) const;
  };

//-------------------------------------------------------
// HashedSet functions

template <class T, class U>
inline int AssocArray<T,U>::has_key(const U &index) const
  {
  int ctr;

  ctr = hashit(index, arrsize);
  while (arr[ctr] != NULL && !data_eq(arr[ctr]->index, index))
        {
        ctr++;
        if (ctr == arrsize)
                ctr = 0;
        }

  if (arr[ctr] == NULL)
        return 0;
  else
        return 1;
  }

template <class T, class U>
inline U *AssocArray<T,U>::new_all_keys(void) const 
  {
  U *all_keys;
  int ctr, dex;

  if (nele == 0)
	return (U*)NULL;

  all_keys = new U[nele];
  assert(all_keys != NULL);

  dex=0;
  for (ctr=0; ctr < arrsize; ctr++)
        {
        if (arr[ctr] == NULL)
                continue;
        all_keys[dex++] = arr[ctr]->index;
        assert(dex <= nele);
        }

  assert(dex == nele);

  return all_keys;
  }

template <class T, class U>
inline T *AssocArray<T,U>::new_all_values(void) const 
  {
  T *all_values;
  int ctr, dex;

  if (nele == 0)
	return NULL;

  all_values = new T[nele];
  assert(all_values != NULL);

  dex=0;
  for (ctr=0; ctr < arrsize; ctr++)
        {
        if (arr[ctr] == NULL)
                continue;
        all_values[dex++] = arr[ctr]->dat;
        assert(dex <= nele);
        }

  assert(dex == nele);

  return all_values;
  }

template <class T, class U>
inline AssocArray<T,U>::AssocArray(T default_value_) 
                : default_value(default_value_),
		   GROWTH_FACTOR(1.531415936535), INIT_SIZE(11)

  {
  int ctr;

  arr = new aa_ref_t<T,U>*[INIT_SIZE];
  for (ctr=0; ctr < INIT_SIZE; ctr++)
        arr[ctr] = (aa_ref_t<T,U>*)NULL;
  assert(arr != NULL);

  arrsize = INIT_SIZE;
  nele = 0;

  return;
  }

template <class T, class U>
inline AssocArray<T,U>::~AssocArray()
  {
  int ctr;

  for (ctr=0; ctr < arrsize; ctr++)
        if (arr[ctr] != NULL)
                delete arr[ctr];

  delete[] arr;
  }

template <class T, class U>
inline void AssocArray<T,U>::forgetAll(int shrinkArray)
  {
  for (int ctr=0; ctr < arrsize; ctr++)
	{
	if (arr[ctr] != NULL)
		{
		delete arr[ctr];
		arr[ctr] = (aa_ref_t<T,U>*)NULL;
		}
	}

  shrinkArray = 0;
  if (shrinkArray && arrsize > INIT_SIZE)
	{
	delete[] arr;
	arr = new aa_ref_t<T,U>*[INIT_SIZE];
	arrsize = INIT_SIZE;
	for (int ctr=0; ctr < arrsize; ctr++)
		arr[ctr] = (aa_ref_t<T,U>*)NULL;
	}

  nele = 0;
  };

template <class T, class U>
inline void AssocArray<T,U>::expand_array(void)
  {
  aa_ref_t<T,U> **old_arr;
  int old_arrsize, old_nele, ctr, dex;

  old_nele = nele;
  old_arr = arr;
  old_arrsize = arrsize;

  // allocate a new array 
  arrsize = (int)round(old_arrsize * GROWTH_FACTOR);
  if (is_even(arrsize))
        arrsize++;
  arr = new aa_ref_t<T,U>*[arrsize];
  assert(arr != NULL);
  for (ctr=0; ctr < arrsize; ctr++)
        arr[ctr] = (aa_ref_t<T,U>*)NULL;
  nele = 0;

  // rehash all the old elements
  for (ctr=0; ctr < old_arrsize; ctr++)
        {
        if (old_arr[ctr] == NULL)
                continue;

//      (*this)[old_arr[ctr]->index] = old_arr[ctr]->dat;
        dex = hashit(old_arr[ctr]->index, arrsize);
        while (arr[dex] != NULL)
                if (++dex == arrsize)
                        dex = 0;
        arr[dex] = old_arr[ctr];
        nele++;
        }

  assert(nele == old_nele);

  delete[] old_arr;

  return;
  }

// like [], but never expands array.
// Precond: index is a key that is already in the hash table.
template <class T, class U>
inline T &AssocArray<T,U>::subscript_helper(U &index)
  {
  int ctr;

  ctr = hashit(index, arrsize);
  while (arr[ctr] != NULL && !data_eq(arr[ctr]->index, index))
        {
        ctr++;
        if (ctr == arrsize)
                ctr = 0;
        }

//  if (arr[ctr] == NULL)
//      {
//      arr[ctr] = new aa_ref_t<T,U>;
//      arr[ctr]->index = index;
//      arr[ctr]->dat = default_value;
//      nele++;
//      }

  assert(arr[ctr] != NULL);

  return arr[ctr]->dat;
  }

template <class T, class U>
inline T &AssocArray<T,U>::operator [](const U &index)
  {
  int ctr;

  ctr = hashit(index, arrsize);
  while (arr[ctr] != NULL && !data_eq(arr[ctr]->index, index))
        {
        ctr++;
        if (ctr == arrsize)
                ctr = 0;
        }

  if (arr[ctr] == NULL)
        {
        arr[ctr] = new aa_ref_t<T,U>;
        arr[ctr]->index = index;
        arr[ctr]->dat = default_value;
        nele++;
        if (nele > arrsize/2)
                {
                expand_array();
                return (*this)[index];
                }
        else
                return arr[ctr]->dat;
        }

  return arr[ctr]->dat;
  }

template <class T, class U>
inline T AssocArray<T,U>::query_val(const U &index) const
  {
  int ctr;

  ctr = hashit(index, arrsize);
  while (arr[ctr] != NULL && !data_eq(arr[ctr]->index, index))
        {
        ctr++;
        if (ctr == arrsize)
                ctr = 0;
        }

  if (arr[ctr] == NULL)
        return default_value;
  else
        return arr[ctr]->dat;
  }

//-------------------------------------------------------
// HashedSet functions

template <class T>
inline T *HashedSet<T>::new_allElements(void) const 
  {
  T *allElements;
  int ctr, dex;

  if (nele == 0)
	return NULL;

  allElements = new T[nele];
  assert(allElements != NULL);

  dex=0;
  for (ctr=0; ctr < arrsize; ctr++)
        {
        if (arr[ctr] == NULL)
                continue;
        allElements[dex++] = *arr[ctr];
        }

  assert(dex == nele);

  return allElements;
  }

template <class T>
inline HashedSet<T>::HashedSet(void) : 
		   GROWTH_FACTOR(1.531415936535), INIT_SIZE(11)
  {
  int ctr;

  arr = new T*[INIT_SIZE];
  assert(arr != NULL);
  for (ctr=0; ctr < INIT_SIZE; ctr++)
        arr[ctr] = NULL;

  arrsize = INIT_SIZE;
  nele = 0;

  return;
  }

template <class T>
inline HashedSet<T>::~HashedSet()
  {
  int ctr;

  for (ctr=0; ctr < arrsize; ctr++)
        if (arr[ctr] != NULL)
                delete arr[ctr];

  delete[] arr;
  }

template <class T>
inline void HashedSet<T>::expand_array(void)
  {
  T **old_arr;
  int old_arrsize, old_nele, ctr, dex;

  old_nele = nele;
  old_arr = arr;
  old_arrsize = arrsize;

  // allocate a new array 
  arrsize = (int)round(old_arrsize * GROWTH_FACTOR);
  if (is_even(arrsize))
        arrsize++;
  arr = new T*[arrsize];
  assert(arr != NULL);
  for (ctr=0; ctr < arrsize; ctr++)
        arr[ctr] = NULL;
  nele = 0;

  // rehash all the old elements
  for (ctr=0; ctr < old_arrsize; ctr++)
        {
        if (old_arr[ctr] == NULL)
                continue;

        dex = hashit(*old_arr[ctr], arrsize);
        while (arr[dex] != NULL)
                if (++dex == arrsize)
                        dex = 0;
        arr[dex] = old_arr[ctr];
        nele++;
        }

  assert(nele == old_nele);

  delete[] old_arr;

  return;
  }

template <class T>
inline void HashedSet<T>::replaceElement(const T &ele)
  {
  T **ptr;
  int ctr;

  ctr = hashit(ele, arrsize);
  while (arr[ctr] != NULL && !data_eq(*arr[ctr], ele))
        {
        ctr++;
        if (ctr == arrsize)
                ctr = 0;
        }

  if (arr[ctr] != NULL)
	{
	// ele is already in the set
	return;
	}
  else
        {
        arr[ctr] = new T(ele);
        nele++;
        if (nele > arrsize/2)
                expand_array();
        }

  return;
  }

template <class T>
inline void HashedSet<T>::insertElement(const T &ele)
  {
  T **ptr;
  int ctr;

  ctr = hashit(ele, arrsize);
  while (arr[ctr] != NULL && !data_eq(*arr[ctr], ele))
        {
        ctr++;
        if (ctr == arrsize)
                ctr = 0;
        }

  if (arr[ctr] != NULL)
	{
	fprintf(stderr,"ERROR: tried to insert duplicate element.\n");
	abort();
	}
  else
        {
        arr[ctr] = new T(ele);
        nele++;
        if (nele > arrsize/2)
                expand_array();
        }

  return;
  }

template <class T>
inline int HashedSet<T>::containsElement(const T &ele)
  {
  int ctr;

  ctr = hashit(ele, arrsize);
  while (arr[ctr] != NULL && !data_eq(*arr[ctr], ele))
        {
        ctr++;
        if (ctr == arrsize)
                ctr = 0;
        }

  return (arr[ctr] == NULL)?0:1;
  }

#endif          // _ASSOCARR_H

