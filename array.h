#ifndef _ARRAY_H
#define _ARRAY_H

#include <assert.h>
#include <string.h>
#include "misc.h"

template<class T>
class Array
  {
  private:
	T *val;
	int n;

  public:
    inline int dim(void) const {return n;}
    inline int n_ele(void) const {return n;}
    inline int length(void) const {return n;}
    inline operator T*() const {return val;}
    inline T *ref(int dex) {return &val[dex];}
    inline void resize(int newSize) // destroys previous data
    { 
        if (n == newSize) return; 
        
        delete[] val; 
        val = new T[newSize]; 
        assert(val != NULL); 
        n=newSize; 
    }
    
    inline void setAllTo(const T& v)
    { 
        for (int dex=0; dex < n; dex++) val[dex] = v; 
    }
    inline Array &operator=(const Array<T> &a) 
    { 
        if (n != a.n) 
            error("ERROR: Array: assignment between arrays of non-equal length\n"); 
        for (int dex=0; dex < n; dex++)
            val[dex] = a.val[dex];
        return *this;
    }
            inline T*allEle(void) { return val; }

            inline Array(int size) : n(size) { assert(size>0); val = new T[size]; assert(val != NULL); 
                bzero(val, size*sizeof(T)); }
                inline Array(const Array<T> &a) { n = a.n; val = new T[n]; assert(val != NULL);
                    for (int dex=0; dex < n; dex++) val[dex] = a.val[dex]; }
                    inline ~Array() { delete[] val; }

                    inline void assert_dim(int n_)
                    { if (n_ != n) { fprintf(stderr, "Array::assert_dim failed. Expected %d, was %d.\n", n_, n); abort(); }}
  };

#endif		// _ARRAY_H
