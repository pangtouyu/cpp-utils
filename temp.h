template <class T>
class HashedSet
  {
  private:
        const double GROWTH_FACTOR = 1.531415936535;
        const int INIT_SIZE = 11;
        T **arr;
        int arrsize, nele;

        inline void expand_array(void);

  public:
        inline int nele_val(void) const { return nele; }
        inline void insertElement(const T &ele);
        inline int containsElement(const T &ele);
        inline HashedSet(void);
        inline ~HashedSet();
        inline T *new_allElements(void) const;
  };

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
        allElements[dex++] = arr[ctr];
        assert(dex <= nele);
        }

  assert(dex == nele);

  return allElements;
  }

template <class T>
inline HashedSet<T>::HashedSet(void)
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

        dex = hashit(old_arr[ctr], arrsize);
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
inline void HasedSet<T>::insertElement(const T &ele)
  {
  T **ptr;
  int ctr;

  ctr = hashit(ele, arrsize);
  while (arr[ctr] != NULL && !data_eq(arr[ctr], ele))
        {
        ctr++;
        if (ctr == arrsize)
                ctr = 0;
        }

  if (arr[ctr] != NULL)
	error("HashedSet: tried to insert duplicate element.");
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
inline void HasedSet<T>::containsElemenent(const T &ele)
  {
  int ctr;

  ctr = hashit(ele, arrsize);
  while (arr[ctr] != NULL && !data_eq(arr[ctr], ele))
        {
        ctr++;
        if (ctr == arrsize)
                ctr = 0;
        }

  return (arr[ctr] == NULL)?0:1;
  }


