#include <stdio.h>
#include <math.h>

#include "sheap.h"

int main(void)
  {
  SHeap h;
  foo_t foo;
  char cmd;
  double d;

  printf("is_empty() = %d\n", h.is_empty());
  for (;;)
	{
	while ((cmd=getchar()) == '\n' || cmd == ' ');
	if (cmd == 'i')
		{
		scanf("%lf%s", &d, foo.str);
		foo.num = (int)floor(d*100);
		h.insert(foo, d);
		}
	else if (cmd == 'D')
		{
		scanf("%lf%s", &d, foo.str);
		foo.num = (int)floor(d*100);
		h.delete_ele(foo, d);
		}
	else if (cmd == 'q')
		{
		foo = h.query_min();
		printf("%s (%d)\n", foo.str, foo.num);
		}
	else if (cmd == 'd')
		{
		foo = h.delete_min();
		printf("%s (%d)\n", foo.str, foo.num);
		}
	else if (cmd == 'e')
		h.empty_heap();
	else 
		printf("Unknown command. (Should be i,q,D,d or e)\n");
	printf("is_empty() = %d\n", h.is_empty());
	}

  }	  	
/*
struct foo_t
  {
  char str[64];
  int num;
  };

  public:
	heap_ele_t Heap::query_min(void);
	int Heap::is_empty(void);
	void Heap::empty_heap(void);
	void Heap::heap_insert(heap_ele_t new_ele, heap_key_t key);
	heap_ele_t Heap::heap_delete_min(void);
	Heap::Heap(int init_max_size);
	Heap::~Heap();

*/
