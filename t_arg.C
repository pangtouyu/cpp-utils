#include <stdio.h>
#include <stdarg.h>

void foo(char *name, ...)
  {
  va_list ap;
  char *p;

  printf("%s: ", name);

  va_start(ap, name);
  do 
	{
	p = va_arg(ap, char*);

	if (p != NULL)
		printf("%s ", p);
	}
  while (p != NULL);
  printf("\n");

  va_end(ap);

  return;
  }

int main(void)
  {
  foo("hello", "bob", "alice", "amy", "peter", 0);

  return 0;
  }

