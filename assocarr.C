/*************************************************

  assocarr.C - Associative arrays

  Andrew Y. Ng, 1996

**************************************************/

#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "assocarr.h"
#include "str.h"
#include "misc.h"

// int:
int hashit(int x, int hashsize)
  {
  return x % hashsize;
  }

int data_eq(int x, int y)
  {
  return x == y;
  }

// double:
int hashit(double d, int hashsize)
  {
  int hash, dex;
  unsigned char *p = (unsigned char*)&d;

  hash = 0;
  for (dex=0; dex < sizeof(double); dex++)
        hash = ((hash<<8)+p[dex]) % hashsize;

  return hash;
  }

int data_eq(double d1, double d2)
  {
  return (d1 == d2);
  }

// const char *:
int hashit(const char *str, int hashsize)
  {
  const unsigned char *p = (const unsigned char*)str;
  int hash, dex;

  hash = 0;
  for (dex=0; p[dex] != 0; dex++)
        hash = ((hash<<8)+p[dex]) % hashsize;

  return hash;
  }

int data_eq(const char *str1, const char *str2)
  {
  return !strcmp(str1, str2);
  }

// String:
int hashit(const String &str, int hashsize)
  {
  return hashit(str.as_char(), hashsize);
  }

int data_eq(const String &str1, const String &str2)
  {
  return (str1 == str2);
  }

//---------------------------------------------------

/*

// test code:
// -1 x     to print out x-th element and arrsize, 
// -2 ?     to print out all elements
// k x      (k != -1,-2) to set x-th element to k
int main(void)
  {
//  AssocArray<int, String> aa(0);
  WordVector aa;
  int *ip, ctr;
  String *sp;
  char buff[4096];
  const char *ptr;
  String str;
  int a, b;

  for (;;)
        {
        scanf("%d%s", &a, buff);
        str = buff;
        if (a == -1)
                printf("aa[%s] = %d (%d,%d)\n", str.as_char(), aa[str], aa.ndifferent_words(), aa.nwords());
//      else if (a == -2)
//              {
//              sp = aa.new_all_keys();
//              ip = aa.new_all_values();
//              for (ctr=0; ctr < aa.nele_val(); ctr++)
//                      printf("%s(%d) ", sp[ctr].as_char(), ip[ctr]);
//              printf("\n");
//              }
        else
                {
                printf("setting aa[%s] = %d\n", str.as_char(), a);
                aa[str] = a;
                }
        }

  // for integral keys:
//  AssocArray<int> aa(0);
//  for (;;)
//      {
//      scanf("%d%d", &a, &b);
//      if (a == -1)
//              printf("aa[%d] = %d\n", b, aa[b]);
//      else
//              {
//              printf("setting aa[%d] = %d\n", b, a);
//              aa[b] = a;
//              }
//      }

  }

*/

