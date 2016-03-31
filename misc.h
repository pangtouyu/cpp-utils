/*************************************************

  misc.h - Header file for misc.C which contains
    miscellaneous, self-intuitive functions.

  Andrew Y. Ng, 1994-96

**************************************************/

#ifndef _MISC_H
#define _MISC_H
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <bits/nan.h>
#include "Random_Number.h"

void bzero(void *s, int n);

inline double nanVal(void) { return NAN; }

void error(const char*msg);
void internal_error(const char *msg);
void prob_internal_error(const char *msg);
void warn(const char*msg);
void internal_warn(const char*msg);
void warn_once(const char *msg);
void not_hang(void);

void *safe_malloc(size_t size);
void *safe_calloc(int n_elems, size_t size);
void my_new_handler(void);
float **alloc_float_arr(int m, int n);
double **alloc_double_arr(int m, int n);
void free_float_arr(float **p, int m);
void free_double_arr(double **p, int m);

FILE *safe_fopen(const char *fn, const char *mode);
int exist_file(const char *fn);
void verify_exist_file(const char *fn);
void verify_not_exist_file(const char *fn);
double *new_doubles_from_file(FILE *fp, int *n_ele);
double *new_doubles_from_file(const char *fn, int *n_ele);
int columns_in_file(FILE *fp, int dont_warn_if_not_at_start=0);
void strip_trailing_white(char *str);
void skip_whitespace(FILE *fp);
double safe_strtod(char *s);
int safe_readInt(FILE *fp);
double safe_readDouble(FILE *fp);
void safe_readStr(FILE *fp, char str[]);
int peekChar(FILE *fp);
inline void safe_ungetc(int c, FILE *fp)
  { if (ungetc(c, fp) == EOF && c != EOF) error("safe_ungetc: failed"); }

inline double round(double x) { return floor(x+0.5); }
inline double fract(double x) { return x - floor(x); }
float float_rand(float lowlim=0.0, float highlim=1.0);
double doub_rand(double lowlim=0.0, double highlim=1.0);
inline int bernoulli(double p) {return next_random_number() < p;}
int int_rand(int lowlim, int highlim);
	 // returns a number between lowlim and highlim INCLUSIVE
inline int rand_int(int lowlim, int highlim) { return int_rand(lowlim, highlim); }
int nonunif_int_dev(double *prob_list, int len);
void select_nbits(int *arr, int nele, int nbits);
void select_frac_bits(int *arr, int nele, double frac);
int *new_select_frac_bits(int nele, double frac);

int get_int(const char *prompt);
double gammln(double xx);
float get_float(const char *prompt);
int factorial(int x);
double float_factorial(double n);
int choose(int m, int n);
double log_float_choose(double m, double n);
double log_double_choose(double m, double n);
double entropy(double f);
inline double sigmoid(double x) { return 1.0/(1.0+exp(-x)); }
inline double logsigmoid(double x) { return -log1p(exp(-x)); }
inline double log2(double x) {return log(x)/log(2.0);}
inline double plog2p(double p) {return (p == 0)?0:p * log2(p);}
double logDiffExp(double a, double b);
double logSumExp(double a, double b);
double logSumExp(double *arr, int nele); 

int min_in_arr(const int *x, int n_points);
int max_in_arr(const int *x, int n_points);
int min_fabs_in_arr(const int *x, int n_points);
int max_fabs_in_arr(const int *x, int n_points);
float min_in_arr(const float *x, int n_points);
float max_in_arr(const float *x, int n_points);
float min_fabs_in_arr(const float *x, int n_points);
float max_fabs_in_arr(const float *x, int n_points);
double min_in_arr(const double *x, int n_points);
double max_in_arr(const double *x, int n_points);
double min_fabs_in_arr(const double *x, int n_points);
double max_fabs_in_arr(const double *x, int n_points);

inline double sqr(double x) {return x*x;}
inline double cube(double x) {return x*x*x;}
inline int sqr(int x) {return x*x;}
inline int is_odd(int x) {return x&1;}
inline int is_even(int x) {return !is_odd(x);}
inline int is_int(double x) {return (floor(x) == x);}
inline int isInt(double x) {return (floor(x) == x);}
inline int isOdd(int x) {return x&1;}
inline int isEven(int x) {return !is_odd(x);}

inline int max(int a, int b) {return a<b?b:a;}
inline int max(int a, int b, int c) {return max(max(a,b),c);}
inline int min(int a, int b) {return a<b?a:b;}
inline int min(int a, int b, int c) {return min(min(a,b),c);}
inline int max(int a, int b, int c, int d) {return max(max(a,b),max(c,d));}
inline int min(int a, int b, int c, int d) {return min(min(a,b),min(c,d));}
inline double max(double a, double b) {return a<b?b:a;}
inline double min(double a, double b) {return a<b?a:b;}
inline double max(double a, double b, double c) {return max(max(a,b),c);}
inline double min(double a, double b, double c) {return min(min(a,b),c);}
inline double max(double a, double b, double c, double d) {return max(max(a,b),max(c,d));}
inline double min(double a, double b, double c, double d) {return min(min(a,b),min(c,d));}
inline double mean(double a, double b) { return (a+b)*0.5; }
int median(int a, int b, int c);
double median(double a, double b, double c);
inline int truncateValue(int lowLim, int val, int highLim)
  { if (val<lowLim) return lowLim; else if (val>highLim) return highLim; else return val; }
inline int saturate(int lowLim, int val, int highLim)
  { if (val<lowLim) return lowLim; else if (val>highLim) return highLim; else return val; }
inline double truncateValue(double lowLim, double val, double highLim)
  { if (val<lowLim) return lowLim; else if (val>highLim) return highLim; else return val; }
inline double saturate(double lowLim, double val, double highLim)
  { if (val<lowLim) return lowLim; else if (val>highLim) return highLim; else return val; }

int sum_arr(const int *arr, int nele);
double sum_arr(const double *arr, int nele);
double mean_in_arr(const int *arr, int nele);
double se_in_arr(const int *arr, int nele);
double mean_in_arr(const double *arr, int nele);
double se_in_arr(const double *arr, int nele);
inline double var_in_arr(const double *arr, int nele) 
		{ return sqr(se_in_arr(arr, nele)); }

inline int sign(double x) {return (x>0)?1:((x<0)?-1:0); }

void qsort_int(int *arr, int n_ele);
void qsort_float(float *arr, int n_ele);
void qsort_double(double *arr, int n_ele);
void reorder_int_arr(int *arr, int n_ele);
void reorder_float_arr(float *arr, int n_ele);
void reorder_double_arr(double *arr, int n_ele);
int is_num(float x);

char *time_str(int strip_cr=0);
const char *get_hostname(void);
const char *get_userid(void);
void bswap(char *b1, char *b2, int length);
inline void swap(double &a, double &b) {double t; t=a;a=b;b=t;}
inline void swap(int &a, int &b) {int t; t=a;a=b;b=t;}
inline void swap(int* &a, int* &b) {int* t; t=a;a=b;b=t;}

inline void zero_array(char *p, int len) {bzero(p, sizeof(char)*len);}
inline void zero_array(int *p, int len) {bzero((char*)p, sizeof(int)*len);}
inline void zero_array(float *p, int len) {bzero((char*)p, sizeof(float)*len);}
inline void zero_array(double *p, int len) {bzero((char*)p, sizeof(double)*len);}
inline void zeroArray(char *p, int len) {bzero(p, sizeof(char)*len);}
inline void zeroArray(int *p, int len) {bzero((char*)p, sizeof(int)*len);}
inline void zeroArray(float *p, int len) {bzero((char*)p, sizeof(float)*len);}
inline void zeroArray(double *p, int len) {bzero((char*)p, sizeof(double)*len);}

int double_eq(double d1, double d2, 
		double epsilon=1.0e-6, double delta=1.0e-10);
int double_leq(double d1, double d2, 
		double epsilon=1.0e-6, double delta=1.0e-10);
int double_geq(double d1, double d2, 
		double epsilon=1.0e-6, double delta=1.0e-10);

int argmax(int *arr, int nele);
int argmax(float *arr, int nele);
int argmax(double *arr, int nele);
int argmin(int *arr, int nele);
int argmin(float *arr, int nele);
int argmin(double *arr, int nele);

inline char yOrN(int x) { return (x!=0)?'y':'n';}
void printCmdLine(int argc, char *argv[], FILE *fp);

inline const char *countingSuffix(int x) { switch (x % 10)
 { case 1: return "st"; case 2: return "nd"; case 3: return "rd"; default: return "th"; }} 
 

template<class T>
inline T **new_array(int m, int n)
  {
  T **ptr = new T*[m];
  assert(ptr != NULL);
  for (int i=0; i < m; i++)
    { ptr[i] = new T[n]; assert(ptr[i] != NULL); }
  return ptr; 
  }

template<class T>
inline void del_array(T** ptr, int m, int n)
  {
  assert(ptr != NULL);
  for (int i=0; i < m; i++)
    { assert(ptr[i] != NULL); delete[] ptr[i]; }
  delete[] ptr;
  }

#ifndef PI
// pi and e, to 100 digits accuracy
const double PI = 
	3.141592653589793238462643383279502884197169399375105820974944592\
307816406286208998628034825342117068;
#endif

const double E = 
        2.7182818284590452353602874713526624977572470936999595749669676277\
24076630353547594571382178525166427;

#endif 		// _MISC_H

