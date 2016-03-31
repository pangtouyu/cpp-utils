#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <values.h>

#include "misc.h"
#include "integrate.h"

#define RHOMBERG_MIN_N (5)
#define RHOMBERG_MAX_N (20)  // the largest value of n for Romberg Integration
			     // must be < 31. (for 32-bit integers)

//#ifndef PI
//const double PI = 3.14159265358979323846;
//#endif

static int g_evals;	// keeps track of number of function evals
			// (used in test code only.)

// the function f() ---------------------------------------------------------

static double sqr(double x, void *p)
  {
  return x*x;
  }

static double small_phi(double x, void *p)
  {
  static double inv_root_2pi = 1.0/sqrt(2.0*PI);

  g_evals++;

  return inv_root_2pi * exp(-0.5 * x * x);
  }

static double studentT(double x, void *n_p)
  {
  int n = *((int*)n_p);
  assert(n > 0);
  g_evals++; if (x < 0) return 0; 
  return exp(gammln((n+1.0)/2) - gammln(n/2.0)) / sqrt(n*PI) 
            * pow(1+sqr(x)/n, -(n+1.0)/2.0);
  }

// Adaptive Simpson's method ------------------------------------------------

double simpson(doub_func_p f, void *parms, 
		double a, double b, double epsilon, int level, 
			int level_max, double *result)
  { 
  double retval;
  double left_simpson, right_simpson;
  double one_simpson, two_simpson, f_a, f_b, f_c;
  double c, d, e;
  double h;

  level++;
  h = b-a;
  c = (a+b)/2.0;

  f_a = (*f)(a, parms);		// to reduce the number of function evals 
  f_b = (*f)(b, parms);
  f_c = (*f)(c, parms);

  one_simpson = h * (f_a + 4.0*f_c + f_b) / 6.0;

  d = (a+c)/2.0;
  e = (c+b)/2.0;

  two_simpson = h * (f_a + 4.0*(*f)(d,parms) + 2.0*f_c+ 4.0*(*f)(e,parms) + f_b)/12.0;

  if (level > level_max)
	{
	retval = two_simpson;
	if (result != NULL)
		*result = retval;
	fprintf(stderr,"Reached level_max. Stopping.\n");
	}
  else if (fabs(two_simpson - one_simpson) < 15.0*epsilon)
	{
	retval = two_simpson;
	if (result != NULL)
		*result = retval;
	}
  else  // continue recusively
	{
	left_simpson = simpson(f,parms,a,c,epsilon/2.0, level, level_max, result);
	right_simpson = simpson(f,parms,c,b,epsilon/2.0, level, level_max, result);
	retval = left_simpson + right_simpson;
	if (result != NULL)
		*result = retval;
	}

  return retval;	
  }

double simpson(doub_func_p f, void *parms, double a, double b, double epsilon)
  {
  return simpson(f, parms, a, b, epsilon, 0, MAXINT, (double*)NULL);
  }

// Romberg Integration ------------------------------------------------------

double romberg(doub_func_p f, void *parms, double a, double b, double epsilon)
  {
  double r[RHOMBERG_MAX_N+1][RHOMBERG_MAX_N+1];
  double sum;
  double h;
  int i, j, k;

  if (!(a<=b))
	fprintf(stderr,"ERROR: romberg: a=%g, b=%g\n", a, b);
  assert(a <= b);

  h = b-a;

  r[0][0] = (h/2.0) * ((*f)(a,parms) + (*f)(b,parms));
  
  for (i = 1; i <= RHOMBERG_MAX_N; i++)
	{
	h = h/2.0;
	sum = 0;
	for (k=1; k <= (1<<(i-1)); k++)
		sum += f(a + (2*k - 1)*h, parms);

	r[i][0] = 0.5 * r[i-1][0] + sum * h;

	for (j = 1; j <= i; j++)
		r[i][j] = r[i][j-1] + (r[i][j-1] - r[i-1][j-1]) /(pow(4,j) - 1);

	if (i >= RHOMBERG_MIN_N
	     && fabs(r[i][i] - r[i-1][i-1]) <= epsilon)	// stopping criteria
		return r[i][i];
	}

  if (fabs(r[i][i] - r[i-1][i-1]) > epsilon)	
	{
	fprintf(stderr,"Romberg Integration reached RHOMBERG_MAX_N. Stopping.\n");
  	return r[i-1][i-1];
	}
  else
	return r[i][i];
  }

double romberg_phi(double x)
  {
  const double EPSILON = 1.0e-6;

  if (x < 0)
	return 1.0 - romberg_phi(-x);

  return 0.5 + romberg(small_phi, NULL, 0.0, x, EPSILON);
  }

double simpson_phi(double x)
  {
  const double EPSILON = 1.0e-6;

  if (x < 0)
	return 1.0 - simpson_phi(-x);

  return 0.5 + simpson(small_phi, NULL, 0.0, x, EPSILON);
  }

double romberg_studentTcdf(double x, int n)
  {
  const double EPSILON = 1.0e-6;

  if (x < 0)
	return 1.0 - romberg_phi(-x);

  return 0.5 + romberg(studentT, (void*)&n, 0.0, x, EPSILON);
  }
/*
// test code
int main(void)
  {
  double simpson_val, romberg_val;
  double f;

  for (;;)
	{
	scanf("%lf", &f);

	g_evals=0;
	simpson_val = simpson_phi(f, NULL);
	printf("Simpson -> %f (%d evals)\n", simpson_val, g_evals);
	g_evals=0;
	romberg_val = romberg_phi(f, NULL);
	printf("Romberg -> %f (%d evals)\n", romberg_val, g_evals);
	}
  }
*/
