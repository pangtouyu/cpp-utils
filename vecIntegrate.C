#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <values.h>

#include "matrix.h"
#include "vecIntegrate.h"

#define RHOMBERG_MIN_N (5)
#define RHOMBERG_MAX_N (20)  // the largest value of n for Romberg Integration
			     // must be < 31. (for 32-bit integers)

static int g_evals;	// keeps track of number of function evals
			// (used in test code only.)


Vector gausPair(double x, void *p)
  {
  Vector retval(2);
  g_evals++;

  static double inv_root_2pi = 1.0/sqrt(2.0*PI);

  double prob = inv_root_2pi * exp(-0.5 * x * x);

  retval[0] = x * prob;
  retval[1] = prob;

  return retval;
  }

// Romberg Integration ------------------------------------------------------

Vector vecRomberg(vec_func_p f, void *parms, double a, double b, Vector &epsilons)
  {
//  epsilons.print();
  int dim = epsilons.dim();

  Vector *r[RHOMBERG_MAX_N+1]; 
  for (int i=0; i <= RHOMBERG_MAX_N; i++) r[i] = new Vector[RHOMBERG_MAX_N+1](dim);
  Vector sum(dim); sum.make_zero();
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
	sum.make_zero();
	for (k=1; k <= (1<<(i-1)); k++)
		sum += f(a + (2*k - 1)*h, parms);

	r[i][0] = 0.5 * r[i-1][0] + sum * h;

	for (j = 1; j <= i; j++)
		r[i][j] = r[i][j-1] + (r[i][j-1] - r[i-1][j-1]) /(pow(4,j) - 1);

	if (i >= RHOMBERG_MIN_N)
	  {
	  int converged=1;
	  for (int d=0; d < dim; d++)
	    if (fabs(r[i][i][d] - r[i-1][i-1][d]) > epsilons[d])	// stopping criteria
	      converged=0;
	  if (converged)
	    {
	    Vector retval = r[i][i];
  	    for (int i=0; i <= RHOMBERG_MAX_N; i++) delete[] r[i];
	    return retval;
	    }
	  }
	}

  int converged=1;
  for (int d=0; d < dim; d++)
    if (fabs(r[i][i][d] - r[i-1][i-1][d]) > epsilons[d])	// stopping criteria
      converged=0;

  if (!converged)
	{
	fprintf(stderr,"Romberg Integration reached RHOMBERG_MAX_N. Stopping.\n");
  	Vector retval = r[i-1][i-1];
  	for (int i=0; i <= RHOMBERG_MAX_N; i++) delete[] r[i];
	return retval;
	}
  else
    	{
  	Vector retval = r[i][i];
  	for (int i=0; i <= RHOMBERG_MAX_N; i++) delete[] r[i];
	return retval;
	}
  }

/*
// test code
int main(void)
  {
  double simpson_val, romberg_val;
  double f;

  Vector epsilons(2);
  epsilons[0] = 1e-4;
  epsilons[1] = 1e-4;

  for (;;)
	{
	scanf("%lf", &f);

	// Calcuate E[X | X >= x] for standard gaussian X
	g_evals=0;
	Vector romberg_val = vecRomberg(gausPair, NULL, f, 5, epsilons);
	printf("Romberg -> %f, %f -> %f (%d evals)\n", romberg_val[0], romberg_val[1], 
	  			romberg_val[0]/romberg_val[1], g_evals);
	}
  }
*/


