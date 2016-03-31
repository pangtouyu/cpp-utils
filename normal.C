#include "Random_Number.h"
#include <math.h>

// from Numerical Recipes 2nd Ed., by W. Press
// pp. 289-290
// (modified to use "Random_Number.h" and doubles instead of floats)
// returns a normally distributed deviate with a zero mean
// and unit variance, using next_random_number() as a source
// of uniform deviates
double gasdev(void)
  {
  static int iset=0;
  static double gset;
  double fac, rsq, v1, v2;

  if (iset == 0)
	{
	// We don't have an extra deviate handy, so...
	do
	  {
	  // pick to uniform numbers in the square extending
	  // from -1 to +1 in each direction, and see if they're
	  // in the unit circle
	  v1 = 2.0 * next_random_number() - 1.0;
	  v2 = 2.0 * next_random_number() - 1.0;
	  rsq = v1 * v1 + v2 * v2;
	  }
	while (rsq >= 1.0 || rsq == 0.0);

	// now make the Box-Muller transformation to get two normal
	// deviates/ Return one and save the other for next time.
	fac = sqrt(-2.0 * log(rsq)/rsq);

	gset = v1 * fac;
	iset = 1;

	return v2 * fac;
	}
  else
	{
	iset = 0;

	return gset;
	}
  }

// returns a normally distributed deviate with the
// specified mean and variance
double normal_dev(double mean, double variance)
  {
  return mean + gasdev()*sqrt(variance);
  }

// returns a normally distributed deviate with the
// specified mean and _standard deviation_
double normal_dev_stddev(double mean, double stddev)
  {
  return mean + gasdev()*stddev;
  }

