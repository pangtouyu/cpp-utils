/* 
 * quasi.c		
 *
 * Implementation of the Quasi-Random Number generator
 * currently hardwired to no more than 52 dimensions
 *
 * modified by Andrew Ng to C++/oop (12/96)
 *
 * See W.H. Press and S.A. Teukolsky, 1989, Quasi- (that is,
 * Sub-) Random Numbers, Computers in Physics V3, No. 6,
 * (Nov/Dec 1989), pp. 76-79
 *
 *
 */

#include <stdio.h>
#include <cstdlib>
#include "quasi.h"

static char rcsid[] = "@(#)quasi.c	1.5 10:15:59 4/18/94   EFC";

/* the primitive polynomial coefficients for up to degree 8  */
static int ip[] = { 0, 1, 1, 2, 1, 4, 2, 4, 7, 11, 13, 14,
			  1, 13, 16, 19, 22, 25,
			  1, 4, 7, 8, 14, 19, 21, 28, 31, 32, 37, 41, 42,
			  50, 55, 56, 59, 62,
			  14, 21, 22, 38, 47, 49, 50, 52, 56, 67, 70, 84,
			  97, 103, 115, 122 };
static int mdeg[] = { 1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5,
			    6, 6, 6, 6, 6, 6,
			    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
			    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8
			     };

static int            maxdim = sizeof(mdeg) / sizeof(int);

static int            maxbit = 30;	    /* must be no more than
					       number of bits in ulong - 1 */

static unsigned long  int*  iv = NULL;
static double               factor = 1.0;
static int instances = 0;

#define INDEX(k,j)	[(k) + (j-1) * maxdim]


static void QuasiInit()	               /* initialize the direction numbers */
{                                      /* only done once */
	int j, k, l, ipp, niv;
	unsigned long int mval, limit;
	unsigned long int i;
        
	niv = maxdim * maxbit;
        iv = new (unsigned long int)[niv];

	for (k = 0; k < niv; k++)
        	iv[k] = 0;

        for (k = 0; k < maxdim; k++)
        	iv[k] = 1;

	mval = 4;
	ipp = 1;

        limit = (1L << (maxbit-1) );

        for (k = maxdim, j = 0; k < niv-1; k += 2)
        {
		iv[k] = ipp;
                if (++j == maxdim)
                {
		        if ( mval < limit )
			          mval *= 2;
                        ipp += 2;
                        j = 0;
                }

		if ( ipp > mval )
                	ipp = 1;
                        
                iv[k+1] = ipp;
                if (++j == maxdim)
                {
		        if ( mval < limit )
			          mval *= 2;
                        ipp += 2;
                        j = 0;
                }
		else
                {
                	ipp += 2;
                        if ( ipp > mval )
                        	ipp = 1;
                }
                

        }

	for (k = 0; k < maxdim; k++)
        {
        	/* normalize the set iv values */
        	for (j = 1; j <= mdeg[k]; j++)
                	iv INDEX(k,j) *= (1L << (maxbit - j));

		/* calcululate the rest of the iv values */
		for (j = mdeg[k] + 1; j <= maxbit; j++)
                {
                	ipp = ip[k];
                        
                        /* calculate Gray code of iv */
                        i = iv INDEX(k, j - mdeg[k]);
                        i ^= i / (1L << mdeg[k]);
                        
                        for (l = mdeg[k] - 1; l >= 1; l--)
                        {
                        	if ( ipp & 1 )
                                	i ^= iv INDEX(k, j-l);
                                ipp /= 2;
                        }

                        iv INDEX(k,j) = i;

        	}
	}

        factor = 1.0 / (1L << maxbit);

                        
}

Quasi::Quasi(int dimension)
{
        int k;

        dim = dimension;
	index = 0;

	if ( dim > maxdim )		// check if dimension is too large   
        {				
        	dim = maxdim;
		fprintf(stderr,"ERROR: QuasiRandomInitialize. dimension "
				"too large.");
		exit(-1);
        }
                                
//	qr->ix = (unsigned long int*)malloc(qr->dim*sizeof(unsigned long int) );
	ix = new (unsigned long int)[dim];

	for (k = 0; k < dim; k++)
        	ix[k] = 0L;

                
	if ( instances++ == 0 )
		QuasiInit();

}

Quasi::~Quasi(void)
{
	if ( --instances == 0 )
	    delete[] iv;

	delete[] ix;
}

void Quasi::QuasiRandomNumber(double x[])
{
	int i, j, k;
	unsigned long int im = index++;

        /* find rightmost zero bit  */
        for (j = 0; j < maxbit; j++, im >>= 1)
        	if ( (im & 1L) == 0 )
                		break;

	i = j * maxdim;

        for (k = 0; k < dim; k++)
        {
        	ix[k] ^= iv[i + k];		/* integer values  */
                x[k]   = (double) ( factor * (double)ix[k] );
        }

}

/* Test code 

int main(void)
  {
  Quasi q(2);
  double x[2];

  for (int ctr=0; ctr < 1024; ctr++)
	{
	q.QuasiRandomNumber(x);
	printf("%f %f\n", x[0], x[1]);
	}

  return 0;
  }
*/
 
