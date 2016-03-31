#include <assert.h>
#include <stdio.h>

#include "lregress.h"
#include "misc.h"
#include "matrix.h"

int g_debugging_flag = 0;

const Vector *g_y;
const Vector *g_pi;
const Vector *g_beta;
const Vector *g_oldBeta;

// debugging code
void bar(void)
  {
  printf("pi: "); g_pi->print();
  printf("y: "); g_y->print();
  printf("beta: "); g_beta->print();
  printf("oldBeta: "); g_oldBeta->print();
  }

void LRegress::old_fitIt(const Matrix &x, const Vector &y, int suppressWarnings)
  {
  const int MINITERS = 2;
  const int MAXITERS = 100;
  const int WARNITERS = 10;

  assert(x.nrows() == y.dim());
  assert(x.ncols() == dim);

  int nsamples = x.nrows();

  Vector w(nsamples); // diag matrix, represeted as vector
  Vector pi(nsamples);
  Vector we(nsamples);

  beta.make_zero();	// starting Beta

  Vector oldPi(nsamples);
  Vector oldBeta(dim);

  int iter;
  int earlyStopping=0, converged = 0;
  for (iter=0; iter < MAXITERS; iter++)
    {
    g_y = &y; g_pi = &pi;
    g_oldBeta  = &oldBeta; g_beta = &beta;

    oldPi = pi;

    for (int i=0; i < nsamples; i++)
	{
	double eta = 0.0;
	for (int j=0; j < dim; j++)
	    eta += x[i][j] * beta[j];
	pi[i] = 1.0/(1.0 + exp(-eta));
	// double fp = exp(eta) / sqr(1.0+exp(eta));
	double fp = 1.0 / ( (1.0+exp(eta)) * (1.0+exp(-eta)) );
	if (isnan(fp))
		{
		printf("%f %f %f \n", eta, fp, x[i][0]);
		beta.print();
		}
	assert(!isnan(fp));
	assert(!isnan(pi[i]));
	we[i] = (y[i]-pi[i]);
	w[i] = fp;
	}

    if (g_debugging_flag) bar();

    Matrix xtwx = x.find_AtransposeWA(w);
    oldBeta = beta;

    int illCond, inverseFailed;
    double add_epsilon = 1.0e-4;
    for (int dex=0; dex < xtwx.nrow(); dex++)
	xtwx[dex][dex] += add_epsilon;		// try to make it non-singular

    Matrix invXtwx = xtwx.find_inverse(&illCond, NULL, &inverseFailed);
    while (inverseFailed || illCond)
      {
      if (inverseFailed & add_epsilon > 1e-2)
	{
	fprintf(stderr,"\nWARNING: lregress: Matrix inversion failed. "
				"(iteration %d; %d dimensional, %d samples)\n", 
				iter+1, x.ncol(), y.dim());
	// fprintf(stderr,"WARNING: oldBeta: "); oldBeta.print(stderr);
	// fprintf(stderr,"WARNING: beta (probably garbage): "); beta.print(stderr);
	}
      if (illCond && add_epsilon > 1e-2)
	{
	fprintf(stderr,"\nWARNING: lregress: (On iteration %d; %d dimensional, %d samples)\n", 
				iter+1, x.ncol(), y.dim());
	fprintf(stderr,"WARNING: lregress: Ill conditioned X^tWX. \n");
	}
      // fprintf(stderr, "WARNING: beta = "); beta.print(stderr);  fprintf(stderr, " \n");

      if (add_epsilon > 1e-2)
        fprintf(stderr, "WARNING: add_epsilon=%g\n", add_epsilon);

      add_epsilon*=2;
      for (int dex=0; dex < xtwx.nrow(); dex++)
	xtwx[dex][dex] += add_epsilon;		// try to make it non-singular

      invXtwx = xtwx.find_inverse(&illCond, NULL, &inverseFailed);
      }


    beta+= invXtwx * (x.find_transpose() * we);

    for (int i=0; i < beta.dim(); i++)
	{
	earlyStopping = 0;
	if (isnan(beta[i]))
	    {
	    fprintf(stderr, "WARNING: lregress: NaN in beta. stopping early. "
				"(iter %d: %d dim, %d samples)", iter+1, x.ncol(), y.dim());
	    beta = oldBeta;
	    earlyStopping = 1; 
	    break;
	    }
	}
    if (earlyStopping)
	break;

    int superConverged = 1;		// "superconverged" if classifying
					// everything correctly.
    for (int i=0; i < nsamples; i++)
	{
	if (fabs(y[i] - pi[i]) > 1.0e-2)
	    {
	    superConverged = 0;
	    break;
	    }
	}
    if (superConverged)
	{
	converged = 1;
	// printf("Super Convergence.\n");
	break;
	}

//    if (iter+1 >= MINITERS && vect_approx_equal(oldPi, pi, 0.0, 1.0e-8))
    if (iter+1 >= MINITERS && vect_approx_equal(oldPi, pi, 0.0, 1.0e-3))
	{
	converged = 1;
	break;
	}

    if (iter+1 % WARNITERS == 0)
	fprintf(stderr, "WARNING: LRegress::fitIt: Newton-Ralphson not "
			"converged at %d iterations.", iter+1);
    }

  if (!converged && !earlyStopping)
	{
	if (!suppressWarnings)
	    {
	    fprintf(stderr,"WARNING: lregress: fitIt did not converge "
			   "(%d dimensional, %d samples)\n", x.ncol(), y.dim());
	    fprintf(stderr, "WARNING: beta: "); beta.print(stderr);
	    }
	}

  // printf("beta = "); beta.print();
  // printf("fitted values: \n");
  // pi.print();

  // printf("lr: %d iterations.\n", iter);

  return;
  }

void LRegress::fitIt(const Matrix &x, const Vector &y, int suppressWarnings,
			double lambda)
  {
  if (lambda < 0) error("LRegress::fitIt: Got negative lambda.");
  const int MINITERS = 2;
  const int MAXITERS = 100;
  const int WARNITERS = 10;

  assert(x.nrows() == y.dim());
  assert(x.ncols() == dim);

  int nsamples = x.nrows();

  Vector w(nsamples); // diag matrix, represeted as vector
  Vector pi(nsamples);
  Vector we(nsamples);

  beta.make_zero();	// starting Beta

  Vector oldPi(nsamples);
  Vector oldBeta(dim);

  int iter;
  int earlyStopping=0, converged = 0;
  for (iter=0; iter < MAXITERS; iter++)
    {
    g_y = &y; g_pi = &pi;
    g_oldBeta  = &oldBeta; g_beta = &beta;

    oldPi = pi;

    for (int i=0; i < nsamples; i++)
	{
	double eta = 0.0;
	for (int j=0; j < dim; j++)
	    eta += x[i][j] * beta[j];
	pi[i] = 1.0/(1.0 + exp(-eta));
	// double fp = exp(eta) / sqr(1.0+exp(eta));
	double fp = 1.0 / ( (1.0+exp(eta)) * (1.0+exp(-eta)) );
	// if (isnan(fp)) { printf("%f %f %f \n", eta, fp, x[i][0]); beta.print(); }
	assert(!isnan(fp)); assert(!isnan(pi[i]));
	we[i] = (y[i]-pi[i]);
	w[i] = fp;
	}

    Matrix xtwx = x.find_AtransposeWA(w);
    if (lambda > 0)
      for (int dex=0; dex < xtwx.nrow(); dex++)
	xtwx[dex][dex] += lambda;	// form the (negative of the) Hessian,
					// for the maximization problem
    oldBeta = beta;

    int illCond, inverseFailed;		// illCond will be ignored
    Matrix invXtwx = xtwx.find_inverse(&illCond, NULL, &inverseFailed);
    if (inverseFailed)
      error("lregress::fitIt: Matrix inversion failed.");

    Vector grad = x.find_transpose() * we;	// grad for the maximization (not min) problem
    if (lambda != 0)
      grad -= lambda * beta; 		// from penalty term  -lambda*beta'*beta

    beta+= invXtwx * grad;
//    beta+= invXtwx * (x.find_transpose() * we);

    for (int i=0; i < beta.dim(); i++)
	{
	earlyStopping = 0;
	if (isnan(beta[i]))
	    {
	    fprintf(stderr, "WARNING: lregress: NaN in beta. stopping early. "
				"(iter %d: %d dim, %d samples)", iter+1, x.ncol(), y.dim());
	    beta = oldBeta;
	    earlyStopping = 1; 
	    break;
	    }
	}
    if (earlyStopping)
	break;

    int superConverged = 1;		// "superconverged" if classifying
					// everything correctly.
    for (int i=0; i < nsamples; i++)
	{
	if (fabs(y[i] - pi[i]) > 1.0e-2)
	    {
	    superConverged = 0;
	    break;
	    }
	}
    if (superConverged)
	{ converged = 1; break; }

//    if (iter+1 >= MINITERS && vect_approx_equal(oldPi, pi, 0.0, 1.0e-8))
    if (iter+1 >= MINITERS && vect_approx_equal(oldPi, pi, 0.0, 1.0e-5))
	{ converged = 1; break; }

    if (iter+1 % WARNITERS == 0)
	fprintf(stderr, "WARNING: LRegress::fitIt: Newton-Ralphson not "
			"converged at %d iterations.", iter+1);
    }

  if (!converged && !earlyStopping)
	{
	if (!suppressWarnings)
	    {
	    fprintf(stderr,"WARNING: lregress: fitIt did not converge "
			   "(%d dimensional, %d samples)\n", x.ncol(), y.dim());
	    fprintf(stderr, "WARNING: beta: "); beta.print(stderr);
	    }
	}

  return;
  }

double LRegress::fittedValue(const Vector &x) const
  {
  x.assert_dim(dim);
  return 1.0/(1.0+exp(-dot_prod(x,beta)));
  }

/*************************************
// test code
int main(void)
  {
  Matrix x(1,1);
  x.initialize_from_file("test.dat");
	// reads in something like:
	//  2 4 5 1
	//  3 1 2 0
	//  2 3 5 0
	// ... last column is desired "output", and will be 
	//  replaced with 1's for constant/intercept term

  Vector y = x.extract_col(x.ncols()-1);
  for (int i=0; i < x.nrows(); i++)
	x[i][x.ncols()-1] = 1.0;	// intercept term

  printf("Training sample (input):\n");
  x.print();
  printf("Training sample (output):\n");
  y.print();

  while (1)
    {
    double lambda; 
    printf("lambda? "); fflush(stdout);
    scanf("%lf", &lambda);
    LRegress lr(x.ncols());
    lr.fitIt(x, y, 0, lambda);

    printf("Parameters converged to: "); lr.betaVal().print();
    printf("Fitted values: \n");
    for (int i=0; i < x.nrows(); i++)
  	printf("%.3f ", lr.fittedValue(x.extract_row(i)));
    printf("\n");
    }

  return 0; 
  }
// g77 dgesvd.o lregress.o matrix.o Random_Number.o misc.o -lm -llapack -lblas
// g77 lregress.C misc.o Random_Number.o matrix.o dgesvd.o -lm -llapack

/*************************************/

