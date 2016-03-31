#include <assert.h>
#include <stdio.h>
#include <math.h>

#include "matrix.h"
#include "misc.h"
#include "normal.h"
#include "randvar.h"
#include "Random_Number.h"

const double PI = 3.14159265358979323846;

// handles the memoized determinant of sigma. If not there, also 
// finds the inverse of sigma, and memoizes that.
double MultiNormal::find_sigmadet(void)
  {
  if (sigmadet == NULL)
	{
	assert(sigmainv == NULL);

	sigmadet = new double;
	sigmainv = new Matrix(sigma.nrow(),sigma.ncol(),0);
	*sigmainv = sigma;
	sigmainv->invert_me(sigmadet);	// note that the abs value of 
					// a covariance matrix's determinant
					// is always just the determinant
	}

  return *sigmadet;
  }

// handles the memoized inverse of sigma. If not there, also 
// finds the determinant of sigma, and memoizes that.
Matrix *MultiNormal::find_sigmainv(void)
  {
  if (sigmainv == NULL)
	{
	assert(sigmadet == NULL);

	sigmadet = new double;
	sigmainv = new Matrix(sigma.nrow(),sigma.ncol(),0);
	*sigmainv = sigma;
	sigmainv->invert_me(sigmadet);	// note that the abs value of 
					// a covariance matrix's determinant
					// is always just the determinant
	}

  return sigmainv;
  }

// Call this whenever anything is changed. Yes, MUST!!
// (Invalidates the memoized stuff.)
void MultiNormal::invalidate_cached_dat(void)
  {
  if (sigmainv != NULL)
	{
	delete sigmainv;
	sigmainv = NULL;
	}

  if (sigmadet != NULL)
	{
	delete sigmadet;
	sigmadet = NULL;
	}

  if (generate_t != NULL)
	{
	delete generate_t;
	generate_t = NULL;
	}

  }

// copy constructor
MultiNormal &MultiNormal::operator=(const MultiNormal &mn)
  {
  if (dim != mn.dim)
	{
	fflush(stdout);
	fprintf(stderr,"MultiNormal ERROR: (Probable internal error.)\n"
                   "MultiNormal ERROR: Tried to assign %dd to %dd\n",
				mn.dim, dim);
	abort();
	}

  mu = mn.mu;
  sigma = mn.sigma;

  if (mn.sigmainv == NULL)
	sigmainv = NULL;
  else
	{
	sigmainv = new Matrix(mn.sigmainv->nrow(), mn.sigmainv->ncol(),0);
	*sigmainv = *mn.sigmainv;
	}

  if (mn.sigmadet == NULL)
	sigmadet = NULL;
  else
	{
	sigmadet = new double;
	*sigmadet = *mn.sigmadet;
	}

  if (mn.generate_t == NULL)
	generate_t = NULL;
  else
	{
	generate_t = new Matrix(mn.generate_t->nrow(), mn.generate_t->ncol(),0);
	*generate_t = *mn.generate_t;
	}

  return *this;
  }

// copy constructor
MultiNormal::MultiNormal(const MultiNormal &mn) 
	: dim(mn.dim), sigma(mn.sigma), mu(mn.mu)
  {
  if (mn.sigmainv == NULL)
	sigmainv = NULL;
  else
	{
	sigmainv = new Matrix(mn.sigmainv->nrow(), mn.sigmainv->ncol(),0);
	*sigmainv = *mn.sigmainv;
	}

  if (mn.sigmadet == NULL)
	sigmadet = NULL;
  else
	{
	sigmadet = new double;
	*sigmadet = *mn.sigmadet;
	}

  if (mn.generate_t == NULL)
	generate_t = NULL;
  else
	{
	generate_t = new Matrix(mn.generate_t->nrow(), mn.generate_t->ncol(),0);
	*generate_t = *mn.generate_t;
	}

  return;
  }

// initializes to dim*iid N(mu0,1)
// ... where the mu0's are Uniform(0,init_range).
// (so, set init_range to 0 to get dim*iid N(0,1). )
//
// Precond: init_range >= 0;
MultiNormal::MultiNormal(int dim_, double init_range) 
		: dim(dim_), sigma(dim_,dim_,1,0), mu(dim_,0)
  {
  int ctr;

  assert(init_range >= 0);

  sigmainv = NULL;
  sigmadet = NULL;
  generate_t = NULL;

  for (ctr=0; ctr < dim_; ctr++)
	sigma[ctr][ctr] = 1;

  for (ctr=0; ctr < dim_; ctr++)
	mu[ctr] = doub_rand(0, init_range);

  return;
  }

Matrix *MultiNormal::find_generate_t(void)
  {
  if (generate_t != NULL)
	return generate_t;

  Matrix u(dim,dim,0);
  Vector diag(dim,0);
  Vector root_diag(dim,0);
  Matrix v(dim,dim,0);
  int ill_conditioned;
  int i, j;

  ill_conditioned = sigma.svd_decomp(u,diag,v);

  // if the covariance matrix is positive definite (as it should be), 
  // then u should be v.transpose. 
  if (!mat_approx_equal(u, v.find_transpose()))
	{
	internal_warn("MultiNormal::find_generate_t(): After SVD on "
			"Covariance Matrix, U != V-transpose. Might "
			"have been given a bad covariance matrix "
			"(not positive (semi)-definite)?? ");
	fprintf(stderr,"WARNING: Was sigma = \n");
	sigma.print(stderr);
	}

  if (ill_conditioned)
	internal_warn("MultiNormal::find_generate_t: Got an ill-conditioned "
			"Covariance Matrix.");

  // root_diag = sqrt(diag) (element-wise sqrt.)
  for (i=0; i < dim; i++)
    {
    assert(diag[i] >= -1.0e-8);	// should be > 0 from SVD's properties
    if (diag[i] < 0)
	diag[i] = 0;
    root_diag[i] = sqrt(diag[i]);
    }

  // find U * sqrt(diag_matrix)
  for (i=0; i < dim; i++)
    for (j=0; j < dim; j++)
	u[i][j] *= root_diag[j];

  generate_t = new Matrix(dim,dim,0);
  (*generate_t) = u;

  return generate_t;
  }

// like the other generate(), but slightly faster. (Saves on 
// some memory allocations/deallocations.)
void MultiNormal::generate(Vector &dev)
  {
  Vector x(dim);
  int ctr;

  dev.assert_dim(dim);

  for (ctr=0; ctr < dim; ctr++)
	x[ctr] = stdnormal_dev();

  dev = ((*(find_generate_t())) * x);
  dev += mu;

  return;
  }

Vector MultiNormal::generate(void)
  {
  Vector x(dim);
  int ctr;

  for (ctr=0; ctr < dim; ctr++)
	x[ctr] = stdnormal_dev();

  return ((*(find_generate_t())) * x) + mu;
  }

void MultiNormal::setparams(const Matrix &newsigma, const Vector &newmu)
  {
  invalidate_cached_dat();

  sigma = newsigma;
  mu = newmu;

  return;
  }

void MultiNormal::setmu(const Vector &newmu)
  {
  invalidate_cached_dat();	// Probably not needed for this, 
				// but nevermind....

  mu = newmu;

  return;
  }

void MultiNormal::setsigma(const Matrix &newsigma)
  {
  invalidate_cached_dat();

  sigma = newsigma;

  return;
  }

double MultiNormal::loglikelihood_at(Vector &x)
  {
  Vector diff  = (x - mu);
  double retval;
 
  assert(find_sigmadet() >= 0);
 
  retval = -log(pow((2.0*PI), (double)dim/2) * sqrt(find_sigmadet()) )
	  - 0.5 * double(Matrix(diff).find_transpose() 
 			* *(find_sigmainv()) 
 			* Matrix(diff));

  // HERE: For debugging only.... remove!
  double temp = likelihood_at(x);
  assert(isnand(temp) || double_eq(exp(retval), temp));

  return retval;
  }

double MultiNormal::likelihood_at(Vector &x)
  {
  Vector diff  = (x - mu);

  assert(find_sigmadet() >= 0);

  return 1.0/(pow((2.0*PI), (double)dim/2) * sqrt(find_sigmadet()) )
	  * exp( - 0.5 * double(Matrix(diff).find_transpose() 
			* *(find_sigmainv()) 
			* Matrix(diff)) );
  }

void MultiNormal::print(FILE *fp)
  {
  if (fp == NULL)
	fp = stdout;

  fprintf(fp, "Means: ");
  mu.print(fp);

  fprintf(fp, "\nCovariance matrix:\n");
  sigma.print(fp);

  return;
  }

//-----------------------------------------------------------------
// DiscreteRV

// Given a vector of probabilities p (that sum to 1), on each
// generat(), DiscreteRV returns i with probability p[i].

int DiscreteRV::generate(void) const
  {
  int index, temp;
  double rand;

  rand = next_random_number();
  temp = (int)floor(rand * n_indices);
  assert(temp >= 0 && temp < n_indices);
  index = indices[temp];

  if (inc_probs[index] < rand)
	{
	while (inc_probs[index] < rand)
		{
		index++;
		assert(rand < inc_probs.length());
		}
	}
  else 
	{
	while (index > 0 && inc_probs[index-1] > rand)
		index--;
	}

  assert(index == 0 || inc_probs[index-1] <= rand);
  assert(inc_probs[index] > rand);

  return index;
  }

void DiscreteRV::init_indices(void)
  {
  double this_index_start;
  int index_ctr, prob_dex;

  assert(indices==NULL);
  n_indices = (int)ceil(INDEX_DENSITY * probs.length());
  n_indices = max(n_indices,1);
  indices = new int[n_indices];

  prob_dex = 0;
  for (index_ctr=0; index_ctr < n_indices; index_ctr++)
	{
	this_index_start = (double)index_ctr/n_indices;
	assert(this_index_start < 1.0);
	while (prob_dex<inc_probs.length()-1 
			&& inc_probs[prob_dex] <= this_index_start)
		prob_dex++;
	assert(this_index_start <= inc_probs[prob_dex]);
	indices[index_ctr] = prob_dex;
	}

  // Check:
  for (index_ctr=0; index_ctr < n_indices; index_ctr++)
	{
	this_index_start = (double)index_ctr/n_indices;
	assert(this_index_start <= inc_probs[indices[index_ctr]]);
	assert(indices[index_ctr] == 0 || 
		     this_index_start >= inc_probs[indices[index_ctr]-1]);
	}

  return;
  }

DiscreteRV::~DiscreteRV()
  {
  delete[] indices;
  }

DiscreteRV::DiscreteRV(Vector &probs_) 
		: probs(probs_), inc_probs(probs.length()), indices(NULL)
  {
  int ctr;

  assert(inc_probs.length() == probs.length());

  inc_probs[0] = probs[0];
  for (ctr=1; ctr < inc_probs.length(); ctr++)
	inc_probs[ctr] = inc_probs[ctr-1] + probs[ctr];
  if (!double_eq(inc_probs[ctr-1], 1.0))
	{
	fprintf(stderr, "ERROR: Saw ");
	probs.print(stderr);
	error("DisceteRV: Given array of probabilities that don't sum to 1.");
	}
  inc_probs[ctr-1] = 1.0;

  init_indices();
  }

