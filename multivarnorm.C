#include <stdio.h>
#include <math.h>

#include "matrix.h"
#include "misc.h"
#include "multivarnorm.h"
#include "normal.h"
#include "Random_Number.h"

//------------------------------------------------------------
// MultivarNormal

MultivarNormal::MultivarNormal(int dim_)
		: dim(dim_), 
		  mu(dim_), sigmainv(dim_,dim_), generateT(dim_,dim_),
		  detGrad(dim_,dim_), detGradValid(0)
  {
  int illCond;

  abs_sigmainvdet = -1.0;	// an invalid value (<0)

  return;
  }

void MultivarNormal::generate(Vector &x)
  {
  x.assert_dim(dim);

  if (abs_sigmainvdet < 0)
	error("MultivarNormal: Apparently tried to generate without "
		"initializing params first.");

  Vector temp(dim);
  int ctr;

  for (ctr=0; ctr < dim; ctr++)
	temp[ctr] = standardNormal_dev();

  x = generateT * temp;
  x += mu;

  return;
  }

void MultivarNormal::setmu(const Vector &mu_)
  {
  mu = mu_;
  }

// returns 1 if ill-conditioned
//
// sets sigmainv, generateT, abs_sigmadet.
int MultivarNormal::setsigmainv(const Matrix &sigmainv_, int warnIfIllCond)
  {
  double TOL = 1.0e-7;
  Matrix u(dim,dim,0), v(dim,dim,0);
  Vector diag(dim,0);
  Vector sqrtDiag(dim,0);
  int illCond;

  detGradValid=0;

  sigmainv = sigmainv_;

  sigmainv.svd_decomp(u, diag, v);
  if (!mat_approx_equal(u.find_transpose(), v))
	{
	if (warnIfIllCond)
		{
		diag.print();
  		u.find_transpose().print();
  		v.print();
		}
	warn("(BAD) MultivarNormal: sigmainv not pos definite. (I think.)");
	}

  double max_diag = max_in_arr(diag, diag.length());
  illCond = 0;
  abs_sigmainvdet = 1.0;
  for (int ctr=0; ctr < dim; ctr++)
	{
        if (diag[ctr] <= TOL*max_diag)
                illCond = 1;
	sqrtDiag[ctr] = sqrt(diag[ctr]);
	abs_sigmainvdet *= diag[ctr];
	}

  assert(abs_sigmainvdet >= 0);

  generateT = v.find_transpose();
  for (int i=0; i < dim; i++)
    for (int j=0; j < dim; j++)
	{
	if (sqrtDiag[j] != 0.0)
		generateT[i][j] /= sqrtDiag[j];
	else
		generateT[i][j] = 0.0;	// VERY BAD if have to do this 
	}

  if (illCond && warnIfIllCond)
	warn("MultivarNormal: Ill-conditioned sigma-inverse");

  if (abs_sigmainvdet == 0 && warnIfIllCond)
	warn("BAD: MultivarNormal: sigmainv is singular!");

  return illCond?1:0;
  }

double MultivarNormal::logLikelihoodAt(const Vector &x) 
  {
  Vector diff  = (x - mu);
  double retval;

  if (abs_sigmainvdet < 0)
    error("MultivarNormal: tried to use without initializing mu and sigmainv?");
  assert(abs_sigmainvdet > 0);
 
  retval = 0.5*log(abs_sigmainvdet) - dim/2.0*log(2.0*PI)
	    - 0.5 * double(Matrix(diff).find_transpose() 
 			* sigmainv
 			* Matrix(diff));

  // HERE: For debugging only.... remove!
  // MultiNormal m(dim, 0.0, 0);
  // m.setmu(mu);
  // m.setsigma(sigmainv.find_inverse());
  // printf("EMultinormal: %f   Multinormal: %f\n", retval,
  //		m.logLikelihoodAt(x));

  return retval;
  }

void MultivarNormal::gradLoglikelihood(const Vector &x, Vector &muGrad, Matrix &sigmainvGrad)
  {
  muGrad.assert_dim(dim); sigmainvGrad.assert_dim(dim,dim);

  if (!detGradValid)
	{
	detGrad = sigmainv.find_det_deriv();
	detGradValid=1;
	}

  Vector xminusmu = x - mu;

  muGrad.make_zero();
  for (int i=0; i < dim; i++)
    for (int j=0; j < dim; j++)
	{
	sigmainvGrad[i][j] = 0.5/abs_sigmainvdet * detGrad[i][j]
				- 0.5 * xminusmu[i] * xminusmu[j];

	if (j == i)
	    muGrad[i] += -0.5 * 2 * (mu[j] - x[j]) * sigmainv[i][j];
	else
	    {
	    muGrad[i] += -0.5 * (mu[j] - x[j]) * sigmainv[i][j];
	    muGrad[j] += -0.5 * (mu[i] - x[i]) * sigmainv[i][j];
	    }
	}

  return;
  }

void MultivarNormal::gradLikelihood(const Vector &x, Vector &muGrad, Matrix &sigmainvGrad)
  {
  gradLoglikelihood(x, muGrad, sigmainvGrad);
  double p = likelihoodAt(x);

  muGrad *= p;
  sigmainvGrad *= p;

  return;
  }

MultivarNormal &MultivarNormal::operator=(const MultivarNormal &mn)
  {
  if (dim != mn.dim)
	{
	fflush(stdout);
	fprintf(stderr,"MultivarNormal ERROR: (Probable internal error.)\n"
                   "MultivarNormal ERROR: Tried to assign %dd distribution "
			"to %dd distribution\n", mn.dim, dim);
	abort();
	}
  sigmainv = mn.sigmainv;
  generateT = mn.generateT;
  abs_sigmainvdet = mn.abs_sigmainvdet;
  mu = mn.mu;
  assert(dim == mn.dim);

  detGrad = mn.detGrad;
  detGradValid = mn.detGradValid;

  return *this;
  }

// copy constructor
MultivarNormal::MultivarNormal(const MultivarNormal &mn)
		: sigmainv(mn.sigmainv), generateT(mn.generateT), 
		   dim(mn.dim), mu(mn.mu), 
		   abs_sigmainvdet(mn.abs_sigmainvdet),
		   detGrad(mn.detGrad), detGradValid(mn.detGradValid)
  {
  return;
  }

//------------------------------------------------------------
// GaussianMix

double GaussianMix::logLikelihoodAt(const Vector &x)
  {
  static Vector temp(ngaussian);	// temporary storage

  if (temp.dim() < ngaussian)
	temp.resize(ngaussian);

  for (int g=0; g < ngaussian; g++)
	temp[g] = log(weights[g]) + gaussians[g].logLikelihoodAt(x);

  // temp[i] = log(wi*loglike_i).
  // we want log(w0*loglike_0 + w1*loglike_1 + ... )

  return logSumExp(temp, ngaussian);
  }

void GaussianMix::setEMReassign(Vector &x, Vector &reassign)
  {
  reassign.assert_dim(ngaussian);

  for (int g=0; g < ngaussian; g++)
	reassign[g] = log(weights[g]) + gaussians[g].logLikelihoodAt(x);
  double logSum = logSumExp(reassign, ngaussian);

  for (int g=0; g < ngaussian; g++)
	reassign[g] = exp(reassign[g]-logSum);

  assert(double_eq(sum_arr(reassign,ngaussian),1.0));

  return;
  }

GaussianMix::GaussianMix(int dim_, int ngaussian_, int initialize)
	: dim(dim_), ngaussian(ngaussian_), weights(ngaussian_,1,1.0/ngaussian_)
  {
  gaussians = new MultivarNormal[ngaussian](dim);

  assert(double_eq(sum_arr(weights, ngaussian),1.0));

  if (initialize)
	{
	Vector mu(dim);
	Matrix sigmainv(dim,dim); 
	sigmainv.make_identity();
	gaussians[0].setmu(mu);
	gaussians[0].setsigmainv(sigmainv);

	for (int ctr=1; ctr < ngaussian; ctr++)
		gaussians[ctr] = gaussians[0];
	}

  return;
  }

GaussianMix::~GaussianMix()
  {
  delete[] gaussians;
  }

GaussianMix &GaussianMix::operator=(const GaussianMix &gm)
  {
  if (dim != gm.dim || ngaussian != gm.ngaussian)
	{
	fflush(stdout);
	fprintf(stderr,"GaussianMix ERROR: (Probable internal error.)\n"
                   "GaussianMix ERROR: Tried to assign %dd %d-gaussian "
		    "distribution to %dd %d-gaussian distribution\n", 
			gm.dim, gm.ngaussian, dim, ngaussian);
	abort();
	}

  weights = gm.weights;;
  for (int g=0; g < ngaussian; g++)
	gaussians[g] = gm.gaussians[g];
  assert(dim==gm.dim);
  assert(ngaussian==gm.ngaussian);

  return *this;
  }

void GaussianMix::generate(Vector &x)
  {
  int gaussDex = nonunif_int_dev(weights, ngaussian);
  gaussians[gaussDex].generate(x);

  return;
  }

/* test code ---------

int main(void)
  {
  Matrix sigma(4,4);
  double temp[] = { 2, 1, 0, 0, 
                    1, 2, 0, 0,
                    0, 0, 3, 0,
                    0, 0, 0, 1 };

  sigma.initialize_from_arr(temp);

  start_random_number(4372,3511);

  MultivarNormal mn(4);

  mn.setsigmainv(sigma.find_inverse());

  Matrix sum(4,4);
  Vector x(4);
  Matrix xt(1,4);
  for (int ctr=0; ctr < 10000; ctr++)
	{
	mn.generate(x);
	xt = Matrix(x).find_transpose();
	sum += Matrix(x) * xt;
	}
  sum /= 10000.0;

  sum.print();

  return 0;
  }

--------*/
