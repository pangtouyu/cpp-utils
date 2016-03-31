#ifndef _MULTIVARNORMAL_H
#define _MULTIVARNORMAL_H

#include <math.h>

#include "matrix.h"

class Distribution
  {
  public:
	virtual int dimension()=0;

	virtual double logLikelihoodAt(const Vector &x)=0;
	virtual double likelihoodAt(const Vector &x) 
			{ return exp(logLikelihoodAt(x)); }

	virtual void generate(Vector &x)=0;
	virtual ~Distribution() {}
  };

class MultivarNormal : public Distribution
  {
  private:
	Matrix sigmainv;
	Matrix generateT;
	double abs_sigmainvdet;
	Vector mu;
	int dim;

  	Matrix detGrad;
	int detGradValid;

  public:
	double logLikelihoodAt(const Vector &x);
	// double likelihoodAt(Vector &x) { return exp(logLikelihoodAt(x)); }
	int dimension() { return dim; }
	void setmu(const Vector &mu_);
	int setsigmainv(const Matrix &sigmainv_, int warnIfIllCond=1);
	void generate(Vector &x);

	const Matrix &sigmainvVal(void) { return sigmainv; }
	const Matrix &siginvVal(void) {return sigmainv; }
	const Vector &muVal(void) { return mu; }
	double abs_siginvdetVal(void) {return abs_sigmainvdet;}

	void gradLoglikelihood(const Vector &x, Vector &muGrad, Matrix &sigmainvGrad);
	void gradLikelihood(const Vector &x, Vector &muGrad, Matrix &sigmainvGrad);

	MultivarNormal &operator=(const MultivarNormal &mn);
	MultivarNormal(int dim_);
	MultivarNormal(const MultivarNormal &mn);
  };

struct GaussianMix : public Distribution
  {
  private:
	Vector weights;
	int dim;
	int ngaussian;
	MultivarNormal *gaussians;

  public:
	void setWeights(Vector &weights_) 
		{ assert(double_eq(sum_arr(weights, ngaussian),1.0)); 
		    weights=weights_; }
	const Vector &weightsVal(void) { return weights; }
	MultivarNormal &gaussianRef(int gaussDex) 
		{ return gaussians[gaussDex]; }
	int dimension(void) { return dim; }
	int ngaussianVal(void) { return ngaussian; } 
	void generate(Vector &x);
	double logLikelihoodAt(const Vector &x);
	void setEMReassign(Vector &x, Vector &reassign);

	GaussianMix &operator=(const GaussianMix &gm);
	GaussianMix(int dim_, int ngaussian_, int initialize=1);
			// initialize initializes to a uniform
			// mixture of 0-mean, unit cov matrices
	~GaussianMix();
  };

#endif _MULTIVARNORMAL_H

