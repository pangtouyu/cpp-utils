#ifndef _DEVIATE_H
#define _DEVIATE_H

#include "matrix.h"

// parent-class for general (multivariate) distributions
class Randvar
  {
  public:
	virtual ~Randvar() {};
	virtual double likelihood_at(Vector &x)=0;
	virtual dimension(void) = 0;	// dimensionality of samples
  };

// multivariate normal
class MultiNormal : public Randvar
  {
  private:
	int dim;	// dim-dimensional multivariate normal 
	Matrix sigma;	// dimxdim matrix
	Vector mu;	// dim-dimensional

	Matrix *sigmainv;	// lazy evaled; NULL if invalid
	double *sigmadet;	// lazy evaled; NULL if invalid

	Matrix *generate_t;	// lazy evaled; NULL if invalid
				// this is the matrix T s.t. if x~N(0,I),
				// then Tx~N(0,sigma)

	double find_sigmadet(void);	// these handle the lazy 
					// evaluations
	Matrix *find_sigmainv(void);
	Matrix *find_generate_t(void);

	void invalidate_cached_dat(void); // to be called whenever
					  // anything (mu,sigma) 
					  // is changed.

  public:
	Matrix *sigma_p(void) {return &sigma;}
	Vector *mu_p(void) {return &mu;}
	int dimension(void) {return dim;}
	void generate(Vector &dev);
	Vector generate(void);
	void setparams(const Matrix &newsigma, const Vector &newmu);
	void setmu(const Vector &newmu);
	void setsigma(const Matrix &newsigma);
	double loglikelihood_at(Vector &x);
	double likelihood_at(Vector &x);
	void print(FILE *fp=NULL);
	MultiNormal(int dim, double init_range);
			// initializes to dim*iid N(mu0,1)
			// where mu0 ~ iid Uniform(0,init_range)
	MultiNormal &operator=(const MultiNormal &mn);
	MultiNormal(const MultiNormal &mn);
  };


class DiscreteRV 
  {
  private:
	const double INDEX_DENSITY = 2.0;
	Vector probs;		// vector of probabilities
	Vector inc_probs;	// will be inc_probs[i] = sum_{j=0}^i probs[j]
	int *indices;
	int n_indices;

  public:
	int generate(void) const;
	void init_indices(void);
	DiscreteRV(Vector &probs_);
	~DiscreteRV();
	double likelihood_at(int x) const 
		{ assert(x>=0&&x<probs.length()); return probs[x]; }
	int dimension(void) const { return probs.length(); }
  };

#endif 		// _DEVIATE_H

