#ifndef _LREGRESS_H
#define _LREGRESS_H

#include "matrix.h"

// Example of using this code for logistic regression:
// 
// To use the logistic regression on, say, training 
// data 2,2,3->1 and 1,0,5->0 set x and y to:
// 
// Matrix x(2,3);
// Vector y(2);
// x[0][0]=2; x[0][1]=2; x[0][2]=3;
// x[1][0]=1; x[1][1]=0; x[1][2]=5;
// and y[0]=0; y[1]=1;
// 
// and...
// 
// LRegress lr(3);         // 3 dimensional inputs
// lr.fitIt(x, y);         // fit a logistic regression model
// 
// To make predictions on, say, 1,2,3:
// 
// Vector q(3);
// q[0]=1; q[1]=2; q[3]=3;
// printf("Predicted value: %f\n", lr.fittedValue(q));
// 
// Incidentally, most "normal" uses of logistic regression dictate 
// adding a constant term, so that 2,2,3 and 1,0,5 above should be 
// expanded to 2,2,3,1 and 1,0,5,1, but this doesn't matter for our 
// problem of ranking states.

class LRegress
  {
  private:
	int dim;	// dimension of input.
	Vector beta;	// length dim 

  public:
	LRegress(int dim_) : beta(dim_), dim(dim_) {}
	const Vector &betaVal(void) const { return beta; }
	void setBeta(Vector &beta_) { beta_.assert_dim(dim); beta = beta_; }

	int dimVal(void) const { return dim; }
	void fitIt(const Matrix &x, const Vector &y, int suppressWarnings=0,
			double lambda=0.0);
		// lambda is the multiplier into a quadratic regularization 
 		//   term -lambda*beta'*beta to encourage weights near 0.
	void old_fitIt(const Matrix &x, const Vector &y, int suppressWarnings=0);
	double fittedValue(const Vector &x) const;
	int prediction(const Vector &x) const
		{ x.assert_dim(dim); return dot_prod(x,beta) > 0; }

	LRegress &operator= (const LRegress &lr)
		{ if (dim != lr.dim) beta.resize(lr.beta.dim()); 
                  beta = lr.beta; dim = lr.dim; return *this; }

  };

#endif		// _LREGRESS_H

