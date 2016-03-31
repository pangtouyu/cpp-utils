#ifndef _RANDVAR_H
#define _RANDVAR_H

#include "matrix.h"

const double DISCRETE_RV_INDEX_DENSITY = 2.0;

class DiscreteRV 
  {
  private:
	Vector probs;		// vector of probabilities
	Vector inc_probs;	// will be inc_probs[i] = sum_{j=0}^i probs[j]
	int *indices;
	int n_indices;

	void init_indices(void);

  public:
	int generate(void) const;
	DiscreteRV(const Vector &probs_);
	~DiscreteRV();
	double likelihood_at(int x) const 
		{ assert(x>=0&&x<probs.length()); return probs[x]; }
	int dimension(void) const { return probs.length(); }
  };

#endif 		// _RANDVAR_H

