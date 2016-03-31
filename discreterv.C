// taken from randvar.[Ch] on ~ang/uai0 on AT&T Labs system

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <values.h>

#include "discreterv.h"
#include "matrix.h"
#include "misc.h"
#include "Random_Number.h"

// Given a vector of probabilities p (that sum to 1), on each
// generate(), DiscreteRV returns i with probability p[i].

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
  n_indices = (int)ceil(DISCRETE_RV_INDEX_DENSITY * probs.length());
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

DiscreteRV::DiscreteRV(const Vector &probs_) 
		: probs(probs_), inc_probs(probs.length()), indices((int*)NULL)
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

