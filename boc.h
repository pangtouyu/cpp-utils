#ifndef _BOC_H
#define _BOC_H

#include "matrix.h"

enum {WEAK_ASSUMPTION, MED_ASSUMPTION, STRONG_ASSUMPTION};
#define BOC_MAX_N_CLASSES 512

// If we have m samples and n-dimentional input, then 
// train_inputs is m by n, and contains one sample's inputs in each row.
// Ditto test_inputs, and the predicted classifications for it
// will be put in predictions
void do_predictions(Matrix &train_inputs, Vector &train_outputs,
		     Matrix &test_inputs, 
			int assumption, Vector &predictions);

#endif	// _BOC_H
