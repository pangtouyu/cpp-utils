/*************************************************

  boc.C - Bayes Optimal Classifier

  Andrew Y. Ng, 1996

**************************************************/
   

#include <assert.h>
#include <stdio.h>
#include <values.h>

#include "boc.h"
#include "matrix.h"
#include "misc.h"

int get_output_vals(Vector &outputs, float vals[])
  {
  int n_seen=0;
  int seen, ctr, dex;

  for (ctr=0; ctr < outputs.length(); ctr++)
	{
	seen=0;
	for (dex=0; dex < n_seen; dex++)
		{
		if (vals[dex] == outputs[ctr])
			{
			seen=1;
			break;
			}
		}
	if (!seen)
		{
		assert(n_seen < BOC_MAX_N_CLASSES);
		vals[n_seen] = outputs[ctr];
		n_seen++;
		}
	}

  return n_seen;
  }

int class_dex(float val, float output_vals[], int n_classes)
  {
  int ctr;

  for (ctr=0; ctr < n_classes; ctr++)
	if (val == output_vals[ctr])
		return ctr;

  error("Bad output value. (Probable internal error.)");
  return 0;	/* to get rid of gcc warning */
  }

Matrix find_mu(Matrix &inputs, Vector &outputs, 
		float output_vals[], int n_classes)
  {
  int n_in_class[BOC_MAX_N_CLASSES] = {0};
  int ctr, dex, this_class;

  assert(inputs.nrow() == outputs.length());

  Matrix mu(n_classes, inputs.ncol());

  for (ctr=0; ctr < inputs.nrow(); ctr++)
	{
	this_class=class_dex(outputs[ctr], output_vals, n_classes);
	assert(this_class >= 0 && this_class < n_classes);
	n_in_class[this_class]++;
	for (dex=0; dex < inputs.ncol(); dex++)
		mu[this_class][dex] += inputs[ctr][dex];
	}

  for (ctr=0; ctr < n_classes; ctr++)
	{
	if (n_in_class[ctr] == 0)
		internal_error("Bad. Class with no elements??");
	for (dex=0; dex < inputs.ncol(); dex++)
		mu[ctr][dex] /= n_in_class[ctr];
	}

  return mu;
  }

Matrix *find_weak_cov(Matrix &inputs, Vector &outputs,
			float output_vals[], int *n_classes_p,
			Matrix *mu_p)
  {
  int n_classes;
  Matrix *inv_cov;
  int ctr, dex, this_class;

  int n_in_class[BOC_MAX_N_CLASSES] = {0};
  int i,j;

  assert(inputs.ncol() > 0);
  n_classes = get_output_vals(outputs, output_vals);
  inv_cov = new Matrix[n_classes](inputs.ncol(),inputs.ncol());
  assert(inv_cov != NULL);
  Matrix mu = find_mu(inputs, outputs, output_vals, n_classes);
		// i,j-element is mean of jth component of input for class i 

  Vector diff_vec(inputs.ncol());
  Vector input_vec(inputs.ncol());
  Vector mu_vec(inputs.ncol());
  Matrix temp_mat(inputs.ncol(), inputs.ncol());

  for (ctr=0; ctr < inputs.nrow(); ctr++)
	{    
	this_class=class_dex(outputs[ctr], output_vals, n_classes);
	mu_vec = mu.extract_row(this_class);
	input_vec = inputs.extract_row(ctr);

	diff_vec = input_vec - mu_vec;	

	for (i=0; i < inputs.ncol(); i++)
	   for (j=0; j < inputs.ncol(); j++)
		temp_mat[i][j] = diff_vec[i] * diff_vec[j];

	inv_cov[this_class] += temp_mat;
	n_in_class[this_class]++;
	}

  for (ctr=0; ctr < n_classes; ctr++)
	{
	assert(n_in_class[this_class] != 0);
	inv_cov[ctr] *= 1.0/n_in_class[ctr];
  	inv_cov[ctr].invert_me();
	}

  *n_classes_p = n_classes;
  mu_p->resize(n_classes, inputs.ncol());
  *mu_p = mu;

  return inv_cov;
  }

Matrix *find_med_cov(Matrix &inputs, Vector &outputs,
			float output_vals[], int *n_classes_p,
			Matrix *mu_p)
  {
  int n_classes;
  int ctr, dex, this_class;
  Matrix *inv_cov;

  int i,j;

  assert(inputs.ncol() > 0);
  n_classes = get_output_vals(outputs, output_vals);
  inv_cov = new Matrix[n_classes](inputs.ncol(),inputs.ncol());
  assert(inv_cov != NULL);
  Matrix mu = find_mu(inputs, outputs, output_vals, n_classes);
		// i,j-element is mean of jth component of input for class i 

  Matrix cov(inputs.ncol(),inputs.ncol());

  Vector diff_vec(inputs.ncol());
  Vector input_vec(inputs.ncol());
  Vector mu_vec(inputs.ncol());
  Matrix temp_mat(inputs.ncol(), inputs.ncol());

  for (ctr=0; ctr < inputs.nrow(); ctr++)
	{    
	this_class=class_dex(outputs[ctr], output_vals, n_classes);
	mu_vec = mu.extract_row(this_class);
	input_vec = inputs.extract_row(ctr);

	diff_vec = input_vec - mu_vec;	

	for (i=0; i < inputs.ncol(); i++)
	   for (j=0; j < inputs.ncol(); j++)
		temp_mat[i][j] = diff_vec[i] * diff_vec[j];

	cov += temp_mat;
	}

  cov *= 1.0/inputs.nrow();

  inv_cov[0] = cov.find_inverse();

  for (ctr=1; ctr < n_classes; ctr++)
	inv_cov[ctr] = inv_cov[0];

  *n_classes_p = n_classes;
  mu_p->resize(n_classes, inputs.ncol());
  *mu_p = mu;

  return inv_cov;
  }

Matrix *find_strong_cov(Matrix &inputs, Vector &outputs,
			float output_vals[], int *n_classes_p,
			Matrix *mu_p)
  {
  int n_classes;
  int ctr, dex, this_class;
  Matrix *inv_cov;

  float sigma_sq;

  assert(inputs.ncol() > 0);
  n_classes = get_output_vals(outputs, output_vals);
  inv_cov = new Matrix[n_classes](inputs.ncol(),inputs.ncol());
  assert(inv_cov != NULL);
  Matrix mu = find_mu(inputs, outputs, output_vals, n_classes);
		// i,j-element is mean of jth component of input for class i 

  sigma_sq = 0;

  for (ctr=0; ctr < inputs.nrow(); ctr++)
     for (dex=0; dex < inputs.ncol(); dex++)
	{
	this_class=class_dex(outputs[ctr], output_vals, n_classes);
	sigma_sq += sqr(inputs[ctr][dex] - mu[this_class][dex]);
	}

  sigma_sq /= inputs.nrow() * n_classes;

  Matrix temp_mat(inputs.ncol(), inputs.ncol());
  for (ctr=0; ctr < inputs.ncol(); ctr++)
	temp_mat[ctr][ctr] = sigma_sq;

  inv_cov[0] = temp_mat.find_inverse();
  for (ctr=1; ctr < n_classes; ctr++)
	inv_cov[ctr] = inv_cov[0];

  *n_classes_p = n_classes;
  mu_p->resize(n_classes, inputs.ncol());
  *mu_p = mu;

  return inv_cov;
  }

float classify_it(Matrix *inv_cov, Vector &test_in, Matrix &mu,
			float output_vals[], int n_classes)
  {
  float best_dist, this_dist;
  int ctr, dex;
  int best_dex;

  Vector vect1(inv_cov[0].nrow());
  Vector vect2(inv_cov[0].nrow());
  Vector mu_vect(inv_cov[0].nrow());

  best_dist = MAXFLOAT;
  for (ctr=0; ctr < n_classes; ctr++)
	{
	mu_vect = mu.extract_row(ctr);
	vect2 = test_in - mu_vect;
	vect1 = inv_cov[ctr]*vect2;

	this_dist = dot_prod(vect1, vect2);

	if (this_dist < best_dist)
		{
		best_dist = this_dist;
		best_dex = ctr;
		}
	}

  return output_vals[best_dex];
  }

// If we have m samples and n-dimentional input, then 
// train_inputs is m by n, and contains one sample's inputs in each row.
// Ditto test_inputs, and the predicted classifications for it
// will be put in predictions
void do_predictions(Matrix &train_inputs, Vector &train_outputs,
		     Matrix &test_inputs, 
			int assumption, Vector &predictions)
  {
  double train_time, test_time;
  float output_values[BOC_MAX_N_CLASSES];
  int n_classes;
  Matrix *inv_cov;	// inverses of the co-variance matrixes 
  int ctr;
  float this_predict;
  Vector test_in(test_inputs.ncol());
  Matrix mu(1,1);	// will be resized later

  assert(predictions.length() == test_inputs.nrow());

  /* training */
  switch (assumption)
	{
	case WEAK_ASSUMPTION:   inv_cov = find_weak_cov(train_inputs, train_outputs,
					output_values, &n_classes, &mu); break;
	case MED_ASSUMPTION:    inv_cov = find_med_cov(train_inputs, train_outputs,
					output_values, &n_classes, &mu); break;
	case STRONG_ASSUMPTION: inv_cov = find_strong_cov(train_inputs, train_outputs,
					output_values, &n_classes, &mu); break;
	default: internal_error("bad assumption type.");
	}

  // make and save the predictions 
  for (ctr=0; ctr < test_inputs.nrow(); ctr++)
	{
	test_in = test_inputs.extract_row(ctr);
	this_predict = classify_it(inv_cov, test_in, mu, 
				output_values, n_classes);
	predictions[ctr] = this_predict;
	}

  delete[] inv_cov;

  return;
  }

// This (test code) is specific to train.ncs and test.ncs,
// which are respectively 600 and 168 lines line, and have
// 8 inputs and 1 output.
/*
int main(void)
  {
  FILE *fp;
  int ctr, dex, err;
  double d;
  int assumption;
  Matrix train_inputs(600,8);
  Matrix test_inputs(168,8);
  Vector train_outputs(600);
  Vector test_outputs(168);
  Vector predictions(test_inputs.nrow());

  // read in data
  fp = safe_fopen("train.ncs","rt");
  for (ctr=0; ctr < 600; ctr++)
    {
    for (dex=0; dex < 8; dex++)
	{
	err=fscanf(fp, "%lf", &train_inputs[ctr][dex]); 
	assert(err==1);
	}
    err=fscanf(fp, "%lf", &train_outputs[ctr]); 
    assert(err==1);
    }
  err=fscanf(fp, "%lf", &d);
  assert(err!=1);
  fclose(fp);

  fp = safe_fopen("test.ncs","rt");
  for (ctr=0; ctr < 168; ctr++)
    {
    for (dex=0; dex < 8; dex++)
	{
	err=fscanf(fp, "%lf", &test_inputs[ctr][dex]); 
	assert(err==1);
	}
    err=fscanf(fp, "%lf", &test_outputs[ctr]); 
    assert(err==1);
    }
  err=fscanf(fp, "%lf", &d);
  assert(err!=1);
  fclose(fp);

  // run the algorithm
  assumption = STRONG_ASSUMPTION;
  do_predictions(train_inputs, train_outputs, 
		test_inputs, assumption, predictions);

  fp = safe_fopen("treslike.txt", "wt");
  for (ctr=0; ctr < 168; ctr++)
	fprintf(fp, "%d %f %f\n", ctr, test_outputs[ctr], predictions[ctr]);
  fclose(fp);

  return 0;
  }
*/
