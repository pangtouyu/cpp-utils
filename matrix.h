/*********************************************

  matrix.h - header for matrix.C, which is code
	  for doing Matrices and Vectors. 

  Andrew Y. Ng, 1996

*********************************************/

#ifndef _MATRIX_H
#define _MATRIX_H

#include <assert.h>
#include "matrix.h"
#include "misc.h"

class Matrix;
void mat_error(const char*msg);

class Vector
  {
  private:
	double *val;
	int n;

  public:
	inline double l2Norm(void) const { return sqrt(l2NormSquared()); }
	inline double l2NormSquared(void) const 
	    { double sum=0; for (int i=0;i<n;i++) sum += sqr(val[i]); return sum; }
	inline int dim(void) const {return n;}
	inline int n_ele(void) const {return n;}
	inline int length(void) const {return n;}
	inline double min(void) const { return min_in_arr(val,n); }
	inline double max(void) const { return max_in_arr(val,n); }
	inline double minFabsValue(void) const { return min_fabs_in_arr(val,n); }
	inline double maxFabsValue(void) const { return max_fabs_in_arr(val,n); }
	inline int argmin(void) const { return ::argmin(val,n); }
	inline int argmax(void) const { return ::argmax(val,n); }
	inline operator double*() const {return val;}
	inline operator double() const 
		{if (n!=1) mat_error("Cannot convert non-1d vector to double");
		 return *val;}
	inline double *ref(int dex) {return &val[dex];}
	inline void dangerousCoerceLength(int newlen) { n = newlen; }
	void resize(int n_, int initialize=1, double init_val=0);
			// destroys previous data
	void assert_dim(int n_) const;
	operator Matrix() const;
	void initialize_from_arr(const double *p);
	void initialize_from_file(FILE *fp);
	void initialize_from_file(const char *fn);
	Vector subVector(int first, int lastp1) const;  // sub-vector including
						  // val[first]..val[lastp1-1]
	void setSubVector(const Vector &v, int first, int lastp1);  
	double sumSubVector(int first, int lastp1) const;  // sums val[first]..val[lastp1-1]
	inline double sum(void) const { return sum_arr(val, n); }
	inline double mean(void) const { return sum()/n; }
	inline double sumSqr(void) const { double sum = 0.0; 
					   for (int i=0; i < n; i++) sum += sqr(val[i]);
					   return sum; }
	void normalize_sum_to(double desiredSum);

	void enlarge(int new_size, int init_new_ele=1, double init_val=0.0);

	Matrix find_transpose(void) const;
	Matrix T(void) const;
	Vector &operator+=(const Vector &vec);
	Vector &operator-=(const Vector &vec);
	Vector &operator*=(double d);
	Vector &operator/=(double d);

	void make_zero(void) {for(int ctr=0; ctr < n; ctr++)
				val[ctr]=0;}
	void randomize(double lowLim, double highLim);
	inline int isZero(void) const {
			for (int i=0; i < n; i++) if (val[i]!=0) return 0;
			return 1; }
	void make_const_vec(double x) {for(int ctr=0; ctr < n; ctr++)
						val[ctr]=x;}
	Vector &operator=(const Vector &vec);
	void print(FILE *fp=(FILE*)NULL, const char *separator=(const char *)NULL, 
				int print_cr=1, const char *printStr=NULL) const ;
	Vector(const Vector &vec);
	Vector(int n_, int initialize=1, double init_val=0);
	Vector(const char *fn);
	~Vector();
  };

// mat[i][j] is the element in the i-th row and j-th column,
//  using zero-indexing
class Matrix	
  {
  private:
	double **val;	// elements are all stored in one 
			// contiguous linear array
	int m, n;

	void init_structure(int m_, int n_);
	double **new_one_indices(void) const;

  public:
	inline int nrow(void) const {return m;} 
	inline int nrows(void) const {return m;} 
	inline int ncol(void) const {return n;}
	inline int ncols(void) const {return n;}
	inline double *ref(int row, int col) {return &val[row][col];}
	double setval(int row, int col, double d);
	void assert_dim(int m_, int n_) const;
	inline int is_square(void) const {return m==n;}
	int is_symmetric(double epsilon=1.0e-6, double delta=1.0e-10) const;
	inline operator double**() const {return val;}
	inline operator double() const 
		{if (m!=1||n!=1) 
		     mat_error("Cannot convert non-1x1 matrix to double");
		 return **val;}
	inline double sum(void) const
	  { double s=0; 
	    for (int i=0; i < m; i++) for (int j=0; j < n; j++) s+=val[i][j];
	    return s; }
	inline double trace(void) const 
		{ assert(is_square()); double sum=0.0; 
		   for (int i=0; i < m; i++) sum += val[i][i]; 
		   return sum; }
	inline double min(void) const 
	  	{ double minval=val[0][0]; for (int i=0; i < m; i++) for (int j=0; j < n; j++)
			minval = ::min(minval, val[i][j]); return minval; }
	inline double max(void) const
	  	{ double maxval=val[0][0]; for (int i=0; i < m; i++) for (int j=0; j < n; j++)
			maxval = ::max(maxval, val[i][j]); return maxval; }
	inline double maxfabs(void) const 
	  	{ double maxval=fabs(val[0][0]); for (int i=0; i < m; i++) for (int j=0; j < n; j++)
			maxval = ::max(maxval, fabs(val[i][j])); return maxval; }

	Matrix &operator+=(const Matrix &mat);
	Matrix &operator-=(const Matrix &mat);
	Matrix &operator*=(double d);
	Matrix &operator/=(double d);

	void transpose_me(void);

	void svd(Matrix &u, Vector &diag, Matrix &v, int *failed_p=(int*)NULL) const;
	Vector solveLeastSquares(const Vector &y, int *illCond=(int*)NULL) const;

	int svd_decomp(Matrix &u, Vector &diag, Matrix &v, int *failed_p=(int*)NULL) const;
	Vector linSolve(const Vector &y, int canDestroyMe, int printDebug=1);
	int invert_me(double *abs_detp=(double*)NULL, int *failed_p=(int*)NULL); // (Uses SVD)
			// returns 1 if ill-conditioned or singular. 
			// if abs_detp is non-NULL, then puts the 
			// ABSOLUTE VALUE of the determinent of the 
			// ORIGINAL matrix there.
	int old_invert_me(double *abs_detp=(double*)NULL, int *failed_p=(int*)NULL); // (Uses SVD)
	Matrix find_inverse(int *ill_conditioned=(int*)NULL, 
				double *abs_detp=(double*)NULL, 
				int *failed_p=(int*)NULL) const;
			// similar behavior to invert_me(..)
	Matrix find_transpose(void) const;
	inline Matrix T(void) const { return find_transpose(); };
	Matrix find_AtransposeWA(const Vector &w) const; // w=diag matrix
	Matrix find_AWAtranspose(const Vector &w) const; // w=diag matrix
	Matrix find_AtransposeA(void) const;
	Matrix find_AAtranspose(void) const;
	Vector find_AtransposeWy(const Vector &w, const Vector &y) const; // w=diag matrix
	Matrix find_mat_sqrt(void) const; // precond: *this is SYMMETRIC
			// returns A, such that A-transpose * A == *this

	void make_zero(void) {for(int i=0; i < nrow(); i++)
				for (int j=0; j < ncol(); j++)
					val[i][j] = 0; }
	void make_identity(void) {assert(nrow()==ncol()); make_zero();
				      for(int i=0; i < nrow(); i++)
					val[i][i] = 1.0; }
	void randomize(double lowLim, double highLim);
	void resize(int m_, int n_, int initialize=1, double init_val=0);
			// destroys all old elements,
	inline void resize_and_assign_to(const Matrix &m) 
		{resize(m.nrow(),m.ncol(),0); *this = m; }

	void initialize_from_arr(const double *p);
	void initialize_from_file(FILE *fp);	
	void initialize_from_file(const char *fn);
	void print(FILE *fp=(FILE*)NULL, const char *separator=(const char*)NULL, 
				const char *printStr=(const char*)NULL) const;
	Vector extract_row(int row) const ;
	Vector extract_col(int col) const ;
	void extract_row(Vector &v, int row) const ;
	void extract_col(Vector &v, int col) const ;
	void set_row(const double dat[], int row);
	void set_col(const double dat[], int col);
	void set_row(const Vector &v, int row);
	void set_col(const Vector &v, int col);
	Matrix &operator=(const Matrix &mat);

	Matrix find_det_deriv(void) const;	// i,j-th element is the
					// deriv of the determinant wrt
					// i,j element of the orig matrix
	double determinant(void);

	// void loadFromFile(const char *fn, int resizeOkay=1);

	Matrix(const Matrix &mat);
	Matrix(int m_, int n_, int initialize=1, double init_val=0);
	Matrix(const char *fn);
	~Matrix();
  };

Vector operator*(double d, const Vector &vec);
Vector operator*(const Vector &vec, double d);
Matrix operator*(const Vector &vec, const Matrix &mat);
Vector operator-(const Vector &vec);
Vector operator-(const Vector &vec1, const Vector &vec2);
Vector operator+(const Vector &vec1, const Vector &vec2);
Vector operator/(const Vector &vec, double d);

Matrix operator*(const Matrix &mat1, const Matrix &mat2);
Matrix operator*(double d, const Matrix &mat);
Matrix operator*(const Matrix &mat, double d);
Matrix operator/(const Matrix &mat, double d);
Matrix operator-(const Matrix &mat);
Matrix operator-(const Matrix &mat1, const Matrix &mat2);

Matrix operator+(const Matrix &mat1, const Matrix &mat2);

Vector operator*(const Matrix &mat, const Vector &vec);
Matrix operator*(const Matrix &mat1, const Matrix &mat2);

double dot_prod(const Vector &vec1, const Vector &vec2);

// these return true if all elements are withing 
// delta of each other in absolute value, or within 1+epsilon
// in relative value
int mat_approx_equal(const Matrix &m1, const Matrix &m2, 
		double epsilon=1.0e-6, double delta=1.0e-10);
int vect_approx_equal(const Vector &v1, const Vector &v2, 
		double epsilon=1.0e-6, double delta=1.0e-10);

double euclideanDistance(const Vector &v1, const Vector &v2);
double squaredEuclideanDistance(const Vector &v1, const Vector &v2);

Vector asVector(const Matrix &m);
inline Vector constantVector(int dim, double val) { Vector v(dim,1,val); return v; }

#endif		// _MATRIX_H
