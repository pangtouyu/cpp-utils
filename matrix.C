/*********************************************

  matrix.C - code for doing Matrices and 
	  Vectors, with nice overloaded C++
	  syntax.

  Andrew Y. Ng, 1996

*********************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "matrix.h"
#include "misc.h"

extern "C" {
void dgesvd_(char *jobu, char *jobvt, 
            int *m, int *n, double *a, 
            int *lda, double *s, double *u, int *ldu, 
            double *vt, int *ldvt, double *work, int *lwork, 
            int *info);
}


#define SHOW_CREATIONS (0)	// if true, will print a message
				// each time a matrix or a vector
				// is created or destroyed.
				// (For debugging.)

void mat_error(const char*msg);
void mat_warn(const char*msg);

// ----------------------------------------------------------------
// Code from Numerical Recipes, to do Singular Value Decomposition
// and LU-decomposition.
// (Changed all floats converted to doubles)

#define SVD_MAXIT 200

static double *vector(int nl, int nh)
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
//	v= new double[nh-nl+1];
	if (!v) mat_error("allocation failure in vector()");
	return v-nl;
}

static void free_vector(double *v, int nl, int nh)
{
	free((char*) (v+nl));
//	delete[] (v+nl);
}

// Numerical Recipes File SVDCMP.C

// if failed_p is non-NULL, the SVD failing to converge will return
// with *failed_p=1. (Otherwise, the program terminates)

static double at,bt,ct;
#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? \
(ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))

static double maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
	(maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void svdcmp(double **a, int m, int n, double *w, double **v, int *failed_p = (int*)NULL)
{
	int flag,i,its,j,jj,k,l,nm;
	double c,f,h,s,x,y,z;
	double anorm=0.0,g=0.0,scale=0.0;
	double *rv1;

 	if (failed_p != NULL)
		*failed_p = 0;

	if (m < n) mat_error("SVDCMP: (Probable internal error) "
				"You must augment A with extra zero rows");
	rv1=vector(1,n);
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				if (i != n) {
					for (j=l;j<=n;j++) {
						for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
						f=s/h;
						for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
					}
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				if (i != m) {
					for (j=l;j<=m;j++) {
						for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
						for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
					}
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=n;i>=1;i--) {
		l=i+1;
		g=w[i];
		if (i < n)
			for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			if (i != n) {
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
					f=(s/a[i][i])*g;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else {
			for (j=i;j<=m;j++) a[j][i]=0.0;
		}
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=SVD_MAXIT;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=PYTHAG(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s=(-f*h);
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k]=(-v[j][k]);
				}
				break;
			}
//			if (its == 60) 
//				{
//				fprintf(stderr, "Matrix/Vector WARNING: %dx%d matrix, "
//				    "No convergence in 60 SVDCMP iterations\n", m, n);
//				}
			if (its == 200) 
				{
				fprintf(stderr, "Matrix/Vector WARNING: %dx%d matrix, "
				    "No convergence in 200 SVDCMP iterations\n", m, n);
				}
			if (its == SVD_MAXIT) 
			    {
			    // for (int ctr1=1; ctr1 <=n; ctr1++)
			    // 	{
			    //	for (int ctr2=1; ctr2 <= n; ctr2++)
			    //	   printf("%.3f ", a[ctr1][ctr2]);
			    //	printf("\n");
			    //	}
			    // foo->print();
			    printf("-----\n");
			    // void bar(void);
			    // bar();
			    if (failed_p != NULL)
				{
				fprintf(stderr,"Matrix/Vector WARNING: No convergence in SVD_MAXIT SVDCMP iterations. "
							"(%dx%d)\n", m,n);
				*failed_p = 1;
				free_vector(rv1,1,n);
				return;
				}
			    fprintf(stderr, "ERROR: (%dx%d matrix)\n", m, n);
			    mat_error("No convergence in SVD_MAXIT SVDCMP iterations");
			    }
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=PYTHAG(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=PYTHAG(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y=y*c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=PYTHAG(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=(c*g)+(s*y);
				x=(c*y)-(s*g);
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_vector(rv1,1,n);
}

#undef SIGN
#undef MAX
#undef PYTHAG

// Numerical Recipes File LUDCMP.c

// if !errorIfSingular, then returns with *d =0 on singularity/failure

#define TINY 1.0e-20;
void ludcmp(double **a, int n, int *indx, double *d, int errorIfSingular=1)
{
	int i,imax,j,k;
	float big,dum,sum,temp;
	double *vv;

//	vv=vector(1,n);
//for (i=1; i <=n; i++) for (j=1; j <=n; j++) a[i][j]++;
//	free_vector(vv,1,n);
//for (i=1; i <=n; i++) for (j=1; j <=n; j++) a[i][j]--;
	vv=vector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) 
			{
			if (errorIfSingular)
				error("Singular matrix in routine LUDCMP");
			else
				*d = 0.0;
			return;
			}
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_vector(vv,1,n);
}

#undef TINY

// End of Numerical Recipes Code
// ----------------------------------------------------------------

void mat_error(const char*msg)
  {
  fflush(stdin);
  fprintf(stderr,"\nMatrix/Vector ERROR: %s\n\n", msg);
  fflush(stderr);

  abort();
//  exit(-1);
  }

void mat_warn(const char*msg)
  {
  fflush(stdin);
  fprintf(stderr,"\nMatrix/Vector WARNING: %s\n\n", msg);
  fflush(stderr);

  return;
  }

int Matrix::is_symmetric(double epsilon, double delta) const
  {
  assert(is_square());

  for (int i=0; i < n; i++)
    for (int j=0; j < i; j++)
	if (!double_eq(val[i][j], val[j][i], epsilon, delta))
		return 0;

  return 1;
  }

// precond: *this is SYMMETRIC and POSITIVE DEFINITE
// finds A, such that A-transpose * A == *this
Matrix Matrix::find_mat_sqrt(void) const 
  {
  assert(is_symmetric());

  Matrix mat = *this;
  Vector diag(m); 
  Matrix v(m,m);

  svd_decomp(mat, diag, v);

  for (int i=0; i < diag.dim(); i++)
	{
	if (diag[i] < 0)
	    mat_error("find_mat_sqrt: matrix not pos definite");
	diag[i] = sqrt(diag[i]);
	}

  for (int i=0; i < m; i++)
    for (int j=0; j < n; j++)
	v[i][j] *= diag[i];

  // check solution and verify *this was pos definite
  // (Hmm... some of this might be redundant/impossible because
  //  of the "diag[i]<0" test earlier...)
  Matrix testMatrix = v.find_AtransposeA();
  if (!mat_approx_equal(testMatrix, *this))
     {
     if (mat_approx_equal(-testMatrix, *this))
	mat_error("find_mat_sqrt: matrix was neg definite");
     else
	mat_error("find_mat_sqrt: matrix not pos definite");
				// at least, I think so....
     }

  return v;
  }

// returns u, diag, and vt, so that
// (*this) = u * (diag(diag)) * vt.
//
// If (*this) is m by n, 
//
// failed_p can be NULL, but if it is, then this will abort execution
// on a failed SVD.
void Matrix::svd(Matrix &u, Vector &diag, Matrix &vt, int *failed_p) const
  {
  int dex=0;
  double *fortran_a = new double[nrows()*ncol()];

  u.assert_dim(m,m);
  diag.assert_dim(::min(m,n));
  vt.assert_dim(n,n);
  
  for (int j=0; j < ncol(); j++)
    for (int i=0; i < nrow(); i++)
      fortran_a[dex++] = val[i][j];

  char *jobu = "A";
  char *jobvt = "A";
  int thisM = nrow(); assert(thisM==m);
  int thisN = ncol(); assert(thisN==n);
  int lda = m;
  double *fortran_u = new double[m*m];
  int ldu = m;
  double *fortran_vt = new double[n*n];
  int ldvt = n;
  const int MIN_WORKSIZE = ::max(3*::min(m,n)+::max(m,n), 5*::min(m,n));
  int lwork = 5*MIN_WORKSIZE;
  double work[lwork];
  int info; 

  dgesvd_(jobu, jobvt, &thisM, &thisN, fortran_a, &lda, &(diag[0]), 
		fortran_u, &ldu, fortran_vt, &ldvt, work, &lwork, &info);
  
  dex=0; for (int j=0; j < m; j++) for (int i=0; i < m; i++) u[i][j] = fortran_u[dex++]; 
  dex=0; for (int j=0; j < n; j++) for (int i=0; i < n; i++) vt[i][j] = fortran_vt[dex++]; 

  if (info > 0)
    {
    if (failed_p == NULL)
      error("Matrix::svd failed.");
    else
      *failed_p = 1;
    }
  else if (info < 0)
    internal_error("Matrix::svd: info < 0");
  else
    {
    if (failed_p != NULL)
      *failed_p = 0;
    }

  delete[] fortran_a; delete[] fortran_u; delete[] fortran_vt; 

  return;
  }

// NOTE: As of 6/5/01, this is deprecated.  Use Matrix::svd instead.
// 
// returns 1 if ill-conditioned.
// SVD will make (*this) = U * diag-matrix * Vt,
// (where diag-matrix is the matrix with diag as its diagonal elements.)
// Note that this is DIFFERENT FROM NUMERICAL RECIPES's call, in that 
// it returns Vt instead of V (I think).
// ayn 060501: Note, changed some of the variable names from v to vt. 
//
// if failed_p is non-NULL, the SVD failing to converge will return
// with *failed_p=1. (Otherwise, the program terminates)
int Matrix::svd_decomp(Matrix &u, Vector &diag, Matrix &vt, int *failed_p) const
  {
  warn_once("This program uses svd_decomp, which is deprecated.");

  const double TOL = 1.0e-7;
//  const double TOL = 1.0e-10;
  int ctr;
  double max_diag;
  double **u_p, **vt_p;	// one-indexed versions of this and Vt
  int ill_conditioned;
  // double abs_det;

  if (failed_p != NULL)
	*failed_p = 0;

  assert(m >= n);	// needed for SVD. (I think... else must check
			// V's dimensions)

  // check dimensions of passed in u, diag, and vt
  u.assert_dim(m,n);	
  diag.assert_dim(n);
  vt.assert_dim(n,n);
  
  u_p = u.new_one_indices();
  vt_p = vt.new_one_indices();

  u = *this;

  svdcmp(u_p, m, n, ((double*)diag)-1, vt_p, failed_p);
		// the -1 is because NR expects a 1-indexed vector

  if (failed_p != NULL && *failed_p)
	return 1;

  vt.transpose_me();	// NRIC returns V, We want Vt.   [ayn: 060501]

  // check if ill-conditioned 
  max_diag = max_in_arr(diag, diag.length());
  ill_conditioned = 0;
  for (ctr=0; ctr < diag.length(); ctr++)
	if (diag[ctr] <= TOL*max_diag)
		ill_conditioned = 1;

  // for debugging
//  printf("U: \n");
//  u.print();
//  printf("-----------------------\n");
//  printf("diag: \n");
//  diag.print();
//  printf("-----------------------\n");
//  printf("v: \n");
//  v.print();

  delete[] u_p;
  delete[] vt_p;

  return ill_conditioned;
  }

// returns 1 if either ill-conditioned or singular.
// notation from Numerical Recipes/
//
// if abs_detp is non-NULL, then puts the ABSOLUTE VALUE
// of the determinent of the ORIGINAL matrix there.
// (Since if you've done an SVD already, then it's cheap
// to find that.)
//
// if failed_p is non-NULL, the SVD failing to converge will return
// with *failed_p=1. (Otherwise, the program terminates)
int Matrix::invert_me(double *abs_detp, int *failed_p)
  {
  const double TOL = 1e-8;
//  const double TOL = 1e-6;

  if (!is_square())
	mat_error("(Possible internal error.) Tried to invert non-square matrix.");

  Matrix u(m,m);
  Vector diag(m);
  Matrix vt(m,m);

  svd(u, diag, vt, failed_p);

  // absolute value of the determinant is the product of the singular values 
  double abs_det=1.0;
  for (int dex=0; dex < m; dex++)
    abs_det *= diag[dex];

  // invert all of diag's elements, unless zero, in which case make 0
  double max_diag = max_in_arr(diag, diag.length());
  int ill_conditioned=0;
  for (int ctr=0; ctr < diag.length(); ctr++)
    if (diag[ctr] < TOL*max_diag || diag[ctr] == 0)
      { diag[ctr] = 0.0; ill_conditioned = 1; }
    else
      diag[ctr] = 1.0/diag[ctr];

  vt.transpose_me();
  // do V * diagonal-matrix
  for (int i=0; i < m; i++)
    for (int j=0; j < m; j++)
	vt[i][j] *= diag[j];

  u.transpose_me();
  // finally, multiply by U-transpose, and we're done!
  *this = vt * u;

  if (abs_detp != NULL)
	*abs_detp = abs_det;

  return ill_conditioned;
  }

// returns 1 if either ill-conditioned or singular.
// notation from Numerical Recipes/
//
// if abs_detp is non-NULL, then puts the ABSOLUTE VALUE
// of the determinent of the ORIGINAL matrix there.
// (Since if you've done an SVD already, then it's cheap
// to find that.)
//
// if failed_p is non-NULL, the SVD failing to converge will return
// with *failed_p=1. (Otherwise, the program terminates)
int Matrix::old_invert_me(double *abs_detp, int *failed_p)
  {
  if (failed_p != NULL)
	*failed_p = 0;

  if (!is_square())
	mat_error("(Possible internal error.) Tried to invert non-square matrix.");

  const double TOL = 1.0e-7;
  int ctr, i, j;
  double max_diag;
  double **u_p, **v_p;	// one-indexed versions of this and V
  Matrix v(n,n,0);
  Vector diag(n,0);	// the non-zero elements of a diagonal matrix
  int ill_conditioned;
  double abs_det;

  u_p = new_one_indices();
  v_p = v.new_one_indices();

  svdcmp(u_p, m, n, ((double*)diag)-1, v_p, failed_p);
		// the -1 is because NR expects a 1-indexed vector
		// Note this replaces the current matrix with U

  if (failed_p != NULL && *failed_p)
	{
	delete[] u_p;
	delete[] v_p;
	return 1;
	}

  // the result of svdcmp is that
  // (original-matrix) = U * diag-matrix * V-transpose
  // Since U and V are (column) orthogonal, the desired inverse is
  // (V-tranpose)^{-1} * 1/diag-matrix *  U^{-1}
  //  = V * 1/diag-matrix * U-tranpose
  // (where diag-matrix is the matrix whose diagonal has diag's elements,
  //  and 1/foo above is the element-wise inverse.)

  // absolute value of the determinant is the product of the
  // singular values 
  abs_det = 1.0;
  for (ctr=0; ctr < n; ctr++)
	abs_det *= 
		diag[ctr];

// for debugging:
//  printf("U: \n");
//  print();
//  printf("diag: \n");
//  diag.print();
//  printf("V: \n");
//  v.print();

  transpose_me();	// replace U with U-transpose

  // invert all of diag's elements, unless zero, in which case make 0
  max_diag = max_in_arr(diag, diag.length());
  ill_conditioned=0;
  for (ctr=0; ctr < diag.length(); ctr++)
	if (diag[ctr] <= TOL*max_diag)
		{
//fprintf(stderr, "Ill conditioned: %f %f\n", diag[ctr], max_diag);
		diag[ctr] = 0.0;
		ill_conditioned = 1;
		}
	else
		diag[ctr] = 1.0/diag[ctr];

  // do V * diagonal-matrix
  for (i=0; i < v.nrow(); i++)
    for (j=0; j < v.ncol(); j++)
	v[i][j] *= diag[j];

  // finally, multiply by U-transpose, and we're done!
  *this = v * (*this);

  internal_error("ayn: double-check this routine.");
  // 060501: Shouldn't we have transposed V at some point?

  delete[] u_p;
  delete[] v_p;

  if (abs_detp != NULL)
	*abs_detp = abs_det;

  return ill_conditioned;
  }

// Will solve for x in "Ax=y" in least-squares style.
Vector Matrix::solveLeastSquares(const Vector &y, int *illCond) const 
  {
  const double TOL = 1e-10;

  y.assert_dim(m);

  Matrix u(m,m);
  Vector diag(::min(m,n));
  Matrix vt(n,n);

  svd(u, diag, vt);

  if (illCond != NULL)
    *illCond = 0;

  Vector temp = u.find_transpose() * y;
  double maxDiag = diag.maxFabsValue();
  for (int dex=0; dex < m; dex++)
    {
    if (dex < ::min(m,n))
      {
      double d = diag[dex];
      if (maxDiag > 0 && d/maxDiag < TOL)
        { 
	if (illCond != NULL)
	  *illCond = 1;
	temp[dex] = 0;
	}
      else
        temp[dex] /= d;
      }
    else
      temp[dex] = 0;
    }
  Vector temp2(n);
  for (int dex=0; dex < n; dex++)
    temp2[dex] = temp[dex];

  return vt.find_transpose() * temp2;
  }

// if canDestroyMe, will muck up *this.
//
// Solves, least-squares+SVD-style, for (Vector)beta in (*this)*beta = y, 
// and returns beta.
Vector Matrix::linSolve(const Vector &y, int canDestroyMe, int printDebug) 
  {
  warn_once("Matrix::linSolve is deprecated (because of svd_decomp)");

  if (!canDestroyMe)
	{
	Matrix tempMat = *this;
	return tempMat.linSolve(y,1,printDebug);
	}

  const double TOL = 1.0e-7;
  int ctr;
  double max_diag;
  double **this_p, **v_p;	// one-indexed versions of this and V
  // double abs_det;
  Matrix v(n,n);
  Vector diag(n); 
  Vector beta(n);

  y.assert_dim(m);

  if (m < n)
	mat_error("Matrix::linSolve: (possible internal error) m < n.");

//  this_p = this.new_one_indices();
  this_p = new_one_indices();
  v_p = v.new_one_indices();

  svdcmp(this_p, m, n, ((double*)diag)-1, v_p);
		// the -1 is because NR expects a 1-indexed vector

  // solve it! (NRIC p.676-7)
  beta.make_zero();
  max_diag = max_in_arr(diag, diag.length());
  if (printDebug)
	printf("condition number: %f\n", max_diag/min_in_arr(diag, diag.length()));
  for (ctr=0; ctr < n; ctr++)
	{
	if (diag[ctr] <= TOL*max_diag)
		continue;	// ill-conditioned

	double ui_b=0;
	for (int dex=0; dex < m; dex++)
		ui_b += (*this)[dex][ctr] * y[dex];

	ui_b /= diag[ctr];
	for (int dex=0; dex < n; dex++)
		beta[dex] += v[dex][ctr] * ui_b;
	}

  delete[] this_p;
  delete[] v_p;

  return beta;
  }
  
void Matrix::transpose_me(void)
  {
  int oldm, oldn;
  double **oldval;
  double *temp;
  int i,j;

  oldm=m;
  oldn=n;
  oldval=val;

  m = oldn;
  n = oldm;
  // allocate matrix
  temp = new double[m*n];
  val = new double*[m+1];
  assert(temp != NULL && val != NULL);
  for (i=0; i < m; i++)
	val[i] = temp+(n*i);
  val[m] = (double*)NULL;	// to try to catch some run-off-array bugs

  for (i=0; i < m; i++)
    for (j=0; j < n; j++)
	val[i][j] = oldval[j][i];

  delete[] *oldval;
  delete[] oldval;

  return;
  }

// terribly inefficient implementation at the moment, but nevermind....
double Matrix::determinant(void)
  {
  assert(is_square());
  int dim = m;

// optimizations:
//  if (dim==1)
//	return val[0][0];
// else if (dim==2)
//	return val[0][0]*val[1][1]-val[0][1]*val[1][0];

  Matrix mat = *this;

  double **a = mat.new_one_indices();
  int *indx = new int[dim];
  assert(a != NULL && indx != NULL);

  double d;

  ludcmp(a, dim, indx-1, &d, 0);	// d will be 0 if singular.
					// indx, not indx, since ludcmp 
					// is 1..n, not 0..n-1 indexed
		

  if (d != 0)
    for (int ctr=1; ctr <= dim; ctr++)	// 1..dim, not 0..dim-1
	d *= a[ctr][ctr];

  delete[] a;
  delete[] indx;

  return d;
  }

Matrix Matrix::find_det_deriv(void) const 
  {
  assert(is_square());
  int dim = m;

  assert(dim >= 1);
  if (dim == 1)
	{
	return Matrix(1,1,1,1.0);	// return 1.0
	}

  Matrix mat(dim,dim,0);
  Matrix temp(dim-1,dim-1,0);

  for (int i=0; i < dim; i++)
    for (int j=0; j < dim; j++)
	{
	for (int i0=0; i0 < i; i0++)
	  for (int j0=0; j0 < j; j0++)
		temp[i0][j0] = val[i0][j0];
	for (int i0=i+1; i0 < dim; i0++)
	  for (int j0=0; j0 < j; j0++)
		temp[i0-1][j0] = val[i0][j0];
	for (int i0=0; i0 < i; i0++)
	  for (int j0=j+1; j0 < dim; j0++)
		temp[i0][j0-1] = val[i0][j0];
	for (int i0=i+1; i0 < dim; i0++)
	  for (int j0=j+1; j0 < dim; j0++)
		temp[i0-1][j0-1] = val[i0][j0];

	// temp is now *this, with row i and col j deleted.

        // debugging: check temp.
	for (int i0=0; i0 < dim; i0++)
	  for (int j0=0; j0 < dim; j0++)
	    {
	    if (i0==i || j0==j) 
		continue;
	    int ii = ((i0<i)?i0:i0-1);
	    int jj = ((j0<j)?j0:j0-1);
	    if (temp[ii][jj] != val[i0][j0])
		{
		printf("---- (%d,%d:%d,%d->%d,%d)\n", i,j,i0,j0,ii,jj);
		print();
		printf("----\n");
		temp.print();
		}
	    assert(temp[ii][jj] == val[i0][j0]);
	    }

	if (is_even(i+j))
		mat[i][j] = temp.determinant();
	else
		mat[i][j] = -temp.determinant();
	}
		
  return mat;
  }

Matrix Matrix::find_transpose(void) const 
  {
  Matrix m = *this;

  m.transpose_me();

  return m;
  }

Matrix Matrix::find_AAtranspose(void)  const 
  {
  Matrix m(nrow(), nrow());
  int dex1, dex2, ctr;
  double sum;

  for (dex1=0; dex1 < nrow(); dex1++)
    for (dex2=dex1; dex2 < nrow(); dex2++)
	{
	sum = 0.0;
	for (ctr=0; ctr < ncol(); ctr++)
		sum += (*this)[dex1][ctr] * (*this)[dex2][ctr];
	m[dex1][dex2] = sum;
	m[dex2][dex1] = sum;
	}

  return m;
  }

// W, a diagonal matrix, is represented as a vector
Matrix Matrix::find_AWAtranspose(const Vector &w) const
  {
  Matrix m(nrow(), nrow());
  w.assert_dim(ncol());
  int dex1, dex2, ctr;
  double sum;

  for (dex1=0; dex1 < nrow(); dex1++)
    for (dex2=dex1; dex2 < nrow(); dex2++)
	{
	sum = 0.0;
	for (ctr=0; ctr < ncol(); ctr++)
		sum += (*this)[dex1][ctr] * (*this)[dex2][ctr] * w[ctr];
	m[dex1][dex2] = sum;
	m[dex2][dex1] = sum;
	}

  return m;
  }

// W, a diagonal matrix, is represented as a vector
Matrix Matrix::find_AtransposeWA(const Vector &w)  const 
  {
  Matrix m(ncol(), ncol());
  int dex1, dex2, ctr;
  double sum;

  for (dex1=0; dex1 < ncol(); dex1++)
    for (dex2=dex1; dex2 < ncol(); dex2++)
	{
	sum = 0.0;
	for (ctr=0; ctr < nrow(); ctr++)
		sum += (*this)[ctr][dex1] * (*this)[ctr][dex2] * w[ctr];
	m[dex1][dex2] = sum;
	m[dex2][dex1] = sum;
	}

  return m;
  }

Matrix Matrix::find_AtransposeA(void)  const 
  {
  Matrix m(ncol(), ncol());
  int dex1, dex2, ctr;
  double sum;

  for (dex1=0; dex1 < ncol(); dex1++)
    for (dex2=dex1; dex2 < ncol(); dex2++)
	{
	sum = 0.0;
	for (ctr=0; ctr < nrow(); ctr++)
		sum += (*this)[ctr][dex1] * (*this)[ctr][dex2];
	m[dex1][dex2] = sum;
	m[dex2][dex1] = sum;
	}

  return m;
  }

Vector Matrix::find_AtransposeWy(const Vector &w, const Vector &y) const
  {
  assert(nrow() == w.dim()); assert(w.dim() == y.dim());
  Vector v(ncol());

  for (int dex1=0; dex1 < ncol(); dex1++)
    for (int dex2=0; dex2 < nrow(); dex2++)
      v[dex1] += val[dex2][dex1] * w[dex2] * y[dex2];

  return v;
  }

// if failed_p is non-NULL, the SVD failing to converge will return
// with *failed_p=1. (Otherwise, the program terminates)
Matrix Matrix::find_inverse(int *ill_conditioned_p, double *abs_detp, int *failed_p) const 
  {
  int ill_conditioned;
  Matrix m = *this;

  if (failed_p != NULL)
	*failed_p = 0;

  ill_conditioned = m.invert_me(abs_detp, failed_p);

  if (ill_conditioned_p != NULL || (failed_p != NULL && *failed_p))
	*ill_conditioned_p = ill_conditioned;

  return m;
  }

void Matrix::assert_dim(int m_, int n_) const
  {
  if (m_ != m || n_ != n)
	{
	fprintf(stderr, "Matrix::assert_dim failed. Expected %dx%d matrix, "
			  " was %dx%d. \n", m_, n_, m, n);
	abort();
	}
  }

// allocates and returns a 1-indexed version of val. 
// (i.e. returns p, so that p[1..m][1..n] are the elements of the
//  matrix.)
double **Matrix::new_one_indices(void) const
  {
  double *temp;
  double **p;
  int i;

  p = new double*[m+2];
  assert(p != NULL);

  temp = val[0];	// so, temp[0..m*n-1] are all the elements
			// of the matrix
  for (i=1; i <= m; i++)
	p[i] = temp+(n*(i-1)) -1;

  p[0] = (double*)NULL;	// to try to catch some run-off-array bugs
  p[m+1] = (double*)NULL;

  return p;
  }

void Vector::initialize_from_arr(const double *p)
  {
  bcopy((char*)p, (char*)val, sizeof(double)*n);
  }

void Vector::initialize_from_file(FILE *fp)
  {
  double *p;
  int n_ele;

  p = new_doubles_from_file(fp, &n_ele);
  assert(n_ele > 0);
  resize(n_ele);
  initialize_from_arr(p);

  delete[] p;

  return;
  }

void Vector::initialize_from_file(const char *fn)
  {
  FILE *fp;

  fp = safe_fopen(fn, "rt");
  initialize_from_file(fp);
  fclose(fp);

  return;
  }

void Matrix::initialize_from_arr(const double *p)
  {
  bcopy((char*)p, (char*)val[0], sizeof(double)*n*m);

  return;
  }

// assumes that the file contains a rectangular array of numbers,
// and reads that directly into a matrix (resizing the matrix
// if necessary).
void Matrix::initialize_from_file(FILE *fp)
  {
  double *p;
  int new_m, new_n, n_ele;

  new_n = columns_in_file(fp);
  p = new_doubles_from_file(fp, &n_ele);
  new_m = n_ele/new_n;

  if (n_ele % new_n != 0)
	mat_error("initialize_from_file: File bad. (Not rectangular array?)");

  resize(new_m, new_n, 0);
  initialize_from_arr(p);

  delete[] p;

  return;
  }

void Matrix::initialize_from_file(const char *fn)
  {
  FILE *fp;

  fp = safe_fopen(fn, "rt");
  initialize_from_file(fp);
  fclose(fp);

  return;
  }

void Matrix::randomize(double lowLim, double highLim)
  {
  for (int i=0; i < m; i++)
    for (int j=0; j < n; j++)
	val[i][j] = doub_rand(lowLim, highLim);
  }

void Matrix::print(FILE *fp, const char *separator, const char *printStr) const 
  {
  const char DEFAULT_SEPARATOR[] = " ";
  const char DEFAULT_PRINTSTR[] = "%.5f";
  int i,j;

  if (printStr == NULL)
    printStr = DEFAULT_PRINTSTR;
  if (separator == NULL)
    separator = DEFAULT_SEPARATOR;
  if (fp == NULL)
	fp = stdout;

  for (i=0; i < m; i++)
    {
    for (j=0; j < n; j++)
        {
    	fprintf(fp, printStr, val[i][j]);
        fprintf(fp, separator);
        }
    fprintf(fp,"\n");
    }

  return;
  }

// allocates memory for an m by n matrix, (and also sets m and n).
void Matrix::init_structure(int m_, int n_)
  {
  int i; // ,j;
  double *temp;

  if (SHOW_CREATIONS)
	fprintf(stderr,"Matrix created\n");

  m = m_;
  n = n_;

  // allocate matrix
  temp = new double[m*n];
  val = new double*[m+1];
  assert(temp != NULL && val != NULL);
  for (i=0; i < m; i++)
	val[i] = temp+(n*i);
  val[m] = (double*)NULL;	// to try to catch some run-off-array bugs

  return;
  }

Vector Matrix::extract_row(int row) const
  {
  int ctr;
  Vector v(n, 0);

  if (row < 0 || row >= m)
	mat_error("extract_row: Out of bounds.");

  for (ctr=0; ctr < n; ctr++)
	v[ctr] = val[row][ctr];

  return v;
  }

Vector Matrix::extract_col(int col) const
  {
  int ctr;
  Vector v(m, 0);

  if (col < 0 || col >= n)
	mat_error("extract_col: Out of bounds.");

  for (ctr=0; ctr < m; ctr++)
	v[ctr] = val[ctr][col];

  return v;
  }

void Matrix::extract_row(Vector &v, int row) const
  {
  int ctr;
  v.assert_dim(n);

  if (row < 0 || row >= m)
	mat_error("extract_row: Out of bounds.");

  for (ctr=0; ctr < n; ctr++)
	v[ctr] = val[row][ctr];

  return;
  }

void Matrix::extract_col(Vector &v, int col) const
  {
  int ctr;
  v.assert_dim(m);

  if (col < 0 || col >= n)
	mat_error("extract_col: Out of bounds.");

  for (ctr=0; ctr < m; ctr++)
	v[ctr] = val[ctr][col];

  return;
  }

void Matrix::set_row(const double dat[], int row)
  {
  int ctr;

  if (row < 0 || row >= m)
	mat_error("set_row: Out of bounds.");

  for (ctr=0; ctr < n; ctr++)
	val[row][ctr] = dat[ctr];

  return;
  }

void Matrix::set_col(const double dat[], int col)
  {
  int ctr;

  if (col < 0 || col >= n)
	mat_error("set_col: Out of bounds.");

  for (ctr=0; ctr < m; ctr++)
	val[ctr][col] = dat[ctr];

  return;
  }

void Matrix::set_row(const Vector &v, int row)
  {
  if (v.dim() != n)
	mat_error("set_row: vector with bad dimensions.");
  set_row((double*)v, row);
  }

void Matrix::set_col(const Vector &v, int col)
  {
  if (v.dim() != m)
	mat_error("set_col: vector with bad dimensions.");
  set_row((double*)v, col);
  }

Matrix &Matrix::operator=(const Matrix &mat)
  {
  if (m != mat.m || n != mat.n)
	{
	fflush(stdout);
	fprintf(stderr,"Matrix ERROR: (Possible internal error.)\n"
                       "Matrix ERROR: Tried to assign %dx%d to %dx%d\n",
				mat.m, mat.n, m,n);
	abort();
//	exit(-1);
	}

  bcopy((char*)mat.val[0],(char*)val[0], m*n*sizeof(double));	

  return *this;
  }

void Matrix::resize(int m_, int n_, int initialize, double init_val)
  {
  int i,j;

  if (m_ != m || n_ != n)
    {
    // free old
    delete[] (*val);
    delete[] val;

    // initialize new
    init_structure(m_,n_);
    }

  if (initialize)
    for (i=0; i<m; i++)
	for (j=0; j<n; j++)
	    val[i][j] = init_val;

  return;
  }

// copy constructor
Matrix::Matrix(const Matrix &mat)
  {
  init_structure(mat.m, mat.n);

  bcopy((char*)mat.val[0],(char*)val[0], m*n*sizeof(double));	
				// recall that all elements are
				// stored in a big contiguous linear array
  }

Matrix::Matrix(int m_, int n_, int initialize, double init_val)
  {
  int i,j;

  if (!initialize && init_val != 0.0)
	internal_warn("Matrix::Matrix(...) : told not to initialize, but "
			"initial value was specified.");

  init_structure(m_,n_);

  if (initialize)
    for (i=0; i<m; i++)
	for (j=0; j<n; j++)
	    val[i][j] = init_val;

  return;
  }


/*
void Matrix::loadFromFile(const char *fn, int resizeOkay)
  {
  FILE *fp = safe_fopen(fn, "rt");
  int temp = columns_in_file(fp, 0);
  if (temp != 2)
    {
    fprintf(stderr, "ERROR: Matrix::loadFromFile: %s\n", fn);
    error("Dimensions not on first line of file.");
    }
  m = safe_readInt(fp);
  n = safe_readInt(fp);

double safe_readDouble(FILE *fp);
  
  dont_warn_if_not_at_start=0);
*/

Matrix::Matrix(const char *fn)
  {
  init_structure(1,1);	// this will be undone in loadFromFile(), but whatever....
  initialize_from_file(fn);
  return;
  }

Matrix::~Matrix()
  {
  if (SHOW_CREATIONS)
	fprintf(stderr,"Matrix destroyed\n");

  delete[] (*val);
  delete[] val;
  }

Matrix operator*(const Matrix &mat1, const Matrix &mat2)
  {
  int i, j, k;
  double temp;

  if (mat1.ncol()!= mat2.nrow())
	mat_error("(Possible internal error.) Tried to multiply "
		"two matrixes with dimensions that don't match up.");

  Matrix mat(mat1.nrow(), mat2.ncol(), 0);

  for (i=0; i < mat.nrow(); i++)
    for (j=0; j < mat.ncol(); j++)
	{
	temp=0.0;
	for (k=0; k < mat1.ncol(); k++)
		temp += mat1[i][k] * mat2[k][j];
	mat[i][j] = temp;
	}

  return mat;
  }

Matrix operator*(double d, const Matrix &mat)
  {
  int i, j;

  Matrix m(mat.nrow(), mat.ncol(), 0);

  for (i=0; i < mat.nrow(); i++)
    for (j=0; j < mat.ncol(); j++)
	m[i][j] = d*mat[i][j];

  return m;
  }

Matrix operator*(const Matrix &mat, double d)
  {
  return d*mat;
  }

Matrix operator-(const Matrix &mat)
  {
  int i,j;

  Matrix m(mat.nrow(), mat.ncol(), 0);

  for (i=0; i < mat.nrow(); i++)
    for (j=0; j < mat.ncol(); j++)
	m[i][j] = -mat[i][j];

  return m;
  }

Matrix operator-(const Matrix &mat1, const Matrix &mat2)
  {
  int i, j;

  if (mat1.nrow() != mat2.nrow()
      || mat1.ncol() != mat2.ncol())
	mat_error("(Possible internal error.) Tried to subtract "
		"two matrices, and dimensions don't match up.");

  Matrix m(mat1.nrow(), mat1.ncol(), 0);

  for (i=0; i < mat1.nrow(); i++)
    for (j=0; j < mat1.ncol(); j++)
	m[i][j] = mat1[i][j] - mat2[i][j];

  return m;
  }

Matrix operator+(const Matrix &mat1, const Matrix &mat2)
  {
  int i, j;

  if (mat1.nrow() != mat2.nrow()
      || mat1.ncol() != mat2.ncol())
	mat_error("(Possible internal error.) Tried to subtract "
		"two matrices, and dimensions don't match up.");

  Matrix m(mat1.nrow(), mat1.ncol(), 0);

  for (i=0; i < mat1.nrow(); i++)
    for (j=0; j < mat1.ncol(); j++)
	m[i][j] = mat1[i][j] + mat2[i][j];

  return m;
  }

Matrix &Matrix::operator*=(double d)
  {
  int i, j;

  for (i=0; i < m; i++)
    for (j=0; j < n; j++)
	val[i][j] *= d;

  return *this;
  }

Matrix &Matrix::operator/=(double d)
  {
  int i, j;

  for (i=0; i < m; i++)
    for (j=0; j < n; j++)
	val[i][j] /= d;

  return *this;
  }

Matrix operator/(const Matrix &mat, double d)
  {
  return mat * (1.0/d);
  }

Matrix &Matrix::operator-=(const Matrix &mat)
  {
  int i, j;

  if (nrow() != mat.nrow() || ncol() != mat.ncol())
	mat_error("(Possible internal error.) Tried to subtract "
		"two matrices, and dimensions don't match up.");

  for (i=0; i < m; i++)
    for (j=0; j < n; j++)
	val[i][j] -= mat[i][j];

  return *this;
  }

Matrix &Matrix::operator+=(const Matrix &mat)
  {
  int i, j;

  if (nrow() != mat.nrow() || ncol() != mat.ncol())
	mat_error("(Possible internal error.) Tried to subtract "
		"two matrices, and dimensions don't match up.");

  for (i=0; i < m; i++)
    for (j=0; j < n; j++)
	val[i][j] += mat[i][j];

  return *this;
  }

void Vector::print(FILE *fp, const char *separator, int print_cr, const char *printStr) const 
  {
  const char DEFAULT_SEPARATOR[] = " ";
  const char DEFAULT_PRINTSTR[] = "%.5f";
  int i;

  if (printStr == NULL)
    printStr = DEFAULT_PRINTSTR;
  if (separator == NULL)
    separator = DEFAULT_SEPARATOR;
  if (fp == NULL)
	fp = stdout;

  for (i=0; i < n; i++)
	{
	fprintf(fp, printStr, val[i]);
	if (i != n-1)
		fprintf(fp, separator);
	}

  if (print_cr)
	fprintf(fp, "\n");

  return;
  }

void Vector::assert_dim(int n_) const
  {
  if (n_ != n)
	{
	fprintf(stderr, "Vector::assert_dim failed. Expected %d dimensional, "
			  " was %d dimensional. \n", n_, n);
	abort();
	}
  }

void Vector::enlarge(int new_size, int init_new_ele, double init_val)
  {
  assert(new_size >= n);

  if (new_size == n)
	return;

  int ctr;
  double *new_val = new double[new_size];
  assert(new_val != NULL);

  for (ctr=0; ctr < n; ctr++)
	new_val[ctr] = val[ctr];
  if (init_new_ele)
	for (; ctr < new_size; ctr++)
		new_val[ctr] = init_val;

  delete[] val;
  val = new_val;
  n = new_size;

  return; 
  }

// copy constructor
Vector &Vector::operator=(const Vector &vec)
  {
  if (n != vec.n)
	{
	fflush(stdout);
	fprintf(stderr,"Vector ERROR: (Probable internal error.)\n"
                   "Vector ERROR: Tried to assign %dd-vector to %dd-vector\n",
				vec.n, n);
	abort();
//	exit(-1);
	}

  bcopy((char*)vec.val, (char*)val, n*sizeof(double));	

  return *this;
  }

Vector::Vector(const Vector &vec)
  {
  n = vec.n;
  val = new double[n];
  assert(val != NULL);

  bcopy((char*)vec.val, (char*)val, n*sizeof(double));

  if (SHOW_CREATIONS)
	fprintf(stderr,"Vector Created\n");

  return;
  }

Vector::Vector(int n_, int initialize, double init_val)
  {
  int i;

  if (SHOW_CREATIONS)
	fprintf(stderr,"Vector Created\n");

  if (!initialize && init_val != 0.0)
	internal_warn("Vector::Vector(...) : told not to initialize, but "
			"initial value was specified.");

  n = n_;

  val = new double[n];
  assert(val != NULL);
  if (initialize)
    for (i=0; i < n; i++)
	val[i] = init_val;

  return;
  }

Vector::Vector(const char *fn)
  {
  FILE *fp = safe_fopen(fn, "rt");
  int nCols = columns_in_file(fp);
  if (nCols <= 0)
    mat_error("Vector::Vector: File bad. (Not vector?)");
  int n_ele; 

  val = new_doubles_from_file(fp, &n_ele);
  assert(val != NULL);
  if (n_ele != nCols && nCols != 1)	// check if looks like either 
    					// row or column vector. 
    mat_error("Vector::Vector: File bad.  (Not vector?)");

  n = n_ele;

  fclose(fp);

  return;
  }

Vector::~Vector()
  {
  if (SHOW_CREATIONS)
	fprintf(stderr, "Vector destroyed\n");

  delete[] val;
  }

Vector operator*(const Matrix &mat, const Vector &vec)
  {
  int i, j;
  double temp;

  if (mat.ncol() != vec.n_ele())
	mat_error("(Possible internal error.) Tried to multiply "
		"matrix with vector, and dimensions don't match up.");

  Vector v(mat.nrow(), 0);

  for (i=0; i < mat.nrow(); i++)
    {
    temp=0.0;
    for (j=0; j < mat.ncol(); j++)
	temp += mat[i][j] * vec[j];
		
    v[i] = temp;
    }

  return v;
  }

//AAA
Matrix operator*(const Vector &vec, const Matrix &mat)
  {
  int i, j;

  if (mat.nrow() != 1)
	mat_error("(Possible internal error.) Tried to multiply "
		"vector with matrix, and dimensions don't match up.");

  Matrix m(vec.dim(), mat.ncol());

  for (i=0; i < vec.dim(); i++)
    for (j=0; j < mat.ncol(); j++)
      m[i][j] = vec[i] * mat[0][j];

  return m;
  }

void Vector::normalize_sum_to(double desiredSum)
  {
  double sum=0.0;
  for (int i=0; i < n; i++)
    {
    if (val[i] < 0)
    	mat_error("normalize_sum_to: called on a non-positive vector.");
    sum += val[i];
    }

  if (sum == 0.0 && desiredSum != 0.0)
    	mat_error("normalize_sum_to: got zero vector.");

  double mult = desiredSum / sum;
  for (int i=0; i < n; i++)
	val[i] *= mult;

  // HERE/debugging only
  // assert(double_eq(sum_arr(val,n), desiredSum));

  return;
  }

void Vector::resize(int n_, int initialize, double init_val)
  {
  if (!initialize && init_val != 0.0)
	internal_warn("Vector::Vector(...) : told not to initialize, but "
			"initial value was specified.");

  if (n != n_)
    {
    delete[] val;
    n = n_;
    val = new double[n];
    assert(val != NULL);
    }

  if (initialize)
    for (int i=0; i < n; i++)
	val[i] = init_val;

  return;
  }

Vector operator*(double d, const Vector &vec)
  {
  int ctr;

  Vector v(vec.n_ele(), 0);

  for (ctr=0; ctr < vec.n_ele(); ctr++)
	v[ctr] = d*vec[ctr];

  return v;
  }

Vector operator*(const Vector &vec, double d)
  {
  return d*vec;
  }

Vector operator/(const Vector &vec, double d)
  {
  return vec*(1.0/d);
  }

Vector operator-(const Vector &vec)
  {
  int ctr;

  Vector v(vec.n_ele(), 0);

  for (ctr=0; ctr < vec.n_ele(); ctr++)
	v[ctr] = -vec[ctr];

  return v;
  }

Vector operator-(const Vector &vec1, const Vector &vec2)
  {
  int ctr;

  if (vec1.n_ele() != vec2.n_ele())
	mat_error("(Possible internal error.) Tried to subtract "
		"two vectors, and dimensions don't match up.");

  Vector v(vec1.n_ele(), 0);

  for (ctr=0; ctr < vec1.n_ele(); ctr++)
	v[ctr] = vec1[ctr] - vec2[ctr];

  return v;
  }

Vector operator+(const Vector &vec1, const Vector &vec2)
  {
  int ctr;

  if (vec1.n_ele() != vec2.n_ele())
	{
	printf("Matrix/Vector ERROR: Dimensions %d, %d\n", vec1.n_ele(), vec2.n_ele());
	mat_error("(Possible internal error.) Tried to add "
		"two vectors, and dimensions don't match up.");
	}

  Vector v(vec1.n_ele(), 0);

  for (ctr=0; ctr < vec1.n_ele(); ctr++)
	v[ctr] = vec1[ctr] + vec2[ctr];

  return v;
  }

Vector &Vector::operator*=(double d)
  {
  int ctr;

  for (ctr=0; ctr < n; ctr++)
	val[ctr] *= d;

  return *this;
  }

Vector &Vector::operator/=(double d)
  {
  int ctr;

  for (ctr=0; ctr < n; ctr++)
	val[ctr] /= d;

  return *this;
  }

Vector &Vector::operator+=(const Vector &vec)
  {
  int ctr;

  if (n != vec.n_ele())
	{
	printf("Matrix/Vector ERROR: (+=) Dimensions %d, %d\n", n, vec.n_ele());
	mat_error("(Possible internal error.) Tried to add "
		"two vectors, and dimensions don't match up.");
	}

  for (ctr=0; ctr < n; ctr++)
	val[ctr] += vec[ctr];

  return *this;
  }

Vector &Vector::operator-=(const Vector &vec)
  {
  int ctr;

  if (n != vec.n_ele())
	mat_error("(Possible internal error.) Tried to subtract "
		"two vectors, and dimensions don't match up.");

  for (ctr=0; ctr < n; ctr++)
	val[ctr] -= vec[ctr];

  return *this;
  }

Vector::operator Matrix() const
  {
  Matrix m(n,1,0); 
  m.initialize_from_arr(val);
  return m;
  }

Matrix Vector::find_transpose(void) const
  { 
  Matrix m(1,n); 
  m.initialize_from_arr(val);
  return m; 
  }

Matrix Vector::T(void) const
  { 
  Matrix m(1,n); 
  m.initialize_from_arr(val);
  return m; 
  }

Vector asVector(const Matrix &m)
  {
  assert(m.ncols()==1);
  Vector v(m.nrows());
  for (int i=0; i < m.nrows(); i++)
	v[i] = m[0][i];
  return v;
  }

double dot_prod(const Vector &vec1, const Vector &vec2)
  {
  double temp;
  int ctr;

  if (vec1.length() != vec2.length())
	{
	fprintf(stderr,"\nMatrix/Vector ERROR: lens: %d, %d\n", 
			vec1.length(), vec2.length());
	mat_error("(Possible internal error.) dot_prod(): Tried to "
		"take dot product of vectors of unequal length");
	}

  temp=0;
  for (ctr=0; ctr < vec1.length(); ctr++)
	temp += vec1[ctr] * vec2[ctr];

  return temp;
  }

// tests if 2 vectors are "approximately equal". See comments
// for double_eq(...) in misc.C for details.
int vect_approx_equal(const Vector &v1, const Vector &v2, 
			double epsilon, double delta)
  {
  int ctr, len;

  if (v1.length() != v2.length())
	mat_error("vect_approx_equal(): Given vectors "
			"with different lengths.");

  len = v1.length();
  for (ctr=0; ctr < len; ctr++)
    if (!double_eq(v1[ctr], v2[ctr], epsilon, delta))
	return 0;

  return 1;
  }

// tests if 2 vectors are "approximately equal". See comments
// for double_eq(...) in misc.C for details.
int mat_approx_equal(const Matrix &m1, const Matrix &m2, 
			double epsilon, double delta)
  {
  int i, j, ncol, nrow;

  if (m1.nrow() != m2.nrow() || m1.ncol() != m2.ncol())
	mat_error("mat_approx_equal(): Given matrices "
			"with different dimensions.");

  nrow = m1.nrow();
  ncol = m1.ncol();

  for (i=0; i < nrow; i++)
    for (j=0; j < ncol; j++)
	if (!double_eq(m1[i][j], m2[i][j], epsilon, delta))
		return 0;

  return 1;
  }

double euclideanDistance(const Vector &v1, const Vector &v2)
  {
  int ctr;
  double sum=0;

  assert(v1.length() == v2.length());
  for (ctr=0; ctr < v1.length(); ctr++)
	sum += sqr(v1[ctr]-v2[ctr]);

  return sqrt(sum);
  }

double squaredEuclideanDistance(const Vector &v1, const Vector &v2)
  {
  int ctr;
  double sum=0;

  assert(v1.length() == v2.length());
  for (ctr=0; ctr < v1.length(); ctr++)
	sum += sqr(v1[ctr]-v2[ctr]);

  return sum;
  }

// sub-vector including val[first]..val[lastp1-1]
Vector Vector::subVector(int first, int lastp1) const 
  {
  assert(first>=0 && lastp1 > first && lastp1 <= n);

  Vector v(lastp1-first,0);

  int i, dex;
  for (i=first, dex=0; i < lastp1; i++, dex++)
	v[dex] = val[i];

  assert(dex==v.length());

  return v;
  }

// sets elements first..lastp1-1 with data from v
void Vector::setSubVector(const Vector &v, int first, int lastp1)
  {
  v.assert_dim(lastp1-first);

  bcopy(v.val, val+first, sizeof(double)*(lastp1-first));

  // debugging. (Check I did the bcopy right....)
  // int i, dex;
  // for (i=first, dex=0; i < lastp1; i++, dex++)
  //	assert(val[i] == v[dex]);
  // assert(dex==v.length());

  return;
  }

// sum's sub-vector, including val[first]..val[lastp1-1]
double Vector::sumSubVector(int first, int lastp1) const 
  {
  assert(first>=0 && lastp1 >= first && lastp1 <= n);

  double sum =0.0;
  for (int i=first; i < lastp1; i++)
	sum += val[i];

  return sum;
  }

void Vector::randomize(double lowLim, double highLim)
  {
  for (int i=0; i < n; i++)
	val[i] = doub_rand(lowLim, highLim);

  return;
  }


/*
// Test code 
int main(void)
  {
  FILE *fp;
  Matrix m(3,3);
  Vector y(3);
  double d[] = {1, 1, 1,		// 3, 4, 5
                1, 1, -1,
		0, -1, 1 };
  double foo;

  y[0] = 12; 
  y[1] = 2;
  y[2] = 1;

  m.initialize_from_arr(d);
  m.print();
  y.print();

  printf("---\n");
  (m.linSolve(y, 1)).print();

  //m.find_inverse(NULL, &foo);
  //printf("abs(det) = %f\n", foo);

  return 0;
  }
*/


/**** To solve a linear of system of equations:

say you have training samples 1,2,3,4->10    5,6,7,8->13,    9,10,11,12->11.
(i.e. V(s) - gamma*E[V(s`)] = [1,2,3,4], when you recieved reinforcement 
  10 going from s to s'.)

Let x = [  1 2  3  4   1]	// note 1 added for constant intercept term.
        [  5 6  7  8   1]
        [ 9 10 11 12   1]

// One thing: the "terminal" state should be represented as [0 0 0 0  0] ALWAYS. (note no intercept.
//	this way, V(terminal) = beta * 0 = 0 always).

y = [10] 
    [13]
    [11]

// to find beta to minize (X*beta-y)^2, do:

Vector beta = x.linSolve(y,1);

****/


