#include <math.h>

#include "dgolden.h"

double dneg_func(double x, void *parms_with_func)
  {
  dparms_w_func_t *p;

  p = (dparms_w_func_t*)parms_with_func;

  return -(*(p->func))(x, p->parms);
  }

// utils from nrutil.h
static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
		(maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define GOLD 1.618034
#define GLIMIT 10.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

// from Numerical Recipes 2nd Ed., by W. Press
// pp. 400-402
// modified to pass a (void*)parms to the function being
// optimized as well.
void dmnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
			double *fc, double (*func)(double,void*), void *parms)
  {
  double ulim, u, r, q, fu, dum;

  *fa = (*func)(*ax, parms);
  *fb = (*func)(*bx, parms);
  if (*fb > *fa)
		{
		SHFT(dum, *ax, *bx, dum)
		SHFT(dum, *fb, *fa, dum)
		}

  *cx=(*bx) + GOLD*(*bx-*ax);
  *fc=(*func)(*cx, parms);

  while (*fb > *fc)
		{
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q - (*bx-*ax)*r)/
				(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
		ulim=(*bx) + GLIMIT*(*cx-*bx);

		if ((*bx-u)*(u-*cx) > 0.0)
			{
			fu=(*func)(u,parms);
			if (fu < *fc)
				{
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;

				return;
				}
			else if (fu > *fb)
				{
				*cx = u;
				*fc = fu;

				return;
				}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u, parms);
			}
		else if ((*cx-u)*(u-ulim) > 0.0)
			{
			fu = (*func)(u,parms);
			if (fu < *fc)
				{
				SHFT(*bx, *cx, u, *cx+GOLD*(*cx-*bx))
				SHFT(*fb, *fc, fu, (*func)(u,parms))
				}
			}
		else if ((u-ulim)*(ulim-*cx) > 0.0)
			{
			u=ulim;
			fu=(*func)(u,parms);
			}
		else
			{
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u,parms);
			}

		SHFT(*ax, *bx, *cx, u)
		SHFT(*fa, *fb, *fc, fu)
		}
  }

void dmnbrak_max(double *ax, double *bx, double *cx, double *fa, double *fb,
			double *fc, double (*func)(double,void*), void *parms)
  {
  dparms_w_func_t p;

  p.func = func;
  p.parms = parms;

  dmnbrak(ax,bx,cx,fa,fb,fc,dneg_func,&p);
  }

#define R 0.61803399
#define C (1.0-R)
#define SHFT2(a,b,c) (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

// assumes that ax, bx and cx bracket a minimum
double dgolden(double ax, double bx, double cx, double (*f)(double, void*),
			 double tol, double *xmin, void *parms)
  {
  double f1, f2, x0, x1, x2, x3;

  x0=ax;
  x3=cx;

  if (fabs(cx-bx) > fabs(bx-ax))
		{
		x1=bx;
		x2=bx+C*(cx-bx);
		}
  else
		{
		x2=bx;
		x1=bx-C*(bx-ax);
		}

  f1 = (*f)(x1, parms);
  f2 = (*f)(x2, parms);

  while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2)) && (fabs(x3-x0) > ABS_TOLERANCE))
		{
		if (f2 < f1)
			{
			SHFT3(x0,x1,x2,R*x1+C*x3)
			SHFT2(f1,f2,(*f)(x2,parms))
			}
		else
			{
			SHFT3(x3,x2,x1,R*x2+C*x0)
			SHFT2(f2,f1,(*f)(x1,parms))
			}
		}

  if (f1 < f2)
	{
	*xmin=x1;
	return f1;
	}
  else
	{
	*xmin=x2;
	return f2;
	}
  }

// other stuff --------------------------------------

// just like golden(), but maximizes instead of minimizes
double dgolden_max(double ax, double bx, double cx, double (*f)(double, void*),
			 double tol, double *xmin, void *parms)
  {
  dparms_w_func_t p;

  p.func = f;
  p.parms = parms;

  return dgolden(ax, bx, cx, dneg_func, tol, xmin, &p);
  }

double doptimize(double *ax, double *bx, double *cx, double *fa, double *fb,
			double *fc, double (*func)(double,void*), void *jnq,
			double accuracy)
  {
  double xmin;

  dmnbrak(ax, bx, cx, fa, fb, fc, func, jnq);
  dgolden(*ax, *bx, *cx, func, accuracy, &xmin, jnq);

  return xmin;
  }

double doptimize_max(double *ax, double *bx, double *cx, double *fa, double *fb,
			double *fc, double (*func)(double,void*), void *jnq,
			double accuracy)
  {
  double xmax;

  dmnbrak_max(ax, bx, cx, fa, fb, fc, func, jnq);
  dgolden_max(*ax, *bx, *cx, func, accuracy, &xmax, jnq);

  return xmax;
  }

double doptimize(double ax, double bx, double (*func)(double,void*), 
		void *parms, double accuracy)
  {
  double ax2, bx2, cx;
  double fa, fb, fc; 

  ax2=ax;
  bx2=bx;
  return doptimize(&ax2, &bx2, &cx, &fa, &fb, &fc, func, parms, accuracy);
  }

double doptimize_max(double ax, double bx, double (*func)(double,void*), 
		void *parms, double accuracy)
  {
  double ax2, bx2, cx;
  double fa, fb, fc; 

  ax2=ax;
  bx2=bx;
  return doptimize_max(&ax2, &bx2, &cx, &fa, &fb, &fc, func, parms, accuracy);
  }

