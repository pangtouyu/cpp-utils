#ifndef _DGOLDEN_H
#define _DGOLDEN_H

const double ABS_TOLERANCE = 1.0e-20;

struct dparms_w_func_t
  {
  void *parms;
  double (*func)(double,void*);
  };

double dneg_func(double x, void *dparms_with_func);

void dmnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
			double *fc, double (*func)(double,void*), void *parms);

double dgolden(double ax, double bx, double cx, double (*f)(double, void*),
			 double tol, double *xmin, void *parms);

double doptimize(double *ax, double *bx, double *cx, double *fa, double *fb,
			double *fc, double (*func)(double,void*), void *jnq,
			double accuracy);

// returns the value of x that minimizes func
double doptimize(double ax, double bx, double (*func)(double,void*), 
		void *parms, double accuracy);
double doptimize_max(double ax, double bx, double (*func)(double,void*), 
		void *parms, double accuracy);

#endif		// _DGOLDEN_H
