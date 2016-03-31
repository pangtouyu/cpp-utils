#ifndef _INTEGRATE_H
#define _INTEGRATE_H

typedef double (*doub_func_p)(double, void *p);

double simpson(doub_func_p f, void *parms, double a, double b, double epsilon);
double romberg(doub_func_p f, void *parms, double a, double b, double epsilon);

#endif 		// _INTEGRATE_H
