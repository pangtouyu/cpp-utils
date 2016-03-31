#ifndef _INTEGRATE_H
#define _INTEGRATE_H

#include "matrix.h"

typedef Vector (*vec_func_p)(double, void *p);

Vector vecRomberg(vec_func_p f, void *parms, double a, double b, Vector &epsilon);

#endif 		// _INTEGRATE_H
