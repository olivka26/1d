#include <stdio.h>
#include "math.h"
#include "chebyshev.hpp"
#include "hermite.hpp"
#include "splines.hpp"

double min(double a, double b);
double max(double a, double b);
void fillpoints(double *x, double *y, double *dy, int n, double a, double b, double (*f)(double), double (*df)(double));
void chebyshevextrema(double *c, double *cx, int n, double a, double b, double *extr);
void hermiteextrema(double *h, double *x, int n, double *extr);
void splineextrema(double *s, double *x, int n, double *extr);
void chebysheverrorextrema(double *c, double *cx, int n, double a, double b, double *extr, double (*f)(double));
void hermiteerrorextrema(double *h, double *x, double *y, int n, double *extr, double (*f)(double));
void splineerrorextrema(double *sp, double *x, int n, double *extr, double (*f)(double));
