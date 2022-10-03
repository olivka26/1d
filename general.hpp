#include <stdio.h>
#include "math.h"
#include "chebyshev.hpp"
#include "hermite.hpp"
#include "splines.hpp"

double f0(double x);
double f1(double x);
double f2(double x);
double f3(double x);
double f4(double x);
double f5(double x);
double f6(double x);
double df0(double x);
double df1(double x);
double df2(double x);
double df3(double x);
double df4(double x);
double df5(double x);
double df6 (double x);
double d2f0(double x);
double d2f1(double x);
double d2f2(double x);
double d2f3(double x);
double d2f4(double x);
double d2f5(double x);
double d2f6 (double x);
double min(double a, double b);
double max(double a, double b);
void fillpoints(double *x, int n, double a, double b);
void fillvalues(double *x, double *y, double *dy, int n, double (*f)(double), double (*df)(double));
/*void chebyshevextrema(double *c, double *cx, int n, double a, double b, double *extr);
void hermiteextrema(double *h, double *x, int n, double *extr);
void splineextrema(double *s, double *x, int n, double *extr);
void chebysheverrorextrema(double *c, double *cx, int n, double a, double b, double *extr, double (*f)(double));
void hermiteerrorextrema(double *h, double *x, double *y, int n, double *extr, double (*f)(double));
void splineerrorextrema(double *sp, double *x, int n, double *extr, double (*f)(double));*/
