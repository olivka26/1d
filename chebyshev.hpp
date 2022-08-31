#include <stdio.h>
#include "math.h"

void chebyshevpoints(double *cx, double *cy, int n, double a, double b, double (*f)(double));
double chebyshevtrans(double a, double b, double x);
void chebyshevfunc(double *c, double (*f)(double), int n, double a, double b);
void chebyshevfunc(double *c, double *cy, int n);
double chebyshevvalue(double x, double a, double b, double *c, int n);
