#include <stdio.h>
#include "math.h"

void splinefunc(int n, double *x, double *y, double *s, double *dy, double a, double b);
double splinevalue(double p, double a, double b, int n, double *s, double *x);
