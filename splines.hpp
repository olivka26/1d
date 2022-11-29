#include <stdio.h>
#include "math.h"

void splinefunc(int n, double *x, double *y, double *s, double *dy);
void splinefunc(int n, double *y, double *s, double *dy, double a, double b);
double splinevaluen(double p, double *s, double * x, int i);
double splinevaluen(double p, double *s, int i);
