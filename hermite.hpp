#include <stdio.h>
#include "math.h"

//void fillpoints(double *x, double *y, double *dy, int n, double a, double b, double (*f)(double), double (*df)(double));
void hermitefunc(int n, double *x, double *y, double *dy, double *h);
void hermitefunc(int n, double *y, double *dy, double *h, double a, double b);
double hermitevaluen(double p, double *x, double *h, int i);
double hermitevaluen(double p, double *h, int i);
