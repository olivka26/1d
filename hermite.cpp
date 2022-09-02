#include "hermite.hpp"

void hermitefunc(int n, double *x, double *y, double *dy, double *h){
    for(int i=0; i<n-1; ++i){
        double dify=y[i+1]-y[i], difx=x[i+1]-x[i], divdif=dify/difx;
        h[4*i]=y[i];
        h[4*i+1]=dy[i];
        h[4*i+2]=(3*divdif-2*dy[i]-dy[i+1])/difx;
        h[4*i+3]=(dy[i]+dy[i+1]-2*divdif)/(difx*difx);
        h[4*i]-=h[4*i+1]*x[i];
        h[4*i]+=h[4*i+2]*x[i]*x[i];
        h[4*i]-=h[4*i+3]*x[i]*x[i]*x[i];
        h[4*i+1]-=2.0*h[4*i+2]*x[i];
        h[4*i+1]+=3.0*h[4*i+3]*x[i]*x[i];
        h[4*i+2]-=3.0*h[4*i+3]*x[i];
    }
}

void hermitefunc(int n, double *y, double *dy, double *h, double a, double b){
    double difx=(b-a)/(n-1);
    for(int i=0; i<n-1; ++i){
        double dify=y[i+1]-y[i], divdif=dify/difx;
        h[4*i]=y[i];
	h[4*i+1]=dy[i];
	h[4*i+2]=(3*divdif-2*dy[i]-dy[i+1])/difx;
	h[4*i+3]=(dy[i]+dy[i+1]-2*divdif)/(difx*difx);
    }
}

double hermitevaluen(double p, double *x, double *h, int i){
    return h[4*i]+h[4*i+1]*(p-x[i])+h[4*i+2]*(p-x[i])*(p-x[i])+h[4*i+3]*(p-x[i])*(p-x[i])*(p-x[i]);
}

double hermitevaluen(double p, double *h, int i){
    return h[4*i]+h[4*i+1]*p+h[4*i+2]*p*p+h[4*i+3]*p*p*p;
}
