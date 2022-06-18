#include "hermite.hpp"

void hermitefunc(int n, double *x, double *y, double *dy, double *h, double a, double b){
    for(int i=0; i<n-1; ++i){
        double dify=y[i+1]-y[i], difx=(b-a)/(n-1), divdif=dify/difx;
        h[4*i]=y[i];
	h[4*i+1]=dy[i];
	h[4*i+2]=(3*divdif-2*dy[i]-dy[i+1])/difx;
	h[4*i+3]=(dy[i]+dy[i+1]-2*divdif)/(difx*difx);
    }
}

double hermitevalue(double p, double a, double b, int n, double *x, double *h){
    int i;
    for(i=0; i<n-1; ++i){
        if(p>=x[i] && p<=x[i+1]){
            break;
        }
    }
    if(i==n-1)
        --i;
    //printf("%d %lf %lf\n", i, x[i], x[i+1]);
    return h[4*i]+h[4*i+1]*(p-x[i])+h[4*i+2]*(p-x[i])*(p-x[i])+h[4*i+3]*(p-x[i])*(p-x[i])*(p-x[i]);
}
