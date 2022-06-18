#include "splines.hpp"
              
void splinefunc(int n, double *x, double *y, double *s, double *dy, double a, double b){
    s[0]=y[0];
    s[1]=dy[0];
    s[3]=4.0*(y[1]+y[0])+(dy[0]*(b-a)/(n-1));
    s[4]=7.0;
    s[5]=1.0;
    for(int i=2;i<n-1;++i){
        s[3*i+1]=6.0-(1.0/s[3*i-2]);
        //s[3*i+2]=1.0;
        s[3*i]=4.0*(y[i]+y[i-1])-(s[3*i-3]/s[3*i-2]);
    }
    s[3*n-2]=7.0-(1.0/s[3*n-5]);
    s[3*n-3]=4.0*(y[n-1]+y[n-2])-(dy[n-1]*(b-a)/(n-1))-(s[3*n-6]/s[3*n-5]);
    s[3*n-3]/=s[3*n-2];
    s[3*n-1]=(s[3*n-3]-y[n-1]);
   // s[3*n-3]/=s[3*n-2];
    s[3*n-1]*=2.0;
    s[3*n-1]*=(n-1);
    s[3*n-1]/=(b-a);
    s[3*n-1]+=dy[n-1];
    s[3*n-1]*=2.0;
    s[3*n-1]*(n-1);
    s[3*n-1]/=(b-a);
    s[3*n-2]=(n-1)*4.0*(y[n-1]-s[3*n-3]);
    s[3*n-2]/=(b-a);
    s[3*n-2]-=dy[n-1];
    for(int i=n-2;i>0;--i){
        s[3*i]-=s[3*(i+1)]; //muly
        s[3*i]/=s[3*i+1];
    }
    s[2]=s[3]-y[0];
    s[2]*=2.0;
    s[2]*=(n-1);
    s[2]/=(b-a);
    s[2]-=dy[0];
    s[2]*=2.0;
    s[2]*=(n-1);
    s[2]/=(b-a);
    for(int i=n-2; i>0; --i){
        s[3*i+1]=4.0*y[i]-3.0*s[3*i]-s[3*(i-1)];
        s[3*i+1]*=(n-1);
        s[3*i+1]/=(b-a);
        s[3*i+2]=s[3*(i+1)]+s[3*i]-2.0*y[i];
        s[3*i+2]*=2.0;
        s[3*i+2]*=(n-1);
        s[3*i+2]*=(n-1);
        s[3*i+2]/=(b-a);
        s[3*i+2]/=(b-a);
    }
}

double splinevalue(double p, double a, double b, int n, double *s, double *x){
    int i;
    double left=a, step=(b-a)/(n-1), right=a+step/2;
    for(i=0; i<n; ++i){
        if(p>=left && p<=right){
            break;
        }
        left=right;
        right+=step;
    }
    //printf("splines: %d %lf %lf\n",i, left, right);
    //printf("ceof: %lf %lf %lf\n",s[3*i], s[3*i+1], s[3*i+2]);
    if(i==n)
        --i;
    return s[3*i]+s[3*i+1]*(p-left)+s[3*i+2]*(p-left)*(p-left);
}
