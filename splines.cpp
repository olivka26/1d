#include "splines.hpp"

void splinefunc(int n, double *x, double *y, double *s, double *dy){
    double dif0=x[1]-x[0], dif1=x[2]-x[1], dif2=x[3]-x[2], coef;
    s[0]=y[0];
    s[1]=dy[0];
    s[3]=y[1]/dif0+y[1]/dif1+2.0*y[0]/dif0+dy[0]/2.0;
    s[4]=3.0/dif0+1.0/(dif0+dif1);
    s[5]=1.0/dif1-1.0/(dif0+dif1);
    for(int i=2;i<n-1;++i){
        dif0=x[i-1]-x[i-2];
        dif1=x[i]-x[i-1];
        dif2=x[i+1]-x[i];
        coef=1.0/dif0-1.0/(dif0+dif1);
        coef/=s[3*i-2];
        s[3*i+1]=1.0/(dif0+dif1)+2.0/dif1+1.0/(dif2+dif1);
        s[3*i+1]-=coef*s[3*i-1];
        s[3*i+2]=1.0/dif2-1.0/(dif1+dif2);
        s[3*i]=y[i]/dif2+y[i]/dif1+y[i-1]/dif1+y[i-1]/dif0;
        s[3*i]-=coef*s[3*(i-1)];
    }
    coef=1.0/dif1-1.0/(dif1+dif2);
    coef/=s[3*n-5];
    s[3*n-3]=y[n-2]/dif2+y[n-2]/dif1+2.0*y[n-1]/dif2-dy[n-1]/2.0;
    s[3*n-3]-=coef*s[3*n-6];
    s[3*n-2]=1.0/(dif2+dif1)+3.0/dif2;
    s[3*n-2]-=coef*s[3*n-4];
    s[3*n-3]/=s[3*n-2];
    s[3*n-1]=(s[3*n-3]-y[n-1]);
    s[3*n-1]*=2.0;
    s[3*n-1]/=dif2;
    s[3*n-1]+=dy[n-1];
    s[3*n-1]*=2.0;
    s[3*n-1]/=dif2;
    s[3*n-2]=y[n-1]-s[3*n-3];
    s[3*n-2]*=4.0;
    s[3*n-2]/=dif2;
    s[3*n-2]-=dy[n-1];
    for(int i=n-2;i>0;--i){
        s[3*i]-=s[3*(i+1)]*s[3*i+2];
        s[3*i]/=s[3*i+1];
    }
    s[2]=s[3]-y[0];
    s[2]*=2.0;
    s[2]/=(x[1]-x[0]);
    s[2]-=dy[0];
    s[2]*=2.0;
    s[2]/=(x[1]-x[0]);
    for(int i=1; i<n-1; ++i){
        dif0=x[i]-x[i-1];
        dif1=x[i+1]-x[i];
        dif2=x[i+1]-x[i-1];
        s[3*i+1]=y[i]/dif0;
        s[3*i+1]+=y[i]/dif1;
        s[3*i+1]-=s[3*i]/dif2;
        s[3*i+1]-=s[3*i]/dif0;
        s[3*i+1]+=s[3*(i+1)]/dif2;
        s[3*i+1]-=s[3*(i+1)]*dif1;
        s[3*i+1]*=2.0;
        s[3*i+2]=s[3*(i+1)]/(dif1*dif2);
        s[3*i+2]-=y[i]/(dif0*dif1);
        s[3*i+2]+=s[3*i]/(dif0*dif2);
        s[3*i+2]*=4.0;
    }
}

void splinefunc(int n, double *y, double *s, double *dy, double a, double b){
    double step=(b-a)/(n-1);
    s[0]=y[0];
    s[1]=dy[0];
    s[3]=4.0*(y[1]+y[0])+(dy[0]*step);
    s[4]=7.0;
    s[5]=1.0;
    for(int i=2;i<n-1;++i){
        s[3*i+1]=6.0-(1.0/s[3*i-2]);
        s[3*i+2]=1.0;
        s[3*i]=4.0*(y[i]+y[i-1])-(s[3*i-3]/s[3*i-2]);
    }
    s[3*n-2]=7.0-(1.0/s[3*n-5]);
    s[3*n-3]=4.0*(y[n-1]+y[n-2])-(dy[n-1]*step)-(s[3*n-6]/s[3*n-5]);
    s[3*n-3]/=s[3*n-2];
    s[3*n-1]=(s[3*n-3]-y[n-1]);
    s[3*n-1]*=2.0;
    s[3*n-1]/=step;
    s[3*n-1]+=dy[n-1];
    s[3*n-1]*=2.0;
    s[3*n-1]/=step;
    s[3*n-2]=4.0*(y[n-1]-s[3*n-3]);
    s[3*n-2]/=step;
    s[3*n-2]-=dy[n-1];
    for(int i=n-2;i>0;--i){
        s[3*i]-=s[3*(i+1)];
        s[3*i]/=s[3*i+1];
    }
    s[2]=s[3]-y[0];
    s[2]*=2.0;
    s[2]/=step;
    s[2]-=dy[0];
    s[2]*=2.0;
    s[2]/=step;
    for(int i=1; i<n-1; ++i){
        s[3*i+1]=4.0*y[i]-3.0*s[3*i]-s[3*(i+1)];
        s[3*i+1]/=step;
        s[3*i+2]=s[3*(i+1)]+s[3*i]-2.0*y[i];
        s[3*i+2]*=2.0;
        s[3*i+2]/=step;
        s[3*i+2]/=step;
    }
    s[0]-=s[1]*a;
    s[0]+=(s[2]*a*a);
    s[1]-=(2.0*s[2]*a);
    double point=a+step/2.0;
    for(int i=1; i<n; ++i){
        s[3*i]-=(s[3*i+1]*point);
        s[3*i]+=(s[3*i+2]*point*point);
        s[3*i+1]-=(2.0*s[3*i+2]*point);
        point+=step;
    }
}

double splinevaluen(double p, double *s, double * x, int i){
    double left=x[0];
    if(i)
        left=(x[i-1]+x[i])/2.0;
    return s[3*i]+s[3*i+1]*(p-left)+s[3*i+2]*(p-left)*(p-left);
}

double splinevaluen(double p, double *s, int i){
    return s[3*i]+s[3*i+1]*p+s[3*i+2]*p*p;
}
