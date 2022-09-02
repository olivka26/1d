#include "chebyshev.hpp"

void chebyshevpoints(double *cx, double *cy, int n, double a, double b, double (*f)(double)){
    double step=3.1415926535897932384626433832795/n;
    double angle=step/2.0;
    double semisum=(a+b)/2;
    double semidif=(b-a)/2;
    for(int i=n-1;i>=0;--i){
        cx[i]=cos(angle);
        cx[i]*=semidif;
        cx[i]+=semisum;
        cy[i]=f(cx[i]);
        angle+=step;
    }
}

double chebyshevtrans(double a, double b, double x){
    return (2*x-(b+a))/(b-a);
}

void chebyshevfunc(double *c, double (*f)(double), int n, double a, double b){
    double step=3.1415926535897932384626433832795/n;
    double angle=step/2.0;
    double semisum=(a+b)/2;
    double semidif=(b-a)/2;
    for(int i=0;i<n;++i){
        c[i]=0;
    }
    for(int i=n-1;i>=0;--i){
        double xx=cos(angle);
        double res1=f(xx*semidif+semisum);
        double res2=res1*xx;
        double res3=2*xx*res2-res1;
        c[0]+=res1;
        c[1]+=res2;
        c[2]+=res3;
        for(int j=3;j<n;++j){
            res1=res2;
            res2=res3;
            res3=2*xx*res2-res1;
            c[j]+=(res3);
        }
        angle+=step;
    }
    for(int i=0;i<n;++i){
        c[i]/=n;
        if(i)
            c[i]*=2;
    }
}

void chebyshevfunc(double *c, double *cy, int n){
    double step=3.1415926535897932384626433832795/n;
    double angle=step/2.0;
    for(int i=0;i<n;++i){
        c[i]=0;
    }
    for(int i=n-1;i>=0;--i){
        double xx=cos(angle);
        double res1=cy[i];
        double res2=res1*xx;
        double res3=2*xx*res2-res1;
        c[0]+=res1;
        c[1]+=res2;
        c[2]+=res3;
        for(int j=3;j<n;++j){
            res1=res2;
            res2=res3;
            res3=2*xx*res2-res1;
            c[j]+=(res3);
        }
        angle+=step;
    }
    for(int i=0;i<n;++i){
        c[i]/=n;
        if(i)
            c[i]*=2;
    }
}


double chebyshevvalue(double x, double a, double b, double *c, int n){
    double t0=1, t1=chebyshevtrans(a,b,x), cheb=2*t1, t2=t1*cheb-t0, res=c[0];
    res+=(c[1]*t1);
    res+=(c[2]*t2);
    t0=t2*cheb-t1;
    t1=t0*cheb-t2;
    t2=t1*cheb-t0;
    for(int i=3;i<n;++i){
        if(i%3==0){
            res+=c[i]*t0;
        }else if(i%3==1){
            res+=c[i]*t1;
        }else{
            res+=c[i]*t2;
            t0=t2*cheb-t1;
            t1=t0*cheb-t2;
            t2=t1*cheb-t0;
        }
    }
    return res;
}

/*double chebyshevvalue(double x, double a, double b, double *c, int n){
    double t0=1,cheb=chebyshevtrans(a,b,x),t1=cheb,t2=2*t1*t1-t0, res=c[0];
    res+=(c[1]*t1);
    res+=(c[2]*t2);
    t0=2*t2*cheb-t1;
    t1=2*t0*cheb-t2;
    t2=2*t1*cheb-t0;
    for(int i=3;i<n;++i){
        if(i%3==0){
            res+=c[i]*t0;
        }else if(i%3==1){
            res+=c[i]*t1;
        }else{
            res+=c[i]*t2;
            t0=2*t2*cheb-t1;
            t1=2*t0*cheb-t2;
            t2=2*t1*cheb-t0;
        }
    }
    return res;
}*/
