#include "general.hpp"

double f0(double x){
    if(x>0)
        return 1;
    return 1;
}

double f1(double x){
    return x;
}

double f2(double x){
    return x*x;
}

double f3(double x){
    return x*x*x;
}

double f4(double x){
    return x*x*x*x;
}

double f5(double x){
    return exp(x);
}

double f6(double x){
    return 1/(25*x*x+1);
}

double df0(double x){
    if(x>0)
        return 0;
    return 0;
}

double df1(double x){
    if(x>0)
        return 1;
    return 1;
}

double df2(double x){
    return 2*x;
}

double df3(double x){
    return 3*x*x;
}

double df4(double x){
    return 4*x*x*x;
}

double df5(double x){
    return exp(x);
}

double df6 (double x){
    return (-50*x)/((25*x*x+1)*(25*x*x+1));
}

double d2f0(double x){
    if(x>0)
        return 0;
    return 0;
}

double d2f1(double x){
    if(x>0)
        return 0;
    return 0;
}

double d2f2(double x){
    if(x>0)
        return 2;
    return 2;
}

double d2f3(double x){
    return 6*x;
}

double d2f4(double x){
    return 12*x*x;
}

double d2f5(double x){
    return exp(x);
}

double d2f6 (double x){
    return (3750*x*x-50)/((25*x*x+1)*(25*x*x+1)*(25*x*x+1));
}

double min(double a, double b){
    return (b < a) ? b : a;
}

double max(double a, double b){
    return (b > a) ? b : a;
}

void fillpoints(double *x, int n, double a, double b){
    x[0]=a;
    double step=(b-a)/(n-1);
    for(int i=1; i<n; ++i){
        x[i]=x[i-1]+step;
    }
}

void fillvalues(double *x, double *y, double *dy, int n, double (*f)(double), double (*df)(double)){
    for(int i=0; i<n; ++i){
        y[i]=f(x[i]);
        dy[i]=df(x[i]);
    }
}

void chebyshevextrema(double *c, double *cx, int n, double a, double b, double *extr){
    double x1=a, y1=chebyshevvalue(x1, a, b, c, n), delta_x=(b-a)/(n-1);
    if(y1<extr[0])
        extr[0]=y1;
    if(y1>extr[1])
        extr[1]=y1;
    for(int i=1; i<=n; ++i){
        double x2, y2;
        if(i<n)
            x2=cx[i];
        else
            x2=b;
        for(double x3=x1+delta_x*0.0001; x3-x2<1e-6; x3+=delta_x*0.0001){
            double y3=chebyshevvalue(x3, a, b, c, n);
            if(y3<extr[0])
                extr[0]=y3;
            if(y3>extr[1])
                extr[1]=y3;
            x1=x3;
            y1=y3;
        }
        y2=chebyshevvalue(x2, a, b, c, n);
        if(y2<extr[0])
             extr[0]=y2;
        if(y2>extr[1])
             extr[1]=y2;
        x1=x2;
        y1=y2;
    }
}

void hermiteextrema(double *h, double *x, int n, double *extr){
    double delta_x=(x[n-1]-x[0])/(n-1);
    for(int i=0;i<n-1;++i){
        double x1=x[i], x2=x[i+1], y1=hermitevaluen(x1, h, i);
        if(y1<extr[0])
            extr[0]=y1;
        if(y1>extr[1])
            extr[1]=y1;
        for(double x3=x1+delta_x*1e-4; x3-x2<1e-6; x3+=delta_x*0.0001){
            double y3=hermitevaluen(x3, h, i);
            if(y3<extr[0])
                extr[0]=y3;
            if(y3>extr[1])
                extr[1]=y3;
            x1=x3;
            y1=y3;
        }
        double y2=hermitevaluen(x2, h, i);
        if(y2<extr[0])
            extr[0]=y2;
        if(y2>extr[1])
            extr[1]=y2;
    }
}

void splineextrema(double *sp, double *x, int n, double *extr){
    double delta_x=(x[n-1]-x[0])/(n-1);
    for(int i=0;i<n;++i){
        double x1=x[i], x2;
        if(i)
            x1=(x[i]+x[i-1])/2.0;
        if(i!=n-1)
            x2=(x[i]+x[i+1])/2.0;
        else
            x2=x[i];
        double y1=splinevaluen(x1, sp, i);
        if(y1<extr[0])
            extr[0]=y1;
        if(y1>extr[1])
            extr[1]=y1;
        for(double x3=x1+delta_x*1e-4; x3-x2<1e-6; x3+=delta_x*0.0001){
            double y3=splinevaluen(x3, sp, i);
            if(y3<extr[0])
                extr[0]=y3;
            if(y3>extr[1])
                extr[1]=y3;
            x1=x3;
            y1=y3;
        }
        double y2=splinevaluen(x2, sp, i);
        if(y2<extr[0])
            extr[0]=y2;
        if(y2>extr[1])
            extr[1]=y2;
    }
}

void chebysheverrorextrema(double *c, double *cx, int n, double a, double b, double *extr, double (*f)(double)){
    double x1=a, y1=f(x1)-chebyshevvalue(x1, a, b, c, n), delta_x=(b-a)/(n-1);
    if(y1<extr[0])
        extr[0]=y1;
    if(y1>extr[1])
        extr[1]=y1;
    for(int i=1; i<=n; ++i){
        double x2, y2;
        if(i<n)
            x2=cx[i];
        else
            x2=b;
        for(double x3=x1+delta_x*0.0001; x3-x2<1e-6; x3+=delta_x*0.0001){
            double y3=f(x3)-chebyshevvalue(x3, a, b, c, n);
            if(y3<extr[0])
                extr[0]=y3;
            if(y3>extr[1])
                extr[1]=y3;
            x1=x3;
            y1=y3;
        }
        y2=f(x2)-chebyshevvalue(x2, a, b, c, n);
        if(y2<extr[0])
             extr[0]=y2;
        if(y2>extr[1])
             extr[1]=y2;
        x1=x2;
        y1=y2;
    }
}

void hermiteerrorextrema(double *h, double *x, double *y, int n, double *extr, double (*f)(double)){
    double delta_x=(x[n-1]-x[0])/(n-1);
    for(int i=0;i<n-1;++i){
        double x1=x[i], x2=x[i+1],y1=y[i]-hermitevaluen(x1, h, i);
        if(y1<extr[0])
            extr[0]=y1;
        if(y1>extr[1])
            extr[1]=y1;
        for(double x3=x1+delta_x*1e-4; x3-x2<1e-6; x3+=delta_x*0.0001){
            double y3=f(x3)-hermitevaluen(x3, h, i);
            if(y3<extr[0])
                extr[0]=y3;
            if(y3>extr[1])
                extr[1]=y3;
            x1=x3;
            y1=y3;
        }
        double y2=y[i+1]-hermitevaluen(x2, h, i);
        if(y2<extr[0])
            extr[0]=y2;
        if(y2>extr[1])
            extr[1]=y2;
    }
}

void splineerrorextrema(double *sp, double *x, int n, double *extr, double (*f)(double)){
    double delta_x=(x[n-1]-x[0])/(n-1);
    for(int i=0;i<n;++i){
        double x1=x[i], x2;
        if(i)
            x1=(x[i]+x[i-1])/2.0;
        if(i!=n-1)
            x2=(x[i]+x[i+1])/2.0;
        else
            x2=x[i];
        double y1=f(x1)-splinevaluen(x1, sp, i);
        if(y1<extr[0])
            extr[0]=y1;
        if(y1>extr[1])
            extr[1]=y1;
        for(double x3=x1+delta_x*1e-4; x3-x2<1e-6; x3+=delta_x*0.0001){
            double y3=f(x3)-splinevaluen(x3, sp, i);
            if(y3<extr[0])
                extr[0]=y3;
            if(y3>extr[1])
                extr[1]=y3;
            x1=x3;
            y1=y3;
        }
        double y2=f(x2)-splinevaluen(x2, sp, i);
        if(y2<extr[0])
            extr[0]=y2;
        if(y2>extr[1])
            extr[1]=y2;
    }
}
