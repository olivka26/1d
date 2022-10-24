#include <QPainter>
#include <stdio.h>
#include <stdlib.h>
#include "window.hpp"
#include "chebyshev.hpp"
#include "hermite.hpp"
#include "splines.hpp"
#include "math.h"
#include "general.hpp"

#define DEFAULT_A -10
#define DEFAULT_B 10
#define DEFAULT_N 30
#define DEFAULT_K 0

using namespace std;

Window::Window(QWidget *parent)
  : QWidget(parent){
      a=DEFAULT_A;
      b=DEFAULT_B;
      n=DEFAULT_N;
      k=DEFAULT_K;
}

Window::~Window(){
    free(x);
    free(y);
    if(c)
        free(c);
    if(cx)
        free(cx);
    if(cy)
        free(cy);
    free(dy);
    free(sp);
    free(h);
}

QSize Window::minimumSizeHint()const{
    return QSize(100, 100);
}

QSize Window::sizeHint()const{
    return QSize(1000, 1000);
}

int Window::parse_command_line(int argc, char *argv[]){
    if(argc!=5)
        return -1;
    if(sscanf(argv[1],"%lf",&a)!=1
       || sscanf(argv[2],"%lf",&b)!=1
       || sscanf(argv[3],"%d",&n)!=1
       || b-a<1.e-6
       || (argc>4 && sscanf(argv[4],"%d",&k)!=1)
       || k<0
       || k>6
       || n<=0)
        return -2;
    k=(k-1)%7;
    allocate();
    fillpoints(x, n, a, b);
    if(n<=50)
        chebyshevpoints(cx, n, a, b);
    delta_x=(b-a)/width();
    change_func();
    return 0;
}

void Window::allocate(){
    if(x)
        free(x);
    x=(double*)malloc(n*sizeof(double));
    if(n<=50){
        if(c)
            free(c);
        c=(double*)malloc(n*sizeof(double));
        if(cx)
            free(cx);
        cx=(double*)malloc(n*sizeof(double));
        if(cy)
            free(cy);
        cy=(double*)malloc(n*sizeof(double));
    }
    if(y)
        free(y);
    y=(double*)malloc(n*sizeof(double));
    if(h)
        free(h);
    h=(double*)malloc(4*(n-1)*sizeof(double));
    if(sp)
        free(sp);
    sp=(double*)malloc(3*n*sizeof(double));
    if(dy)
        free(dy);
    dy=(double*)malloc(n*sizeof(double));
}

void Window::chebyshevextrema(){
    double y1=chebyshevvalue(a, a, b, c, n);
    if(y1<extr[0])
        extr[0]=y1;
    if(y1>extr[1])
        extr[1]=y1;
    y1=chebyshevvalue(b, a, b, c, n);
    if(y1<extr[0])
        extr[0]=y1;
    if(y1>extr[1])
        extr[1]=y1;
    for(double x1=a+delta_x; x1-b<1e-6; x1+=delta_x){
        y1=chebyshevvalue(x1, a, b, c, n);
        if(y1<extr[0])
            extr[0]=y1;
        if(y1>extr[1])
            extr[1]=y1;
    }
}
    /*for(int i=1; i<=n; ++i){
        double x2, y2;
        if(i<n)
            x2=cx[i];
        else
            x2=b;
        for(double x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
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
    }*/

void Window::hermiteextrema(){
    for(int i=0;i<n-1;++i){
        double x1=x[i], x2=x[i+1], y1=hermitevaluen(x1, h, i);
        if(y1<extr[0])
            extr[0]=y1;
        if(y1>extr[1])
            extr[1]=y1;
        y1=hermitevaluen(x2, h, i);
        if(y1<extr[0])
            extr[0]=y1;
        if(y1>extr[1])
            extr[1]=y1;
        if(fabs(h[4*i+3])>1e-6){
            double discr=h[4*i+2]*h[4*i+2]-3*h[4*i+3]*h[4*i+1], extr_p;
            if(discr>1e-6){
                extr_=-h[4*i+2]-sqrt(discr);
                extr_p/=3;
                extr_p/=h[4*i+3];
                if(extr_p>=x1 && extr_p<=x2){
                    y1=hermitevaluen(extr_p, h, i);
                    if(y1<extr[0])
                        extr[0]=y1;
                    if(y1>extr[1])
                        extr[1]=y1;
                }
                extr_=-h[4*i+2]+sqrt(discr);
                extr_p/=3;
                extr_p/=h[4*i+3];
                if(extr_p>=x1 && extr_p<=x2){
                    y1=hermitevaluen(extr_p, h, i);
                    if(y1<extr[0])
                        extr[0]=y1;
                    if(y1>extr[1])
                        extr[1]=y1;
                }
            }
        }
    }
}
        /*for(double x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
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
            extr[1]=y2;*/

void Window::splineextrema(){
    for(int i=0;i<n;++i){
        double x1=x[i], x2, y1;
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
        y1=splinevaluen(x2, sp, i);
        if(y1<extr[0])
            extr[0]=y1;
        if(y1>extr[1])
            extr[1]=y1;
        if(fabs(sp[3*i+2])>1e-6){
            double extr_p=(-sp[3*i+1])/(2*sp[3*i+2]);
            if(extr_p>==x1 && extr_p<=x2){
                y1=splinevaluen(extr_p, sp, i);
                if(y1<extr[0])
                    extr[0]=y1;
                if(y1>extr[1])
                    extr[1]=y1;
            }
        }
    }
}
        /*for(double x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
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
            extr[1]=y2;*/

void Window::chebysheverrorextrema(){
    double y1=chebyshevvalue(a, a, b, c, n)-f(a);
    if(y1<extr[0])
        extr[0]=y1;
    if(y1>extr[1])
        extr[1]=y1;
    y1=chebyshevvalue(x[n/2], a, b, c, n)-y[n/2];
    if(y1<extr[0])
        extr[0]=y1;
    if(y1>extr[1])
        extr[1]=y1;
    y1=chebyshevvalue(b, a, b, c, n)-f(b);
    if(y1<extr[0])
        extr[0]=y1;
    if(y1>extr[1])
        extr[1]=y1;
    for(double x1=a+delta_x; x1-x[n/2]<1e-6; x1+=delta_x){
        y1=chebyshevvalue(x1, a, b, c, n)-f(x1);
        if(y1<extr[0])
            extr[0]=y1;
        if(y1>extr[1])
            extr[1]=y1;
    }
    for(double x1=x[n/2]+delta_x; x1-b<1e-6; x1+=delta_x){
        y1=chebyshevvalue(x1, a, b, c, n)-f(x1);
        if(y1<extr[0])
            extr[0]=y1;
        if(y1>extr[1])
            extr[1]=y1;
    }
}
    /*double x1=a, y1=f(x1)-chebyshevvalue(x1, a, b, c, n);
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
        for(double x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
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
    }*/

void Window::hermiteerrorextrema(){
    if(k=0){
        for(int i=0;i<n-1;++i){
            double x1=x[i], x2=x[i+1], y1=hermitevaluen(x1, h, i)-y[i];
            if(y1<extr[0])
                extr[0]=y1;
            if(y1>extr[1])
                extr[1]=y1;
            y1=hermitevaluen(x2, h, i)-y[i+1];
            if(y1<extr[0])
                extr[0]=y1;
            if(y1>extr[1])
                extr[1]=y1;
            if(fabs(h[4*i+3])>1e-6){
                double discr=h[4*i+2]*h[4*i+2]-3*h[4*i+3]*h[4*i+1], extr_p;
                if(discr>1e-6){
                    extr_=-h[4*i+2]-sqrt(discr);
                    extr_p/=3;
                    extr_p/=h[4*i+3];
                    if(extr_p>=x1 && extr_p<=x2){
                        y1=hermitevaluen(extr_p, h, i)-f(extr_p);
                        if(y1<extr[0])
                            extr[0]=y1;
                        if(y1>extr[1])
                            extr[1]=y1;
                    }
                    extr_=-h[4*i+2]+sqrt(discr);
                    extr_p/=3;
                    extr_p/=h[4*i+3];
                    if(extr_p>=x1 && extr_p<=x2){
                        y1=hermitevaluen(extr_p, h, i)-f(extr_p);
                        if(y1<extr[0])
                            extr[0]=y1;
                        if(y1>extr[1])
                            extr[1]=y1;
                    }
                }
            }
        }
        return;
    }
    if(k==1){
        for(int i=0;i<n-1;++i){
            double x1=x[i], x2=x[i+1], y1=hermitevaluen(x1, h, i)-y[i];
            if(y1<extr[0])
                extr[0]=y1;
            if(y1>extr[1])
                extr[1]=y1;
            y1=hermitevaluen(x2, h, i)-y[i+1];
            if(y1<extr[0])
                extr[0]=y1;
            if(y1>extr[1])
                extr[1]=y1;
            if(fabs(h[4*i+3])>1e-6){
                double discr=h[4*i+2]*h[4*i+2]-3*h[4*i+3]*(h[4*i+1]-1), extr_p;
                if(discr>1e-6){
                    extr_=-h[4*i+2]-sqrt(discr);
                    extr_p/=3;
                    extr_p/=h[4*i+3];
                    if(extr_p>=x1 && extr_p<=x2){
                        y1=hermitevaluen(extr_p, h, i)-f(extr_p);
                        if(y1<extr[0])
                            extr[0]=y1;
                        if(y1>extr[1])
                            extr[1]=y1;
                    }
                    extr_=-h[4*i+2]+sqrt(discr);
                    extr_p/=3;
                    extr_p/=h[4*i+3];
                    if(extr_p>=x1 && extr_p<=x2){
                        y1=hermitevaluen(extr_p, h, i)-f(extr_p);
                        if(y1<extr[0])
                            extr[0]=y1;
                        if(y1>extr[1])
                            extr[1]=y1;
                    }
                }
            }
        }
        return;
    }
    if(k==2){
        for(int i=0;i<n-1;++i){
            double x1=x[i], x2=x[i+1], y1=hermitevaluen(x1, h, i)-y[i];
            if(y1<extr[0])
                extr[0]=y1;
            if(y1>extr[1])
                extr[1]=y1;
            y1=hermitevaluen(x2, h, i)-y[i+1];
            if(y1<extr[0])
                extr[0]=y1;
            if(y1>extr[1])
                extr[1]=y1;
            if(fabs(h[4*i+3])>1e-6){
                double discr=(h[4*i+2]-1)*(h[4*i+2]-1)-3*h[4*i+3]*h[4*i+1], extr_p;
                if(discr>1e-6){
                    extr_=-h[4*i+2]+1-sqrt(discr);
                    extr_p/=3;
                    extr_p/=h[4*i+3];
                    if(extr_p>=x1 && extr_p<=x2){
                        y1=hermitevaluen(extr_p, h, i)-f(extr_p);
                        if(y1<extr[0])
                            extr[0]=y1;
                        if(y1>extr[1])
                            extr[1]=y1;
                    }
                    extr_=-h[4*i+2]+1+sqrt(discr);
                    extr_p/=3;
                    extr_p/=h[4*i+3];
                    if(extr_p>=x1 && extr_p<=x2){
                        y1=hermitevaluen(extr_p, h, i)-f(extr_p);
                        if(y1<extr[0])
                            extr[0]=y1;
                        if(y1>extr[1])
                            extr[1]=y1;
                    }
                }
            }
        }
        return;
    }
    if(k==3){
        for(int i=0;i<n-1;++i){
            double x1=x[i], x2=x[i+1], y1=hermitevaluen(x1, h, i)-y[i];
            if(y1<extr[0])
                extr[0]=y1;
            if(y1>extr[1])
                extr[1]=y1;
            y1=hermitevaluen(x2, h, i)-y[i+1];
            if(y1<extr[0])
                extr[0]=y1;
            if(y1>extr[1])
                extr[1]=y1;
            if(fabs(h[4*i+3]-1)>1e-6){
                double discr=h[4*i+2]*h[4*i+2]-3*(h[4*i+3]-1)*h[4*i+1], extr_p;
                if(discr>1e-6){
                    extr_=-h[4*i+2]-sqrt(discr);
                    extr_p/=3;
                    extr_p/=(h[4*i+3]-1);
                    if(extr_p>=x1 && extr_p<=x2){
                        y1=hermitevaluen(extr_p, h, i)-f(extr_p);
                        if(y1<extr[0])
                            extr[0]=y1;
                        if(y1>extr[1])
                            extr[1]=y1;
                    }
                    extr_=-h[4*i+2]+sqrt(discr);
                    extr_p/=3;
                    extr_p/=(h[4*i+3]-1);
                    if(extr_p>=x1 && extr_p<=x2){
                        y1=hermitevaluen(extr_p, h, i)-f(extr_p);
                        if(y1<extr[0])
                            extr[0]=y1;
                        if(y1>extr[1])
                            extr[1]=y1;
                    }
                }
            }
        }
        return;
    }
    for(int i=0;i<n-1;++i){
        double x1=x[i], x2=x[i+1],y1=hermitevaluen(x1, h, i)-y[i];
        if(y1<extr[0])
            extr[0]=y1;
        if(y1>extr[1])
            extr[1]=y1;
        for(double x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
            double y3=hermitevaluen(x3, h, i)-f(x3);
            if(y3<extr[0])
                extr[0]=y3;
            if(y3>extr[1])
                extr[1]=y3;
            x1=x3;
            y1=y3;
        }
        double y2=hermitevaluen(x2, h, i)-y[i];
        if(y2<extr[0])
            extr[0]=y2;
        if(y2>extr[1])
            extr[1]=y2;
    }
}
    /*for(int i=0;i<n-1;++i){
        double x1=x[i], x2=x[i+1],y1=y[i]-hermitevaluen(x1, h, i);
        if(y1<extr[0])
            extr[0]=y1;
        if(y1>extr[1])
            extr[1]=y1;
        for(double x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
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
    }*/

void Window::splineerrorextrema(){
    if(k==0){
        for(int i=0;i<n;++i){
            double x1=x[i], x2, y1;
            if(i)
                x1=(x[i]+x[i-1])/2.0;
            if(i!=n-1)
                x2=(x[i]+x[i+1])/2.0;
            else
                x2=x[i];
            double y1=splinevaluen(x1, sp, i)-f(x1);
            if(y1<extr[0])
                extr[0]=y1;
            if(y1>extr[1])
                extr[1]=y1;
            y1=splinevaluen(x2, sp, i)-f(x2);
            if(y1<extr[0])
                extr[0]=y1;
            if(y1>extr[1])
                extr[1]=y1;
            if(i==n/2+1){
                y1=splinevaluen(x[n/2], sp, i)-y[n/2];
                if(y1<extr[0])
                    extr[0]=y1;
                if(y1>extr[1])
                    extr[1]=y1;
            }
            if(fabs(sp[3*i+2])>1e-6){
                double extr_p=(-sp[3*i+1])/(2*sp[3*i+2]);
                if(extr_p>==x1 && extr_p<=x2){
                    y1=splinevaluen(extr_p, sp, i)-f(extr_p);
                    if(y1<extr[0])
                        extr[0]=y1;
                    if(y1>extr[1])
                        extr[1]=y1;
                }
            }
        }
        return;
    }
    if(k==1){
        for(int i=0;i<n;++i){
            double x1=x[i], x2, y1;
            if(i)
                x1=(x[i]+x[i-1])/2.0;
            if(i!=n-1)
                x2=(x[i]+x[i+1])/2.0;
            else
                x2=x[i];
            double y1=splinevaluen(x1, sp, i)-f(x1);
            if(y1<extr[0])
                extr[0]=y1;
            if(y1>extr[1])
                extr[1]=y1;
            y1=splinevaluen(x2, sp, i)-f(x2);
            if(y1<extr[0])
                extr[0]=y1;
            if(y1>extr[1])
                extr[1]=y1;
            if(i==n/2+1){
                y1=splinevaluen(x[n/2], sp, i)-y[n/2];
                if(y1<extr[0])
                    extr[0]=y1;
                if(y1>extr[1])
                    extr[1]=y1;
            }
            if(fabs(sp[3*i+2])>1e-6){
                double extr_p=(-sp[3*i+1]+1)/(2*sp[3*i+2]);
                if(extr_p>==x1 && extr_p<=x2){
                    y1=splinevaluen(extr_p, sp, i)-f(extr_p);
                    if(y1<extr[0])
                        extr[0]=y1;
                    if(y1>extr[1])
                        extr[1]=y1;
                }
            }
        }
        return;
    }
    if(k==2){
        for(int i=0;i<n;++i){
            double x1=x[i], x2, y1;
            if(i)
                x1=(x[i]+x[i-1])/2.0;
            if(i!=n-1)
                x2=(x[i]+x[i+1])/2.0;
            else
                x2=x[i];
            double y1=splinevaluen(x1, sp, i)-f(x1);
            if(y1<extr[0])
                extr[0]=y1;
            if(y1>extr[1])
                extr[1]=y1;
            y1=splinevaluen(x2, sp, i)-f(x2);
            if(y1<extr[0])
                extr[0]=y1;
            if(y1>extr[1])
                extr[1]=y1;
            if(i==n/2+1){
                y1=splinevaluen(x[n/2], sp, i)-y[n/2];
                if(y1<extr[0])
                    extr[0]=y1;
                if(y1>extr[1])
                    extr[1]=y1;
            }
            if(fabs(sp[3*i+2]-1)>1e-6){
                double extr_p=(-sp[3*i+1])/(2*sp[3*i+2]-2);
                if(extr_p>==x1 && extr_p<=x2){
                    y1=splinevaluen(extr_p, sp, i)-f(extr_p);
                    if(y1<extr[0])
                        extr[0]=y1;
                    if(y1>extr[1])
                        extr[1]=y1;
                }
            }
        }
        return;
    }
    if(k==3){
        for(int i=0;i<n;++i){
            double x1=x[i], x2, y1;
            if(i)
                x1=(x[i]+x[i-1])/2.0;
            if(i!=n-1)
                x2=(x[i]+x[i+1])/2.0;
            else
                x2=x[i];
            double y1=splinevaluen(x1, sp, i)-f(x1);
            if(y1<extr[0])
                extr[0]=y1;
            if(y1>extr[1])
                extr[1]=y1;
            y1=splinevaluen(x2, sp, i)-f(x2);
            if(y1<extr[0])
                extr[0]=y1;
            if(y1>extr[1])
                extr[1]=y1;
            if(i==n/2+1){
                y1=splinevaluen(x[n/2], sp, i)-y[n/2];
                if(y1<extr[0])
                    extr[0]=y1;
                if(y1>extr[1])
                    extr[1]=y1;
            }
            double discr=sp[3*i+2]*sp[3*i+2]+3*sp[3*i+1];
            if(fabs(discr)>=1e-6){
                double extr_p=(sp[3*i+2]-sqrt(discr))/3;
                if(extr_p>==x1 && extr_p<=x2){
                    y1=splinevaluen(extr_p, sp, i)-f(extr_p);
                    if(y1<extr[0])
                        extr[0]=y1;
                    if(y1>extr[1])
                        extr[1]=y1;
                }
                extr_p=(sp[3*i+2]+sqrt(discr))/3;
                if(extr_p>==x1 && extr_p<=x2){
                    y1=splinevaluen(extr_p, sp, i)-f(extr_p);
                    if(y1<extr[0])
                        extr[0]=y1;
                    if(y1>extr[1])
                        extr[1]=y1;
                }
            }
        }
        return;
    }
    for(int i=0;i<n;++i){
        double x1=x[i], x2;
        if(i)
            x1=(x[i]+x[i-1])/2.0;
        if(i!=n-1)
            x2=(x[i]+x[i+1])/2.0;
        else
            x2=x[i];
        double y1;
        if(i==n/2+1){
            y1=splinevaluen(x[n/2], sp, i)-y[n/2];
            if(y1<extr[0])
                extr[0]=y1;
            if(y1>extr[1])
                extr[1]=y1;
        }
        y1=splinevaluen(x1, sp, i)-f(x1);
        if(y1<extr[0])
            extr[0]=y1;
        if(y1>extr[1])
            extr[1]=y1;
        for(double x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
            double y3=splinevaluen(x3, sp, i)-f(x3);
            if(y3<extr[0])
                extr[0]=y3;
            if(y3>extr[1])
                extr[1]=y3;
            x1=x3;
            y1=y3;
        }
        double y2=splinevaluen(x2, sp, i)-f(x2);
        if(y2<extr[0])
            extr[0]=y2;
        if(y2>extr[1])
            extr[1]=y2;
    }
}

void Window::extrema_hunt(){
    if(view_id!=4){
        if(p>0){
            extr[1]=y[n/2];
            if(extr[1]>max_y)
                extr[1]=max_y;
            extr[0]=min_y;
        }else if(p<0){
            extr[0]=y[n/2];
            if(extr[0]<min_y)
                extr[0]=min_y;
            extr[1]=max_y;
        }else{
            extr[0]=min_y;
            extr[1]=max_y;
        }
    }else{
        extr[0]=0;
        extr[1]=0;
    }
    if((view_id==0 || view_id==3) && n<=50)
        chebyshevextrema();
    if(view_id==1 || view_id==3)
        hermiteextrema();
    if(view_id==2 || view_id==3)
        splineextrema();
    if(view_id==4 && n<=50)
        chebysheverrorextrema();
    if(view_id==4){
        hermiteerrorextrema();
        splineerrorextrema();
    }
}

void Window::print_console(){
    printf("format: %d\n", view_id);
    printf("segment: [%lf;%lf]\n", a, b);
    printf("squeeze-stretch: %d\n", s);
    printf("points: %d\n", n);
    printf("disturbance: %d\n", p);
    printf("function absmax: %10.3e\n", absmax);
    if(k==0 && p==0 && view_id!=4){
        printf("factic absmax: %10.3e\n", absmax);
        printf("extrema: %10.3e %10.3e\n\n", absmax, absmax);
    }else{
        printf("factic absmax: %10.3e\n", max(fabs(extr[0]), fabs(extr[1])));
        printf("extrema: %10.3e %10.3e\n\n", extr[0], extr[1]);
    }
}

void Window::change_func(){
    k=(k+1)%7;
    switch(k){
        case 0:
            f_name="k=0 f(x)=1";
            f=f0;
            df=df0;
            d2f=d2f0;
            min_y=1;
            max_y=1;
            absmax=1;
            break;
        case 1:
            f_name="k=1 f(x)=x";
            f=f1;
            df=df1;
            d2f=d2f1;
            min_y=a;
            max_y=b;
            absmax=max(fabs(min_y), fabs(max_y));
            break;
        case 2:
            f_name="k=2 f(x)=x^2";
            f=f2;
            df=df2;
            d2f=d2f2;
            min_y=min(a*a, b*b);
            max_y=max(b*b, a*a);
            if(a*b<0){
                min_y=0;
            }
            absmax=max_y;
            break;
        case 3:
            f_name="k=3 f(x)=x^3";
            f=f3;
            df=df3;
            d2f=d2f3;
            min_y=a*a*a;
            max_y=b*b*b;
            absmax=max(fabs(min_y), fabs(max_y));
            break;
        case 4:
            f_name="k=4 f(x)=x^4";
            f=f4;
            df=df4;
            d2f=d2f4;
            min_y=min(a*a*a*a, b*b*b*b);
            max_y=max(b*b*b*b, a*a*a*a);
            if(a*b<0){
                min_y=0;
            }
            absmax=max_y;
            break;
        case 5:
            f_name="k=5 f(x)=e^x";
            f=f5;
            df=df5;
            d2f=d2f5;
            min_y=exp(a);
            max_y=exp(b);
            absmax=max_y;
            break;
        case 6:
            f_name="k=6 f(x)=1/(25x^2+1)";
            f=f6;
            df=df6;
            d2f=d2f6;
            if(a>=0 && a<=0){
                 min_y=1/(25*b*b+1);
                 max_y=1;
                 absmax=1;
             }else if(b>=0 && b<=0){
                 min_y=1/(25*a*a+1);
                 max_y=1;
                 absmax=1;
             }else if(a*b>0){
                 min_y=min(1/(25*a*a+1), 1/(25*b*b+1));
                 max_y=max(1/(25*a*a+1), 1/(25*b*b+1));
                 absmax=max_y;
             }else if(a*b<0){
                 min_y=min(1/(25*a*a+1), 1/(25*b*b+1));
                 max_y=1;
                 absmax=1;
             }
            break;
    }
    fillvalues(x, y, dy, n, f, df);
    if(p)
        y[n/2]+=(p*0.1*absmax);
    if(n<=50){
        chebyshevf(cx, cy, n, f);
        for(int i=0; i<n; ++i){
            if(fabs(cx[i]-x[n/2])<1e-6){
                cy[i]+=(p*0.1*absmax);
                break;
            }
        }
        chebyshevfunc(c, cy, n);
    }
    hermitefunc(n, x, y, dy, h);
    //splinefunc(n, x, y, sp, dy);
    splinefunc(n, y, sp, dy, a, b);
    extrema_hunt();
    print_console();
    update();
}

void Window::change_view(){
    view_id=(view_id+1)%5;
    extrema_hunt();
    print_console();
    update();
}

void Window::enhance(){
    double t=(b-a)/2;
    a=a-t;
    b=b+t;
    ++s;
    delta_x=(b-a)/width();
    fillpoints(x, n, a, b);
    fillvalues(x, y, dy, n, f, df);
    if(p)
        y[n/2]+=(p*0.1*absmax);
    if(n<=50){
        chebyshevpoints(cx, n, a, b);
        chebyshevf(cx,cy,n,f);
        for(int i=0; i<n; ++i){
            if(fabs(cx[i]-x[n/2])<1e-6){
                cy[i]+=(p*0.1*absmax);
                break;
            }
        }
        chebyshevfunc(c, cy, n);
    }
    hermitefunc(n, x, y, dy, h);
    //splinefunc(n, x, y, sp, dy);
    splinefunc(n, y, sp, dy, a, b);
    extrema_hunt();
    print_console();
    update();
}

void Window::squeeze(){
    double t=(b-a)/4;
    a=a+t;
    b=b-t;
    --s;
    delta_x=(b-a)/width();
    fillpoints(x, n, a, b);
    fillvalues(x, y, dy, n, f, df);
    if(p)
        y[n/2]+=(p*0.1*absmax);
    if(n<=50){
        chebyshevpoints(cx, n, a, b);
        chebyshevf(cx,cy,n,f);
        for(int i=0; i<n; ++i){
            if(fabs(cx[i]-x[n/2])<1e-6){
                cy[i]+=(p*0.1*absmax);
                break;
            }
        }
        chebyshevfunc(c, cy, n);
    }
    hermitefunc(n, x, y, dy, h);
    //splinefunc(n, x, y, sp, dy);
    splinefunc(n, y, sp, dy, a, b);
    extrema_hunt();
    print_console();
    update();
}

void Window::more_points(){
    n*=2;
    allocate();
    fillpoints(x, n, a, b);
    fillvalues(x, y, dy, n, f, df);
    if(p)
        y[n/2]+=(p*0.1*absmax);
    if(n<=50){
        chebyshevpoints(cx, n, a, b);
        chebyshevf(cx,cy,n,f);
        for(int i=0; i<n; ++i){
            if(fabs(cx[i]-x[n/2])<1e-6){
                cy[i]+=(p*0.1*absmax);
                break;
            }
        }
        chebyshevfunc(c, cy, n);
    }
    hermitefunc(n, x, y, dy, h);
    //splinefunc(n, x, y, sp, dy);
    splinefunc(n, y, sp, dy, a, b);
    extrema_hunt();
    print_console();
    update();
}

void Window::less_points(){
    n/=2;
    allocate();
    fillpoints(x, n, a, b);
    fillvalues(x, y, dy, n, f, df);
    if(p)
        y[n/2]+=(p*0.1*absmax);
    if(n<=50){
        chebyshevpoints(cx, n, a, b);
        chebyshevf(cx,cy,n,f);
        for(int i=0; i<n; ++i){
            if(fabs(cx[i]-x[n/2])<1e-6){
                cy[i]+=(p*0.1*absmax);
                break;
            }
        }
        chebyshevfunc(c, cy, n);
    }
    hermitefunc(n, x, y, dy, h);
    //splinefunc(n, x, y, sp, dy);
    splinefunc(n, y, sp, dy, a, b);
    extrema_hunt();
    print_console();
    update();
}

void Window::delta_up(){
    ++p;
    y[n/2]+=(0.1*absmax);
    if(n<=50){
        for(int i=0; i<n; ++i){
            if(fabs(cx[i]-x[n/2])<1e-6){
                cy[i]+=(0.1*absmax);
                break;
            }
        }
        chebyshevfunc(c, cy, n);
    }
    hermitefunc(n, x, y, dy, h);
    //splinefunc(n, x, y, sp, dy);
    splinefunc(n, y, sp, dy, a, b);
    extrema_hunt();
    print_console();
    update();
}

void Window::delta_down(){
    --p;
    y[n/2]-=(0.1*absmax);
    if(n<=50){
        for(int i=0; i<n; ++i){
            if(fabs(cx[i]-x[n/2])<1e-6){
                cy[i]-=(0.1*absmax);
                break;
            }
        }
        chebyshevfunc(c, cy, n);
    }
    hermitefunc(n, x, y, dy, h);
    //splinefunc(n, x, y, sp, dy);
    splinefunc(n, y, sp, dy, a, b);
    extrema_hunt();
    print_console();
    update();
}

void Window::paintEvent(QPaintEvent* /*event*/){
    QPainter painter(this);
    double x1, x2, y1, y2, x3, y3;
    double delta_x=(b-a)/width();
    QPen pen_black(Qt::black, 0, Qt::SolidLine);
    QPen pen_red(Qt::red, 0, Qt::SolidLine);
    QPen pen_green(Qt::green, 0, Qt::SolidLine);
    QPen pen_blue(Qt::blue, 0, Qt::SolidLine);
    QPen pen_cyan(Qt::cyan, 0, Qt::SolidLine);
    painter.save();
    // make Coordinate Transformations
    painter.translate(0.5*width(), 0.5*height());
    if(k==0 && p==0 && view_id!=4){
        painter.scale(width()/(b-a), -height()/2);
        painter.translate(-0.5*(a+b), -1);
    }else{
        painter.scale(width()/(b-a), -height()/(extr[1]-extr[0]));
        painter.translate(-0.5*(a+b), -0.5*(extr[1]+extr[0]));
    }
    painter.setPen(pen_red);
    painter.drawLine(QPointF(a, 0), QPointF(b, 0));
    if(k==0 && p==0 && view_id!=4)
        painter.drawLine(QPointF(0, 0), QPointF(0, 2));
    else
        painter.drawLine(QPointF(0, extr[0]), QPointF(0, extr[1]));
    if(view_id!=4){
        painter.setPen(pen_black);
        for(int i=0;i<n-1;++i){
            x1=x[i];
            x2=x[i+1];
            y1=y[i];
            for(x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
                y3=f(x3);
                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=y[i+1];
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
        }
    }
    if((view_id==0 || view_id==3) && n<=50){
        painter.setPen(pen_cyan);
        x1=a;
        y1=chebyshevvalue(x1, a, b, c, n);
        for(int i=1; i<=n; ++i){
            if(i<n)
                x2=cx[i];
            else
                x2=b;
            for(x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
                y3=chebyshevvalue(x3, a, b, c, n);
                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=chebyshevvalue(x2, a, b, c, n);
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
            x1=x2;
            y1=y2;
        }
    }
    if(view_id==1 || view_id==3){
        painter.setPen(pen_green);
        for(int i=0;i<n-1;++i){
            x1=x[i];
            x2=x[i+1];
            y1=hermitevaluen(x1, h, i);
            for(x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
                y3=hermitevaluen(x3, h, i);
                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=hermitevaluen(x2, h, i);
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
        }
    }
    if(view_id==2 || view_id==3){
        painter.setPen(pen_blue);
        for(int i=0;i<n;++i){
            x1=x[i];
            if(i)
                x1=(x[i]+x[i-1])/2.0;
            if(i!=n-1)
                x2=(x[i]+x[i+1])/2.0;
            else
                x2=x[i];
            y1=splinevaluen(x1, sp, i);
            for(x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
                y3=splinevaluen(x3, sp, i);
                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=splinevaluen(x2, sp, i);
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
        }
    }
    /*if((view_id==0 || view_id==3) && n<=50){
        painter.setPen(pen_cyan);
        x1=a;
        y1=chebyshevvalue(x1, a, b, c, n);
        for(int i=1; i<=n; ++i){
            if(i<n)
                x2=cx[i];
            else
                x2=b;
            for(x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
                y3=chebyshevvalue(x3, a, b, c, n);
                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=chebyshevvalue(x2, a, b, c, n);
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
            x1=x2;
            y1=y2;
        }
    }*/
    if(view_id==4 && n<=50){
        painter.setPen(pen_cyan);
        x1=a;
        y1=chebyshevvalue(x1, a, b, c, n)-f(x1);
        for(int i=1; i<=n; ++i){
            if(i<n)
                x2=cx[i];
            else
                x2=b;
            for(x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
                y3=f(x3)-chebyshevvalue(x3, a, b, c, n);
                if(fabs(x3-x[n/2])<1e-6)
                    y3+=(p*0.1*absmax);
                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=chebyshevvalue(x2, a, b, c, n)-f(x2);
            if(fabs(x3-y[n/2])<1e-6)
                y2+=(p*0.1*absmax);
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
            x1=x2;
            y1=y2;
        }
    }
    if(view_id==4){
        painter.setPen(pen_green);
        for(int i=0;i<n-1;++i){
            x1=x[i];
            x2=x[i+1];
            y1=hermitevaluen(x1, h, i)-y[i];
            for(x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
                y3=hermitevaluen(x3, h, i)-f(x3);
                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=hermitevaluen(x2, h, i)-y[i+1];
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
        }
        painter.setPen(pen_blue);
        for(int i=0;i<n;++i){
            x1=x[i];
            if(i)
                x1=(x[i]+x[i-1])/2.0;
            if(i!=n-1)
                x2=(x[i]+x[i+1])/2.0;
            else
                x2=x[i];
            y1=splinevaluen(x1, sp, i)-f(x1);
            for(x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
                y3=splinevaluen(x3, sp, i)-f(x3);
                if(fabs(x3-x[n/2])<1e-6)
                    y3+=(p*0.1*absmax);
                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=splinevaluen(x2, sp, i)-f(x2);
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
        }
    }
    /*if(view_id==4 && n<=50){
        painter.setPen(pen_cyan);
        x1=a;
        y1=f(x1)-chebyshevvalue(x1, a, b, c, n);
        for(int i=1; i<=n; ++i){
            if(i<n)
                x2=cx[i];
            else
                x2=b;
            for(x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
                y3=f(x3)-chebyshevvalue(x3, a, b, c, n);
                if(fabs(x3-x[n/2])<1e-6)
                    y3+=(p*0.1*absmax);
                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=f(x2)-chebyshevvalue(x2, a, b, c, n);
            if(fabs(x3-y[n/2])<1e-6)
                y2+=(p*0.1*absmax);
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
            x1=x2;
            y1=y2;
        }
    }*/
    painter.restore();
    painter.setPen("black");
    painter.drawText(0, 20, f_name);
    painter.drawText(10, 30, QString("format: %1").arg(view_id));
    painter.drawText(10, 45, QString("scale: %1 %2 %3").arg(s).arg(a).arg(b));
    painter.drawText(10, 60, QString("points: %1").arg(n));
    painter.drawText(10, 75, QString("p: %1").arg(p));
    if(k==0 && p==0 && view_id!=4)
        painter.drawText(10, 90, QString("absmax(fact): %1( %2)").arg(absmax).arg(absmax));
    else
        painter.drawText(10, 90, QString("absmax(fact): %1( %2)").arg(absmax).arg(max(fabs(extr[0]/*+delta_y*/), fabs(extr[1]/*-delta_y*/))));
}
