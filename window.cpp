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

static double f0(double x){
    return 1;
}

static double f1(double x){
    return x;
}

static double f2(double x){
    return x*x;
}

static double f3(double x){
    return x*x*x;
}

static double f4(double x){
    return x*x*x*x;
}

static double f5(double x){
    return exp(x);
}

static double f6(double x){
    return 1/(25*x*x+1);
}

static double df0(double x){
    return 0;
}

static double df1(double x){
    return 1;
}

static double df2(double x){
    return 2*x;
}

static double df3(double x){
    return 3*x*x;
}

static double df4(double x){
    return 4*x*x*x;
}

static double df5(double x){
    return exp(x);
}

static double df6 (double x){
    return (-50*x)/((25*x*x+1)*(25*x*x+1));
}

static double d2f0(double x){
    return 0;
}

static double d2f1(double x){
    return 0;
}

static double d2f2(double x){
    return 2;
}

static double d2f3(double x){
    return 6*x;
}

static double d2f4(double x){
    return 12*x*x;
}

static double d2f5(double x){
    return exp(x);
}

static double d2f6 (double x){
    return (3750*x*x-50)/((25*x*x+1)*(25*x*x+1)*(25*x*x+1));
}

Window::Window(QWidget *parent)
  : QWidget(parent){
      a=DEFAULT_A;
      b=DEFAULT_B;
      n=DEFAULT_N;
      k=DEFAULT_K;
      change_func();
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
    change_func();
    return 0;
}

void Window::allocate(){
    x=(double*)malloc(n*sizeof(double));
    if(n<=50){
        c=(double*)malloc(n*sizeof(double));
        cx=(double*)malloc(n*sizeof(double));
        cy=(double*)malloc(n*sizeof(double));
    }
    y=(double*)malloc(n*sizeof(double));
    if(view_id==1 || view_id==3 || view_id==4)
        h=(double*)malloc(4*(n-1)*sizeof(double));
    if(view_id==2 || view_id==3 || view_id==4)
        sp=(double*)malloc(3*n*sizeof(double));
    dy=(double*)malloc(n*sizeof(double));
    extr=(double*)malloc(2*sizeof(double));
    extr[0]=min_y;
    extr[1]=max_y;
}

void Window::print_console(double delta_y){
    printf("format: %d\n", view_id);
    printf("segment: [%lf;%lf]\n", a, b);
    printf("squeeze-stretch: %d\n", s);
    printf("points: %d\n", n);
    printf("disturbance: %d\n", p);
    printf("function absmax: %lf\n", absmax);
    if(k==0 && p==0 && view_id!=4){
        printf("factic absmax: %lf\n", absmax);
        printf("extrema: %lf %lf\n\n", absmax);
    }else{
        printf("factic absmax: %lf\n", max(fabs(extr[0]+delta_y), fabs(extr[1]-delta_y)));
        printf("extrema: %lf %lf\n\n", extr[0]+delta_y, extr[1]-delta_y);
    }
}

void Window::destroy(){
    free(x);
    free(y);
    if(n<=50 && (view_id==0 || view_id==3 || view_id==4)){
        free(c);
        free(cx);
        free(cy);
    }
    free(dy);
    if(view_id==2 || view_id==3 || view_id==4)
       free(sp);
    if(view_id==1 || view_id==3 || view_id==4)
       free(h);
    free(extr);
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
    update();
}

void Window::change_view(){
    view_id=(view_id+1)%5;
    update();
}

void Window::enhance(){
    double t=(b-a)/2;
    a=a-t;
    b=b+t;
    ++s;
    update();
}

void Window::squeeze(){
    double t=(b-a)/4;
    a=a+t;
    b=b-t;
    --s;
    update();
}

void Window::more_points(){
    n*=2;
    update();
}

void Window::less_points(){
    n/=2;
    update();
}

void Window::delta_up(){
    ++p;
    update();
}

void Window::delta_down(){
    --p;
    update();
}

void Window::paintEvent(QPaintEvent *event){
    QPainter painter(this);
    double x1, x2, y1, y2, x3, y3;
    double delta_y, delta_x=(b-a)/(n-1);
    delta_x/=10000;
    QPen pen_black(Qt::black, 0, Qt::SolidLine);
    QPen pen_red(Qt::red, 0, Qt::SolidLine);
    QPen pen_green(Qt::green, 0, Qt::SolidLine);
    QPen pen_blue(Qt::blue, 0, Qt::SolidLine);
    QPen pen_cyan(Qt::cyan, 0, Qt::SolidLine);
    allocate();
    fillpoints(x, y, dy, n, a, b, f, df);
    if(p!=0){
        y[n/2]+=(p*0.1*absmax);
        if(p>0 && y[n/2]>max_y){
            max_y=y[n/2];
            extr[1]=y[n/2];
        }
        if(p<0 && y[n/2]<min_y){
            min_y=y[n/2];
            extr[0]=y[n/2];
        }
        absmax=max(fabs(min_y), fabs(max_y));
    }
    if((view_id==0 || view_id==3 || view_id==4) && n<=50){
        chebyshevpoints(cx, cy, n, a, b, f);
        for(int i=0; i<n; ++i){
            if(cx[i]<=x[n/2] && cx[i]>=x[n/2]){
                cy[i]+=(p*0.1*max_y);
                if(cy[i]>extr[1])
                    extr[1]=cy[i];
                if(cy[i]<extr[0])
                    extr[0]=cy[i];
                break;
            }
        }
        chebyshevfunc(c, cy, n);
    }
    if(view_id==1 || view_id==3 || view_id==4)
       hermitefunc(n, x, y, dy, h);
    //splinefunc(n, x, y, sp, dy);
    if(view_id==2 || view_id==3 || view_id==4)
        splinefunc(n, y, sp, dy, a, b);
    if((view_id==0 || view_id==3) && n<=50)
        chebyshevextrema(c, cx, n, a, b, extr);
    if(view_id==1 || view_id==3)
        hermiteextrema(h, x, n, extr);
    if(view_id==2 || view_id==3)
        splineextrema(sp, x, n, extr);
    if(view_id==4){
        extr[0]=0;
        extr[1]=0;
        if(n<=50)
            chebysheverrorextrema(c, cx, n, a, b, extr, f);
        hermiteerrorextrema(h, x, y, n, extr, f);
        splineerrorextrema(sp, x, n, extr, f);
    }
    delta_y=0.01*(extr[1]-extr[0]);
    extr[0]-=delta_y;
    extr[1]+=delta_y;
    print_console(delta_y);
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
    if(view_id==1 || view_id==3){
        painter.setPen(pen_green);
        for(int i=0;i<n-1;++i){
            x1=x[i];
            x2=x[i+1];
            y1=hermitevaluen(x1, h, i);
            printf("i= %d , y[i] = %lf , hermite = %lf\n", i, y[i], y1);
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
    if((view_id==0 || view_id==3) && n<=50){
        painter.setPen(pen_cyan);
        x1=a;
        y1=chebyshevvalue(x1, a, b, c, n);
        for(int i=0; i<=n; ++i){
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
    if(view_id==4){
        painter.setPen(pen_green);
        for(int i=0;i<n-1;++i){
            x1=x[i];
            x2=x[i+1];
            y1=y[i]-hermitevaluen(x1, h, i);
            for(x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
                y3=f(x3)-hermitevaluen(x3, h, i);
                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=y[i+1]-hermitevaluen(x2, h, i);
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
            y1=f(x1)-splinevaluen(x1, sp, i);
            for(x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
                y3=f(x3)-splinevaluen(x3, sp, i);
                if(fabs(x3-y[n/2])<1e-6)
                    y3+=(p*0.1*absmax);
                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=f(x2)-splinevaluen(x2, sp, i);
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
        }
    }
    if(view_id==4 && n<=50){
        painter.setPen(pen_cyan);
        x1=a;
        y1=f(x1)-chebyshevvalue(x1, a, b, c, n);
        for(int i=0; i<=n; ++i){
            if(i<n)
                x2=cx[i];
            else
                x2=b;
            for(x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
                y3=f(x3)-chebyshevvalue(x3, a, b, c, n);
                if(fabs(x3-y[n/2])<1e-6)
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
    }
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
        painter.drawText(10, 90, QString("absmax(fact): %1( %2)").arg(absmax).arg(max(fabs(extr[0]+delta_y), fabs(extr[1]-delta_y))));
    destroy();
}
