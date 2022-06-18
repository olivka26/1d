#include <QPainter>
#include <stdio.h>
#include <stdlib.h>
#include "window.hpp"
#include "hermite.hpp"
#include "splines.hpp"
#include "math.h"

#define DEFAULT_A -10
#define DEFAULT_B 10
#define DEFAULT_N 30
#define DEFAULT_K 0

using namespace std;

double chebyshevtrans(double a, double b, double x){
    return (2*x-(b+a))/(b-a);
}

double chebyshevfunc(double x, double a, double b, double *c, int n){
    double t0=1,cheb=chebyshevtrans(a,b,x),t1=cheb,t2=2*t1*t1-t0, res=0;
    for(int i=0;i<n-2;i+=3){
        res+=c[i]*t0;
        res+=c[i+1]*t1;
        res+=c[i+2]*t2;
        t0=2*t2*cheb-t1;
        t1=2*t0*cheb-t2;
        t2=2*t1*cheb-t0;
    }
    return res;
}

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
    change_func();
    return 0;
}

void Window::change_func(){
    k=(k+1)%7;
    //double dif=y[n/2];
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
    QPen pen_black(Qt::black, 0, Qt::SolidLine);
    QPen pen_red(Qt::red, 0, Qt::SolidLine);
    QPen pen_green(Qt::green, 0, Qt::SolidLine);
    QPen pen_blue(Qt::blue, 0, Qt::SolidLine);
    QPen pen_cyan(Qt::cyan, 0, Qt::SolidLine);
    x=(double*)malloc(n*sizeof(double));
    y=(double*)malloc(n*sizeof(double));
    h=(double*)malloc(4*n*sizeof(double));
    sp=(double*)malloc(3*n*sizeof(double));
    dy=(double*)malloc(n*sizeof(double));
    c=(double*)malloc(n*sizeof(double));
    printf("%d\n", view_id);
    printf("[%lf, %lf]\n",a,b);
    x[0]=a;
    y[0]=f(a);
    dy[0]=df(a);
    for(int i=1; i<n; ++i){
        x[i]=x[i-1]+(b-a)/(n-1);
        y[i]=f(x[i]);
        if(i==n/2){
            y[i]+=p*0.1*max_y;
            if(p>0 && y[i]>max_y){
                max_y=y[i];
            }
            if(p<0 && y[i]<min_y){
                min_y=y[i];
            }
            absmax=max(fabs(min_y), fabs(max_y));
        }
        dy[i]=df(x[i]);
        c[i]=0;
    }
    
    delta_y=0.01*(max_y-min_y);
    min_y-=delta_y;
    max_y+=delta_y;
    painter.save();
    // make Coordinate Transformations
    painter.translate(0.5*width(), 0.5*height());
    painter.scale(width()/(b-a), -height()/(max_y-min_y));
    painter.translate(-0.5*(a+b), -0.5*(min_y+max_y));
   // hermitefunc(n, x, y, dy, h);
    for(int i=0; i<n-1; ++i){
        double dify=y[i+1]-y[i], divdif=dify/delta_x;
        h[4*i]=y[i];
        h[4*i+1]=dy[i];
        h[4*i+2]=(3*divdif-2*dy[i]-dy[i+1])/delta_x;
        h[4*i+3]=(dy[i]+dy[i+1]-2*divdif)/(delta_x*delta_x);
    }
    splinefunc(n, x, y, sp, dy, a, b);
    for(int i=0;i<n;++i){
        double xx=cos(3.1415926535897932384626433832795*(2*i+1)/(2*n));
        double res1=f((a+b)/2+(b-a)*xx/2);
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
    }
    for(int i=0;i<n;++i){
        c[i]/=n;
        if(i)
            c[i]*=2;
    }
    if(view_id!=4){
        painter.setPen(pen_black);
        x1=a;
        y1=f(x1);
        for(x2=x1+delta_x; x2-b<1e-6; x2+=delta_x){
            for(x3=x1+delta_x*1e-4; x3-x2<1e-6; x3+=delta_x*0.0001){
                y3=f(x3);
                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=f(x2);
            if(x2==a+delta_x*n/2){
                y2+=p*0.1*absmax;
            }
            x1=x2;
            y1=y2;
        }
        y2=f(b);
        //printf("Graph: %lf %lf %lf %lf\n",x1, y1, b, y2);
        painter.drawLine(QPointF(x1, y1), QPointF(b, y2));
    }
    if(view_id==0 || view_id==3){
        painter.setPen(pen_green);
        x1=a;
        y1=hermitevalue(x1, a, b, n, x, h);
        for(x2=x1+delta_x; x2-b<1e-6; x2+=delta_x){
            for(x3=x1+delta_x*1e-4; x3-x2<1e-6; x3+=delta_x*0.0001){
                y3=hermitevalue(x3, a, b, n, x, h);
                //printf("ZW: %lf %lf\n",y3, f(x3));
                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=hermitevalue(x2, a, b, n, x, h);
            //printf("ZW: %lf %lf\n",y2, f(x2));
            x1=x2;
            y1=y2;
        }
        y2=hermitevalue(b, a, b, n, x, h);
        painter.drawLine(QPointF(x1, y1), QPointF(b, y2));
        //printf("Hermite: %lf %lf %lf %lf\n",x1, y1, b, y2);
    }
    if(view_id==1 || view_id==3){
        painter.setPen(pen_blue);
        x1=a;
        y1=splinevalue(a, a, b, n, sp, x);
        if(y1<min_y)
            min_y=y1;
        if(y1>max_y)
            max_y=y1;
        for(x2=x1+delta_x/2; x2-b<1e-6; x2+=delta_x){
            for(x3=x1+delta_x*1e-4; x3-x2<1e-6; x3+=delta_x*0.0001){
                y3=splinevalue(x3, a, b, n, sp, x);
                if(y3<min_y)
            min_y=y3;
        if(y3>max_y)
            max_y=y3;
                //printf("ZW: %lf %lf\n",y3, f(x3));
                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=splinevalue(x2, a, b, n, sp, x);
            if(y2<min_y)
            min_y=y2;
        if(y2>max_y)
            max_y=y2;
            //printf("ZW: %lf %lf\n",y2, f(x2));
            x1=x2;
            y1=y2;
        }
        y2=splinevalue(b, a, b, n, sp, x);
        painter.drawLine(QPointF(x1, y1), QPointF(b, y2));
    }
    if((view_id==2 || view_id==3) && n<=50){
        painter.setPen(pen_cyan);
        x1=a;
        y1=chebyshevfunc(x1, a, b, c, n);
        if(y1<min_y)
            min_y=y1;
        if(y1>max_y)
            max_y=y1;
        x2=x1+0.0001*delta_x;
        y2=chebyshevfunc(x2, a, b, c, n);
        if(y2<min_y)
            min_y=y2; 
        if(y2>max_y)
            max_y=y2;
        while(x2<b-1e-6){
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
            //printf("Chebyshev: %lf, %lf, %lf, %lf\n", x1, y1, x2, y2);
                //painter.drawLine(x1, y1, x3, y3);
            x1=x2;
            y1=y2;
            x2+=0.0001*delta_x;
            y2=chebyshevfunc(x2, a, b, c, n);
            if(y2<min_y)
                min_y=y2; 
            if(y2>max_y)
                max_y=y2;
        }
        x2=b;
        y2=chebyshevfunc(x2, a, b, c, n);
        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
        //painter.drawLine(x1, y1, x2, y2);
       // update();
    }

    if(view_id==4){
        painter.setPen(pen_green);
        x1=a;
        y1=f(x1)-hermitevalue(x1, a, b, n, x, h);
        for(x2=x1+delta_x; x2<b-1e-6; x2+=delta_x){
            y2=f(x2)-hermitevalue(x2, a, b, n, x, h);
/*
            if(i==n/2)
                y2+=delta*0.1*(abs(max_y)>abs(min_y) ? abs(max_y) : abs(min_y));
*/            
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
           // painter.drawLine(x1, y1, x2, y2);
            x1=x2;
            y1=y2;
        }
        x2=b;
        y2=f(x2)-hermitevalue(x2, a, b, n, x, h);
        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
        //painter.drawLine(x1, y1, x2, y2);
        painter.setPen(pen_blue);
        x1=a;
        y1=f(x1)-splinevalue(x1, a, b, n, sp, x);
        for(int i=0, x2=x1+delta_x; x2<b-1e-6; x2+=delta_x, ++i){
            if(i==0 || i==n-2)
                x2-=(delta_x/2);
            y2=f(x2)-splinevalue(x2, a, b, n, sp, x);
  /*          if(i==n/2)
                y2+=delta*0.1*(abs(max_y)>abs(min_y) ? abs(max_y) : abs(min_y));*/
            painter.drawLine(QPointF(x1, y1), QPointF (x2, y2));
            //painter.drawLine(x1, y1, x2, y2);
            x1=x2;
            y1=y2;
        }
        x2=b;
        y2=f(x2)-hermitevalue(x2, a, b, n, x, h);
        painter.drawLine (QPointF(x1, y1), QPointF(x2, y2));
        //painter.drawLine(x1, y1, x2, y2);
       // update();
        if(n<=50){
        painter.setPen(pen_cyan);
        x1=a;
        y1=f(x1)-chebyshevfunc(x1, a, b, c, n);
        x2=x1+0.0001*delta_x;
        y2=f(x2)-chebyshevfunc(x2, a, b, c, n);
        while(x2<b-1e-6){
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
            //printf("Chebyshev: %lf, %lf, %lf, %lf\n", x1, y1, x2, y2);
                //painter.drawLine(x1, y1, x3, y3);
            x1=x2;
            y1=y2;
            x2+=0.0001*delta_x;
            y2=f(x2)-chebyshevfunc(x2, a, b, c, n);
        }
        x2=b;
        y2=f(x2)-chebyshevfunc(x2, a, b, c, n); 
        painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
        }
    }
    painter.setPen(pen_red);
    painter.drawLine(a, 0, b, 0);
    painter.drawLine(0, max_y, 0, min_y);
    delta_y=0.01*(max_y-min_y);
    painter.restore();
    painter.setPen("black");
    painter.drawText(0, 20, f_name);
    painter.drawText(10, 30, QString("Format = %1").arg(view_id + 1));
    painter.drawText(10, 45, QString("Scale = %1 %2 %3").arg(s).arg(a).arg(b));
    painter.drawText(10, 60, QString("n = %1").arg(n));
    painter.drawText(10, 75, QString("p = %1").arg(p));
    printf("points: %d\n", n);
    printf("max_y: %lf\n", max_y);
    free(x);
    free(y);
    free(c);
    free(dy);
    free(sp);
    free(h);
}

/*    switch(k){
        case 0:
             min_y=1;
             max_y=1;
             absmax=1;
             break;
        case 1:
             min_y=a;
             max_y=b;
             absmax=max(fabs(min_y), fabs(max_y));
             break;
        case 2:
             min_y=min(a*a, b*b);
             max_y=max(b*b, a*a);
             absmax=max_y;
             break;
        case 3:
             min_y=a*a*a;
             max_y=b*b*b;
             absmax=max(fabs(min_y), fabs(max_y));
             break;
        case 4:
             min_y=min(a*a*a*a, b*b*b*b);
             max_y=max(b*b*b*b, a*a*a*a);
             absmax=max_y;
             break;
        case 5:
             min_y=exp(a);
             max_y=exp(b);
             absmax=max_y;
             break;
        case 6:
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
    }*/
