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
       || b-a<1e-6
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

void Window::funcextrema(){
    switch(k){
        case 0:
            min_y=1.0;
            max_y=1.0;
            absmax=1.0;
            break;
        case 1:
            min_y=a;
            max_y=b;
            absmax=max(fabs(min_y), fabs(max_y));
            break;
        case 2:
            min_y=min(a*a, b*b);
            max_y=max(b*b, a*a);
            if(a*b<0){
                min_y=0.0;
            }
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
            if(a*b<0){
                min_y=0.0;
            }
            absmax=max_y;
            break;
        case 5:
            min_y=exp(a);
            max_y=exp(b);
            absmax=max_y;
            break;
        case 6:
            d2f=d2f6;
            if(a>=0 && a<=0){
                min_y=1.0/(25.0*b*b+1.0);
                max_y=1.0;
                absmax=1;
            }else if(b>=0 && b<=0){
                min_y=1.0/(25.0*a*a+1.0);
                max_y=1.0;
                absmax=1.0;
            }else if(a*b>0){
                min_y=min(1.0/(25.0*a*a+1.0), 1.0/(25.0*b*b+1.0));
                max_y=max(1.0/(25.0*a*a+1.0), 1.0/(25.0*b*b+1.0));
                absmax=max_y;
            }else if(a*b<0){
                min_y=min(1.0/(25.0*a*a+1.0), 1.0/(25.0*b*b+1.0));
                max_y=1.0;
                absmax=1.0;
            }
            break;
    }
}

void Window::chebyshevextrema(){
    extrc[0]=10000000;
    extrc[1]=-10000000;
    double y1=chebyshevvalue(a, a, b, c, n);
    if(y1<extrc[0])
        extrc[0]=y1;
    if(y1>extrc[1])
        extrc[1]=y1;
    for(double x1=a+delta_x; x1-b<delta_x; x1+=delta_x){
        y1=chebyshevvalue(x1, a, b, c, n);
        if(y1<extrc[0])
            extrc[0]=y1;
        if(y1>extrc[1])
            extrc[1]=y1;
    }
    y1=chebyshevvalue(b, a, b, c, n);
    if(y1<extrc[0])
        extrc[0]=y1;
    if(y1>extrc[1])
        extrc[1]=y1;
}

void Window::chebysheverrorextrema(){
    extrc[0]=10000000;
    extrc[1]=-10000000;
    double y1=y[0]-chebyshevvalue(a, a, b, c, n);
    if(y1<extrc[0])
        extrc[0]=y1;
    if(y1>extrc[1])
        extrc[1]=y1;
    for(double x1=a+delta_x; x1-b<delta_x; x1+=delta_x){
        y1=f(x1)-chebyshevvalue(x1, a, b, c, n);
        if(fabs(x1-x[n/2])<delta_x && p)
            y1+=(0.1*absmax*p);
        if(y1<extrc[0])
            extrc[0]=y1;
        if(y1>extrc[1])
            extrc[1]=y1;
    }
    y1=y[n-1]-chebyshevvalue(b, a, b, c, n);
    if(y1<extrc[0])
        extrc[0]=y1;
    if(y1>extrc[1])
        extrc[1]=y1;
}

void Window::hermiteextrema(){
    extrh[0]=10000000;
    extrh[1]=-10000000;
    double x1=a+delta_x,y1,x2;
    y1=hermitevaluen(a,h,0);
    if(y1<extrh[0])
        extrh[0]=y1;
    if(y1>extrh[1])
        extrh[1]=y1;
    for(int i=0;i<n-1;++i){
        x2=x[i+1];
        while(x1-x2<delta_x){
            y1=hermitevaluen(x1,h,i);
            if(y1<extrh[0])
                extrh[0]=y1;
            if(y1>extrh[1])
                extrh[1]=y1;
            x1+=delta_x;
        }
    }
    y1=hermitevaluen(b,h,n-2);
    if(y1<extrh[0])
        extrh[0]=y1;
    if(y1>extr[1])
        extrh[1]=y1;
}

void Window::hermiteerrorextrema(){
    extrh[0]=10000000;
    extrh[1]=-10000000;
    double x1=a+delta_x,y1,x2;
    y1=y[0]-hermitevaluen(a,h,0);
    if(y1<extrh[0])
        extrh[0]=y1;
    if(y1>extrh[1])
        extrh[1]=y1;
    for(int i=0;i<n-1;++i){
        x2=x[i+1];
        while(x1-x2<delta_x){
            y1=f(x1)-hermitevaluen(x1,h,i);
            if(fabs(x1-x[n/2])<delta_x && p)
                y1+=(p*0.1*absmax);
            if(y1<extrh[0])
                extrh[0]=y1;
            if(y1>extrh[1])
                extrh[1]=y1;
            x1+=delta_x;
        }
    }
    y1=y[n-1]-hermitevaluen(b,h,n-2);
    if(y1<extrh[0])
        extrh[0]=y1;
    if(y1>extrh[1])
        extrh[1]=y1;
}

void Window::splineextrema(){
    extrs[0]=10000000;
    extrs[1]=-10000000;
    double x1=a+delta_x,y1,x2;
    y1=splinevaluen(a, sp, 0.0);
    if(y1<extrs[0])
        extrs[0]=y1;
    if(y1>extrs[1])
        extrs[1]=y1;
    for(int i=0;i<n;++i){
        if(i!=n-1)
            x2=(x[i]+x[i+1])/2.0;
        else
            x2=x[i];
        while(x1-x2<delta_x){
            y1=splinevaluen(x1,sp,i);
            if(y1<extrs[0])
                extrs[0]=y1;
            if(y1>extrs[1])
                extrs[1]=y1;
            x1+=delta_x;
        }
    }
    y1=splinevaluen(b, sp, n-1);
    if(y1<extrs[0])
        extrs[0]=y1;
    if(y1>extrs[1])
        extrs[1]=y1;
}

void Window::splineerrorextrema(){
    extrs[0]=10000000;
    extrs[1]=-10000000;
    double x1=a+delta_x,y1,x2;
    y1=y[0]-splinevaluen(a, sp, 0.0);
    if(y1<extrs[0])
        extrs[0]=y1;
    if(y1>extrs[1])
        extrs[1]=y1;
    for(int i=0;i<n;++i){
        if(i!=n-1)
            x2=(x[i]+x[i+1])/2.0;
        else
            x2=x[i];
        while(x1-x2<delta_x){
            y1=f(x1)-splinevaluen(x1,sp,i);
            if(fabs(x1-x[n/2])<delta_x && p)
                y1+=(p*0.1*absmax);
            if(y1<extrs[0])
                extrs[0]=y1;
            if(y1>extrs[1])
                extrs[1]=y1;
            x1+=delta_x;
        }
    }
    y1=y[n-1]-splinevaluen(b, sp, n-1);
    if(y1<extrs[0])
        extrs[0]=y1;
    if(y1>extrs[1])
        extrs[1]=y1;
}

void Window::extrema_hunt(){
    /*if(view_id!=4){
       if(p>0){
            extr[1]=y[n/2];
            if(extr[1]<max_y)
                extr[1]=max_y;
            extr[0]=min_y;
        }else if(p<0){
            extr[0]=y[n/2];
            if(extr[0]>min_y)
                extr[0]=min_y;
            extr[1]=max_y;
        }else{
            extr[0]=min_y;
            extr[1]=max_y;
       }
    }
    if(view_id==4){
        extr[0]=0;
        extr[1]=0;
    }*/
    extr[0]=1000000;
    extr[1]=-1000000;
    if(view_id==2 || view_id==3){
        splineextrema();
        if(extrs[0]<extr[0])
            extr[0]=extrs[0];
        if(extrs[1]>extr[1])
            extr[1]=extrs[1];
        //printf("S: %10.3e;%10.3e\n",extr[0],extr[1]);
    }
    if(view_id==1 || view_id==3){
        hermiteextrema();
        if(extrh[0]<extr[0])
            extr[0]=extrh[0];
        if(extrh[1]>extr[1])
            extr[1]=extrh[1];
        //printf("H: %10.3e;%10.3e\n",extr[0],extr[1]);
    }
    if((view_id==0 || view_id==3) && n<=50){
        chebyshevextrema();
        if(extrc[0]<extr[0])
            extr[0]=extrc[0];
        if(extrc[1]>extr[1])
            extr[1]=extrc[1];
        //printf("C: %10.3e;%10.3e\n",extr[0],extr[1]);
    }
    if(view_id==4){
        splineerrorextrema();
        if(extrs[0]<extr[0])
            extr[0]=extrs[0];
        if(extrs[1]>extr[1])
            extr[1]=extrs[1];
        //printf("SE: %10.3e;%10.3e\n",extr[0],extr[1]);
        hermiteerrorextrema();
        if(extrh[0]<extr[0])
            extr[0]=extrh[0];
        if(extrh[1]>extr[1])
            extr[1]=extrh[1];
        //printf("HE: %10.3e;%10.3e\n",extr[0],extr[1]);
        if(n<=50){
            chebysheverrorextrema();
            if(extrc[0]<extr[0])
                extr[0]=extrc[0];
            if(extrc[1]>extr[1])
                extr[1]=extrc[1];
            //printf("CE: %10.3e;%10.3e\n",extr[0],extr[1]);
        }
    }
    if(view_id!=4){
        if(extr[1]<max_y)
            extr[1]=max_y;
        if(extr[0]>min_y)
            extr[0]=min_y;
        if(p>0 && y[n/2]>extr[1]){
            extr[1]=y[n/2];
        }
        if(p<0 && y[n/2]<extr[0]){
            extr[0]=y[n/2];
        }
    }
}

void Window::print_console(){
    printf("format: %d\n", view_id);
    printf("segment: [%lf;%lf]\n", a, b);
    printf("squeeze-stretch: %d\n", s);
    printf("points: %d\n", n);
    printf("disturbance: %d\n", p);
    printf("function absmax: %10.3e\n", absmax);
    /*if(k==0 && p==0 && view_id!=4){
        printf("factic absmax: %10.3e\n", absmax);
        printf("extrema: %10.3e %10.3e\n\n", absmax, extr[0],extr[1]);
    }else{*/
        printf("factic absmax: %10.3e\n", max(fabs(extr[0]), fabs(extr[1])));
        printf("extrema: %10.3e %10.3e\n\n", extr[0], extr[1]);
        if(view_id==4 || view_id==3 || view_id==1)
            printf("Hermite: %10.3e %10.3e\n\n", extrh[0], extrh[1]);
        if(view_id==4 || view_id==3 || view_id==2)
            printf("Splines: %10.3e %10.3e\n\n", extrs[0], extrs[1]);
        if(n<=50 && (view_id==4 || view_id==3 || view_id==0))
            printf("Chebyshev: %10.3e %10.3e\n\n", extrc[0], extrc[1]);
    //}
}

void Window::change_func(){
    k=(k+1)%7;
    switch(k){
        case 0:
            f_name="k=0 f(x)=1";
            f=f0;
            df=df0;
            d2f=d2f0;
            break;
        case 1:
            f_name="k=1 f(x)=x";
            f=f1;
            df=df1;
            d2f=d2f1;
            break;
        case 2:
            f_name="k=2 f(x)=x^2";
            f=f2;
            df=df2;
            d2f=d2f2;
            break;
        case 3:
            f_name="k=3 f(x)=x^3";
            f=f3;
            df=df3;
            d2f=d2f3;
            break;
        case 4:
            f_name="k=4 f(x)=x^4";
            f=f4;
            df=df4;
            d2f=d2f4;
            break;
        case 5:
            f_name="k=5 f(x)=e^x";
            f=f5;
            df=df5;
            d2f=d2f5;
            break;
        case 6:
            f_name="k=6 f(x)=1/(25x^2+1)";
            f=f6;
            df=df6;
            d2f=d2f6;
            break;
    }
    funcextrema();
    fillvalues(x, y, dy, n, f, df);
    if(p){
        y[n/2]+=(p*0.1*absmax);
    }
    if(n<=50){
        chebyshevf(cx, cy, n, f);
        for(int i=0; i<n; ++i){
            if(fabs(cx[i]-x[n/2])<delta_x){
                cy[i]+=(p*0.1*absmax);
                break;
            }
        }
        chebyshevfunc(c, cx,cy, n,a,b);
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
    funcextrema();
    extrema_hunt();
    print_console();
    update();
}

void Window::enhance(){
    double t=(b-a)/2.0;
    a=a-t;
    b=b+t;
    ++s;
    delta_x=(b-a)/width();
    funcextrema();
    fillpoints(x, n, a, b);
    fillvalues(x, y, dy, n, f, df);
    if(p)
        y[n/2]+=(p*0.1*absmax);
    if(n<=50){
        chebyshevpoints(cx, n, a, b);
        chebyshevf(cx,cy,n,f);
        for(int i=0; i<n; ++i){
            if(fabs(cx[i]-x[n/2])<delta_x){
                cy[i]+=(p*0.1*absmax);
                break;
            }
        }
        chebyshevfunc(c, cx,cy, n,a,b);
    }
    hermitefunc(n, x, y, dy, h);
    //splinefunc(n, x, y, sp, dy);
    splinefunc(n, y, sp, dy, a, b);
    extrema_hunt();
    print_console();
    update();
}

void Window::squeeze(){
    double t=(b-a)/4.0;
    a=a+t;
    b=b-t;
    --s;
    delta_x=(b-a)/width();
    funcextrema();
    fillpoints(x, n, a, b);
    fillvalues(x, y, dy, n, f, df);
    if(p)
        y[n/2]+=(p*0.1*absmax);
    if(n<=50){
        chebyshevpoints(cx, n, a, b);
        chebyshevf(cx,cy,n,f);
        for(int i=0; i<n; ++i){
            if(fabs(cx[i]-x[n/2])<delta_x){
                cy[i]+=(p*0.1*absmax);
                break;
            }
        }
        chebyshevfunc(c, cx,cy, n,a,b);
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
    funcextrema();
    allocate();
    fillpoints(x, n, a, b);
    fillvalues(x, y, dy, n, f, df);
    if(p)
        y[n/2]+=(p*0.1*absmax);
    if(n<=50){
        chebyshevpoints(cx, n, a, b);
        chebyshevf(cx,cy,n,f);
        for(int i=0; i<n; ++i){
            if(fabs(cx[i]-x[n/2])<delta_x){
                cy[i]+=(p*0.1*absmax);
                break;
            }
        }
        chebyshevfunc(c, cx,cy, n,a,b);
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
    funcextrema();
    allocate();
    fillpoints(x, n, a, b);
    fillvalues(x, y, dy, n, f, df);
    if(p)
        y[n/2]+=(p*0.1*absmax);
    if(n<=50){
        chebyshevpoints(cx, n, a, b);
        chebyshevf(cx,cy,n,f);
        for(int i=0; i<n; ++i){
            if(fabs(cx[i]-x[n/2])<delta_x){
                cy[i]+=(p*0.1*absmax);
                break;
            }
        }
        chebyshevfunc(c, cx,cy, n,a,b);
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
    funcextrema();
    y[n/2]+=(0.1*absmax);
    if(n<=50){
        for(int i=0; i<n; ++i){
            if(fabs(cx[i]-x[n/2])<delta_x){
                cy[i]+=(0.1*absmax);
                break;
            }
        }
        chebyshevfunc(c, cx,cy, n,a,b);
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
    funcextrema();
    y[n/2]-=(0.1*absmax);
    if(n<=50){
        for(int i=0; i<n; ++i){
            if(fabs(cx[i]-x[n/2])<delta_x){
                cy[i]-=(0.1*absmax);
                break;
            }
        }
        chebyshevfunc(c, cx,cy, n,a,b);
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
    double x1, x2, y1, x3, y3;
    delta_x=(b-a)/width();
    QPen pen_black(Qt::black, 0, Qt::SolidLine);
    QPen pen_red(Qt::red, 0, Qt::SolidLine);
    QPen pen_green(Qt::green, 0, Qt::SolidLine);
    QPen pen_blue(Qt::blue, 0, Qt::SolidLine);
    QPen pen_cyan(Qt::cyan, 0, Qt::SolidLine);
    painter.save();
    // make Coordinate Transformations
    painter.translate(0.5*width(), 0.5*height());
    double delta_y;
    //=(extr[1]-extr[0])/height();
    if(k==0 && p==0 && view_id!=4){
        delta_y=0.05;
        painter.scale(width()/(b-a+2.0*delta_x), -height()/(2.0+2.0*delta_y));
        painter.translate(-0.5*(a+b), -1.0);
    }else{
        delta_y=(extr[1]-extr[0])/100.0;
        painter.scale(width()/(b-a+2.0*delta_x), -height()/(extr[1]-extr[0]+2*delta_y));
        painter.translate(-0.5*(a+b), -0.5*(extr[1]+extr[0]));
    }
    painter.setPen(pen_red);
    painter.drawLine(QPointF(a-delta_x, 0), QPointF(b+delta_x, 0));
    if(k==0 && p==0 && view_id!=4)
        painter.drawLine(QPointF(0, -delta_y), QPointF(0, 2.0+delta_y));
    else
        painter.drawLine(QPointF(0, extr[0]-delta_y), QPointF(0, extr[1]+delta_y));
    if(view_id!=4){
        painter.setPen(pen_black);
        /*for(int i=0;i<n-1;++i){
            x1=x[i];
            x2=x[i+1];
            y1=y[i];
            for(x3=x1+delta_x; x3-x2<delta_x; x3+=delta_x){
                y3=f(x3);
                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=y[i+1];
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
            if(fabs(x3-x[n/2])<delta_x && p)
                y3+=(0.1*absmax*p);
        }*/
        x1=a;
        y1=y[0];
        for(x3=a+delta_x; x3-b<delta_x; x3+=delta_x){
            y3=f(x3);
            if(fabs(x3-x[n/2])<delta_x && p)
                y3+=(0.1*absmax*p);
            painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
            x1=x3;
            y1=y3;
        }
        y3=y[n-1];
        painter.drawLine(QPointF(x1, y1), QPointF(b, y3));
    }
    if((view_id==0 || view_id==3) && n<=50){
        painter.setPen(pen_cyan);
        x1=a;
        y1=chebyshevvalue(a, a, b, c, n);
        for(x3=a+delta_x; x3-b<delta_x; x3+=delta_x){
            y3=chebyshevvalue(x3, a, b, c, n);
            /*if(y3>extr[1])
                printf("ExceedC:%10.3e,%10.3e\n",x3,y3);
            if(y3<extr[0])
                printf("ExceedC:%10.3e,%10.3e\n",x3,y3);*/
            painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
            x1=x3;
            y1=y3;
        }
        y3=chebyshevvalue(b, a, b, c, n);
        painter.drawLine(QPointF(x1, y1), QPointF(b, y3));
    }
    if(view_id==1 || view_id==3){
        painter.setPen(pen_green);
        x1=a;
        x3=a+delta_x;
        y1=hermitevaluen(a,h,0);
        for(int i=0;i<n-1;++i){
            x2=x[i+1];
            while(x3-x2<delta_x){
                y3=hermitevaluen(x3,h,i);
                /*if(y3>extr[1])
                    printf("ExceedH %d:%10.3e,%10.3e\n",i,x3,y3);
                if(y3<extr[0])
                    printf("ExceedH %d:%10.3e,%10.3e\n",i,x3,y3);*/
                painter.drawLine(QPointF(x1,y1),QPointF(x3,y3));
                x1=x3;
                y1=y3;
                x3+=delta_x;
            }
        }
        y3=hermitevaluen(b,h,n-2);
        painter.drawLine(QPointF(x1,y1),QPointF(b,y3));
    }
    if(view_id==2 || view_id==3){
        painter.setPen(pen_blue);
        x1=a;
        x3=a+delta_x;
        y1=splinevaluen(a, sp, 0);
        for(int i=0;i<n;++i){
            if(i!=n-1)
                x2=(x[i]+x[i+1])/2.0;
            else
                x2=x[i];
            while(x3-x2<delta_x){
                y3=splinevaluen(x3,sp,i);
                /*if(y3>extr[1])
                    printf("ExceedS %d:%10.3e,%10.3e\n",i,x3,y3);
                if(y3<extr[0])
                    printf("ExceedS %d:%10.3e,%10.3e\n",i,x3,y3);*/
                painter.drawLine(QPointF(x1,y1),QPointF(x3,y3));
                x1=x3;
                y1=y3;
                x3+=delta_x;
            }
        }
        y3=splinevaluen(b, sp, n-1);
        painter.drawLine(QPointF(x1,y1),QPointF(b,y3));
    }
    if(view_id==4 && n<=50){
        painter.setPen(pen_cyan);
        x1=a;
        y1=y[0]-chebyshevvalue(a, a, b, c, n);
        for(x3=a+delta_x; x3-b<delta_x; x3+=delta_x){
            y3=f(x3)-chebyshevvalue(x3, a, b, c, n);
            if(fabs(x3-x[n/2])<delta_x && p)
                y3+=(0.1*absmax*p);
            /*if(y3>extr[1])
                printf("ExceedCE:%10.3e,%10.3e\n",x3,y3);
            if(y3<extr[0])
                printf("ExceedCE:%10.3e,%10.3e\n",x3,y3);*/
            painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
            x1=x3;
            y1=y3;
        }
        y3=y[n-1]-chebyshevvalue(b, a, b, c, n);
        painter.drawLine(QPointF(x1, y1), QPointF(b, y3));
    }
    if(view_id==4){
        painter.setPen(pen_green);
        x1=a;
        x3=a+delta_x;
        y1=y[0]-hermitevaluen(a,h,0);
        for(int i=0;i<n-1;++i){
            x2=x[i+1];
            while(x3-x2<delta_x){
                y3=f(x3)-hermitevaluen(x3,h,i);
                if(fabs(x3-x[n/2])<delta_x && p)
                    y3+=(0.1*p*absmax);
                /*if(y3>extr[1])
                    printf("ExceedHE %d:%10.3e,%10.3e\n",i,x3,y3);
                if(y3<extr[0])
                    printf("ExceedHE %d:%10.3e,%10.3e\n",i,x3,y3);*/
                painter.drawLine(QPointF(x1,y1),QPointF(x3,y3));
                x1=x3;
                y1=y3;
                x3+=delta_x;
            }
        }
        y3=y[n-1]-hermitevaluen(b,h,n-2);
        painter.drawLine(QPointF(x1,y1),QPointF(b,y3));
        painter.setPen(pen_blue);
        x1=a;
        x3=a+delta_x;
        y1=y[0]-splinevaluen(a, sp, 0);
        for(int i=0;i<n;++i){
            if(i!=n-1)
                x2=(x[i]+x[i+1])/2.0;
            else
                x2=x[i];
            while(x3-x2<delta_x){
                y3=f(x3)-splinevaluen(x3,sp,i);
                if(fabs(x3-x[n/2])<delta_x && p)
                    y3+=(p*0.1*absmax);
                /*if(y3>extr[1])
                    printf("ExceedSE %d:%10.3e,%10.3e\n",i,x3,y3);
                if(y3<extr[0])
                    printf("ExceedSE %d:%10.3e,%10.3e\n",i,x3,y3);*/
                painter.drawLine(QPointF(x1,y1),QPointF(x3,y3));
                x1=x3;
                y1=y3;
                x3+=delta_x;
            }
        }
        y3=y[n-1]-splinevaluen(b, sp, n-1);
        painter.drawLine(QPointF(x1,y1),QPointF(b,y3));
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
        painter.drawText(10, 90, QString("absmax(fact): %1( %2)").arg(absmax).arg(max(fabs(extr[0]/*+delta_y*/), fabs(extr[1]/*-delta_y*/))));
}
