
#include <QDebug>
#include <QString>
#include "calc.h"
#include "chebyshev.hpp"
#include "hermite.hpp"
#include "splines.hpp"
#include "math.h"
#include "general.hpp"

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

int calcAndDraw::Draw(mglGraph *gr)
{
          gr->Title(f_name);
          gr->SetOrigin(0,0);
          gr->SetRanges(a,b,f(a),f(b));
          gr->Axis();
          gr->Label('x', "x", 1);
          gr->Label('y', "f(x)", 1);
          gr->Grid();

//        gr->SubPlot(2,2,1);	gr->Title("2 axis");
//        gr->SetRanges(-1,1,-1,1);	gr->SetOrigin(-1,-1,-1);	// first axis
//        gr->Axis();	gr->Label('y',"axis 1",0);	gr->FPlot("sin(pi*x)","r2");
//        gr->SetRanges(0,1,0,1);		gr->SetOrigin(1,1,1);		// second axis
//        gr->Axis();	gr->Label('y',"axis 2",0);	gr->FPlot("cos(pi*x)");

//        gr->SubPlot(2,2,3);	gr->Title("More axis");	gr->SetOrigin(NAN,NAN);	gr->SetRange('x',-1,1);
//        gr->Axis();	gr->Label('x',"x",0);	gr->Label('y',"y_1",0);	gr->FPlot("x^2","k");
//        gr->SetRanges(-1,1,-1,1);	gr->SetOrigin(-1.3,-1);	// second axis
//        gr->Axis("y","r");	gr->Label('y',"#r{y_2}",0.2);	gr->FPlot("x^3","r");

//        gr->SubPlot(2,2,2);	gr->Title("4 segments, inverted axis");		gr->SetOrigin(0,0);
//        gr->InPlot(0.5,1,0.5,1);	gr->SetRanges(0,10,0,2);	gr->Axis();
//        gr->FPlot("sqrt(x/2)");		gr->Label('x',"W",1);	gr->Label('y',"U",1);
//        gr->InPlot(0,0.5,0.5,1);	gr->SetRanges(1,0,0,2);	gr->Axis("x");
//        gr->FPlot("sqrt(x)+x^3");	gr->Label('x',"\\tau",-1);
//        gr->InPlot(0.5,1,0,0.5);	gr->SetRanges(0,10,4,0);	gr->Axis("y");
//        gr->FPlot("x/4");	gr->Label('y',"L",-1);
//        gr->InPlot(0,0.5,0,0.5);	gr->SetRanges(1,0,4,0);	gr->FPlot("4*x^2");

    double x1, x2, y1, y2, x3, y3;
    double delta_y, delta_x=(b-a)/(n-1);
    delta_x /= 10000.0;

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

    if(view_id!=4){
        for(int i=0;i<n-1;++i){
            x1=x[i];
            x2=x[i+1];
            y1=y[i];
            for(x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
                y3=f(x3);
//                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                gr->Line(mglPoint(x1, y1),mglPoint(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=y[i+1];
//            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
            gr->Line(mglPoint(x1, y1),mglPoint(x2, y2));

        }
    }
    if(view_id==1 || view_id==3){
//        painter.setPen(pen_green);
        for(int i=0;i<n-1;++i){
            x1=x[i];
            x2=x[i+1];
            y1=hermitevaluen(x1, h, i);
            printf("i= %d , y[i] = %lf , hermite = %lf\n", i, y[i], y1);
            for(x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
                y3=hermitevaluen(x3, h, i);
//                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                gr->Line(mglPoint(x1, y1),mglPoint(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=hermitevaluen(x2, h, i);
//            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
            gr->Line(mglPoint(x1, y1),mglPoint(x2, y2));
        }
    }
    if(view_id==2 || view_id==3){
//        painter.setPen(pen_blue);
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
//                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                gr->Line(mglPoint(x1, y1),mglPoint(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=splinevaluen(x2, sp, i);
//            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
            gr->Line(mglPoint(x1, y1),mglPoint(x2, y2));
        }
    }
    if((view_id==0 || view_id==3) && n<=50){
//        painter.setPen(pen_cyan);
        x1=a;
        y1=chebyshevvalue(x1, a, b, c, n);
        for(int i=0; i<=n; ++i){
            if(i<n)
                x2=cx[i];
            else
                x2=b;
            for(x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
                y3=chebyshevvalue(x3, a, b, c, n);
//                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                gr->Line(mglPoint(x1, y1),mglPoint(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=chebyshevvalue(x2, a, b, c, n);
//            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
            gr->Line(mglPoint(x1, y1),mglPoint(x2, y2));
            x1=x2;
            y1=y2;
        }
    }
    if(view_id==4){
//        painter.setPen(pen_green);
        for(int i=0;i<n-1;++i){
            x1=x[i];
            x2=x[i+1];
            y1=y[i]-hermitevaluen(x1, h, i);
            for(x3=x1+delta_x; x3-x2<1e-6; x3+=delta_x){
                y3=f(x3)-hermitevaluen(x3, h, i);
//                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                gr->Line(mglPoint(x1, y1),mglPoint(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=y[i+1]-hermitevaluen(x2, h, i);
//            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
            gr->Line(mglPoint(x1, y1),mglPoint(x2, y2));
        }
//        painter.setPen(pen_blue);
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
//                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                gr->Line(mglPoint(x1, y1),mglPoint(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=f(x2)-splinevaluen(x2, sp, i);
//            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
            gr->Line(mglPoint(x1, y1),mglPoint(x2, y2));
        }
    }
    if(view_id==4 && n<=50){
//        painter.setPen(pen_cyan);
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
//                painter.drawLine(QPointF(x1, y1), QPointF(x3, y3));
                gr->Line(mglPoint(x1, y1),mglPoint(x3, y3));
                x1=x3;
                y1=y3;
            }
            y2=f(x2)-chebyshevvalue(x2, a, b, c, n);
            if(fabs(x3-y[n/2])<1e-6)
                y2+=(p*0.1*absmax);
//            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2));
            gr->Line(mglPoint(x1, y1),mglPoint(x2, y2));
            x1=x2;
            y1=y2;
        }
    }

    return 0;
}

void calcAndDraw::allocate(){
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

void calcAndDraw::print_console(double delta_y){
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

void calcAndDraw::destroy(){
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

void calcAndDraw::change_func(){
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

void calcAndDraw::change_view(){
    view_id=(view_id+1)%5;
    if (qmgl) qmgl->update();
}

void calcAndDraw::enhance(){
    double t=(b-a)/2;
    a=a-t;
    b=b+t;
    ++s;
    if (qmgl) qmgl->update();
}

void calcAndDraw::squeeze(){
    double t=(b-a)/4;
    a=a+t;
    b=b-t;
    --s;
    if (qmgl) qmgl->update();
}

void calcAndDraw::more_points(){
    n*=2;
    if (qmgl) qmgl->update();
}

void calcAndDraw::less_points(){
    n/=2;
    if (qmgl) qmgl->update();
}

void calcAndDraw::delta_up(){
    ++p;
    if (qmgl) qmgl->update();
}

void calcAndDraw::delta_down(){
    --p;
    if (qmgl) qmgl->update();
}
