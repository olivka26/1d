#ifndef CALC_H
#define CALC_H

#include <QMainWindow>
#include <QObject>
#include <mgl2/qmathgl.h>

#define DEFAULT_A -10.0
#define DEFAULT_B 10.0
#define DEFAULT_N 30
#define DEFAULT_K 0

class calcAndDraw : public QMainWindow, public mglDraw
{
    Q_OBJECT

private:
    const char *f_name;
    double a; //left end
    double b; //right end
    int n; //amount of points
    int k; //id of the approxiamted function
    double(*f)(double); //function
    double(*df)(double);//derivative
    double(*d2f)(double);//second derivative
    int view_id=0;//what to be viewed
    int s=0;//squeeze-stretch
    //int delta=0;
    int p=0;//desturbance;
    double *x;//points array
    double *cx;//Chebyshev points array
    double *cy;//Chebyshev values array
    double *y;//function values array
    double *dy;//derivative values array
    double *sp;//splines coefficients array
    double *h;//Hermite coefficients array
    double *c;//Chebyshev coefficients array
    double *extr;
    double max_y, min_y, absmax;

    QMathGL *qmgl;
public:
    calcAndDraw(double startA = DEFAULT_A,
                double startB = DEFAULT_B,
                int startN = DEFAULT_N,
                int startK = DEFAULT_K,
                QMainWindow *parent = 0) : QMainWindow(parent), mglDraw() {a = startA;
                                                                           b = startB;
                                                                           n = startN;
                                                                           k = startK;
                                                                           change_func();}
    void SetMGL(QMathGL *mgl) { qmgl = mgl; }
    int Draw(mglGraph *gr);
    void allocate();
    void print_console(double delta_y);
    void destroy();

    private slots:
    void change_func();
    void change_view();
    void enhance();
    void squeeze();
    void more_points();
    void less_points();
    void delta_up();
    void delta_down();
};

#endif // CALC_H
