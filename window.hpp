#ifndef WINDOW_H
#define WINDOW_H

#include <QtWidgets/QtWidgets>

class Window : public QWidget{
    Q_OBJECT
private:
    int func_id; //approximating fucntion
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
    double *y;//function values array
    double *dy;//derivative values array
    double *sp;//splines coefficients array
    double *h;//Hermite coefficients array
    double *c;//Chebyshev coefficients array
    double max_y, min_y, absmax;
public:
    Window(QWidget *parent);
//    ~Window();
    QSize minimumSizeHint()const;
    QSize sizeHint()const;
    int parse_command_line(int argc, char *argv[]);
public slots:
    void change_func();
    void change_view();
    void enhance();
    void squeeze();
    void more_points();
    void less_points();
    void delta_up();
    void delta_down();
protected:
    void paintEvent (QPaintEvent *event);
};

#endif
