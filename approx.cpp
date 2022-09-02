#include <QDebug>
#include <QCommandLineOption>
#include <QCommandLineParser>
#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QScrollArea>
#include <QtWidgets/QMenuBar>
#include "calc.h"

enum CommandLineParseResult
{
    CommandLineOk,
    CommandLineError,
    CommandLineVersionRequested,
    CommandLineHelpRequested
};

double a, b;
int n,k;

CommandLineParseResult parseCommandLine(QCommandLineParser &parser, QString *errorMessage)
{
    parser.setApplicationDescription(QCoreApplication::translate("main",
                                     "Construction of approximations of functions on a segment AB."));
//    QCommandLineOption dimension(QStringList() << "d" << "dimension",
//                                 QCoreApplication::translate("main", "Function dimension select (2,3)."),
//                                 QCoreApplication::translate("main", "dim"));
//    parser.addOption(dimension);
    QCommandLineOption seqmentA(QStringList() << "a" << "coordA",
                                QCoreApplication::translate("main", "The coordinate of the end A of the segment."),
                                QCoreApplication::translate("main", "endSegmentA"));
    parser.addOption(seqmentA);
    QCommandLineOption seqmentB(QStringList() << "b" << "coordB",
                                QCoreApplication::translate("main", "The coordinate of the end B of the segment."),
                                QCoreApplication::translate("main", "endSegmentB"));
    parser.addOption(seqmentB);

    QCommandLineOption countP(QStringList() << "n" << "count",
                              QCoreApplication::translate("main", "Points count."),
                              QCoreApplication::translate("main", "segmentPointsCount"));
    parser.addOption(countP);
    QCommandLineOption funcType(QStringList() << "k" << "type",
                                QCoreApplication::translate("main", "Type of function.\n"
                                                                    "k=0: f(x)=1\nk=1: f(x)=x\nk=2: f(x)=x^2\nk=3: f(x)=x^3\nk=4: f(x)=x^4\nk=5: f(x)=exp(x)\nk=6: f(x)=1/(25x+1)"),
                                QCoreApplication::translate("main", "typeOfFunc"));
    parser.addOption(funcType);
    QCommandLineOption helpOption = parser.addHelpOption();
    QCommandLineOption versionOption = parser.addVersionOption();

    if (!parser.parse(QCoreApplication::arguments())) {
        *errorMessage = parser.errorText();
        return CommandLineError;
    }

    if (parser.isSet(versionOption))
        return CommandLineVersionRequested;

    if (parser.isSet(helpOption))
        return CommandLineHelpRequested;

    a = parser.value(seqmentA).toDouble();
    b = parser.value(seqmentB).toDouble();
    n = parser.value(countP).toInt();
    k = parser.value(funcType).toInt();

    return CommandLineOk;
}

int main(int argc,char **argv)
{
  QApplication app(argc,argv);
  QCommandLineParser parser;

  QCoreApplication::setApplicationName("approx");
  QCoreApplication::setApplicationVersion("Version: 1.0");

  QString errorMessage;
  switch (parseCommandLine(parser, &errorMessage)) {
  case CommandLineOk:
      break;
  case CommandLineError:
      fputs(qPrintable(errorMessage), stderr);
      fputs("\n\n", stderr);
      fputs(qPrintable(parser.helpText()), stderr);
      return 1;
  case CommandLineVersionRequested:
      printf("%s %s\n", qPrintable(QCoreApplication::applicationName()),
             qPrintable(QCoreApplication::applicationVersion()));
      return 0;
  case CommandLineHelpRequested:
      parser.showHelp();
      Q_UNREACHABLE();
  }

  QMainWindow *wnd = new QMainWindow;
  QMenuBar *toolbar = new QMenuBar(wnd);

  wnd->resize(810,610);  // for fill up the QMGL, menu and toolbars
  wnd->setWindowTitle("Approximation");
  // here I allow to scroll QMathGL -- the case
  // then user want to prepare huge picture
  calcAndDraw *draw = new calcAndDraw(a, b, n, k, wnd);
  QScrollArea *scroll = new QScrollArea(draw);

  QAction *actionFunc = toolbar->addAction("Function", draw, SLOT(change_func()));
  actionFunc->setShortcut(QString("0"));
  QAction *actionMeth = toolbar->addAction ("Methods", draw, SLOT(change_view()));
  actionMeth->setShortcut(QString("1"));
  QAction *actionEnch = toolbar->addAction("Enhance", draw, SLOT(enhance()));
  actionEnch->setShortcut(QString("2"));
  QAction *actionSq = toolbar->addAction("Squeeze", draw, SLOT(squeeze()));
  actionSq->setShortcut(QString("3"));
  QAction *actionMore = toolbar->addAction("More_points", draw, SLOT(more_points()));
  actionMore->setShortcut(QString("4"));
  QAction *actionLess = toolbar->addAction("Less_points", draw, SLOT(less_points()));
  actionLess->setShortcut(QString("5"));
  QAction *actionDeltaPl = toolbar->addAction("Delta+", draw, SLOT(delta_up()));
  actionDeltaPl->setShortcut(QString("6"));
  QAction *actionDeltaMin = toolbar->addAction("Delta-", draw, SLOT(delta_down()));
  actionDeltaMin->setShortcut(QString("7"));

  // Create and setup QMathGL
  QMathGL *QMGL = new QMathGL(draw);
  draw->SetMGL(QMGL);
  QMGL->setDraw(draw);
  QMGL->update();  

  // continue other setup (menu, toolbar and so on)
  QAction *actionExit = toolbar->addAction("Exit", draw, SLOT(close()));
  actionExit->setShortcut(QString("Ctrl+X"));
  toolbar->setMaximumHeight(30);
  draw->setMenuBar(toolbar);
  scroll->setWidget(QMGL);
  draw->setCentralWidget(scroll);
  draw->show();

  return app.exec();
}
