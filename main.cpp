#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>
#include "window.hpp"

int main (int argc, char *argv[]){
    QApplication app(argc, argv);
    QMainWindow *window=new QMainWindow;
    QMenuBar *tool_bar=new QMenuBar(window);
    Window *graph_area=new Window(window);
    QAction *action;
    if(graph_area->parse_command_line(argc, argv)){
        QMessageBox::warning(0, "Bad input!",
                            "Bad input!");
        return -1;
    }
    action=tool_bar->addAction("function", graph_area, SLOT(change_func()));
    action->setShortcut(QString("0"));
    action=tool_bar->addAction ("methods", graph_area, SLOT(change_view()));
    action->setShortcut(QString("1"));
    action=tool_bar->addAction("enhance", graph_area, SLOT(enhance()));
    action->setShortcut(QString("2"));
    action=tool_bar->addAction("squeeze", graph_area, SLOT(squeeze()));
    action->setShortcut(QString("3"));
    action=tool_bar->addAction("more_points", graph_area, SLOT(more_points()));
    action->setShortcut(QString("4"));
    action=tool_bar->addAction("less_points", graph_area, SLOT(less_points()));
    action->setShortcut(QString("5"));
    action=tool_bar->addAction("delta+", graph_area, SLOT(delta_up()));
    action->setShortcut(QString("6"));
    action=tool_bar->addAction("delta-", graph_area, SLOT(delta_down()));
    action->setShortcut(QString("7"));
    action=tool_bar->addAction("exit", window, SLOT(close()));
    action->setShortcut(QString("Ctrl+X"));
    tool_bar->setMaximumHeight(30);
    window->setMenuBar(tool_bar);
    window->setCentralWidget(graph_area);
    window->show();
    app.exec();
    delete window;
    return 0;
}
