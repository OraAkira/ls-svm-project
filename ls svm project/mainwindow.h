#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "LS_SVM.h"
#include <QMainWindow>
#include <QApplication>
#include <QMessageBox>
#include <QDesktopWidget>
#include <QPoint>
#include <QMainWindow>
#include <QStringList>
#include <QString>
#include <QTime>
#include <vector>
#include <map>
#include <string>
#include <utility>
#include <cfloat>
using namespace std;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    
private slots:
    void on_pushButton_clicked();

    void on_radioButton_clicked();

    void on_radioButton_3_clicked();

    void on_radioButton_2_clicked();

    void on_checkBox_9_clicked();

    void on_checkBox_5_clicked();

    void on_checkBox_3_clicked();

    void on_checkBox_4_clicked();


    void on_spinBox_valueChanged(int arg1);

    void on_comboBox_currentIndexChanged(int index);

    void on_doubleSpinBox_9_valueChanged(double arg1);

    void on_doubleSpinBox_8_valueChanged(double arg1);

    void on_doubleSpinBox_12_valueChanged(double arg1);

    void on_comboBox_3_currentIndexChanged(const QString &arg1);

    void on_comboBox_2_currentIndexChanged(const QString &arg1);

    void on_checkBox_12_clicked(bool checked);

    void on_checkBox_13_clicked(bool checked);

    void on_checkBox_14_clicked(bool checked);

    void on_checkBox_15_clicked(bool checked);

private:
    Ui::MainWindow *ui;

    void adjust_zoom(double x_min, double x_max);
};

#endif // MAINWINDOW_H
