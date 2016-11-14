#ifndef GRAPH2D_H
#define GRAPH2D_H

#include <QGLWidget>
#include <QString>
#include <GL/glu.h>
#include <GL/glext.h>
#include <vector>
#include <QMutex>
#include "mainwindow.h"

using namespace std;

namespace QGraph2dFigures
{
    const int CIRCLE  = 0;
    const int CROSS   = 1;
    const int RHOMBUS = 2;
}

class QGraph2dPoint
{
public:
    QGraph2dPoint();
    QGraph2dPoint(double x, double y);
    float x;
    float y;
    bool operator < (const QGraph2dPoint & t) const;
};

class QGraph2d : public QGLWidget
{
    Q_OBJECT
protected:
    void timerEvent(QTimerEvent *);
public:
    vector<QGraph2dPoint> data;
    void initializeGL();
    void resizeGL(int nWidth, int nHeight);
    void paintGL();
    QGraph2d(QWidget* parent = 0);
    void precond();
    void setDraw(bool status);
    bool draw();
    void setFigures(bool status);
    bool figures();
    void clear();
    void setMaximumX(double max_x);
    void setMaximumY(double max_y);
    void setMinimumX(double min_x);
    void setMinimumY(double min_y);
    void setColorFigures(float r, float g, float b);
    void setColorLine(float r, float g, float b);
    void setLineFunc(float(*line_func)(float x));
    void setFigure(int figure_num);
    void setLinePointsNum(int line_points_num);
    int linePointsNum();
private:
    void adjustAxis(float & min, float & max, int & numTicks);
    QString axis_x, axis_y;
    float max_x, max_y;
    float min_x, min_y;
    float size_x, size_y;
    QMutex mtx;
    bool need_draw;
    int num_ticks_x, num_ticks_y;
    float diff;
    GLfloat color_figures[3];
    GLfloat color_line[3];
    float(*line_func)(float x);
    int figure_num;
    int line_points_num;
    bool need_figures;
};

#endif // GRAPH2D_H

