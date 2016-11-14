#include "graph2d.h"

#include <QtGui>
#include <QFont>
#include <cmath>
#include <cfloat>
#include <QString>
#include <QTextCodec>

const GLfloat SIZE_CIRCLE  = 0.014f; // размер круга
const GLfloat SIZE_CROSS   = 0.014f; // размер креста
const GLfloat SIZE_RHOMBUS = 0.015f; // размер ромба

const GLfloat WIDTH_FIGURES = 2.0f; // толщина линий фигурок
const GLfloat WIDTH_LINES   = 2.0f; // толщина линий графика

bool QGraph2dPoint::operator < (const QGraph2dPoint & t) const
{
    return x < t.x;
}

QGraph2dPoint::QGraph2dPoint() { }

QGraph2dPoint::QGraph2dPoint(double x, double y)
{
    this->x = (float)x;
    this->y = (float)y;
}

void draw_circle(GLfloat x, GLfloat y, GLfloat diff)
{
    const GLfloat radius = SIZE_CIRCLE;
    const int num_segments = 20;
    glBegin(GL_POLYGON);
    for(int i = 0; i < num_segments; i++)
    {
        float theta = 2.0f * 3.1415926f * float(i) / float(num_segments);
        float x_s = radius * cos(theta) * diff;
        float y_s = radius * sin(theta) / diff;
        glVertex2f(x + x_s, y + y_s);
    }
    glEnd();
}

void draw_cross(GLfloat x, GLfloat y, GLfloat diff)
{
    const GLfloat length = SIZE_CROSS;
    glBegin(GL_LINES);
    glVertex2d(x - length * diff, y - length / diff);
    glVertex2d(x + length * diff, y + length / diff);
    glVertex2d(x - length * diff, y + length / diff);
    glVertex2d(x + length * diff, y - length / diff);
    glEnd();
}

void draw_rhombus(GLfloat x, GLfloat y, float diff)
{
    const GLfloat length = SIZE_RHOMBUS;
    glBegin(GL_LINE_STRIP);
    glVertex2d(x - length * diff, y);
    glVertex2d(x, y + length / diff);
    glVertex2d(x + length * diff, y);
    glVertex2d(x, y - length / diff);
    glVertex2d(x - length * diff, y);
    glEnd();
}

void(*figures_func[3])(GLfloat x, GLfloat y, GLfloat diff) =
{
    draw_circle,
    draw_cross,
    draw_rhombus
};

void QGraph2d::setMaximumX(double max_x)
{
    this->max_x = (float)max_x;
}

void QGraph2d::setMaximumY(double max_y)
{
    this->max_y = (float)max_y;
}

void QGraph2d::setMinimumX(double min_x)
{
    this->min_x = (float)min_x;
}

void QGraph2d::setMinimumY(double min_y)
{
    this->min_y = (float)min_y;
}

void QGraph2d::setColorFigures(float r, float g, float b)
{
    color_figures[0] = r;
    color_figures[1] = g;
    color_figures[2] = b;
}

void QGraph2d::setColorLine(float r, float g, float b)
{
    color_line[0] = r;
    color_line[1] = g;
    color_line[2] = b;
}

void QGraph2d::setLineFunc(float(*line_func)(float x))
{
    this->line_func = line_func;
}

void QGraph2d::setFigure(int figure_num)
{
    this->figure_num = figure_num;
}

void QGraph2d::setLinePointsNum(int line_points_num)
{
    this->line_points_num = line_points_num;
}

int QGraph2d::linePointsNum()
{
    return line_points_num;
}

QGraph2d::QGraph2d(QWidget* parent) : QGLWidget(parent)
{
    startTimer(500);
    need_draw = false;
    axis_x = trUtf8("x");
    axis_y = trUtf8("y");
    need_figures = false;
}

void QGraph2d::timerEvent(QTimerEvent *)
{
    updateGL();
}

void QGraph2d::initializeGL()
{
    qglClearColor(Qt::white);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_FLAT);
    glEnable(GL_CULL_FACE);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
}

void QGraph2d::resizeGL(int nWidth, int nHeight)
{
    float ortho_x_beg = -0.12, ortho_x_end = 1.05;
    float ortho_y_beg = /*-0.1*/-0.2, ortho_y_end = /*1.13*/1.08;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(ortho_x_beg, ortho_x_end, ortho_y_beg, ortho_y_end, -10.0, 1.0);
    glViewport(0, 0,(GLint)nWidth, (GLint)nHeight);
    diff = (float)nHeight / (float)nWidth * (ortho_x_end - ortho_x_beg) / (ortho_y_end - ortho_y_beg);
    diff *= 1.5;
}

void QGraph2d::setDraw(bool status)
{
    mtx.lock();
    need_draw = status;
    mtx.unlock();
}

bool QGraph2d::draw()
{
    bool result;
    mtx.lock();
    result = need_draw;
    mtx.unlock();
    return result;
}

void QGraph2d::setFigures(bool status)
{
    need_figures = status;
}

bool QGraph2d::figures()
{
    return need_figures;
}

void QGraph2d::clear()
{
    mtx.lock();
    need_draw = false;
    data.clear();
    mtx.unlock();
}

void QGraph2d::adjustAxis(float & min, float & max, int & numTicks)
{
    const double axis_epsilon = 1.0 / 10000.0;
    if(max - min < axis_epsilon)
    {
        min -= 2.0 * axis_epsilon;
        max += 2.0 * axis_epsilon;
    }

    const int MinTicks = 4;
    double grossStep = (max - min) / MinTicks;
    double step = pow(10, floor(log10(grossStep)));

    if (5 * step < grossStep)
        step *= 5;
    else if (2 * step < grossStep)
        step *= 2;

    numTicks = (int)(ceil(max / step) - floor(min / step));
    min = floor(min / step) * step;
    max = ceil(max / step) * step;
}

void QGraph2d::precond()
{
    // Поправляем значения мин / макс чтобы влазило в сетку
    adjustAxis(min_x, max_x, num_ticks_x);
    adjustAxis(min_y, max_y, num_ticks_y);
    size_x = max_x - min_x;
    size_y = max_y - min_y;
}

void QGraph2d::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if(need_draw)
    {
        // координатные оси
        glColor3f(0.0f, 0.0f, 0.0f);
        glLineWidth(2.0f);

        glBegin(GL_LINES);
        glVertex2d(0.0, -0.005);
        glVertex2d(0.0, 1.02);

        glVertex2d(-0.005, 0.0);
        glVertex2d(1.02, 0.0);
        glEnd();

        // подписи осей
        QFont font("Courier", 8);
        QFont font2("Times", 10);
        font.setLetterSpacing(QFont::PercentageSpacing, 75.0);
        renderText(1.01, -0.06, 0.0f, axis_x, font2);
        renderText(-0.05, 1.01f, 0.0f, axis_y, font2);

        // координатная сетка
        glColor3f(0.7f,0.7f,0.7f);
        glLineWidth(1.0f);
        for (int i = 0; i <= num_ticks_x; ++i)
        {
            float x = (float)i / (float)num_ticks_x;
            glBegin(GL_LINES);
            glVertex2f(x, -0.01f);
            glVertex2f(x, 1.0f);
            glEnd();
        }
        for (int i = 0; i <= num_ticks_y; ++i)
        {
            float y = (float)i / (float)num_ticks_y;
            glBegin(GL_LINES);
            glVertex2f(-0.01f, y);
            glVertex2f(1.0f, y);
            glEnd();
        }

        // отрисовка шкалы
        glColor3f(0.0f,0.0f,0.0f);
        glLineWidth(2.0f);
        for (int i = 0; i < num_ticks_x; ++i)
        {
            float x = (float)i / (float)num_ticks_x;
            float x_real = (float)(floor((x * size_x + min_x) * 10000.0 + 0.5)) / 10000.0;
            QString st = QString::number(x_real);
            renderText(x - 0.01f, -0.06f, 0.001f, st, font);
        }
        for (int i = 0; i < num_ticks_y; ++i)
        {
            float y = (float)i / (float)num_ticks_y;
            float y_real = (float)(floor((y * size_y + min_y) * 10000.0 + 0.5)) / 10000.0;
            QString st = QString("%1").arg(y_real, 5, 'g', -1, ' ');
            renderText(-0.1f, y - 0.01f, 0.001f, st, font);
        }

        // рисуем фигурки и график
        mtx.lock();
        if(need_draw)
        {
            if(need_figures)
            {
                // рисуем фигурки
                glColor3f(color_figures[0], color_figures[1], color_figures[2]);
                glLineWidth(WIDTH_FIGURES);
                for(int i = 0; i < (int)data.size(); i++)
                {
                    float x_loc = (data[i].x - min_x) / size_x;
                    float y_loc = (data[i].y - min_y) / size_y;
                    figures_func[figure_num](x_loc, y_loc, diff);
                }
            }

            // отрисовка легенды
            float x_leg = 0.15f;
            float y_leg = /*1.07f*/-0.15f;
            glColor3f(0.0f,0.0f,0.0f);
            renderText(x_leg, y_leg - 0.015, 0.001f, trUtf8("Ŷ"), font2);
            // Рисуем точки графика
            glColor3f(color_line[0], color_line[1], color_line[2]);
            glLineWidth(WIDTH_FIGURES);
            glBegin(GL_LINES);
            glVertex2f(x_leg - 0.02f, y_leg);
            glVertex2f(x_leg - 0.12f, y_leg);
            glEnd();
            x_leg += 0.32f;
            if(need_figures)
            {
                glColor3f(0.0f,0.0f,0.0f);
                renderText(x_leg, y_leg - 0.015, 0.001f, trUtf8("Y"), font2);
                // рисуем фигурки
                glColor3f(color_figures[0], color_figures[1], color_figures[2]);
                glLineWidth(WIDTH_FIGURES);
                figures_func[figure_num](x_leg - 0.037f, y_leg, diff);
            }



            // отрисовка графика
            glColor3f(color_line[0], color_line[1], color_line[2]);
            glLineWidth(WIDTH_LINES);
            glBegin(GL_LINE_STRIP);
            float x = min_x;
            float dx = size_x / (float)line_points_num;
            for(int i = 0; i <= line_points_num + 1/* && x <= max_x*/; i++)
            {
                float y = line_func(x);
                float x_loc = (x - min_x) / size_x;
                float y_loc = (y - min_y) / size_y;
                glVertex3f(x_loc, y_loc, 0.0001);
                x = min_x + dx * (float)i;
            }
            glEnd();
        }
        mtx.unlock();
    }
}


