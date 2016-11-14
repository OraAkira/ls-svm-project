#include "mainwindow.h"
#include "ui_mainwindow.h"

static vector<pair<double, double> > variable;
static map<int, vector<pair<QString, int> > > stack_func;
static vector<double> slice_vals; // параметры сечений

static realiz prime;

float get_y(float x);

static int slice_curr = 0;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // Перемещение в центр экрана
    QPoint center = QApplication::desktop()->availableGeometry().center();
    QPoint corner = QApplication::desktop()->availableGeometry().topLeft();
    center.setX(center.x() - this->width() / 2);
    center.setY(center.y() - this->height() / 2);
    if(center.x() <= corner.x() || center.y() <= corner.y())
        this->move(corner);
    else
        this->move(center);

    // Заголовок формы
    this->setWindowTitle(trUtf8("LS SVM"));

    // Запрет изменения размеров формы
    this->setMinimumSize(this->size());
    this->setMaximumSize(this->size());

    // Запрет изменения полей вывода
    ui->lineEdit->setReadOnly(true);
    ui->lineEdit_2->setReadOnly(true);
    ui->lineEdit_3->setReadOnly(true);

    // Установим активной вкладку с исходными данными
    ui->tabWidget->setCurrentIndex(0);

    // ограничения на занесение данных в поля
    ui->spinBox->setMinimum(1); // ограничения на ко-во факторов
    ui->spinBox->setMaximum(5);

    ui->spinBox_2->setMinimum(1); // ограничения на ко-во экспериментов
    ui->spinBox_2->setMaximum(500);

    ui->doubleSpinBox_9->setMinimum(-20.0); // ограничения на уровни варьирования факторов
    ui->doubleSpinBox_9->setMaximum(20.0);
    ui->doubleSpinBox_8->setMinimum(-20.0);
    ui->doubleSpinBox_8->setMaximum(20.0);

    // ограничения на значение шума
    ui->doubleSpinBox_10->setMinimum(0.0);
    ui->doubleSpinBox_10->setMaximum(100.0);

    // ограничения на значение гаммы
    ui->doubleSpinBox_11->setMinimum(0.0001);
    ui->doubleSpinBox_11->setMaximum(99999.9);
    ui->doubleSpinBox_11->setDecimals(4);

    // ограничения на значение границ поиска гаммы
    // левая граница
    ui->doubleSpinBox_6->setMinimum(0.0001);
    ui->doubleSpinBox_6->setMaximum(99999.9999);
    ui->doubleSpinBox_6->setDecimals(4);
    // правая граница
    ui->doubleSpinBox_7->setMinimum(0.0001);
    ui->doubleSpinBox_7->setMaximum(99999.9999);
    ui->doubleSpinBox_7->setDecimals(4);

    // ограничения для параметров линейного ядра
    ui->doubleSpinBox_2->setMinimum(0.01);
    ui->doubleSpinBox_2->setMaximum(20);

    // ограничения для параметров полиномиального ядра
    // констанста с
    ui->doubleSpinBox_3->setMinimum(0.01);
    ui->doubleSpinBox_3->setMaximum(120);
    // константа a
    ui->doubleSpinBox_4->setMinimum(0.01);
    ui->doubleSpinBox_4->setMaximum(20);
    // степень полинома d
    ui->doubleSpinBox_5->setMinimum(0.1);
    ui->doubleSpinBox_5->setMaximum(2.8);

    // Ограничения для параметров RBF - ядра
    ui->doubleSpinBox->setMinimum(0.0003);
    ui->doubleSpinBox->setMaximum(3162.228);
    ui->doubleSpinBox->setDecimals(4);

    // Задаем умолчательные значения количества факторов
    ui->spinBox_2->setValue(10);

    // Задаем умолчательные значения варьирования факторов
    ui->doubleSpinBox_9->setValue(-10.0);
    ui->doubleSpinBox_8->setValue(10.0);
    int max_factor_number = ui->spinBox->maximum();
    variable.resize(max_factor_number);
    QStringList factor_text;
    for(int i = 0; i < max_factor_number; i++)
    {
        factor_text << QString::number(i + 1);
        variable[i] = make_pair(-10.0, 10.0);
    }
    ui->comboBox->clear();
    ui->comboBox->addItems(factor_text);
    ui->comboBox->setCurrentIndex(0);

    // Задаем умолчательные значения границ поиска гаммы
    ui->doubleSpinBox_6->setValue(0.0001);
    ui->doubleSpinBox_7->setValue(99999.9999);

    //Изменение шага накрутки счетчика для d полиномиального ядра
    ui->doubleSpinBox_5->setSingleStep(0.1);

    // Заполнение хранилища функций
    // Функции с числом факторов 1
    stack_func[1].push_back(make_pair(trUtf8("sin(x1)/x1"), 1));
    stack_func[1].push_back(make_pair(trUtf8("1 + 2*x1 + x1^x1"), 2));
    stack_func[1].push_back(make_pair(trUtf8("sin(x1) * cos(x1)"), 3));
    // Функции с числом факторов 2
    stack_func[2].push_back(make_pair(trUtf8("sin(x1)/x2"), 4));
    stack_func[2].push_back(make_pair(trUtf8("1 + 2*x1 + x2^2"), 5));
    stack_func[2].push_back(make_pair(trUtf8("sin(x1) * cos(x2)"), 6));
    // Функции с числом факторов 3
    stack_func[3].push_back(make_pair(trUtf8("(x1 * x1 + 2 * x2)/ x3"), 7));

    // Занесение списка на форму
    QStringList funcs_text;
    for(int i = 0; i < (int)stack_func[ui->spinBox->value()].size(); i++)
        funcs_text << stack_func[ui->spinBox->value()][i].first;
    ui->listWidget_2->clear();
    ui->listWidget_2->insertItems(0, funcs_text);
    ui->listWidget_2->setCurrentRow(0);

    ui->spinBox->setValue(1); // Число факторов
    on_spinBox_valueChanged(ui->spinBox->value());

    // активируем кнопки
    ui-> radioButton ->click();
    ui-> checkBox_9 ->click();
    ui-> checkBox_3 ->click();
    ui-> checkBox_5 ->click();
    ui-> checkBox_4 ->click();

    ui-> checkBox_6 ->click();
    ui-> checkBox_7 ->click();
    ui-> checkBox_10 ->click();
    ui-> checkBox_11 ->click();

    // Задание начальных данных сечений
    ui->doubleSpinBox_12->setMinimum(ui->doubleSpinBox_9->minimum());
    ui->doubleSpinBox_12->setMaximum(ui->doubleSpinBox_9->maximum());
    ui->doubleSpinBox_12->setValue(0.0);
    slice_vals.resize(ui->spinBox->maximum());
    for(int i = 0; i < (int)slice_vals.size(); i++)
        slice_vals[i] = 0.0;
    slice_curr = 0;

    // Инициализация виджета с графиком
    // Пока данных нет, не рисуем
    ui->Graph2dWidget->setDraw(false);
    // Цвета
    //ui->Graph2dWidget->setColorFigures(0.82f, 0.31f, 0.27f);
    ui->Graph2dWidget->setColorLine(0.28f, 0.48f, 0.6f);
    //ui->Graph2dWidget->setColorLine(0.45f, 0.69f, 0.08f);
    ui->Graph2dWidget->setColorFigures(0.45f, 0.69f, 0.08f);

    // Фигурки
    ui->Graph2dWidget->setFigure(QGraph2dFigures::CROSS); // CIRCLE, CROSS, RHOMBUS
    // Функция
    ui->Graph2dWidget->setLineFunc(get_y);
    // Число точек на линии
    ui->Graph2dWidget->setLinePointsNum(100);

    // Отключаем все элементы управления графиком
    ui->label_25->setVisible(false);
    ui->comboBox_2->setVisible(false);
    ui->label_26->setVisible(false);
    ui->comboBox_3->setVisible(false);
    ui->label_27->setVisible(false);
    ui->doubleSpinBox_12->setVisible(false);
    ui->checkBox_12->setVisible(false);
}

MainWindow::~MainWindow()
{
    delete ui;
}

// Функция для отправки в построитель графиков
float get_y(float x)
{
    if(slice_curr < 0 || prime.input.amt_factor < 1) return 0.0;
    if(prime.input.amt_factor > 1)
    {
        int index = slice_curr;
        double * x_a = new double [prime.input.amt_factor];
        for(int i = 0; i < prime.input.amt_factor; i++)
            x_a[i] = slice_vals[i];
        x_a[index] = x;
        double y = prime.estimate_get(x_a, prime.c + 1, prime.c[0]);
        delete [] x_a;
        return y;
    }
    double x_a = x;
    return prime.estimate_get(&x_a, prime.c + 1, prime.c[0]);
}

void MainWindow::adjust_zoom(double x_min, double x_max)
{
    float max_y = DBL_MIN, min_y = DBL_MAX;

    float x = x_min;
    int line_points_num = ui->Graph2dWidget->linePointsNum();
    float dx = (x_max - x_min) / (float)line_points_num;
    for(int i = 0; i <= line_points_num + 1; i++)
    {
        float y = get_y(x);
        if(y > max_y) max_y = y;
        if(y < min_y) min_y = y;
        x = x_min + dx * (float)i;
    }

  //  if(prime.input.amt_factor == 1 && ui->Graph2dWidget->figures())
    if(prime.input.amt_factor == 1)
    {
        for(int i = 0; i < prime.input.amt_exp; i++)
        {
            if(prime.y[i] > max_y) max_y = prime.y[i];
            if(prime.y[i] < min_y) min_y = prime.y[i];
        }
    }

    ui->Graph2dWidget->setMaximumX(x_max);
    ui->Graph2dWidget->setMinimumX(x_min);
    ui->Graph2dWidget->setMaximumY(max_y);
    ui->Graph2dWidget->setMinimumY(min_y);
    ui->Graph2dWidget->precond();
}

void MainWindow::on_pushButton_clicked()
{
    // Отключаем рисование
    ui->Graph2dWidget->setDraw(false);
    ui->Graph2dWidget->clear();

    // Очищаем старые данные, если они есть
    prime.clear_memory();

    // присвоение входных значений параметрам
    prime.input.amt_factor = ui->spinBox->value();          // количество факторов
    prime.input.amt_exp = ui->spinBox_2->value();           // количество экспериментов
    prime.input.num_func = stack_func
            [prime.input.amt_factor]
            [ui->listWidget_2->currentRow()]
            .second;                                        // номер функции
    prime.input.noise = ui->doubleSpinBox_10->value();      // процент зашумленности
    prime.input.gamma = ui->doubleSpinBox_11->value();      // параметр регуляризации
    prime.input.auto_gamma =ui->checkBox_9->isChecked();    // наличие автоподбора гаммы
    prime.input.left_border = ui->doubleSpinBox_6->value(); // левая граница поиска гаммы
    prime.input.right_border = ui->doubleSpinBox_7->value();// правая граница поиска гаммы

    // Если границы заданы неверно, меняем местами
    if(prime.input.left_border > prime.input.right_border)
        swap(prime.input.left_border, prime.input.right_border);

    // Выбор вида ядра
    if (ui->radioButton->isChecked())    prime.curr_kernel = 1; // линейное ядро
    if (ui->radioButton_2->isChecked())  prime.curr_kernel = 2; // RBF - ядро
    if (ui->radioButton_3->isChecked())  prime.curr_kernel = 3; // полиномиальное ядро

    // Заполнение параметров линейного ядра
    prime.kernel1.coef = ui->doubleSpinBox_2->value();
    // Заполнение параметров RBF - ядра
    prime.kernel2.sigma =  ui->doubleSpinBox->value();
    // Заполнение параметров полиномиального ядра
    prime.kernel3.degree = ui->doubleSpinBox_5->value();
    prime.kernel3.coef = ui->doubleSpinBox_3->value();
    prime.kernel3.a = ui->doubleSpinBox_4->value();

    // Наличие автоподбора параметров ядра
    if (prime.curr_kernel == 1) prime.input.auto_param = false;
    if (prime.curr_kernel == 2) prime.input.auto_param = ui->checkBox->isChecked();
    if (prime.curr_kernel == 3) prime.input.auto_param = ui->checkBox_2->isChecked();

    // Границы варьирования факторов
    prime.input.var_level.clear();
    for(int i = 0; i < prime.input.amt_factor; i++)
    {
        if(variable[i].first < variable[i].second)
            prime.input.var_level.push_back(pair<double, double>(variable[i].first, variable[i].second));
        else
            prime.input.var_level.push_back(pair<double, double>(variable[i].second, variable[i].first));
    }

    // Будем обновлять выборку или нет
    if(ui->checkBox_8->isChecked())
        srand(time(NULL));
    else
        srand(0x12344321U); //магическое число
    //  три константы, скорее всего одна из них отвечает за генератор
    // 0x12344321U, 0xEE11DD22U, 0xFEDCBA98

    // Вычисления
    prime.init();
    double fast_loo_val;
    QTime time_start, time_stop;
    try
    {
        time_start = QTime::currentTime();
        prime.ls_svm_slae();
        prime.auto_fit();
        fast_loo_val = prime.fast_loo_cv();
        prime.ls_svm_slae();
        prime.calc_y_estimate();
        time_stop = QTime::currentTime();
    }

    catch(int errnum)
    {
        QString error_descr;
        if(errnum == 0)
            error_descr = trUtf8("Неверно задана функция ядра");
        else if(errnum == 1)
            error_descr = trUtf8("Обратная матрица не существует");
        else if(errnum == 2)
            error_descr = trUtf8("Неверно задана функция отклика");
        else
            error_descr = trUtf8("Неизвестная ошибка");
        QMessageBox::critical(this, trUtf8("Ошибка"), error_descr, QMessageBox::Ok, QMessageBox::Ok);
        return;
    }

    // Вывод результата MSE в форму
    if (ui->checkBox_10->isChecked()) ui->lineEdit->setText(QString::number(prime.calc_MSE(), 'e', 7));
    else ui->lineEdit->setText(trUtf8("не задано"));

    // Вывод результата LOO в форму
    if (ui->checkBox_7->isChecked()) ui->lineEdit_2->setText(QString::number(fast_loo_val, 'e', 7));
    else ui->lineEdit_2->setText(trUtf8("не задано"));

    // Вывод отклика и оценок на форму в табличку
    ui->tableWidget->clear();
    if(ui->checkBox_6->isChecked())
    {
        // Зададим три колонки
        ui->tableWidget->setColumnCount(3);
        // Зададим заголовок
        QStringList header_text;
        header_text << trUtf8("Y") << trUtf8("Ŷ_learn") << trUtf8("Ŷ");
        ui->tableWidget->setHorizontalHeaderLabels(header_text);
        // Зададим стоки таблицы
        ui->tableWidget->setRowCount(prime.input.amt_exp);
        // Заполним таблицу
         cout << "yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy" << endl;
        for(int i = 0; i < prime.input.amt_exp; i++)
        {
            QTableWidgetItem * newItem;
            newItem = new QTableWidgetItem(QString::number(prime.y[i]));
            ui->tableWidget->setItem(i, 0, newItem);
            newItem = new QTableWidgetItem(QString::number(prime.y_new[i]));
            ui->tableWidget->setItem(i, 1, newItem);
            newItem = new QTableWidgetItem(QString::number(prime.y2_new[i])); // Тут новое название игрека
            ui->tableWidget->setItem(i, 2, newItem);

            cout << prime.y[i] << endl;
        }
        cout << "ooooooooooooooooooooooooooooooooooooooooo" << endl;
        for(int i = 0; i <  prime.input.amt_exp; i++)
           cout << prime.y2_new[i] << endl;
        // Урежем ширину колонок
        ui->tableWidget->setColumnWidth(0, 77);
        ui->tableWidget->setColumnWidth(1, 77);
        ui->tableWidget->setColumnWidth(2, 77);
    }

    // Вывод времени в форму
    if(ui->checkBox_11->isChecked())
    {
        double time = (double)time_start.msecsTo(time_stop) / 1000.0;
        ui->lineEdit_3->setText(QString::number(time));
    }
    else ui->lineEdit_3->setText(trUtf8("не задано"));

    // Передаем данные в рисовалку графика
    if(prime.input.amt_factor == 1)
    {
        for(int i = 0; i < prime.input.amt_exp; i++)
            ui->Graph2dWidget->data.push_back(QGraph2dPoint(prime.x[i][0], prime.y[i]));
        sort(ui->Graph2dWidget->data.begin(), ui->Graph2dWidget->data.end());

        // Отключим виджеты управления сечением
        ui->label_25->setVisible(false);
        ui->comboBox_2->setVisible(false);
        ui->label_26->setVisible(false);
        ui->comboBox_3->setVisible(false);
        ui->label_27->setVisible(false);
        ui->doubleSpinBox_12->setVisible(false);
        // Включим флажок отклика
        ui->checkBox_12->setVisible(true);
        ui->checkBox_12->setChecked(true);
        ui->Graph2dWidget->setFigures(true);

        adjust_zoom(prime.input.var_level[0].first, prime.input.var_level[0].second);
    }
    else
    {
        // Занесение факторов в сечение
        QStringList factor_text3;
        for(int i = 0; i < prime.input.amt_factor; i++)
            factor_text3 << QString::number(i + 1);
        ui->comboBox_2->clear();
        ui->comboBox_2->addItems(factor_text3);
        ui->comboBox_2->setCurrentIndex(0);
        on_comboBox_2_currentIndexChanged("1");

        // Включим виджеты управления сечением
        ui->label_25->setVisible(true);
        ui->comboBox_2->setVisible(true);
        ui->label_26->setVisible(true);
        ui->comboBox_3->setVisible(true);
        ui->label_27->setVisible(true);
        ui->doubleSpinBox_12->setVisible(true);
        // Отключим флажок отклика
        ui->checkBox_12->setVisible(false);
        ui->Graph2dWidget->setFigures(false);
    }
    ui->Graph2dWidget->setDraw(true);

    // Установим активной вкладку с результатами
    ui->tabWidget->setCurrentIndex(1);
}

// При выборе радиокнопки линейного ядра
void MainWindow::on_radioButton_clicked()
{
    ui->groupBox_3->setVisible(true);
    ui->groupBox_2->setVisible(false);
    ui->groupBox_4->setVisible(false);
}

// При выборе радиокнопки полиномиального ядра
void MainWindow::on_radioButton_3_clicked()
{
    ui->groupBox_3->setVisible(false);
    ui->groupBox_2->setVisible(false);
    ui->groupBox_4->setVisible(true);
}

// При выборе радиокнопки RBF - ядра
void MainWindow::on_radioButton_2_clicked()
{
    ui->groupBox_3->setVisible(false);
    ui->groupBox_2->setVisible(true);
    ui->groupBox_4->setVisible(false);
}

// При активном автоподборе гаммы блок с вводом границ скрывается и наоборот
void MainWindow::on_checkBox_9_clicked()
{
    if (ui->checkBox_9->isChecked())ui->groupBox_6->setEnabled(true);
    else ui->groupBox_6->setEnabled(false);
}


// Установка значений "по умолчанию" для линейного ядра
void MainWindow::on_checkBox_5_clicked()
{
    if (ui->checkBox_5->isChecked())
    {
        ui->doubleSpinBox_2->setValue(1.0);
        ui->doubleSpinBox_2->setEnabled(false);
    }
    else  ui->doubleSpinBox_2->setEnabled(true);
}

// Установка значений "по умолчанию" для полиномиального ядра
void MainWindow::on_checkBox_3_clicked()
{
    if (ui->checkBox_3->isChecked())
    {
        ui->doubleSpinBox_3->setValue(1.0);
        ui->doubleSpinBox_4->setValue(1.0);
        ui->doubleSpinBox_5->setValue(2.0);
        ui->doubleSpinBox_3->setEnabled(false); // поле неактивно
        ui->doubleSpinBox_4->setEnabled(false);
        ui->doubleSpinBox_5->setEnabled(false);
    }
    else
    {
        ui->doubleSpinBox_3->setEnabled(true);// поле активно
        ui->doubleSpinBox_4->setEnabled(true);
        ui->doubleSpinBox_5->setEnabled(true);
    }
}

// Установка значений "по умолчанию" для RBF - ядра

void MainWindow::on_checkBox_4_clicked()
{
    if (ui->checkBox_4->isChecked())
    {
        ui->doubleSpinBox->setValue(prime.kernel2.sigma_select[0]);
        ui->doubleSpinBox->setEnabled(false);
    }
    else  ui->doubleSpinBox->setEnabled(true);
}

// Действие при изменении числа факторов
void MainWindow::on_spinBox_valueChanged(int arg1)
{
    // Занесение списка на форму
    QStringList funcs_text;
    for(int i = 0; i < (int)stack_func[arg1].size(); i++)
        funcs_text << stack_func[arg1][i].first;
    ui->listWidget_2->clear();
    ui->listWidget_2->insertItems(0, funcs_text);
    ui->listWidget_2->setCurrentRow(0);

    // Занесение номеров факторов на форму
    QStringList factor_text;
    for(int i = 0; i < arg1; i++)
    {
        factor_text << QString::number(i + 1);
    }
    ui->comboBox->clear();
    ui->comboBox->addItems(factor_text);
}

// Действие при выборе номера фактора
void MainWindow::on_comboBox_currentIndexChanged(int index)
{
    if(index < 0 || !variable.size()) return;
    ui->doubleSpinBox_9->setValue(variable[index].first);
    ui->doubleSpinBox_8->setValue(variable[index].second);
}

// Изменение левой границы варьирования факторов
void MainWindow::on_doubleSpinBox_9_valueChanged(double arg1)
{
    int index = ui->comboBox->currentIndex();
    if(index < 0 || !variable.size()) return;
    variable[index].first = arg1;
}

// Изменение правой границы варьирования факторов
void MainWindow::on_doubleSpinBox_8_valueChanged(double arg1)
{
    int index = ui->comboBox->currentIndex();
    if(index < 0 || !variable.size()) return;
    variable[index].second = arg1;
}

// Изменение значения переменной сечения
void MainWindow::on_doubleSpinBox_12_valueChanged(double arg1)
{
    bool draw = ui->Graph2dWidget->draw();
    ui->Graph2dWidget->setDraw(false);

    int index = ui->comboBox_3->currentText().toInt() - 1;
    if(index < 0 || !slice_vals.size()) return;
    slice_vals[index] = arg1;

    if(prime.input.var_level.size() && slice_curr >= 0)
        adjust_zoom(prime.input.var_level[slice_curr].first, prime.input.var_level[slice_curr].second);

    ui->Graph2dWidget->setDraw(draw);
}

// Изменение переменной сечения
void MainWindow::on_comboBox_3_currentIndexChanged(const QString &arg1)
{
    bool draw = ui->Graph2dWidget->draw();
    ui->Graph2dWidget->setDraw(false);

    int index = arg1.toInt() - 1;
    if(index < 0 || !slice_vals.size()) return;
    ui->doubleSpinBox_12->setValue(slice_vals[index]);
    if(prime.input.var_level.size())
    {
        ui->doubleSpinBox_12->setMinimum(prime.input.var_level[index].first);
        ui->doubleSpinBox_12->setMaximum(prime.input.var_level[index].second);
    }

    if(prime.input.var_level.size() && slice_curr >= 0)
        adjust_zoom(prime.input.var_level[slice_curr].first, prime.input.var_level[slice_curr].second);

    ui->Graph2dWidget->setDraw(draw);
}

// Изменение переменной по которой сечение
void MainWindow::on_comboBox_2_currentIndexChanged(const QString &arg1)
{
    bool draw = ui->Graph2dWidget->draw();
    ui->Graph2dWidget->setDraw(false);

    QStringList factor_text2;
    int index_slice = arg1.toInt() - 1;
    slice_curr = index_slice;
    for(int i = 0; i < ui->spinBox->value(); i++)
    {
        if(i != index_slice)
            factor_text2 << QString::number(i + 1);
    }
    ui->comboBox_3->clear();
    ui->comboBox_3->addItems(factor_text2);

    if(prime.input.var_level.size() && slice_curr >= 0)
        adjust_zoom(prime.input.var_level[slice_curr].first, prime.input.var_level[slice_curr].second);

    ui->Graph2dWidget->setDraw(draw);
}

// Включение отображения отклика
void MainWindow::on_checkBox_12_clicked(bool checked)
{
    if(checked)
        ui->Graph2dWidget->setFigures(true);
    else
        ui->Graph2dWidget->setFigures(false);
    adjust_zoom(prime.input.var_level[0].first, prime.input.var_level[0].second);
}
