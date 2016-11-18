#include "mainwindow.h"
#include "ui_mainwindow.h"

static vector<pair<double, double> > variable;
static map<int, vector<pair<QString, int> > > stack_func;
static vector<double> slice_vals; // параметры сечений

static loo_cv_realiz prime_loo_cv;
static regularity_realiz prime_reg;
static weight_loo_cv_realiz prime_weight_loo_cv;

float get_y(float x);
float get_y2(float x);
float get_y3(float x);

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

    // Задаем умолчательные значения количества наблюдений
    ui->spinBox_2->setValue(100);

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
    ui-> radioButton_2 ->click();
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
    ui->Graph2dWidget->setColorLine(0, 0.28f, 0.48f, 0.6f);
    ui->Graph2dWidget->setColorLine(1, 0.82f, 0.31f, 0.27f);
    ui->Graph2dWidget->setColorLine(0, 0.28f, 0.48f, 0.6f);
    //ui->Graph2dWidget->setColorLine(2, 0.06f, 0.05f, 0.04f);
    ui->Graph2dWidget->setColorFigures(0.45f, 0.69f, 0.08f);

    // Фигурки
    ui->Graph2dWidget->setFigure(QGraph2dFigures::CROSS); // CIRCLE, CROSS, RHOMBUS
    // Функция
    ui->Graph2dWidget->setLineFunc(0, get_y);
    ui->Graph2dWidget->setLineFunc(1, get_y2);
    ui->Graph2dWidget->setLineFunc(2, get_y3);

    // Число точек на линии
    ui->Graph2dWidget->setLinePointsNum(100);

    // Отключаем все элементы управления графиком
    ui->label_25->setVisible(false);
    ui->comboBox_2->setVisible(false);
    ui->label_26->setVisible(false);
    ui->comboBox_3->setVisible(false);
    ui->checkBox_12->setVisible(false);

    // Всегда видимы, по умолчанию активны
    ui->checkBox_13->setChecked(true);
    ui->label_27->setVisible(false);
    ui->doubleSpinBox_12->setVisible(false);
    ui->checkBox_14->setChecked(true);
    ui->checkBox_15->setChecked(true);
    ui->Graph2dWidget->setLines(0, true);
    ui->Graph2dWidget->setLines(1, true);
    ui->Graph2dWidget->setLines(2, true);
}

MainWindow::~MainWindow()
{
    delete ui;
}

// Функция для отправки в построитель графиков
float get_y(float x)
{
    prime_loo_cv.input.weighted = 0;
    if(slice_curr < 0 || prime_loo_cv.input.amt_factor < 1) return 0.0;
    if(prime_loo_cv.input.amt_factor > 1)
    {
        int index = slice_curr;
        double * x_a = new double [prime_loo_cv.input.amt_factor];
        for(int i = 0; i < prime_loo_cv.input.amt_factor; i++)
            x_a[i] = slice_vals[i];
        x_a[index] = x;
        double y = prime_loo_cv.estimate_get(x_a, prime_loo_cv.c_normal + 1, prime_loo_cv.c_normal[0]);
        delete [] x_a;
        return y;
    }
    double x_a = x;
    return prime_loo_cv.estimate_get(&x_a, prime_loo_cv.c_normal + 1, prime_loo_cv.c_normal[0]);
}

float get_y2(float x)
{
    prime_loo_cv.input.weighted = 0;
    if(slice_curr < 0 || prime_weight_loo_cv.input.amt_factor < 1) return 0.0;
    if(prime_weight_loo_cv.input.amt_factor > 1)
    {
        int index = slice_curr;
        double * x_a = new double [prime_weight_loo_cv.input.amt_factor];
        for(int i = 0; i < prime_weight_loo_cv.input.amt_factor; i++)
            x_a[i] = slice_vals[i];
        x_a[index] = x;
        double y = prime_weight_loo_cv.estimate_get(x_a, prime_weight_loo_cv.c_WCV + 1, prime_weight_loo_cv.c_WCV[0]);
        delete [] x_a;
        return y;
    }
    double x_a = x;
    return prime_weight_loo_cv.estimate_get(&x_a, prime_weight_loo_cv.c_WCV + 1, prime_weight_loo_cv.c_WCV[0]);
    }
float get_y3(float x)
{
    prime_loo_cv.input.weighted = 1;
    if(slice_curr < 0 || prime_loo_cv.input.amt_factor < 1) return 0.0;
    if(prime_loo_cv.input.amt_factor > 1)
    {
        int index = slice_curr;
        double * x_a = new double [prime_loo_cv.input.amt_factor];
        for(int i = 0; i < prime_loo_cv.input.amt_factor; i++)
            x_a[i] = slice_vals[i];
        x_a[index] = x;
        double y = prime_loo_cv.estimate_get(x_a, prime_loo_cv.c_weight + 1, prime_loo_cv.c_weight[0]);
        delete [] x_a;
        return y;
    }
    double x_a = x;
    return prime_loo_cv.estimate_get(&x_a, prime_loo_cv.c_weight + 1, prime_loo_cv.c_weight[0]);
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
        y = get_y2(x);
        if(y > max_y) max_y = y;
        if(y < min_y) min_y = y;
        y = get_y3(x);
        if(y > max_y) max_y = y;
        if(y < min_y) min_y = y;
        x = x_min + dx * (float)i;
    }

  //  if(prime_loo_cv.input.amt_factor == 1 && ui->Graph2dWidget->figures())
    if(prime_loo_cv.input.amt_factor == 1)
    {
        for(int i = 0; i < prime_loo_cv.input.amt_exp; i++)
        {
            if(prime_loo_cv.y[i] > max_y) max_y = prime_loo_cv.y[i];
            if(prime_loo_cv.y[i] < min_y) min_y = prime_loo_cv.y[i];
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
    prime_loo_cv.clear_memory();
    prime_reg.clear_memory();
    prime_weight_loo_cv.clear_memory();

    // присвоение входных значений параметрам
    prime_loo_cv.input.amt_factor = ui->spinBox->value();          // количество факторов
    prime_loo_cv.input.amt_exp = ui->spinBox_2->value();           // количество экспериментов
    prime_loo_cv.input.num_func = stack_func
            [prime_loo_cv.input.amt_factor]
            [ui->listWidget_2->currentRow()]
            .second;                                               // номер функции
    prime_loo_cv.input.noise = ui->doubleSpinBox_10->value();      // процент зашумленности
    prime_loo_cv.input.gamma = ui->doubleSpinBox_11->value();      // параметр регуляризации
    prime_loo_cv.input.auto_gamma =ui->checkBox_9->isChecked();    // наличие автоподбора гаммы
    prime_loo_cv.input.left_border = ui->doubleSpinBox_6->value(); // левая граница поиска гаммы
    prime_loo_cv.input.right_border = ui->doubleSpinBox_7->value();// правая граница поиска гаммы

    // Если границы заданы неверно, меняем местами
    if(prime_loo_cv.input.left_border > prime_loo_cv.input.right_border)
        swap(prime_loo_cv.input.left_border, prime_loo_cv.input.right_border);

    // Выбор вида ядра
    if (ui->radioButton->isChecked())    prime_loo_cv.curr_kernel = 1; // линейное ядро
    if (ui->radioButton_2->isChecked())  prime_loo_cv.curr_kernel = 2; // RBF - ядро
    if (ui->radioButton_3->isChecked())  prime_loo_cv.curr_kernel = 3; // полиномиальное ядро

    // Заполнение параметров линейного ядра
    prime_loo_cv.kernel1.coef = ui->doubleSpinBox_2->value();
    // Заполнение параметров RBF - ядра
    prime_loo_cv.kernel2.sigma =  ui->doubleSpinBox->value();
    // Заполнение параметров полиномиального ядра
    prime_loo_cv.kernel3.degree = ui->doubleSpinBox_5->value();
    prime_loo_cv.kernel3.coef = ui->doubleSpinBox_3->value();
    prime_loo_cv.kernel3.a = ui->doubleSpinBox_4->value();

    // Наличие автоподбора параметров ядра
    if (prime_loo_cv.curr_kernel == 1) prime_loo_cv.input.auto_param = false;
    if (prime_loo_cv.curr_kernel == 2) prime_loo_cv.input.auto_param = ui->checkBox->isChecked();
    if (prime_loo_cv.curr_kernel == 3) prime_loo_cv.input.auto_param = ui->checkBox_2->isChecked();

    // Границы варьирования факторов
    prime_loo_cv.input.var_level.clear();
    for(int i = 0; i < prime_loo_cv.input.amt_factor; i++)
    {
        if(variable[i].first < variable[i].second)
            prime_loo_cv.input.var_level.push_back(pair<double, double>(variable[i].first, variable[i].second));
        else
            prime_loo_cv.input.var_level.push_back(pair<double, double>(variable[i].second, variable[i].first));
    }

    // Будем обновлять выборку или нет
    if(ui->checkBox_8->isChecked())
        srand(time(NULL));
    else
        srand(0x12344321U); //магическое число
    //  три константы, скорее всего одна из них отвечает за генератор
    // 0x12344321U, 0xEE11DD22U, 0xFEDCBA98

    /* -------- Активация модификаций -----------*/

    prime_loo_cv.input.new_gen = 0;                // 0 - отсутствие/ 1 - наличие генерации
    prime_weight_loo_cv.input.weight_cv_exist = 0; // 0 - отсутствие/ 1 - наличие взвешенной процедуры скользящего контроля
    prime_reg.input.cr_reg_exist = 0;              // 0 - отсутствие/ 1 - наличие критерия регулярности
    prime_loo_cv.input.auto_param_c = 0;           // 0 - отсутствие/ 1 - наличие подбора константы с
    prime_loo_cv.input.weighted = 0;               // 0 - отсутствие/ 1 - наличие взвешенного алгоритма робастности
    prime_loo_cv.input.outliers = 0;               // 0 - отсутствие/ 1 - наличие загрязнений (выбросов)/ 2 - один выброс
    prime_loo_cv.input.c_setting = 1;              // 0 - отсутствие/ 1 - наличие подборки границ взвешенного метода
    prime_loo_cv.input.alpha = 0.3;                // определяет процент выбросов (значение от 0 до 1) при input.outliers == 1
    prime_loo_cv.input.c_1 = 2.5;                  // стандартные значения граничных констант с1 = 2.5
    prime_loo_cv.input.c_2 = 3;                    // с2 = 3 (предположение о нормальном распределении помехи)

    /* -------- Блок вычислений -----------------*/

    // копирование данных из класса loo_cv в оставшиеся классы критериев
    prime_reg.copy_input(dynamic_cast<realiz *>(&prime_loo_cv));
    prime_weight_loo_cv.copy_input(dynamic_cast<realiz *>(&prime_loo_cv));
    prime_loo_cv.init();
    prime_reg.copy_init(dynamic_cast<realiz *>(&prime_loo_cv));
    prime_weight_loo_cv.copy_init(dynamic_cast<realiz *>(&prime_loo_cv));

    // Переменные значений критериев соответствующих классов
    double fast_loo_val,cr_regularity_val, weight_loo_val;
    double loo_MSE_val, reg_MSE_val, weight_loo_MSE_val;

    // Начало работы приложения
    QTime time_start, time_stop;
    try
    {
        time_start = QTime::currentTime();

    /* 1 - критерий регулярности для стандартного LS SVM */

        if ( prime_reg.input.cr_reg_exist == 1)
            {
              prime_reg.input.num_interval = 25; // количество интервалов разбиения тестовой выборки
              // выборка делится на количество интервалов, в каждом из которых последовательно исключается наблюдение

              // получение параметров с использованием критерия регулярности
              prime_reg.sort();
              prime_reg.auto_fit();
              cr_regularity_val = prime_reg.cr_regularity();

              // получение оценок по текущим параметрам, определенным с использованием критерия регулярности
              prime_reg.ls_svm_slae();
              prime_reg.solve_slau();
              prime_reg.calc_y_estimate();
              reg_MSE_val = prime_reg.calc_MSE();

            /*cout << "parametrs regular cr." << endl;
              cout << "gamma_r" << endl;
              cout << prime_reg.input.gamma << endl;
              cout << "sigma_r" << endl;
              cout << prime_reg.kernel2.sigma << endl;
              cout << "regular" << endl;
              cout <<  cr_regularity_val << " ";
              cout << "mse_reg" << endl;
              cout << reg_MSE_val << endl;*/
            }

   /* 2 - критерий LOO CV для стандартного LS SVM */

        // получение параметров с использованием кросс-проверки
        prime_loo_cv.auto_fit(); // параметры определяются выбранным способом
        fast_loo_val = prime_loo_cv.fast_loo_cv(); // подсчет значения критерия fst loo cv по полученным параметрам

        cout << "parametrs Loo fast" << endl;
        cout << "gamma" << endl;
        cout << prime_loo_cv.input.gamma << endl;
        cout << "sigma2" << endl;
        cout << prime_loo_cv.kernel2.sigma << endl;

        // получение оценок по текущим параметрам , определенным с использованием критерия Loo_fast
        prime_loo_cv.ls_svm_slae();                // построение системы
        prime_loo_cv.solve_slau();                 // решение системы
        prime_loo_cv.calc_y_estimate();            // получение оценок отклика
        loo_MSE_val = prime_loo_cv.calc_MSE();     // пдсчет MSE

        cout << "loo" << endl;
        cout <<  fast_loo_val << " ";
        cout << "mse_loo" << endl;
        cout << loo_MSE_val << endl;

        /*переприсваивание вектора параметров, полученных из решения системы,
        для отрисовки в связи с перезаписью вектора для взвешенного метода */
        for(int i = 0; i <  prime_loo_cv.input.amt_exp+1; i++)
            prime_loo_cv.c_normal[i] = prime_loo_cv.c[i];

  /* 3 - критерий LOO CV для взвешенного LS SVM */

       prime_loo_cv.input.weighted = 1;
       if (prime_loo_cv.input.weighted == 1)
         {
             // проверка необходимости тюнинга параметра "с" взвешенной процедуры
             if (prime_loo_cv.input.auto_param_c == 1) prime_loo_cv.robast_constant_param_get();
             prime_loo_cv.ls_svm_weight();                    // решение системы с помщью взвешенного ls svm
             prime_loo_cv.calc_y_estimate();                  // получения оценок отклика
             loo_MSE_val = prime_loo_cv.calc_MSE();           // подсчет значения MSE для взвешенного ls svm
             fast_loo_val = prime_loo_cv.fast_loo_cv();       // подсчет значения fast loo cv для взвешенного ls svm

             cout << "loo for weight" << endl;
             cout <<  fast_loo_val << endl;
             cout << "mse_loo  for weight" << endl;
             cout << loo_MSE_val<< endl;

            //переприсваивание вектора параметров для отрисовки
            for(int i = 0; i <  prime_loo_cv.input.amt_exp+1; i++)
                prime_loo_cv.c_weight[i] = prime_loo_cv.c[i];
         }

         // переключение флагов
         prime_weight_loo_cv.input.weighted = 0;
         prime_weight_loo_cv.input.weight_cv_exist = 1;
         prime_weight_loo_cv.input.auto_param_c = 1;

   /* 4 - взвешенный (робастный)критерий LOO CV (loo RCV) для стандартного и взвешенного LS SVM */

         if (prime_weight_loo_cv.input.weight_cv_exist == 1)

           {
             // loo RCV для стандартного  LS SVM
             prime_weight_loo_cv.auto_fit();            // осуществляется подбор параметров с использованием loo RCV
             prime_weight_loo_cv.ls_svm_slae();         // строится и решается слау LS SVM
             prime_weight_loo_cv.solve_slau();

             // loo RCV для взвешенного LS SVM (процедура взвешенного алгоритма следует за шагом стандартного алгоритма)
             prime_weight_loo_cv.input.weighted = 1;    // включение взвешенного LS SVM

             // проверка необходимости тюнинга параметра "с" взвешенной процедуры
             if (prime_weight_loo_cv.input.auto_param_c == 1) prime_weight_loo_cv.robast_constant_param_get();

             // Итерация взешенного LS SVM с параметрами, полученными при использовании loo RCV
             if (prime_weight_loo_cv.input.weighted == 1) prime_weight_loo_cv.ls_svm_weight();

             // получение оценок и вычисление критериев
             prime_weight_loo_cv.calc_y_estimate();
             weight_loo_val = prime_weight_loo_cv.weight_fast_loo_cv();
             weight_loo_MSE_val = prime_weight_loo_cv.calc_MSE();

             cout << "parametrs с for RCV ----- posle" << endl;
             cout <<  prime_weight_loo_cv.input.c_setting << endl;
             cout << "loo for weight and RCV" << endl;
             cout <<  weight_loo_val << endl;
             cout << "mse_loo for weight and RCV" << endl;
             cout << weight_loo_MSE_val << endl;

             cout << "gamma_robust CV" << endl;
             cout << prime_weight_loo_cv.input.gamma << endl;
             cout << "sigma_robust CV" << endl;
             cout << prime_weight_loo_cv.kernel2.sigma << endl;

             //переприсваивание вектора параметров для отрисовки
             for(int i = 0; i <  prime_weight_loo_cv.input.amt_exp+1; i++)
                 prime_weight_loo_cv.c_WCV[i] = prime_weight_loo_cv.c[i];
             }

       time_stop = QTime::currentTime();
    }


    catch(string error_descr)
    {
        QMessageBox::critical(this, trUtf8("Ошибка"), trUtf8(error_descr.c_str()), QMessageBox::Ok, QMessageBox::Ok);
        return;
    }

    // Вывод результата MSE в форму
    if (ui->checkBox_10->isChecked()) ui->lineEdit->setText(QString::number(prime_loo_cv.calc_MSE(), 'e', 7));
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
        // Зададим строки таблицы
        ui->tableWidget->setRowCount(prime_loo_cv.input.amt_exp);
        // Заполним таблицу
        for(int i = 0; i < prime_loo_cv.input.amt_exp; i++)
        {
            QTableWidgetItem * newItem;
            newItem = new QTableWidgetItem(QString::number(prime_loo_cv.y[i]));
            ui->tableWidget->setItem(i, 0, newItem);
            newItem = new QTableWidgetItem(QString::number(prime_loo_cv.y_new[i]));
            ui->tableWidget->setItem(i, 1, newItem);
            newItem = new QTableWidgetItem(QString::number(prime_loo_cv.y2_new[i])); // Тут новое название игрека
            ui->tableWidget->setItem(i, 2, newItem);
        }

        for(int i = 0; i <  prime_loo_cv.input.amt_exp; i++)
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
    if(prime_loo_cv.input.amt_factor == 1)
    {
        for(int i = 0; i < prime_loo_cv.input.amt_exp; i++)
            ui->Graph2dWidget->data.push_back(QGraph2dPoint(prime_loo_cv.x[i][0], prime_loo_cv.y[i]));
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

        adjust_zoom(prime_loo_cv.input.var_level[0].first, prime_loo_cv.input.var_level[0].second);
    }
    else
    {
        // Занесение факторов в сечение
        QStringList factor_text3;
        for(int i = 0; i < prime_loo_cv.input.amt_factor; i++)
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
        ui->doubleSpinBox->setValue(prime_loo_cv.kernel2.sigma_select[6]);
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

    if(prime_loo_cv.input.var_level.size() && slice_curr >= 0)
        adjust_zoom(prime_loo_cv.input.var_level[slice_curr].first, prime_loo_cv.input.var_level[slice_curr].second);

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
    if(prime_loo_cv.input.var_level.size())
    {
        ui->doubleSpinBox_12->setMinimum(prime_loo_cv.input.var_level[index].first);
        ui->doubleSpinBox_12->setMaximum(prime_loo_cv.input.var_level[index].second);
    }

    if(prime_loo_cv.input.var_level.size() && slice_curr >= 0)
        adjust_zoom(prime_loo_cv.input.var_level[slice_curr].first, prime_loo_cv.input.var_level[slice_curr].second);

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

    if(prime_loo_cv.input.var_level.size() && slice_curr >= 0)
        adjust_zoom(prime_loo_cv.input.var_level[slice_curr].first, prime_loo_cv.input.var_level[slice_curr].second);

    ui->Graph2dWidget->setDraw(draw);
}

// Включение отображения отклика
void MainWindow::on_checkBox_12_clicked(bool checked)
{
    if(checked)
        ui->Graph2dWidget->setFigures(true);
    else
        ui->Graph2dWidget->setFigures(false);
    adjust_zoom(prime_loo_cv.input.var_level[0].first, prime_loo_cv.input.var_level[0].second);
}

void MainWindow::on_checkBox_13_clicked(bool checked)
{
    if(checked)
        ui->Graph2dWidget->setLines(0, true);
    else
        ui->Graph2dWidget->setLines(0, false);
    adjust_zoom(prime_loo_cv.input.var_level[0].first, prime_loo_cv.input.var_level[0].second);
}

void MainWindow::on_checkBox_14_clicked(bool checked)
{
    if(checked)
        ui->Graph2dWidget->setLines(1, true);
    else
        ui->Graph2dWidget->setLines(1, false);
    adjust_zoom(prime_loo_cv.input.var_level[0].first, prime_loo_cv.input.var_level[0].second);

}

void MainWindow::on_checkBox_15_clicked(bool checked)
{
    if(checked)
        ui->Graph2dWidget->setLines(2, true);
    else
        ui->Graph2dWidget->setLines(2, false);
    adjust_zoom(prime_loo_cv.input.var_level[0].first, prime_loo_cv.input.var_level[0].second);
}
