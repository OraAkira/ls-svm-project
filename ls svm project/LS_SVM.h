#pragma once
#include <cstdio>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <utility>
#include <time.h>
#include <map>

using namespace std;

bool rule_sorting (double, double);                     // правило для сортировки вектора
double dot_product(double *, double *, int m);          // скалярное произведение
void solve_gauss(double **, double *, double *, int n); // Решение системы по методу Гаусса

/* входные данные  */
class data_in
{
 public:
    double noise;		  // процент зашумленности
    double gamma;		  // параметр регуляризации
    double left_border;   // левая граница поиска гаммы
    double right_border;  // правая граница поиска гаммы
    int cr_reg_exist;     // наличие/отсутствие ОКР
    int num_func;		  // номер выбранной функции
    int amt_factor;		  // количество факторов эксперимента
    int amt_exp;		  // количество экспериментов
    int num_interval;     // количество интервалов разбиения выборки
    int kernel_type;      // флаг задает тип ядра: 1 - линеное ядро, 2 - полиномиальное, 3 - RBF-ядро
    bool auto_gamma;	  // флаг 0 - наличия, 1 - отсутствия автоподбора гаммы
    bool auto_param;	  // флаг 0 - наличия, 1 - отсутствия автоподбора параметров ядра
    int criterion_type;   // флаг задает тип критерия: 0 - Loo CV, 1 - weight Loo CV, 2 - ОКР
    int new_gen;          // флаг обновление генератора: 0 - все фиксируется, 1 - обновляем шум, 2 - обновляем все

 	vector < pair<double,double> >  var_level; // ограничения на уровни варьирования факторов

    /* Взвешенный LS SVM  */
    int weighted;         // задает отсутствите(0) / наличие(1) взвешенного LS SVM
    int outliers;         // флаг 0 - отсутствие / 1 - наличие загрязнений (выбросов) / 2 - один выброс на конце отрезка
    int weight_cv_exist;  // наличие/отсутствие взвешенного loo CV
    double alpha;         // определяет процент выбросов (значение от 0 до 1)
    double c_1, c_2;      // граничные константы обычно выбираются в качестве с1=2,5 и с2=3
    double c_setting;     // переменная, отвечающая за автоподбор граничных констант
    bool auto_param_c;	  // флаг 0 - отсутствия, 1 - наличия автоподбора параметра с робастного метода
};

/* класс "линейное ядро"  */
class kernel_linear 
{	
 public:
    int m;                                 // количество факторов
    double coef;                           // свободный коэффициент
    double kernel_func(double *, double *);// функция вычисления значения для линейного  ядра
};

/* класс "RBF -ядро"  */
class kernel_RBF
{	
 public:
    int m;	                                // количество факторов
    double sigma;                           // параметр RBF - ядра
    kernel_RBF();                           // конструктор для генерации набора допустимых значений параметра sigma (сетка)
    double sigma_select[100];               // набор значений параметра сигма
    double kernel_func(double *, double *); // функция вычисления значения для RBF - ядра
};

/* класс "полиномиальное ядро"  */
class kernel_polynomial
{	
 public:
    int m;	                                // количество факторов
    double a;                               // коэффициент 1 (при факторе)
    double coef;                            // коэффициент 2 (свободный)
    double degree;                          // степень полинома
    double kernel_func(double *, double *); // функция вычисления значения для полимиального ядра
};

/* абстрактный класс  */
class realiz
{
    protected:

       double golden(); // метод золотого сечения

       /* СЛАУ по методу LS SVM A*c=b */
       double **A;
       double *b;

       /* Функции, используемые для обращения матрицы */
       double **InversionMatrix(double **, int n, int m);				   // функция обращения матрицы
       bool SLAU_solution_Gauss(double **, double **, double **, int n);   // решение слау АХ=I
       bool trianglematrix1(double **, double **, int n);				   // приведение матрицы к верхнему треугольному виду
       void transf1(double **, double **, int i, int n, double &A_d);	   // перестановка строк

       virtual double criterion_choice() = 0;   // выбор типа критерия качества работы алгоритма

       /* Генераторы случайных чисел */
       double rav_num_gen(double left, double right);         // генератор равномерно-распределенных случайных чисел
       double normal_num_gen(double m, double sigma22);       // генератор нормально-распределенных случайных чисел
       double cauchy_gen(double m, double sig);               // генератор распределенных по Коши случайных чисел

       double function_choice(int i);                         // выбор типа функции
       double kernel_func_choice(double *, double *);         // выбор соответствующей функции ядра
       double euclidean_distance(double *, double *, int m);  // Евклидово расстояние
       void gamma_get();								      // функция автоподбора параметра gamma
       void kernel_param_get();						          // функция автоподбора параметров ядра


      public:

          /* СЛАУ по методу LS SVM A*c=b */
          double *c;        // переменная для хранения параметров по методу LS SVM при построении системы
          double *c_normal; // вектор параметров для отрисовки результатов стандартного алгоритма LS SVM
          double *c_weight; // вектор параметров для отрисовки результатов взвешенного алгоритма LS SVM с стандартным loo CV
          double *c_WCV;    // вектор параметров для отрисовки результатов взвешенного алгоритма LS SVM с взвешенным loo CV (loo WCV)

          int curr_kernel;
          data_in input;             // объект класса входных данных
          kernel_linear kernel1;     // объект класса линейного ядра
          kernel_RBF kernel2;        // объект класса RBF-ядра
          kernel_polynomial kernel3; // объект класса полиномиального ядра

          double **x;         // ссылка на матрицу наблюдений
          double **x_tmp;     // вспомогательная матрица
          double *y_ist;      // ссылка на массив незашумленных откликов
          double *y;          // ссылка на массив зашумленных откликов
          double *y_new;      // ссылка на массив оценок отклика полученных с обучением
          double *y2_new;     // ссылка на массив оценок отклика
          double *y_tmp;      // вспомогательный вектор

          /* Взвешенный LS SVM  */
          double **A_weight;                        // СЛАУ по методу weight LS SVM A*c=b
          double *weight;                           // ссылка на массив весовых коэффициентов v_k
          double *ek;                               // значения из алгоритма ek=альфа/гамма
          double *rk;                               // оcтатки для взвешенного loo_cv
          double *y_weight;                         // ссылка на массив оценок отклика для взвешенного ls svm
          double mediana(double *, int);            // функция вычисления медианы
          double mad_calc(double *, int);           // функция вычисления медианы абсолютного отклонения
          void ls_svm_weight();                     // функция робастного ls-svm
          void weight_calc(double *, int, double);  // функция подсчета весовых коэффициентов
          void robast_constant_param_get();         // функция автоподбора параметра "с" робастного метода

          realiz();  // конструтор класса
          ~realiz(); // деструктор класса

          /* функции копирования данных для использования нескольких критериев */
          void copy_input(class realiz *);  // копирование входных данных из формы
          void copy_init(class realiz *);   // копирование генерируемых данных и выделение памяти

          void init();                                       // генерация матрицы наблюдений и получение отклика
          void auto_fit();								     // функция выбора и автоподбора параметров
          void clear_memory();							     // функция очистки памяти
          void solve_slau();                                 // функция вызова решателя СЛАУ
          void ls_svm_slae();                                // функция составления СЛАУ по методу LS SVM
          void calc_y_estimate();                            // функция получения вектора оценок отклика по текущим параметрам
          double calc_MSE();                                 // функция вычисление MSE
          double estimate_get(double *, double *, double b); // функция получения оценки
};

 /* класс loo CV критерия качества работы алгоритма */
class loo_cv_realiz: public realiz
{
 protected:
   double criterion_choice();   // выбор типа критерия качества работы алгоритма
 public:   
   double fast_loo_cv();        // функция подсчета значения по быстрому алгоритму LOO_SV
};

 /* класс взвешенного loo CV критерия качества работы алгоритма */
class weight_loo_cv_realiz: public realiz
{
 protected:
   double criterion_choice();   // выбор типа критерия качества работы алгоритма
 public:
   double weight_fast_loo_cv(); // функция подсчета значения по быстрому взвешенному алгоритму LOO_SV
};

 /* класс ОКР (обобщенный критерий регулярности) критерия качества работы алгоритма */
class regularity_realiz: public realiz
{
 protected:
   double criterion_choice(); // выбор типа критерия качества работы алгоритма
 public:
   void sort();                            // сортировка выборки c помощью контейнера
   void exclusion(int l, int n, int s);    // исключение наблюдений не участвующих в обучении
   double cr_regularity();                 // функция подсчета значения по ОКР
};




