#pragma once
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <utility>
#include <time.h>

using namespace std;

double dot_product(double *, double *, int n); // скалярное произведение

//входные данные
class data_in
{
 public:
	int num_func;		  // номер выбранной функции
	int amt_factor;		  // количество факторов эксперимента
	int amt_exp;		  // количество экспериментов
	double noise;		  // процент зашумленности 
	double gamma;		  // параметр регуляризации
	bool auto_gamma;	  // флаг 0 - наличия, 1- отсутствия автоподбора гаммы
	bool auto_param;	  // флаг 0 - наличия, 1- отсутствия автоподбора параметров ядра
	double left_border;   // левая граница поиска гаммы
	double right_border;  // правая граница поиска гаммы
    int kernel_type;      // тип ядра: 1 - линеное ядро, 2 - полиномиальное, 3 - RBF-ядро
 	vector < pair<double,double> >  var_level; // ограничения на уровни варьирования факторов
};

// класс "линейное ядро"
class kernel_linear 
{	
 public:
	int m; // количество факторов
	double coef;
    double kernel_func(double *, double *);
};

// класrnel_RBFс "RBF -ядро"
class kernel_RBF
{	
 public:
	int m;	// количество факторов
	double sigma;
	kernel_RBF();
	double sigma_select[10];
    double kernel_func(double *, double *);
};

// класс "полиномиальное ядро"
class kernel_polynomial
{	
 public:
	int m;	// количество факторов
	double a;
	double coef;
	double degree;
    double kernel_func(double *, double *);
};

// основной класс
class realiz
{
 private:

    double golden(); // метод золотого сечения

    //СЛАУ по методу LS SVM A*c=b
    double **A;

    /* Функции, используемые для обращения матрицы*/
    double **InversionMatrix(double **, int n, int m);				   // функция обращения матрицы
    bool SLAU_solution_Gauss(double **, double **, double **, int n);  // решение слау АХ=I
    bool trianglematrix1(double **, double **, int n);				   // приведение матрицы к верхнему треугольному виду
    void transf1(double **, double **, int i, int n, double &A_d);	   // перестановка строк

    double function_choice(int i);                     // выбор типа функции
    double rav_num_gen(double left, double right);     // генератор равномерно-распределенных случайных чисел
    double normal_num_gen(double m, double sigma22);   // генератор нормально-распределенных случайных чисел
    void gamma_get();								   // функция автоподбора параметра gamma
    void kernel_param_get();						   // функция автоподбора параметров ядра
    double kernel_func_choice(double *, double *);     // выбор соответствующей функции ядра

 public:
    int curr_kernel;
    data_in input;             // объект класса входных данных
    kernel_linear kernel1;     // объект класса линейного ядра
    kernel_RBF kernel2;        // объект класса RBF ядра
    kernel_polynomial kernel3; // объект класса полиномиального ядра

    double **x;    // ссылка на матрицу наблюдений
    double *y;     // ссылка на массив откликов
    double *y_new; // ссылка на массив оценок отклика полученных с обучением
    double *y2_new; // ссылка на массив оценок отклика
    double *y_ist; // ссылка на массив незашумленных откликов

    //СЛАУ по методу LS SVM A*c=b
    double *c;

    realiz();  // конструтор класса
    ~realiz(); // деструктор класса

    void init();                                       // генерация матрицы наблюдений и получение отклика
    void auto_fit();								   // функция выбора и автоподбора параметров
    void clear_memory();							   // функция очистки памяти
    void ls_svm_slae();								   // функция составления и решения СЛАУ по методу LS SVM
    void calc_y_estimate();                            // функция получения вектора оценок отклика по текущим параметрам
    double fast_loo_cv();							   // функция подсчета значения по быстрому алгоритму LOO_SV
    double calc_MSE();                                 // функция вычисление MSE
    double estimate_get(double *, double *, double b); // функция получения оценки
};



