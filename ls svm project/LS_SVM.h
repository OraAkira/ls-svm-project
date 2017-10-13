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

double dot_product(double *, double *, int n); // scalar product ---- скалярное произведение

// input data ---- входные данные
class data_in
{
 public:
	int num_func;		  // number of the selected function ---- номер выбранной функции
	int amt_factor;		  // number of factors in the experiment ---- количество факторов эксперимента
	int amt_exp;		  // number of experiments ---- количество экспериментов
	double noise;		  // percent of noisiness ---- процент зашумленности 
	double gamma;		  // regularization parameter ---- параметр регуляризации
	bool auto_gamma;	  // false: presence; true: auto-selection of gamma ---- флаг 0 - наличия, 1- отсутствия автоподбора гаммы
	bool auto_param;	  // false: presence; true: lack of auto-selection of kernel parameters ---- флаг 0 - наличия, 1- отсутствия автоподбора параметров ядра
	double left_border;   // left border of gamma search ---- левая граница поиска гаммы
	double right_border;  // rightmost gamma search boundary ---- правая граница поиска гаммы
    int kernel_type;      // type of kernel: 1. the lined core, 2. polynomial, 3. RBF-core ---- тип ядра: 1 - линеное ядро, 2 - полиномиальное, 3 - RBF-ядро
 	vector < pair<double,double> >  var_level; // restrictions on the levels of variation of factors ---- ограничения на уровни варьирования факторов
};

// linear core ---- класс "линейное ядро"
class kernel_linear 
{	
 public:
	int m; // number of factors ---- количество факторов
	double coef;
    double kernel_func(double *, double *);
};

// RBF-core ---- класrnel_RBFс "RBF -ядро"
class kernel_RBF
{	
 public:
	int m;	// number of factors ---- количество факторов
	double sigma;
	kernel_RBF();
	double sigma_select[10];
    double kernel_func(double *, double *);
};

// polynomial kernel ---- класс "полиномиальное ядро"
class kernel_polynomial
{	
 public:
	int m;	// number of factors ---- количество факторов
	double a;
	double coef;
	double degree;
    double kernel_func(double *, double *);
};

// main class ---- основной класс
class realiz
{
 private:

    double golden(); // golden section method ---- метод золотого сечения

    // SLAE the method A*c=b ---- СЛАУ по методу LS SVM A*c=b
    double **A;

    /* Functions used to refer to a matrix ---- Функции, используемые для обращения матрицы */
    double **InversionMatrix(double **, int n, int m);				   // matrix inversion function ---- функция обращения матрицы
    bool SLAU_solution_Gauss(double **, double **, double **, int n);  // solution slack AX = I ---- решение слау АХ=I
    bool trianglematrix1(double **, double **, int n);				   // reduction of the matrix to the upper triangular form ---- приведение матрицы к верхнему треугольному виду
    void transf1(double **, double **, int i, int n, double &A_d);	   // permutation of strings ---- перестановка строк

    double function_choice(int i);                     // function type selection --- выбор типа функции
    double rav_num_gen(double left, double right);     // uniformly distributed random number generator ---- генератор равномерно-распределенных случайных чисел
    double normal_num_gen(double m, double sigma22);   // normal-distributed random number generator ---- генератор нормально-распределенных случайных чисел
    void gamma_get();								   // auto-parameter function gamma ---- функция автоподбора параметра gamma
    void kernel_param_get();						   // automatic kernel parameter selection function ---- функция автоподбора параметров ядра
    double kernel_func_choice(double *, double *);     // Select the appropriate kernel function ---- выбор соответствующей функции ядра

 public:
    int curr_kernel;
    data_in input;             // object of the input data class ---- объект класса входных данных
    kernel_linear kernel1;     // object of the kernel_linear class ---- объект класса линейного ядра
    kernel_RBF kernel2;        // object of the kernel_RBF class ---- объект класса RBF ядра
    kernel_polynomial kernel3; // object of the kernel_polynomial class ---- объект класса полиномиального ядра

    double **x;        // reference to an observation matrix ---- ссылка на матрицу наблюдений
    double *y;         // reference to an array of responses ---- ссылка на массив откликов
    double *y_new;     // reference to an array of feedback estimates received with training ---- ссылка на массив оценок отклика полученных с обучением
    double *y2_new;    // reference to an array of response estimates ---- ссылка на массив оценок отклика
    double *y_ist;     // reference to an array of noiseless responses ---- ссылка на массив незашумленных откликов

    // SLAE the method LS-SVM A*c=b ---- СЛАУ по методу LS SVM A*c=b
    double *c;

    realiz();  // class constructor ---- конструтор класса
    ~realiz(); // class destructor ---- деструктор класса

    void init();                                       // generation of an observation matrix and obtaining a response ---- генерация матрицы наблюдений и получение отклика
    void auto_fit();								   // function of selection and auto-selection of parameters ---- функция выбора и автоподбора параметров
    void clear_memory();							   // memory clear function ---- функция очистки памяти
    void ls_svm_slae();								   // function of composing and solving SLAU by LS LSM method ---- функция составления и решения СЛАУ по методу LS SVM
    void calc_y_estimate();                            // function for obtaining a vector of response estimates by current parameters ---- функция получения вектора оценок отклика по текущим параметрам
    double fast_loo_cv();							   // function of calculating the value of the fast algorithm LOO_SV ---- функция подсчета значения по быстрому алгоритму LOO_SV
    double calc_MSE();                                 // function calculation of MSE ---- функция вычисление MSE
    double estimate_get(double *, double *, double b); // evaluation function ---- функция получения оценки
};



