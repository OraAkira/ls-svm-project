#include "LS_SVM.h"

// конструктор
realiz::realiz()
{
    rk =NULL;
    ek = NULL;
    weight = NULL;
	y_ist = NULL;
    y_tmp = NULL;
	y_new = NULL;
    y2_new = NULL;
    y_weight = NULL;
	y = NULL;
	x = NULL;
    x_tmp = NULL;
	A = NULL;
    b = NULL;
	c = NULL;
    c_normal = NULL;
    c_weight = NULL;
    c_WCV = NULL;
 }

// деструктор
realiz::~realiz()
{
	clear_memory();
}

// правило для сортировки вектора
bool rule_sorting (double a,double b) { return (a<b);}

// выбор типа критерия качества работы алгоритма
double loo_cv_realiz::criterion_choice()
{
    return fast_loo_cv();
}
double regularity_realiz::criterion_choice()
{
    return cr_regularity();
}

double  weight_loo_cv_realiz::criterion_choice()
{
    return weight_fast_loo_cv();
}

// функция генерации равномерно-распределенных чисел
double realiz::rav_num_gen(double left, double right)
{
//    if ((input.new_gen == 0)||(input.new_gen == 1)) srand(0x12344321U);
	double alpha1 = rand() / (double)RAND_MAX;
	double rand_num =  left + alpha1 * (right - left);

	return rand_num;
}
/*
// функция генерации случайных чисел распределенных по Коши (ТЯЖЕЛЫЕ ХВОСТЫ)
// вариант #1
double realiz::cauchy_gen(double m, double sig)
{
    double x, z;
    //x = rand() / (double)RAND_MAX;
    x = rav_num_gen(0,1);
    //cout
    z = 0.5 + atan((x - m)/sig)/M_PI;
    //return (m + sig*(tan(M_PI*(z - 0.5))));
    return z;
}

// вариант #2
double realiz::cauchy_gen(double x0, double gamma) {
    double x, y;
    do {
        x = rav_num_gen(-1,1);
        y = rav_num_gen(-1,1);
    } while (x * x + y * y > 1.0 || y == 0.0);
    return x0 + gamma * x / y;
}*/

// вариант #3
double realiz::cauchy_gen(double m, double sig)
{
    double x;
   // x  = rav_num_gen(0,1);
   x = double(rand())/RAND_MAX;
   // return 0.5 + (1.0/M_PI) * atan2(x - m, sig);

    return m + (tan(M_PI*(x - 0.5)))*sig;

}

// функция генерации нормально-распределенных чисел (СИММЕТРИЧНОЕ ЗАГРЯЗНЕНИЕ)
double realiz::normal_num_gen(double m, double sigma22)
{
    if (input.new_gen == 0) srand(0x12344321U);

    double alpha1 = rand() / (double)RAND_MAX;
    double alpha2 = rand() / (double)RAND_MAX;
    double ksi = sqrt(-2.0 * log(alpha1)) * cos(2.0 * 3.14159 * alpha2);

    return m + sqrt(sigma22) * ksi;
}

// функция вычисления Евклидова расстояния
 double realiz::euclidean_distance(double *vect1, double *vect2, int m)
{
    double d = 0.0;

    for(int i = 0; i < m; i++)
        d += (vect1[i] - vect2[i])*(vect1[i] - vect2[i]);

    return sqrt(d);
}

 // функция вычисления медианы массива данных
 double realiz::mediana(double *array, int arraySize)
 {
     int m;

     vector <double> vecEk(input.amt_exp);  // сортировка
     for(int i = 0; i < input.amt_exp; i++)
         vecEk[i] = fabs(array[i]);
     sort(vecEk.begin(), vecEk.end(), rule_sorting);                      
     // определяем, четное ли число элементов массива
     m = arraySize%2;

     return(m ? vecEk[arraySize/2]: (vecEk[arraySize/2 - 1] + vecEk[arraySize/2])/2); // тернарная условная операция
 }

 // функция вычисления медианы абсолютного отклонения (MAD)
 double realiz::mad_calc(double *array, int arraySize)
 {
     double median = 0.0;
     double *ek_temp;
     ek_temp = new double[arraySize];

     for(int i = 0; i < arraySize; i++)
       ek_temp[i] = array[i];
     median = mediana(ek_temp, arraySize);
     for(int i = 0; i < arraySize; i++)
       ek_temp[i] = ek_temp[i] - median;

     return(mediana(ek_temp, arraySize));
 }

 // функция подсчета весовых коэффициентов v_k
  void realiz::weight_calc(double *array, int arraySize, double s)
 {
    double temp_calc_c1 = 0.0;
    double temp_calc_c1_c2 = 0.0;
    double temp_calc_c2 = 0.0;
    double temp = 0.0;
    for(int i = 0; i < arraySize; i++)
        {   temp = fabs(array[i]/(s*input.c_setting));
            if(temp <= input.c_1) {weight[i] = 1.0; temp_calc_c1++;}
            else {
                    if(temp <= input.c_2) {weight[i] = (input.c_2 - temp)/(input.c_2 - input.c_1); temp_calc_c1_c2++;}
                    else {weight[i] = 0.0001; temp_calc_c2++;}
                 }
        }
 }

 // функция автоподбора оптимального значения параметра "c" для робастного метода
  void realiz::robast_constant_param_get()
  {
      input.c_setting = 0.1;
      double min = 0.0;
      double tmp = 0.0;
      double c_tmp = input.c_setting;
      double max_value = 1.5;
      double step = 0.01;
          ls_svm_weight();
          min = criterion_choice();
          while(input.c_setting <= max_value)
           {
              ls_svm_weight();
              tmp = criterion_choice();
              if (tmp < min)
                  {
                      min = tmp;
                      c_tmp = input.c_setting;
                  }
              input.c_setting += step;
           }
       input.c_setting = c_tmp;
   }

 // функция робастного ls-svm
 void realiz::ls_svm_weight()
 {
     // получение параметров для СЛАУ НЕвзвешенного метода
     input.weighted = 0;
     ls_svm_slae();
     solve_slau();
     // вычисление весов и получение параметров уже для СЛАУ взвешенного метода
     input.weighted = 1;
     double s = 0.0;

     for(int i = 0; i < input.amt_exp; i++)
     {
         ek[i] = c[i+1]/input.gamma;
     }

     // вычисление робастной оценки стандартного отклонения переменных ошибок eк с предположением об их Гауссовом распределении
     s = 1.483 * mad_calc(ek, input.amt_exp);
     weight_calc(ek, input.amt_exp, s);
     ls_svm_slae();
     solve_slau();
 }

// метод золотого сечения
double realiz::golden()
{
	double eps = 0.0001;
    double a = input.left_border, b = input.right_border;
    double q1 = (3.0 - sqrt(5.0)) / 2.0, q2 = (sqrt(5.0) - 1.0) / 2.0;
    double x1 = a + q1 * (b - a), x2 = a + q2 * (b - a);
	input.gamma = x1;
    double f1 = criterion_choice();
	input.gamma = x2;
    double f2 = criterion_choice();

    while(true)
    {
        if(f1 <= f2)
        {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + q1 * (b - a);
			input.gamma = x1;
            f1 = criterion_choice();
        }
        else
        {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + q2 * (b - a);
			input.gamma = x2;
            f2 = criterion_choice();
        }
        if(fabs(a - b) < eps)

            return (a + b) / 2.0;
    }
}

// функция автоподбора параметра gamma с использованием  метода золотого сечения
void realiz::gamma_get() 
{
	input.gamma = golden();
}

// функция автоподбора параметров ядра
void realiz::kernel_param_get()
{ 
    int index = 0;
	double degree;
    double max_degree = 2.8;
    double min =0.0;
    double step = 0.1;
    double tmp, tmp2 = 0.0;

	if (curr_kernel == 2) 
	{	
		kernel2.sigma = kernel2.sigma_select[0];
        min = criterion_choice();
        for (int i = 1; i < 100; i++)
		 {
			kernel2.sigma = kernel2.sigma_select[i];
            tmp = criterion_choice();
			if (tmp < min) 
				{
					min = tmp;
					index = i;
				}
		 }
	    kernel2.sigma = kernel2.sigma_select[index];
	 }

	if (curr_kernel == 3) 
	{	degree = tmp = 1;
		kernel3.degree = tmp;
        min = criterion_choice();
		while (tmp < max_degree)
        {
            tmp = (double)(floor((tmp + step) * 100.0 + 0.5)) / 100.0;
			kernel3.degree = tmp;
            tmp2 = criterion_choice();
			if (tmp2 < min) 
				{
					min = tmp2;
					degree = tmp;
				}
		 }
	    kernel3.degree = degree;
	}
}

// функция выбора и автоподбора параметра
void realiz::auto_fit()
{
    if ((input.auto_gamma == 1) && (input.auto_param == 1))
    {
        if ((input.auto_gamma == 1) && (input.auto_param == 1))
        {
            double old_gamma, old_kernel_param, eps = 1e-6;
            int iter = 0, maxiter = 1000;

            double gamma_check = -1.0, param_check = -1.0;
            int counter = 10, max_counter = 10;
            double eps_loop = 1e-10;

            if (curr_kernel == 2)
            {
                do
                {
                    old_gamma = input.gamma;
                    old_kernel_param = kernel2.sigma;
                    kernel_param_get();
                    gamma_get();
                    iter++;

                    if(fabs((gamma_check - input.gamma) / gamma_check) < eps_loop && fabs((param_check - kernel2.sigma) / param_check)  < eps_loop)
                    {
                        cerr << "Loop detected" << endl;
                        break;
                    }
                    if(counter >= max_counter)
                    {
                        counter = 0;
                        gamma_check = input.gamma;
                        param_check = kernel2.sigma;
                    }
                    else
                    {
                        counter++;
                    }
                }

                while(!(fabs((old_gamma - input.gamma) / old_gamma) < eps && fabs((old_kernel_param - kernel2.sigma) / old_kernel_param) < eps)
                    && (iter < maxiter));
            }

            if (curr_kernel == 3)
            {
                do
                {
                    old_gamma = input.gamma;
                    old_kernel_param = kernel3.degree;
                    kernel_param_get();
                    gamma_get();
                    iter++;

                    if(fabs((gamma_check - input.gamma) / gamma_check) < eps_loop && fabs((param_check - kernel3.degree) / param_check)  < eps_loop)
                    {
                        cerr << "Loop detected" << endl;
                        break;
                    }
                    if(counter >= max_counter)
                    {
                        counter = 0;
                        gamma_check = input.gamma;
                        param_check = kernel3.degree;
                    }
                    else
                    {
                        counter++;
                    }
                }

                while(!(fabs((old_gamma - input.gamma) / old_gamma) < eps && fabs((old_kernel_param - kernel3.degree) / old_kernel_param) < eps)
                    && iter < maxiter);
            }
        }
    }
    if ((input.auto_gamma == 1) && (input.auto_param == 0))
        {
            gamma_get();
        }
    if ((input.auto_gamma == 0) && (input.auto_param == 1))
        {
            kernel_param_get();
        }
}

// функция инициализации
void realiz::init()
{	
    // инициализация
    kernel1.m = kernel2.m = kernel3.m = input.amt_factor;
    c_normal = new double [input.amt_exp+1]; // вектор полученных параметров для ls svm с использованием loo cv
    c_weight = new double [input.amt_exp+1]; // вектор полученных параметров для взвешенного ls svm с использованием loo cv
    c_WCV = new double [input.amt_exp+1];    // вектор полученных параметров для взвешенного ls svm с использованием loo rcv
    weight = new double [input.amt_exp];     // значения весов
    ek = new double [input.amt_exp];         // переменные ошибки для подсчета остатков
    rk = new double [input.amt_exp];         // переменные остатков
    y_weight = new double [input.amt_exp];   // переменные отклика для взвешенного метода
    y_ist = new double [input.amt_exp];      // истинные значения отклика
    y = new double [input.amt_exp];          // зашумленные значения отклика
    y_tmp = new double [input.amt_exp];      // промежуточные значения отклика
    y_new = new double [input.amt_exp];      // оценки отклика при использовании кросс-проверки
    y2_new = new double [input.amt_exp];     // оценки отклика по текущим значениям параметров
    x = new double *[input.amt_exp];         // матрица значений наблюдений
    x_tmp = new double *[input.amt_exp];     // матрица промежуточных значений наблюдений
    for(int i = 0; i < input.amt_exp; i++)
     {x[i] = new  double [input.amt_factor];
      x_tmp[i] = new double [input.amt_exp];}

    // элементы СЛАУ для метода LS_SVM
    c = new double [input.amt_exp+1];
    b = new double [input.amt_exp+1];
    A = new double *[input.amt_exp+1];
    for(int i = 0; i < input.amt_exp + 1; i++)
     A[i] = new  double [input.amt_exp + 1];

    // заполнение матрицы наблюдений
	for(int i = 0; i < input.amt_exp; i++)
     for(int j = 0; j < input.amt_factor; j++)
       x[i][j] = rav_num_gen(input.var_level[j].first, input.var_level[j].second);

	 // вычисление массива откликов
	for(int i = 0; i < input.amt_exp; i++)
        y[i] = y_ist[i] = function_choice(i);
	
    // вычисление параметра масштаба для зашумления отклика
	double *aid;    
	aid = new double [input.amt_exp];
	double y_mean = 0.0; // среднее значение отклика
	double w2 = 0.0;	 // мощность сигнала
	double sigma2;		 // дисперсия
	double e;			 // ошибка 

	for(int i = 0; i < input.amt_exp; i++)
		y_mean += y[i];
	y_mean = y_mean/input.amt_exp;
	for(int i = 0; i < input.amt_exp; i++)
		aid[i] = y[i] - y_mean;
	for(int i = 0; i < input.amt_exp; i++)
		w2 += aid[i] * aid[i];
	w2 = w2 / (input.amt_exp - 1);
    sigma2 = input.noise * w2 / 100;

    /* Наличие выбросов предусматривается лишь в одномерном случае! */

    // Сортируем наблюдения по возрастанию, чтобы определить последнее наблюдение
    int max_index = 0;
    double max_value = x[0][0]; // Временная переменная для сортировки вектора x
    for(int i = 1; i < input.amt_exp; i++)
        if(x[i][0] > max_value)
        {
            max_index = i;
            max_value = x[i][0];
        }

    // Генерация ошибки, зашумление откликов, добавление выбросов (при необходимости)
    double x_rnd = 0.0;

    if (input.outliers == 0)      // 0 - Обычные наблюдения
     {
        for(int i = 0; i < input.amt_exp; i++)
         { e = normal_num_gen(0.0, sigma2);
           y[i] += e;
         }
     }
    else {
        if (input.outliers == 1)  // 1 - Выбросы распределенные по Коши (тяжелые хвосты)
               for(int i = 0; i < input.amt_exp; i++)
                {   x_rnd = double(rand())/RAND_MAX;
                    if (x_rnd <= input.alpha) e = cauchy_gen(0.0, sqrt(sigma2));
                    else e = normal_num_gen(0.0, sigma2);
                    y[i] += e;
                }
        else for(int i = 0; i < input.amt_exp; i++)                   // 2 - Добавляем выброс к значению отклика на конце отрезка
             {
                e = normal_num_gen(0.0, sigma2);
                if (i == max_index) y[i] += 0.9; // Размер выброса 0.9
                else y[i] += e;
             }

         }
	delete [] aid;
}

// функция копирования входных данных из формы
void realiz::copy_input(class realiz *r)
{
    input.alpha = r-> input.alpha;
    input.outliers = r-> input.outliers;
    input.c_setting = r-> input.c_setting;
    input.c_1 = r -> input.c_1;
    input.c_2 = r -> input.c_2;
    input.weighted = r -> input.weighted;
    input.new_gen = r -> input.new_gen;
    input.criterion_type = r -> input.criterion_type;
    input.num_func = r -> input.num_func;
    input.amt_factor = r -> input.amt_factor;
    input.amt_exp = r -> input.amt_exp;
    input.num_interval = r -> input.num_interval;
    input.noise = r -> input.noise;
    input.gamma = r -> input.gamma;
    input.auto_gamma = r -> input.auto_gamma;
    input.auto_param = r -> input.auto_param;
    input.auto_param_c = r -> input.auto_param_c;
    input.left_border = r -> input.left_border;
    input.right_border = r -> input.right_border;
    input.kernel_type = r -> input.kernel_type;
    input.var_level = r -> input.var_level;

    curr_kernel = r -> curr_kernel;
    kernel1.coef = r -> kernel1.coef;
    kernel2.sigma = r -> kernel2.sigma;
    kernel3.degree = r -> kernel3.degree;
    kernel3.coef = r -> kernel3.coef;
    kernel3.a = r -> kernel3.a;
 }

// функция копирования генерируемых данных и выделение памяти
void realiz::copy_init(class realiz *r)
{
    // инициализация
    kernel1.m = kernel2.m = kernel3.m = input.amt_factor;
    y_weight = new double [input.amt_exp];
    c_WCV = new double [input.amt_exp+1];
    weight = new double [input.amt_exp];
    ek = new double [input.amt_exp];
    rk = new double [input.amt_exp];

    // элементы СЛАУ для метода LS_SVM
    c = new double [input.amt_exp+1];
    b = new double [input.amt_exp+1];
    A = new double *[input.amt_exp+1];
    for(int i = 0; i < input.amt_exp + 1; i++)
     A[i] = new  double [input.amt_exp + 1];

    y_ist = new double [input.amt_exp];      // истинное значение отклика
    y = new double [input.amt_exp];          // зашумленное значение отклика
    y_tmp = new double [input.amt_exp];      // зашумленное значение отклика
    y_new = new double [input.amt_exp];      // оценка отклика при использовании кросс-проверки
    y2_new = new double [input.amt_exp];     // оценка отклика ппо текущим значениям параметров
     x = new double *[input.amt_exp];
     x_tmp = new double *[input.amt_exp];
    for(int i = 0; i < input.amt_exp; i++)
     {x[i] = new  double [input.amt_factor];
      x_tmp[i] = new double [input.amt_exp];}

    for(int i = 0; i < input.amt_exp; i++)
     {for(int j = 0; j < input.amt_factor; j++)
         x[i][j] = r -> x[i][j];
      y_ist[i] = r -> y_ist[i];
      y[i] = r -> y[i];
     }
}

// функция получения оценки в точке
double realiz::estimate_get(double *vect ,double *sup_val, double b_coef)
{
    double tmp;
    tmp = b_coef;
	for(int i = 0; i < input.amt_exp; i++)
        tmp+= kernel_func_choice(vect, x[i])*sup_val[i];

	return tmp;
}

// выбор типа функции и получения вектора отклика
double realiz::function_choice(int i)
{
    switch (input.num_func)
    {
    case 1:
        return sin(x[i][0])/x[i][0];
    case 2:
        return 0.1*x[i][0];
    case 3:
        return sin(x[i][0]) * cos(x[i][0]);
    case 4:
        return sin(x[i][0])/x[i][1];
    case 5:
        return 1 + 2*x[i][0] + x[i][1] * x[i][1];
    case 6:
        return sin(x[i][0]) * cos(x[i][1]);
    case 7:
        return (x[i][0] * x[i][0] * 2.0 * x[i][1]) / x[i][2];
    }
    throw(string("Неверно задана функция отклика"));
 
}

// решение СЛАУ A*c=b Гаусс
void solve_gauss(double ** A, double * b, double *rez, int n)
{
    // выделение памяти под доп. матрицы СЛАУ
	double **AA;
	AA = new double *[n];
	// переприсваивание входной матрицы, с целью избежания ее изменения
	for(int i = 0; i < n; i++)
	{ AA[i] = new  double [n];
		for(int j = 0; j < n; j++)
				AA[i][j] = A[i][j];
	}

    //верхний треугольный вид
    for(int i = 0; i < n; i++)
    {
        if(!AA[i][i])
        {
            bool flag = false;
            for(int j = i + 1; j < n && !flag; j++)
                if(AA[j][i])
                {
                    for(int k = i; k < n; k++)
                    {
                        double tmp = AA[i][k];
                        AA[i][k] = AA[j][k];
                        AA[j][k] = tmp;
                    }
                    double tmp = b[i];
                    b[i] = b[j];
                    b[j] = tmp;
                    flag = true;
                }
        }
        b[i] = b[i] / AA[i][i];
        for(int j = n - 1; j >= i; j--)
            AA[i][j] = AA[i][j] / AA[i][i];
        for(int j = i + 1; j < n; j++)
        {
            b[j] -= b[i] * AA[j][i];
            for(int k = n - 1; k >= i; k--)
                AA[j][k] -= AA[i][k] * AA[j][i];
        }
    }

    //диагональный вид
    for(int i = n - 1; i > 0; i--)
        for(int j = i - 1; j >= 0; j--)
            b[j] -= AA[j][i] * b[i];

    for(int i = 0; i < n; i++)
        rez[i] = b[i];

	for(int i = 0; i < n; i++)
		delete []  AA[i];
	delete []  AA;
}

// функция составления  СЛАУ A*c=b
void realiz::ls_svm_slae()
{
    // заполнение
    A[0][0] = 0.0;
    b[0] = 0.0;

    for(int i = 1; i < input.amt_exp+1; i++)
    {
        A[i][0] = 1.0;
        A[0][i] = 1.0;
        b[i] = y[i - 1];
         for(int j = 1; j < input.amt_exp+1; j++)
         {
            A[i][j] = kernel_func_choice(x[i-1],x[j-1]);
         }
       if (input.weighted == 0) A[i][i] += 1.0 / input.gamma;
         else A[i][i] += 1.0 / (input.gamma * weight[i-1]);
     }
  }

// функция вызова решателя СЛАУ
 void realiz::solve_slau()
 {
    // вызов функции решения СЛАУ по методу Гаусса
    solve_gauss(A, b, c, input.amt_exp+1);
 }

// функция подсчета значения по быстрому алгоритму loo cv
 double loo_cv_realiz::fast_loo_cv()
{
    // получение значений: матрица А и вектор с
    ls_svm_slae();
    // вызов функции решения СЛАУ
    solve_slau();

	// выделение памяти
    double b_coef;
    double *alfa_coef;
	double **A_inv;
	double tmp = 0.0;
    b_coef = c[0];
    alfa_coef = new double [input.amt_exp];

	//обращение матрицы
	A_inv = InversionMatrix(A, input.amt_exp + 1, input.amt_exp + 1);

	// подсчитаем новый вектор неизвестных при условии кросс-проверки
    for(int p = 1; p < input.amt_exp+1; p++) // цикл по итерациям кросс-проверки
        {
            for(int i = 0; i < input.amt_exp; i++) // цикл по наблюдениям
                    alfa_coef[i] = c[i+1] - (A_inv[i+1][p]*c[p])/A_inv[p][p];
            b_coef = c[0] - (A_inv[0][p]*c[p])/A_inv[p][p];
            alfa_coef[p-1] = 0.0;
            y_new[p-1] = estimate_get(x[p-1], alfa_coef, b_coef);
            tmp+= (y[p-1] - y_new[p-1])*(y[p-1] - y_new[p-1]);
        }

	// очистка памяти
    delete [] alfa_coef;
	for(int i = 0; i < input.amt_exp+1; i++)
		delete [] A_inv[i];
	delete [] A_inv;

    return tmp/input.amt_exp;
}

 // функция подсчета значения по быстрому взвешенному алгоритму LOO_SV
  double  weight_loo_cv_realiz:: weight_fast_loo_cv()
 {
     // получение значений: матрица А и вектор с
     ls_svm_slae();
     // вызов функции решения СЛАУ
     solve_slau();

     // выделение памяти
     double b_coef;
     double *alfa_coef;
     double **A_inv;
     double tmp_r = 0.0;
     b_coef = c[0];
     alfa_coef = new double [input.amt_exp];

     //обращение матрицы
     A_inv = InversionMatrix(A, input.amt_exp + 1, input.amt_exp + 1);

     // подсчитаем новый вектор неизвестных при условии кросс-проверки
     for(int p = 1; p < input.amt_exp+1; p++) // цикл по итерациям кросс-проверки
     {
         for(int i = 0; i < input.amt_exp; i++) // цикл по наблюдениям
            alfa_coef[i] = c[i+1] - (A_inv[i+1][p]*c[p])/A_inv[p][p];
         b_coef = c[0] - (A_inv[0][p]*c[p])/A_inv[p][p];
         alfa_coef[p-1] = 0.0;
         y_new[p-1] = estimate_get(x[p-1], alfa_coef, b_coef);
      }

    // очистка памяти
    delete [] alfa_coef;
    for(int i = 0; i < input.amt_exp+1; i++)
        delete [] A_inv[i];
    delete [] A_inv;

    return tmp_r/input.amt_exp;
 }

 // функция сортировки выборки для подсчета ОКР
 void regularity_realiz::sort()
    {
     map <double, int> d_sort;              // контейнер для сортировки
     double *centr;                         // вектор координат центра эксперимента
     double *distance;                      // вектор расстояний от центра эксперимента до точки
     centr = new double [input.amt_factor];
     distance = new double [input.amt_exp];

     // временные переменные для сортировки
     double *y_ist_tmp;
     y_ist_tmp = new double [input.amt_exp];

     //вычисление координат центра эксперимента
     for(int i = 0; i < input.amt_factor; i++)
         centr[i] = (input.var_level[i].second + input.var_level[i].first)/2.0;

     //вычисление расстояний
     for(int i = 0; i < input.amt_exp; i++)
         distance[i] = euclidean_distance(x[i], centr, input.amt_factor);

     // заполнение контейнера
     for(int i = 0; i < input.amt_exp; i++)
         d_sort[distance[i]] = i;

     //сортировка матрицы наблюдений и вектора отклика
     int i = 0;
     for(map<double, int>::const_iterator it = d_sort.begin(); it != d_sort.end(); it++, i++)
        {
            for(int j = 0; j < input.amt_factor; j++)
                x_tmp[i][j] = x[it->second][j];
            y_tmp[i] = y[it->second];
            y_ist_tmp[i] = y_ist[it->second];
        }

     for(int i = 0; i < input.amt_exp; i++)
        {
            for(int j = 0; j < input.amt_factor; j++)
                x[i][j] = x_tmp[i][j];
            y[i] = y_tmp[i];
            y_ist[i] = y_ist_tmp[i];
        }

    // очистка памяти
    delete [] y_ist_tmp;
 }

 // функция исключения наблюдений для ОКР
 void regularity_realiz::exclusion(int p, int n, int size)
 {
    // p - номер тестового наблюдения в интервале
    int tmp = 0;

    //перезаписываем матрицу х и вектор y для тестовой выборки
    for(int i = 0; i < n; i++)
    {
        if (i != p)
        {  for(int j = 0; j < input.amt_factor; j++)
                x[tmp][j] = x_tmp[i][j];
           y[tmp] = y[i];
           tmp++;
        }
        else p+= size;
    }
 }

// функция подсчета значения по ОКР
 double regularity_realiz::cr_regularity()
 {    
    int amt_exp_tmp = input.amt_exp;
    //input.num_interval = 25;         // количество интервалов разбиения
    int size_int;                    // количество наблюдений в одном интервале
    double tmp = 0.0;

    // Подсчет количества наблюдений в зависимости от количества интервалов
    size_int = input.amt_exp / input.num_interval;
    input.amt_exp -= input.num_interval; // изменяем размерность

    // подсчитаем новый вектор неизвестных при условии ОКР
    for(int p = 0; p < size_int; p++) // цикл по итерациям проверки, p - номер тестового наблюдения в интервале
    {
        exclusion(p, amt_exp_tmp, size_int); // исключение тестовых наблюдений
        ls_svm_slae();                       // получение значений: матрица А и вектор b
        solve_slau();                        // вызов функции решения СЛАУ

        // получение оценок тестовой выборки
        for(int i = p; i < amt_exp_tmp; i+=size_int)
        {
           y_new[i] = estimate_get(x_tmp[i], c + 1, c[0]);
            tmp+= (y_tmp[i] - y_new[i])*(y_tmp[i] - y_new[i]);
        }
    }

     //восстанавление размерности
    input.amt_exp = amt_exp_tmp;

    for(int i = 0; i < input.amt_exp; i++)
    {
     for(int j = 0; j < input.amt_factor; j++)
        x[i][j] = x_tmp[i][j];
     y[i] = y_tmp[i];

    }
    return tmp/input.amt_exp;
}

// функция получения оценок отклика по текущим значениям параметров
 void realiz::calc_y_estimate()
 {
    for(int i = 0; i < input.amt_exp; i++)
     y2_new[i] = estimate_get(x[i], c + 1, c[0]);
 }


 // функция подсчета MSE (среднеквадратическая ошибка)
 double realiz::calc_MSE()
 {
    double mse = 0.0;
    double tmp;
    for(int i = 0; i < input.amt_exp; i++)
    {
        tmp = (y_ist[i] - y2_new[i]);
        mse += tmp*tmp;
    }
    return mse/input.amt_exp;

}

 // очистка памяти
 void realiz::clear_memory()
 {
	if(x)
	{	for(int i = 0; i < input.amt_exp; i++)
			delete [] x[i];
        delete [] x;
		x = NULL;
    }

	if(y)
	{
		delete [] y;
		y = NULL;
	}

    if(y_tmp)
    {
        delete [] y_tmp;
        y_tmp = NULL;
    }

	if(y_ist)
	{
		delete [] y_ist;
		y_ist = NULL;
	}

	if(y_new)
	{
		delete [] y_new;
		y_new = NULL;
	}

    if(y2_new)
    {
        delete [] y2_new;
        y2_new = NULL;
    }

	if(A)
	{
		for(int i = 0; i < input.amt_exp+1; i++)
			delete [] A[i];
		delete [] A;
		A = NULL;
	}

	if(c) 
	{	
		delete [] c;
		c = NULL;
	}

    if(b)
    {
        delete [] b;
        b = NULL;
    }
    if (weight)
    {
        delete [] weight;
        weight = NULL;
    }

    if (ek)
    {
        delete [] ek;
        ek = NULL;
    }

    if (rk)
    {
        delete [] rk;
        rk = NULL;
    }

    if (y_weight)
    {
        delete [] y_weight;
        y_weight = NULL;
    }
    if (c_normal)
    {
        delete [] c_normal;
        c_normal = NULL;
    }
    if (c_weight)
    {
        delete [] c_weight;
        c_weight = NULL;
    }
    if (c_WCV)
    {
        delete [] c_WCV;
        c_WCV = NULL;
    }
 }

 


