#include "LS_SVM.h"

// конструктор
realiz::realiz()
{
	y_ist = NULL;
	y_new = NULL;
    y2_new = NULL;
	y = NULL;
	x = NULL;
	A = NULL;
	c = NULL;
}

// деструктор
realiz::~realiz()
{
	clear_memory();
}

// функция генерации равномерно-распределенных чисел
double realiz::rav_num_gen(double left, double right)
{
	double alpha1 = rand() / (double)RAND_MAX;
	double rand_num =  left + alpha1 * (right - left);
	return rand_num;
}

// функция генерации нормально-распределенных чисел
double realiz::normal_num_gen(double m, double sigma22)
{
    double alpha1 = rand() / (double)RAND_MAX;
    double alpha2 = rand() / (double)RAND_MAX;
    double ksi = sqrt(-2.0 * log(alpha1)) * cos(2.0 * 3.14159 * alpha2);
     return m + sqrt(sigma22) * ksi;
}

// метод золотого сечения
double realiz::golden()
{
	double eps = 0.0001;
    double a = input.left_border, b = input.right_border;
    double q1 = (3.0 - sqrt(5.0)) / 2.0, q2 = (sqrt(5.0) - 1.0) / 2.0;
    double x1 = a + q1 * (b - a), x2 = a + q2 * (b - a);
	input.gamma = x1;
    double f1 = fast_loo_cv();
	input.gamma = x2;
	double f2 = fast_loo_cv();

    while(true)
    {
        if(f1 <= f2)
        {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + q1 * (b - a);
			input.gamma = x1;
            f1 = fast_loo_cv();
        }
        else
        {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + q2 * (b - a);
			input.gamma = x2;
            f2 = fast_loo_cv();
        }
        if(fabs(a - b) < eps)
            return (a + b) / 2.0;
    }
}


// функция автоподбора параметров gamma
void realiz::gamma_get() 
{
	input.gamma = golden();
}

void realiz::kernel_param_get() // функция автоподбора параметров ядра
{ 
	double min;
	double tmp;
	double tmp2;
    double max_degree = 2.8;
    double step = 0.1;
	double degree;
	int index = 0;
	if (curr_kernel == 2) 
	{	
		kernel2.sigma = kernel2.sigma_select[0];
		min = fast_loo_cv();
		for (int i = 1; i < 10; i++)
		 {
			kernel2.sigma = kernel2.sigma_select[i];
			tmp = fast_loo_cv();
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
		min = fast_loo_cv();
		while (tmp < max_degree)
		{
            //tmp += step;
            tmp = (double)(floor((tmp + step) * 100.0 + 0.5)) / 100.0;
			kernel3.degree = tmp;
			tmp2 = fast_loo_cv();
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

    c = new double [input.amt_exp+1];
    A = new double *[input.amt_exp+1];
    for(int i = 0; i < input.amt_exp + 1; i++)
     A[i] = new  double [input.amt_exp + 1];

	y_ist = new double [input.amt_exp];
	y_new = new double [input.amt_exp];
    y2_new = new double [input.amt_exp];
	y = new double [input.amt_exp];
	x = new double *[input.amt_exp];

	for(int i = 0; i < input.amt_exp; i++)
     x[i] = new  double [input.amt_factor];

	// заполнение матрицы наблюдений
	for(int i = 0; i < input.amt_exp; i++)
	 for(int j = 0; j <  input.amt_factor; j++)
	 {
		x[i][j] = rav_num_gen (input.var_level[j].first, input.var_level[j].second);
	 }

	 // вычисление массива откликов
	for(int i = 0; i < input.amt_exp; i++)
		y[i] = y_ist[i] = function_choice(i);  
	
	// вычисление дисперсии
	double *aid;    
	aid = new double [input.amt_exp];
	double y_mean = 0.0; // среднее значение отклика
	double w2 = 0.0;	 // мощность сигнала
	double sigma2;		 // дисперсия
	double e;			 // ошибка 
    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<< endl;
	// зашумление отклика
	for(int i = 0; i < input.amt_exp; i++)
		y_mean += y[i];

	y_mean = y_mean/input.amt_exp;

	for(int i = 0; i < input.amt_exp; i++)
		aid[i] = y[i] - y_mean;

	for(int i = 0; i < input.amt_exp; i++)
		w2 += aid[i] * aid[i];

	w2 = w2 / (input.amt_exp - 1);
	sigma2 = input.noise * w2 / 100;

	for(int i = 0; i < input.amt_exp; i++)
	{
		e = normal_num_gen(0.0, sigma2);
		y[i] += e; 
        // распечатка
        cout << x[i][0] << endl;
	}

	delete [] aid;
}

// функция получения оценки
double realiz::estimate_get(double *vect,double *sup_val, double b)
{
	double tmp;
	tmp = b;
	for(int i = 0; i < input.amt_exp; i++)
		tmp+= kernel_func_choice(vect,x[i])*sup_val[i];

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
        return 1 + 2*x[i][0] + x[i][0] * x[i][0];
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
    throw 2;
 
}

// решение СЛАУ A*c=b Гаусс

void solve_gauss(double ** A, double * b, double * rez, int n)
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


// функция составления и решения СЛАУ A*c=b
void realiz::ls_svm_slae()
{
    double *b;
    b = new double [input.amt_exp+1];

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
		 A[i][i] += 1.0 / input.gamma;
	}

	// вызов функции решения СЛАУ
    solve_gauss(A, b, c, input.amt_exp+1);  // Гаусс
	delete [] b;
}


// функция подсчета значения по быстрому алгоритму LOO_SV
 double realiz::fast_loo_cv()
{
	// получение значений: матрица А и вектор с
	ls_svm_slae(); 

	// выделение памяти
	double b; 
	double *alfa; 
	double **A_inv;
	double tmp = 0.0;
	b = c[0];
	alfa = new double [input.amt_exp];

	//обращение матрицы
	A_inv = InversionMatrix(A, input.amt_exp + 1, input.amt_exp + 1);

	// подсчитаем новый вектор неизвестных при условии кросс-проверки
    for(int p = 1; p < input.amt_exp+1; p++) // цикл по итерациям кросс-проверки
    {
        for(int i = 0; i < input.amt_exp; i++) // цикл по наблюдениям
        {
            alfa[i] = c[i+1] - (A_inv[i+1][p]*c[p])/A_inv[p][p];
        }
        b = c[0] - (A_inv[0][p]*c[p])/A_inv[p][p];
        alfa[p-1] = 0.0;
        y_new[p-1] = estimate_get(x[p-1], alfa, b);
        tmp+= (y[p-1] - y_new[p-1])*(y[p-1] - y_new[p-1]);
    }


	// очистка памяти
	delete [] alfa;
	for(int i = 0; i < input.amt_exp+1; i++)
		delete [] A_inv[i];
	delete [] A_inv;

return tmp/input.amt_exp;
}

// функция получения оценок отклика по текущим значениям параметров
 void realiz::calc_y_estimate()
 {
    for(int i = 1; i < input.amt_exp+1; i++)
    y2_new[i-1] = estimate_get(x[i-1], c + 1, c[0]);

}


 // функция подсчета MSE
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
 }

 


