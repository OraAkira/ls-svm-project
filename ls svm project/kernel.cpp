#include "LS_SVM.h"

// функция вычисления скалярного произведения
double dot_product(double *vect1, double *vect2, int m)
{	
	double d_p = 0.0;
	for(int i = 0; i < m; i++)
        d_p += vect1[i] * vect2[i];
    return d_p;
}

// конструктор для создания набора допустимых значений параметра RBF - ядра
kernel_RBF::kernel_RBF()
{
	double p = 3.5;
	for(int i = 0; i < 10; i++)
	{ 
		sigma_select[i] = pow(10,p);
		p -= 0.5;
	}
}

// функция вычисления значения для линейного ядра
double kernel_linear::kernel_func(double *vect1, double *vect2)
{
	return (dot_product(vect1, vect2, m) + coef);
}

// функция вычисления значения  RBF - ядра 
double kernel_RBF::kernel_func(double *vect1, double *vect2)
{
	double tmp = 0.0;
	for(int i = 0; i < m; i++)
	{ 
			tmp += (vect1[i] - vect2[i])*(vect1[i] - vect2[i]); 	
	}
	tmp /= - (2.0 * sigma);
    return exp(tmp);
}

// функция вычисления значения полиномиального ядра
double kernel_polynomial::kernel_func(double *vect1, double *vect2)
{
	return pow( a*(dot_product(vect1, vect2, m)) + coef, degree);
}

// выбор соответствующей функции ядра
double realiz::kernel_func_choice(double *vect1, double *vect2)
{
   switch (curr_kernel)
 {
 case 1:
 return kernel1.kernel_func(vect1, vect2);
 case 2:
 return kernel2.kernel_func(vect1, vect2);
 case 3:
 return kernel3.kernel_func(vect1, vect2);
 }
 throw 0;
}
