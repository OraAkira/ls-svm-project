#include "LS_SVM.h"

/* Обращения матрицы производится с использованием метода решения квадратичных СЛАУ Гаусса - Жордана. 
Решение СЛАУ вида АХ=I, где A - известная матрица, X - неизвестная матрица, I - единичная матрица. */  


//Функция приведение матрицы к верхнему треугольному виду 
bool realiz::trianglematrix1(double **A, double **X, int m)
{ 
	bool flag = true;

	for(int i = 0 ; i < m ; i++)
	{ 
        double A_d = A[i][i];//ведущий элемент

		transf1(A,X,i,m,A_d);//переставим строки местами 

        if(fabs(A_d) < 1E-20)  return flag = false;

        for(int j = 0; j < m; j++)
		{ 
			if (i != j)
            {
			   double A_j = A[j][i] / A[i][i]; 

               for(int k = 0; k < m; k++) 
                   { 
					A[j][k] = A[j][k] - A[i][k]*A_j; 
					X[j][k] = X[j][k] - X[i][k]*A_j; 
					}
            } 
        } 
    } 
    return flag; 
} 

// Функция перестановка строк
void realiz::transf1(double **A, double **X, int i, int m, double &A_d)
{ 
    int line1 = i;

    for(int j = i; j < m; j++) 
        if( fabs(A[j][i]) > fabs(A_d)) 
		{
			line1 = j; 
			A_d = A[j][i];
		}

	double c;   
	for(int j = 0; j < m; j++)
		{ 
			//меняем местами строки в матрице А
            c = A[i][j]; 
            A[i][j] = A[line1][j]; 
            A[line1][j] = c; 
			//меняем местами строки в матрице Х
			c = X[i][j]; 
			X[i][j] = X[line1][j]; 
			X[line1][j] = c; 
		} 
}

// Функция решения СЛАУ
bool realiz::SLAU_solution_Gauss(double **A, double **I, double **X, int m)
{ 
	for(int i = 0 ; i < m; i++) 
		for(int j = 0 ; j < m; j++) 
			X[i][j] = I[i][j]; //предполагаем, что обратная матрица единичная 


    trianglematrix1(A,X,m); //приведем матрицу А к треугольному виду 
	bool flag = trianglematrix1(A,X,m);

    if(flag == false) return false; 

	double sum=0; 

	for (int i = 0; i < m ; i++)
	{
		sum = A[i][i];
		A[i][i] = A[i][i]/sum;

		for(int k = 0; k < m; k++) 
			X[i][k] = X[i][k]/sum;
	} 

    return true; 
} 

//о Функция обращения матрицы A = A^(-1)
double **realiz::InversionMatrix(double **F, int n, int m)
{

	// выделение памяти под доп. матрицы СЛАУ
	double **A;
	double **I;
	double **rez;
	I = new double*[n];
	rez = new double*[n];
	A = new double *[n];
	// переприсваивание входной матрицы, с целью избежания ее изменения
	for(int i = 0; i < n; i++)
	{	I[i]=new double[m];
		rez[i]=new double[m];
		A[i] = new  double [m];
		for(int j = 0; j < m; j++)
				A[i][j] = F[i][j];
	}

 // заполнение единичной матрицы
	for(int i = 0 ; i < n; i++) 
		for(int j = 0 ; j < m; j++) 
		 {
			if (i==j) I[i][j] = 1.0;
			else I[i][j] = 0.0;
		}

	// вызов решалки СЛАУ
	if(!SLAU_solution_Gauss(A,I,rez,n)) throw 1;

	// очистка памяти
	for(int i = 0; i < n; i++)
	{	delete [] I[i];
		delete [] A[i];}
	delete [] I;
	delete [] A;
	return rez;
}
