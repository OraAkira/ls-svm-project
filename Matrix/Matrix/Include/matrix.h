#ifndef __MATRIX_H__

#define __MATRIX_H__

#include <vector>

class Matrix
{
    private:
        int row;
        int column;
        double **data;

        friend Matrix operator + (Matrix &, Matrix &);
        friend Matrix operator - (Matrix &, Matrix &);
        friend Matrix operator * (Matrix &, Matrix &);
        friend Matrix operator * (const double, Matrix &);
        friend Matrix operator / (Matrix &, Matrix &);
        friend Matrix operator ^ (Matrix &, const int32_t);
        friend Matrix operator ^ (Matrix &, const double);
        friend Matrix operator ^ (const double, Matrix &);
    public:
        Matrix();
        Matrix(int, int);
        Matrix(Matrix&);
		Matrix(std::vector<double> &);
        ~Matrix();

        bool SwitchRow(int, int);
        bool SwitchColumn(int, int);
        bool RowTransform(int, int, double);
        bool ColumnTransform(int, int, double);
        
        // basic operations
		int Row();                            // rows of matrix object
        int Column();                         // columns of matrix object
        int Rank();                           // rank of matrix object
        Matrix DotProduct(Matrix &);          // dot product of matrix object
        Matrix DotDivide(Matrix &);           // dot divide of matrix onject
        Matrix Transpose();                   // transpose matrix of matrix object
        Matrix Sum();                         // the sum of every column
        Matrix Mean();                        // the average of every column
        Matrix Max();                         // the maximun of every column
        Matrix Min();                         // the minimun of every column

        // about squrare matrix
        Matrix Inverse();                 // inverse of matrix object
        double Determinant();             // Determinant of a square
        double Trace();                   // trace of matrix object
        
        // the value of matrix
        double & operator () (int i, int j);

        // operators of matrix
        Matrix & operator += (Matrix &);
        Matrix & operator -= (Matrix &);
        Matrix & operator = (Matrix &);
        bool operator == (Matrix &);

		friend std::ostream & operator << (std::ostream &os, Matrix &c);
		friend std::istream & operator >> (std::istream &is, Matrix &c);

        // the method of matrix object
        void Print();

		static void errorlog(int);
        static Matrix eyes(int order);
        static Matrix zeros(int order);
        static Matrix linspace(int left, int right, int num);
};

#endif