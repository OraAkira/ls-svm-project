#include "../Include/matrix.h"
#include <iostream>
#include <vector>

// Switch
bool Matrix::SwitchRow(int x, int y) {
    if( (x < 1 && x > this->row) || (y < 1 && y > this->column) ) {
        Matrix::errorlog(5);
        return false;
    }
    double temp = 0;
    for(int i = 1; i <= this->column; i++) {
        temp = (*this)(x, i);
        (*this)(x, i) = (*this)(y, i);
        (*this)(y, i) = temp;
    }
    return true;
}

bool Matrix::SwitchColumn(int x, int y) {
    if( (x < 1 && x > this->row) || (y < 1 && y > this->column) ) {
        Matrix::errorlog(5);
        return false;
    }
    double temp = 0;
    for(int i = 1; i <= this->row; i++) {
        temp = (*this)(i, x);
        (*this)(i, x) = (*this)(i, y);
        (*this)(i, y) = temp;
    }
    return true;
}

bool Matrix::RowTransform(int x, int y, double k) {
    if( (x < 1 && x > this->row) || (y < 1 && y > this->column) ) {
        Matrix::errorlog(5);
        return false;
    }
    for(int i = 1; i <= this->column; i++) {
        (*this)(x, i) += (*this)(y, i) * k;
    }
    return true;
}

bool Matrix::ColumnTransform(int x, int y, double k) {
    if( (x < 1 && x > this->row) || (y < 1 && y > this->column) ) {
        Matrix::errorlog(5);
        return false;
    }
    for(int i = 1; i <= this->row; i++) {
        (*this)(i, x) += (*this)(i, y) * k;
    }
    return true;
}


// operator override
Matrix operator + (Matrix & A, Matrix & B) {
    int row = A.Row();
    int column = A.Column();
    if(row != B.Row() || column != B.Column()) {
        Matrix::errorlog(1);
        return A;
    }
    Matrix temp(row, column);
    for(int i = 1; i <= row; i++) {
        for(int j = 1; j <= column; j++) {
            temp(i, j) = A(i, j) + B(i, j);
        }
    }
    return temp;
}

Matrix operator - (Matrix & A, Matrix & B) {
    int row = A.Row();
    int column = A.Column();
    if(row != B.Row() || column != B.Column()) {
        Matrix::errorlog(1);
        return A;
    }
    Matrix temp(row, column);
    for(int i = 1; i <= row; i++) {
        for(int j = 1; j <= column; j++) {
            temp(i, j) = A(i, j) - B(i, j);
        }
    }
    return temp;
}

Matrix operator * (Matrix & A, Matrix & B) {
    int row = A.Row();
    int column = B.Column();
    if(A.Column() != B.Row()) {
        Matrix::errorlog(1);
        return A;
    }

    int cnt = A.Column();
    Matrix temp(row, column);
    for(int i = 1; i <= row; i++) {
        for(int j = 1; j <= column; j++) {
            temp(i, j) = 0;
            for(int k = 1; k <= cnt; k++) {
                temp(i, j) += A(i, k) * B(k, j);
            }
        }
    }
    return temp;
}

Matrix operator * (const double k, Matrix & A) {
    int row = A.Row();
    int column = A.Column();
    Matrix temp(row, column);
    for(int i = 1; i <= row; i++) {
        for(int j = 1; j <= column; j++) {
            temp(i, j) = k * A(i, j);
        }
    }
    return temp;
}

Matrix operator / (Matrix & A, Matrix & B) {
    return A * B.Inverse();
}

/*
Matrix operator \ (const Matrix & A, const Matrix & B) {
    return A.Inverse() * B;
}
*/

Matrix operator ^ (Matrix & A, int32_t n) {
	Matrix AA(A);
    Matrix temp = Matrix::eyes(A.Row());
    while(n) {
        if(n & 1) {
            temp = temp * AA;
        }
        n >>= 1;
        AA = AA * AA;
    }
    return temp;
}

// TODO  return x^A
Matrix operator ^ (const double x, Matrix & A) {
	return A;
}

std::ostream & operator << (std::ostream & os, Matrix & c)
{
	int row = c.Row();
	int col = c.Column();
	os << std::endl;
	for (int i = 1; i <= row; i++) {
		os << c(i, 1);
		for (int j = 2; j <= col; j++) {
			os << '\t' << c(i, j);
		}
		os << std::endl;
	}
	return os;
}

std::istream & operator >> (std::istream & is, Matrix & c)
{
	int row = c.Row();
	int col = c.Column();
	for (int i = 1; i <= row; i++) {
		for (int j = 1; j <= col; j++) {
			is >> c(i, j);
		}
	}
	return is;
}

// TODO return A^p
Matrix operator ^ (Matrix & A, const double p) {
	return A;
}
    


bool Matrix::operator == (Matrix & other) {
    if(this->row != other.Row() || this->column != other.Column()) {
        return false;
    }
    for(int i = 1; i <= this->row; i++) {
		for (int j = 1; j <= this->column; j++) {
			if ((*this)(i, j) != other(i, j)) {
				return false;
			}
		}
    }
	return true;
}

double & Matrix::operator () (int i, int j) {
    if(i < 1 || j < 1) {
        Matrix::errorlog(2);
        return this->data[ 0 ][ 0 ];
    }
	return this->data[i - 1][j - 1];
}


Matrix & Matrix::operator += (Matrix & other) {
    if(this->row != other.Row() || this->column != other.Column()) {
        Matrix::errorlog(1);
        return (*this);
    }
    for(int i = 1; i <= this->row; i++) {
        for(int j = 1; j <= this->column; j++) {
            (*this)(i, j) += other(i, j);
        }
    }
    return (*this);
}

Matrix & Matrix::operator -= (Matrix & other) {
    if(this->row != other.Row() || this->column != other.Column()) {
        Matrix::errorlog(1);
        return (*this);
    }
    for(int i = 1; i <= this->row; i++) {
        for(int j = 1; j <= this->column; j++) {
            (*this)(i, j) -= other(i, j);
        }
    }
    return (*this);
}

Matrix & Matrix::operator = (Matrix & other) {
    if(this != &other) {
        for(int i = 0; i < this->row; i++) {
            delete [] this->data[ i ];
        }
        delete [] this->data;

        this->row = other.Row();
        this->column = other.Column();
        this->data = new double * [ this->row ];

        for(int i = 0; i < this->row; i++) {
            this->data[ i ] = new double [ this->column ];
        }
    
        for(int i = 0; i < this->row; i++) {
            for(int j = 0; j < this->column; j++) {
                this->data[ i ][ j ] = other(i+1, j+1);
            }
        }
    }
    return (*this);
}
