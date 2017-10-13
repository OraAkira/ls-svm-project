#include "../Include/matrix.h"
#include <iostream>
#include <vector>
#include <cmath>

// Basic operations
int Matrix::Row() {
    if(this == nullptr) {
		Matrix::errorlog(0);
        return -1;
    }
    return this->row;
}

int Matrix::Column() {
    if(this == nullptr) {
        Matrix::errorlog(0);
        return -1;
    }
    return this->column;
}

int Matrix::Rank() {
    if(this == nullptr) {
        Matrix::errorlog(0);
        return -1;
    }
    Matrix temp(*this);
    if(this->row > this->column) {
        temp = temp.Transpose();
    }

	int rows = temp.Row();
	int cols = temp.Column();
    for(int j = 1; j <= rows; j++) {
        // Select the element with the largest absolute value 
        double maxEle = fabs(temp(j, j));
        int maxRow = j;
        int maxCol = j;
        for(int i = j; i <= rows; i++) {
            for(int k = j; k <= cols; k++) {
                if(maxEle < fabs(temp(i, k))) {
                    maxEle = fabs(temp(i, k));
                    maxRow = i;
                    maxCol = k;
                }
            }
        }
        if(maxEle < 10e-6) {
            return (j - 1);
        }
        // Switch to diagonal
        if(maxRow != j)
            temp.SwitchRow(maxRow, j);
        if(maxCol != j)
            temp.SwitchColumn(maxCol, j);
        // Elimination
        for(int i = j+1; i <= rows; i++) {
			temp.RowTransform(i, j, -(temp(i, j) / temp(j, j)));
        }
    }
    return temp.Row();
}

Matrix Matrix::DotProduct(Matrix &other) {
    if(this == nullptr) {
        Matrix::errorlog(0);
        return (*this);
    }
    if(this->row != other.Row() || this->column != other.Column()) {
        Matrix::errorlog(1);
        return (*this);
    }

    Matrix temp(this->row, this->column);
    for(int i = 0; i < this->row; i++) {
        for(int j = 0; j < this->column; j++) {
            temp(i+1, j+1) = this->data[ i ][ j ] * other(i+1, j+1);
        }
    }
    return temp;
}

Matrix Matrix::DotDivide(Matrix &other) {
    if(this == nullptr) {
        Matrix::errorlog(0);
        return (*this);
    }
    if(this->row != other.Row() || this->column != other.Column()) {
        Matrix::errorlog(1);
        return (*this);
    }

    Matrix temp(this->row, this->column);
    for(int i = 0; i < this->row; i++) {
        for(int j = 0; j < this->column; j++) {
            temp(i+1, j+1) = this->data[ i ][ j ] / other(i+1, j+1);
        }
    }
    return temp;
}

Matrix Matrix::Transpose() {
    if(this == nullptr) {
        Matrix::errorlog(0);
        return (*this);
    }
    
    Matrix temp(this->column, this->row);
    for(int i = 0; i < this->column; i++) {
        for(int j = 0; j < this->row; j++) {
            temp(i+1, j+1) = this->data[ j ][ i ];
        }
    }
    return temp;
}


Matrix Matrix::Sum() {
    if(this == nullptr) {
        Matrix::errorlog(0);
        return (*this);
    }
    Matrix temp(1, this->column);
    for(int i = 1; i <= this->column; i++) {
        temp(1, i) = 0;
        for(int j = 1; j <= this->row; j++) {
            temp(1, i) += this->data[ j-1 ][ i-1 ];
        }
    }
    return temp;
}

Matrix Matrix::Mean() {
    if(this == nullptr) {
        Matrix::errorlog(0);
        return (*this);
    }
    Matrix temp(1, this->column);
    for(int i = 1; i <= this->column; i++) {
        temp(1, i) = this->Sum()(1, i) / this->row;
    }
    return temp;
}

Matrix Matrix::Max() {
    if(this == nullptr) {
        Matrix::errorlog(0);
        return (*this);
    }
    Matrix temp(1, this->column);
    for(int i = 1; i <= this->column; i++) {
        temp(1, i) = this->data[ 0 ][ i-1 ];
        for (int j = 2; j <= this->row; j++) {
            if(temp(1, i) < this->data[ j-1 ][ i-1 ]) {
				temp(1, i) = this->data[j - 1][i - 1];
            }
        }
    }
    return temp;
}

Matrix Matrix::Min() {
    if(this == nullptr) {
        Matrix::errorlog(0);
        return (*this);
    }
    Matrix temp(1, this->column);
    for(int i = 1; i <= this->column; i++) {
        temp(1, i) = this->data[ 0 ][ i-1 ];
        for (int j = 2; j <= this->row; j++) {
            if(temp(1, i) > this->data[ j-1 ][ i-1 ]) {
				temp(1, i) = this->data[j - 1][i - 1];
            }
        }
    }
    return temp;
}

// Square matrix
Matrix Matrix::Inverse() {
    if(this == nullptr) {
        Matrix::errorlog(0);
        return (*this);
    }
    if(this->row != this->column) {
        Matrix::errorlog(6);
        return (*this);
    }
    if(this->Determinant() == 0) {
        Matrix::errorlog(7);
        return (*this);
    }

    int len = this->row;
    Matrix L = Matrix::eyes(len);
    Matrix U = Matrix::eyes(len);

    for(int i = 1; i <= len; i++) {
        L(i, i) = 1;
        U(1, i) = (*this)(1, i);
        L(i, 1) = (*this)(i, 1) / U(1, 1);
    }
    for(int i = 2; i <= len; i++) {
        for(int j = i; j <= len; j++) {
            U(i, j) = (*this)(i, j);
            for(int k = 1; k < i; k++) {
                U(i, j) -= L(i, k) * U(k, j);
            }
        }
        for(int j = 2; ((i < len) && (j <= i)); j++) {
            L(i+1, j) = (*this)(i+1, j);
            for(int k = 1; k < j; k++) {
                L(i+1, j) -= L(i+1, k) * U(k, j);
            }
			L(i + 1, j) /= U(j, j);
        }
    }

    Matrix l = Matrix::zeros(len);
    Matrix u = Matrix::zeros(len);

    for(int i = 1; i <= len; i++) {
		l(i, i) = 1;
		for (int k = i + 1; k <= len; k++) {
			for (int j = i; j < k; j++) {
				l(k, i) = l(k, i) - L(k, j) * l(j, i);
			}
		}
    }
    for(int i = 1; i <= len; i++) {
		u(i, i) = 1.0 / U(i, i);
		for (int k = i - 1; k >= 1; k--) {
			double s = 0;
			for (int j = k + 1; j <= i;j++) {
				s += U(k, j) * u(j, i);
			}
			u(k, i) = -s / U(k, k);
		}
    }

    return (u * l);
}

double Matrix::Determinant() {
    if(this == nullptr) {
        Matrix::errorlog(0);
        return 0;
    }
    if(this->row != this->column) {
        Matrix::errorlog(6);
        return 0;
    }
	Matrix temp(*this);
	double result = 1;
	int count = 1;
	int rows = temp.Row();
	int cols = temp.Column();
	for (int j = 1; j <= rows; j++) {
		// Select the element with the largest absolute value 
		double maxEle = fabs(temp(j, j));
		int maxRow = j;
		int maxCol = j;
		for (int i = j; i <= rows; i++) {
			for (int k = j; k <= cols; k++) {
				if (maxEle < fabs(temp(i, k))) {
					maxEle = fabs(temp(i, k));
					maxRow = i;
					maxCol = k;
				}
			}
		}
		if (maxEle < 10e-6) {
			return 0;
		}
		// Switch to diagonal
		if (maxRow != j) {
			temp.SwitchRow(maxRow, j);
			count *= -1;
		}
		if (maxCol != j) {
			temp.SwitchColumn(maxCol, j);
			count *= -1;
		}
		result *= temp(j, j);
		// Elimination
		for (int i = j + 1; i <= rows; i++) {
			temp.RowTransform(i, j, -(temp(i, j) / temp(j, j)));
		}
	}
    return result * count;
}

double Matrix::Trace() {
    if (this == nullptr) {
        Matrix::errorlog(0);
        return -1;
    }
    if (this->row != this->column) {
        Matrix::errorlog(6);
        return -1;
    }
    double temp = 0;
    for (int i = 0; i < this->row; i ++) {
        temp += this->data[i][i];
    }
    return temp;
}


void Matrix::Print() {
    std::cout << typeid(*this).name() << " = " << std::endl;
    for(int i = 1; i <= this->row; i++) {
		std::cout << (*this)(i, 1);
        for(int j = 2; j <= this->column; j++) {
            std::cout << '\t' << (*this)(i, j);
        }
        std::cout << std::endl;
    }
}

