#include "../Include/matrix.h"
#include <iostream>
#include <vector>

// Constructor
Matrix::Matrix() {
    this->row = 1;
    this->column = 1;

    this->data = new double * [ this->row ];
    for(int i = 0; i < this->row; i++) {
        this->data[ i ] = new double [ this->column ];
    }

    for(int i = 0; i < this->row; i++) {
        for(int j=0;j<this->column;j++) {
            this->data[ i ][ j ] = 0;
        }
    }
}

Matrix::Matrix(Matrix &other) {
    this->row = other.Row();
    this->column = other.Column();
    
    this->data = new double * [ this->row ];
    for(int i = 0; i < this->row; i++) {
        this->data[ i ] = new double [ this->column ];
    }
    
    for(int i = 0; i < this->row; i++) {
        for(int j = 0; j < this->column; j++) {
            this->data[ i ][ j ] = other(i + 1, j + 1);
        }
    }
}

Matrix::Matrix(std::vector<double> &v) {
    this->row = v.size();
    this->column = 1;

    this->data = new double * [ this->row ];
    for(int i = 0; i < this->row; i++) {
        this->data[ i ] = new double [ this->column ];
    }
    
    for(int i = 0; i < this->row; i++) {
        for(int j = 0; j < this->column; j++) {
            this->data[ i ][ j ] = v[ i ];
        }
    }
}

Matrix::Matrix(int rows, int columns) {
    this->row = rows;
    this->column = columns;
    
    this->data = new double * [ this->row ];
    for(int i = 0; i < this->row; i++) {
        this->data[ i ] = new double [ this->column ];
    }

    for(int i = 0; i < this->row; i++) {
        for(int j=0;j<this->column;j++) {
            this->data[ i ][ j ] = 0;
        }
    }
}

// Destructor
Matrix::~Matrix() {
    for(int i = 0; i < this->row; i++) {
        delete [] this->data[ i ];
    }
    delete [] this->data;
}

void Matrix::errorlog(int id) {
    switch(id) {
        case 0:
            std::cout << "error: Empty object!\n" << std::endl;
            break;
        case 1:
            std::cout << "error: Dimension doesn't match!\n" << std::endl;
            break;
        case 2:
            std::cout << "error: Indices must be real positive integers!\n" << std::endl;
            break;
        case 3:
            std::cout << "error: Undefined function or object!\n" << std::endl;
            break;
        case 4:
            std::cout << "error: Matrix dimensions must agree!\n" << std::endl;
            break;
        case 5:
            std::cout << "error: Index exceeds matrix dimensions!\n" << std::endl;
            break;
        case 6:
            std::cout << "error: This maxtrix isn't a square matrix!\n" << std::endl;
            break;
        case 7:
            std::cout << "error: This maxtrix isn't a full rank matrix!\n" << std::endl;
            break;
        default:
            std::cout << "error: Undefined error!\n" << std::endl;
            break;
    }
}

Matrix Matrix::eyes(int order)
{
	Matrix temp(order, order);
	for (int i = 1; i <= order; i++) {
		temp(i, i) = 1;
	}
	return temp;
}

Matrix Matrix::zeros(int order)
{
	return Matrix(order, order);
}

// TODO
Matrix Matrix::linspace(int left, int right, int num)
{
	return Matrix();
}
