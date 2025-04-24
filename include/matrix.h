#ifndef _MATRIX_
#define _MATRIX_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

class Matrix {
public:
    int n_row, n_column;
	double **data;

    // Parameterized constructor
    Matrix();
	Matrix(const int n);
    Matrix(const int n_row, const int n_column);
	
	// Member operators
	double& operator () (const int row, const int column);
	Matrix& operator + (Matrix &m);
	Matrix& operator - (Matrix &m);
	Matrix& operator * (double &d);
	Matrix& operator * (Matrix &m);
	
	// Non-member operators
	friend ostream& operator << (ostream &o, Matrix &m);
};

// Operator overloading
ostream& operator << (ostream &o, Matrix &m);

// Methods
Matrix& zeros(const int n);
Matrix& zeros(const int n_row, const int n_column);
Matrix& eye(const int n);
void assign_row(Matrix &m, double *d, const int i);
void assign_column(Matrix &m, double *d, const int j);
double* extract_row(Matrix &m, const int i);
double* extract_column(Matrix &m, const int j);

#endif