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
	double& operator () (const int n);
	double& operator () (const int row, const int column);
	Matrix& operator + (const double d);
	Matrix& operator + (Matrix &m);
	Matrix& operator - (const double d);
	Matrix& operator - (Matrix &m);
	Matrix& operator * (const double d);
	Matrix& operator * (Matrix &m);
	Matrix& operator / (const double d);
	Matrix operator / (Matrix &m);
	void operator = (Matrix &m);
	
	// Non-member operators
	friend ostream& operator << (ostream &o, Matrix &m);
	friend Matrix& operator - (Matrix &m);
};

// Operator overloading
ostream& operator << (ostream &o, Matrix &m);
Matrix& operator - (Matrix &m);

// Methods
Matrix& zeros(const int n);
Matrix& zeros(const int n_row, const int n_column);
Matrix& eye(const int n);
void assign_row(Matrix &m1, Matrix &m2, const int i);
void assign_column(Matrix &m1, Matrix &m2, const int j);
Matrix& extract_row(Matrix &m, const int i);
Matrix& extract_column(Matrix &m, const int j);
Matrix& inv(Matrix &m);
Matrix& transpose(Matrix &m);
double dot(Matrix &m1, Matrix &m2);
double norm(Matrix &m);
Matrix& extract_vector(Matrix &m, int n, int o);
Matrix& union_vector(Matrix &m1, Matrix &m2);

#endif