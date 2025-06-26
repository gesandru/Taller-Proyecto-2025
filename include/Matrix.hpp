/**
 * @file Matrix.hpp
 * @brief Clase Matrix
 *
 * @author George
 */

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
	/**
	 * @brief constructor de una matríz vacía
	 *
	 */
    Matrix();
	
	/**
	 * @brief constructor de una matríz de una sola fila
	 *
	 * @param n número de columnas de la matríz
	 */
	Matrix(const int n);
	
	/**
	 * @brief constructor de una matríz
	 *
	 * @param n_row número de filas de la matríz
	 * @param n_column número de columnas de la matríz
	 */
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
/**
 * @brief Construye una matríz de ceros de una sola fila
 *
 * @param n número de columnas de la matríz
 */
Matrix& zeros(const int n);

/**
 * @brief Construye una matríz de ceros
 *
 * @param n_row número de filas de la matríz
 * @param n_column número de columnas de la matríz
 */
Matrix& zeros(const int n_row, const int n_column);

/**
 * @brief Construye una matríz identidad
 *
 * @param n tamaño de la matríz identidad
 */
Matrix& eye(const int n);

/**
 * @brief asigna una fila de una matríz a otra
 *
 * @param m1 matríz a la que se le asigna la fila
 * @param m2 matríz de una sola fila que se toma
 * @param i número de fila en la que insertar
 */
void assign_row(Matrix &m1, Matrix &m2, const int i);

/**
 * @brief asigna una columna de una matríz a otra
 *
 * @param m1 matríz a la que se le asigna la columna
 * @param m2 matríz de una sola fila que se toma
 * @param j número de columna en la que insertar
 */
void assign_column(Matrix &m1, Matrix &m2, const int j);

/**
 * @brief extrae una fila de una matríz
 *
 * @param m matríz de la que se le saca la fila
 * @param i número de fila a extraer
 */
Matrix& extract_row(Matrix &m, const int i);

/**
 * @brief extrae una columna de una matríz
 *
 * @param m matríz de la que se le saca la fila
 * @param j número de columna a extraer
 */
Matrix& extract_column(Matrix &m, const int j);

/**
 * @brief devuelve la inversa de una matríz
 *
 * @param m matríz a invertir
 */
Matrix& inv(Matrix &m);

/**
 * @brief devuelve la traspuesta de una matríz
 *
 * @param m matríz a transponer
 */
Matrix& transpose(Matrix &m);

/**
 * @brief devuelve el producto escalar entre dos matrices
 *
 * @param m1 primera matríz
 * @param m2 segunda matríz
 */
double dot(Matrix &m1, Matrix &m2);

/**
 * @brief devuelve la norma de una matriz
 *
 * @param m matríz de la que hacer la norma
 */
double norm(Matrix &m);

/**
 * @brief devuelve una cierta cantidad de elementos de una matríz en forma de otra matríz
 *
 * @param m matríz de la que extraer un vector
 * @param n elemento del cual empezar a tomar
 * @param o elemento donde parar de tomar
 */
Matrix& extract_vector(Matrix &m, int n, int o);

/**
 * @brief devuelve una matríz de una sola fila como la unión de dos matrices de una sola fila
 *
 * @param m1 primera matríz a unir
 * @param m2 segunda matríz a unir
 */
Matrix& union_vector(Matrix &m1, Matrix &m2);

#endif