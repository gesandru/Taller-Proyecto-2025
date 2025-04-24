#include "..\include\matrix.h"

Matrix::Matrix() {
    this->n_row = 0;
	this->n_column = 0;
	this->data = nullptr;
}

Matrix::Matrix(const int n) {
    if (n <= 0) {
		cout << "Matrix create: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	this->n_row = n;
	this->n_column = n;
	this->data = (double **) malloc(n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < n; i++) {
		this->data[i] = (double *) malloc(n_column*sizeof(double));
	}
}

Matrix::Matrix(const int n_row, const int n_column) {
    if (n_row <= 0 || n_column <= 0) {
		cout << "Matrix create: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	this->n_row = n_row;
	this->n_column = n_column;
	this->data = (double **) malloc(n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < n_row; i++) {
		this->data[i] = (double *) malloc(n_column*sizeof(double));
	}
}

double& Matrix::operator () (const int row, const int column) {
	if (row <= 0 || row > this->n_row || column <= 0 || column > this->n_column) {
		cout << "Matrix get: error in row/column\n";
        exit(EXIT_FAILURE);
	}
	
	return this->data[row - 1][column - 1];
}

Matrix& Matrix::operator + (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sum: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) + m(i,j);
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator - (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sub: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) - m(i,j);
		}
	}
	
	return *m_aux;
}

ostream& operator << (ostream &o, Matrix &m) {
	for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++)
			printf("%5.20lf ", m(i,j));
        o << "\n";
    }
	
    return o;
}

Matrix& zeros(const int n) {
	Matrix *m_aux = new Matrix(n);
	
	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= n; j++) {
			(*m_aux)(i,j) = 0;
		}
	}
	
	return (*m_aux);
}

Matrix& zeros(const int n_row, const int n_column) {
	Matrix *m_aux = new Matrix(n_row, n_column);
	
	for(int i = 1; i <= n_row; i++) {
		for(int j = 1; j <= n_column; j++) {
			(*m_aux)(i,j) = 0;
		}
	}
	
	return (*m_aux);
}

Matrix& eye(const int n){
	Matrix *m_aux = new Matrix(n, n);
	
	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= n; j++) {
			if(i==j){
				(*m_aux)(i,j) = 1;
			}else{
				(*m_aux)(i,j) = 0;
			}
		}
	}
	
	return (*m_aux);
}

void assign_row(Matrix &m, double *d, const int i){
	if (sizeof(d)>m.n_row || i>m.n_row) {
		cout << "error in assign_row\n";
        exit(EXIT_FAILURE);
	}
	
	for(int j = 1; j <=m.n_column; j++){
		m(i, j) = d[j-1];
	}
}

void assign_column(Matrix &m, double *d, const int j){
	if (sizeof(d)>m.n_column || j>m.n_column) {
		cout << "error in assign_column\n";
        exit(EXIT_FAILURE);
	}
	
	for(int i = 1; i <=m.n_column; i++){
		m(i, j) = d[i-1];
	}
}

double* extract_row(Matrix &m, const int i){
	if (i>m.n_row) {
		cout << "error in extract_row\n";
        exit(EXIT_FAILURE);
	}
	
	double *d = (double *) malloc(m.n_row*sizeof(double));
	
	for(int j = 1; j <=m.n_column; j++){
		d[j-1] = m(i, j);
	}
	
	return d;
}