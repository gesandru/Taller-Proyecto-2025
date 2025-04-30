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
	
	this->n_row = 1;
	this->n_column = n;
	this->data = (double **) malloc(n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	this->data[0] = (double *) malloc(n_column*sizeof(double));
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

double& Matrix::operator () (const int n) {
	if (n <= 0 || this->n_column <= 0 || n > this->n_column) {
		cout << "Matrix get: error in column\n";
        exit(EXIT_FAILURE);
	}
	
	return this->data[(n-1)/this->n_column][(n-1)%this->n_column];
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

Matrix& Matrix::operator * (const double d){
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) * d;
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator * (Matrix &m) {
	if (this->n_column != m.n_row) {
		cout << "Matrix mult: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, m.n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= m.n_column; j++) {
			(*m_aux)(i,j) = 0.0;
			for (int k = 1; k <= this->n_column; k++){
				
				(*m_aux)(i,j) = (*m_aux)(i,j) + (*this)(i,k) * m(k,j);
			}
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator / (Matrix &m){
	Matrix *m_aux = new Matrix(3,3);
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
	for(int j = 1; j <= n; j++) {
		(*m_aux)(1,j) = 0;
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

void assign_row(Matrix &m1, Matrix &m2, const int i){
	if (m2.n_column!=m1.n_column || i>m1.n_row) {
		cout << "error in assign_row\n";
        exit(EXIT_FAILURE);
	}
	
	for(int j = 1; j <=m1.n_column; j++){
		m1(i, j) = m2(1,j);
	}
}

void assign_column(Matrix &m1, Matrix &m2, const int j){
	if (m2.n_row!=m1.n_row || j>m1.n_column) {
		cout << "error in assign_column\n";
        exit(EXIT_FAILURE);
	}
	
	for(int i = 1; i <=m1.n_column; i++){
		m1(i, j) = m2(1,i);
	}
}

Matrix& extract_row(Matrix &m1, const int i){
	if (i>m1.n_row) {
		cout << "error in extract_row\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m2 = new Matrix(m1.n_column);
	
	for(int j = 1; j <=m1.n_column; j++){
		(*m2)(1,j) = m1(i,j);
	}
	
	return (*m2);
}

//Sacado de aquí: https://www.geeksforgeeks.org/finding-inverse-of-a-matrix-using-gauss-jordan-method/
Matrix& inv(Matrix &m)
{
	if (m.n_row!=m.n_column) {
		cout << "error in inv\n";
        exit(EXIT_FAILURE);
	}
	
	int order = 2*m.n_row;
	
	Matrix *m_aug = new Matrix(order, m.n_column);
	
	Matrix m_eye = eye(m.n_row);
	
	// Create the augmented matrix
    // Add the identity matrix
    // of order at the end of original matrix.
    for (int i = 1; i <= order; i++) {
 
        for (int j = 1; j <= 2 * order; j++) {
 
            // Add '1' at the diagonal places of
            // the matrix to create a identity matrix
            if (j == (i + order))
                (*m_aug)(i,j) = 1;
        }
    }
 
    // Interchange the row of matrix,
    // interchanging of row will start from the last row
    for (int i = order ; i > 0; i--) {
 
        // Swapping each and every element of the two rows
        // if (matrix[i - 1][0] < matrix[i][0])
        // for (int j = 0; j < 2 * order; j++) {
        //
        //        // Swapping of the row, if above
        //        // condition satisfied.
        // temp = matrix[i][j];
        // matrix[i][j] = matrix[i - 1][j];
        // matrix[i - 1][j] = temp;
        //    }
 
        // Directly swapping the rows using pointers saves
        // time
 
        if (matrix[i - 1][0] < matrix[i][0]) {
            float* temp = matrix[i];
            matrix[i] = matrix[i - 1];
            matrix[i - 1] = temp;
        }
    }
 
    // Print matrix after interchange operations.
    printf("\n=== Augmented Matrix ===\n");
    PrintMatrix(matrix, order, order * 2);
 
    // Replace a row by sum of itself and a
    // constant multiple of another row of the matrix
    for (int i = 0; i < order; i++) {
 
        for (int j = 0; j < order; j++) {
 
            if (j != i) {
 
                temp = matrix[j][i] / matrix[i][i];
                for (int k = 0; k < 2 * order; k++) {
 
                    matrix[j][k] -= matrix[i][k] * temp;
                }
            }
        }
    }
 
    // Multiply each row by a nonzero integer.
    // Divide row element by the diagonal element
    for (int i = 0; i < order; i++) {
 
        temp = matrix[i][i];
        for (int j = 0; j < 2 * order; j++) {
 
            matrix[i][j] = matrix[i][j] / temp;
        }
    }
 
    // print the resultant Inverse matrix.
    printf("\n=== Inverse Matrix ===\n");
    PrintInverse(matrix, order, 2 * order);
 
    return;
}

Matrix& transpose(Matrix &m){
	
	Matrix *m2 = new Matrix(m.n_column, m.n_row);
	
	for(int i = 1; i <= m.n_row; i++) {
        for(int j = 1; j <= m.n_column; j++) {
            (*m2)(j,i) = m(i,j);
        }
    }
	
	return (*m2);
}