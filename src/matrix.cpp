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

Matrix& Matrix::operator + (const double d) {
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) + d;
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator - (const double d) {
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) - d;
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

Matrix& Matrix::operator / (const double d){
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) / d;
		}
	}
	
	return *m_aux;
}

Matrix Matrix::operator / (Matrix &m){
	Matrix m_aux = inv(m);
	Matrix m_result = (*this) * m_aux;
	
	return m_result;
}

void Matrix::operator = (Matrix &m) {
	
	this->n_row = m.n_row;
	this->n_column = m.n_column;
	
	this->data = (double **) malloc(m.n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix =: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < n_row; i++) {
		this->data[i] = (double *) malloc(m.n_column*sizeof(double));
	}
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*this)(i,j) = m(i,j);
		}
	}
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
	if (m2.n_column!=m1.n_row || j>m1.n_column) {
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

Matrix& extract_column(Matrix &m1, const int j){
	if (j>m1.n_column) {
		cout << "error in extract_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m2 = new Matrix(m1.n_row);
	for(int i = 1; i <=m1.n_row; i++){
		(*m2)(1,i) = m1(i,j);
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
	
	int order = m.n_row;
	Matrix *m_aug = new Matrix(order, 2*order);
	Matrix temp;
	double d;
	
	// Create the augmented matrix
	// Add the identity matrix
    // of order at the end of original matrix.
    for (int i = 1; i <= order; i++) {
 
        for (int j = 1; j <= 2 * order; j++) {
 
            // Add '1' at the diagonal places of
            // the matrix to create a identity matrix
            if (j == (i + order)){
                (*m_aug)(i,j) = 1;
			} else{
				(*m_aug)(i,j) = 0;
			}
		}
    }

	for (int i = 1; i <= order; i++) {
 
        for (int j = 1; j <= order; j++) {
            (*m_aug)(i,j) = m(i,j);
        }
    }
	

    // Interchange the row of matrix,
    // interchanging of row will start from the last row
    for (int i = order ; i > 1; i--) {
 
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
 
        if ((*m_aug)(i-1,1) < (*m_aug)(i,1)) {
            temp = extract_row((*m_aug),i);
            assign_row((*m_aug), extract_row((*m_aug),i-1),i);
            assign_row((*m_aug), temp,i-1);
        }
    }
 
    // Replace a row by sum of itself and a
    // constant multiple of another row of the matrix
    for (int i = 1; i <= order; i++) {
 
        for (int j = 1; j <= order; j++) {
 
            if (j != i) {
 
                d = (*m_aug)(j,i) / (*m_aug)(i,i);
                for (int k = 1; k <= 2 * order; k++) {
 
                    (*m_aug)(j,k) -= (*m_aug)(i,k) * d;
                }
            }
        }
    }
 
    // Multiply each row by a nonzero integer.
    // Divide row element by the diagonal element
    for (int i = 1; i <= order; i++) {
 
        d = (*m_aug)(i,i);
        for (int j = 1; j <= 2 * order; j++) {
 
            (*m_aug)(i,j) = (*m_aug)(i,j) / d;
        }
    }
	
	Matrix *m_result = new Matrix(order,order);
	for(int i = 1; i <= order; i++){
		assign_column((*m_result), extract_column((*m_aug),i+order),i);
	}
 
    return (*m_result);
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

double dot(Matrix &m1, Matrix &m2) {
	if (m1.n_row!=m2.n_row || m1.n_column!=m2.n_column) {
		cout << "error in dot\n";
        exit(EXIT_FAILURE);
	}
	
    double result = 0.0;

    for (int j = 1; j <= m1.n_column; j++){
            result += (m1(1,j) * m2(1,j));
    }
    return result;
}

double norm(Matrix &m){
	if (m.n_row>1) {
		cout << "error in norm\n";
        exit(EXIT_FAILURE);
	}
	
    double result = 0.0;

    for (int j = 1; j <= m.n_column; j++){
            result += pow(m(1,j),2);
    }
	result = sqrt(result);
    return result;
}