/**
 * @file tests_matrix.cpp
 * @brief Archivo de tests unitarios de las funciones de Matrix
 *
 * @author George
 */

#include "..\include\Matrix.hpp"
#include <cstdio>
#include <cmath>

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

int m_equals(Matrix A, Matrix B, double p) {
	if (A.n_row != B.n_row || A.n_column != B.n_column)
		return 0;
	else
		for(int i = 1; i <= A.n_row; i++)
			for(int j = 1; j <= A.n_column; j++)
				if(fabs(A(i,j)-B(i,j)) > p) {
					printf("%2.20lf %2.20lf\n",A(i,j),B(i,j));
					return 0;
				}
	
	return 1;
}

int m_sum_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = 2; C(1,2) =  2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = 8; C(2,2) = -3; C(2,3) = 1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = -2; C(3,3) = 0; C(3,4) = 7;
	
	Matrix R = A + B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_sum_02() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix C(f, c);
	C(1,1) = 1.5; C(1,2) =  3.5; C(1,3) = 9.5; C(1,4) = 1.5;
	C(2,1) = 2.5; C(2,2) = 0.5; C(2,3) = 1.5; C(2,4) = 1.5;
	C(3,1) = 1.5; C(3,2) = 2.5; C(3,3) = 1.5; C(3,4) = 6.5;
	
	Matrix R = A + 1.5;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_sub_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) = 0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = -2; C(1,2) = 2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = -6; C(2,2) = 1; C(2,3) = -1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = 4; C(3,3) = 0; C(3,4) = 3;
	
	Matrix R = A - B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_sub_02() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix C(f, c);
	C(1,1) = -1.5; C(1,2) =  0.5; C(1,3) = 6.5; C(1,4) = -1.5;
	C(2,1) = -0.5; C(2,2) = -2.5; C(2,3) = -1.5; C(2,4) = -1.5;
	C(3,1) = -1.5; C(3,2) = -0.5; C(3,3) = -1.5; C(3,4) = 3.5;
	
	Matrix R = A - 1.5;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_sub_03() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = -4; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = -9; A(3,4) = 5;
	
	Matrix B = -A;
	
	Matrix C(f, c);
	C(1,1) = 4; C(1,2) =  -2; C(1,3) = -8; C(1,4) = 0;
	C(2,1) = -1; C(2,2) = 1; C(2,3) = 0; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = -1; C(3,3) = 9; C(3,4) = -5;
    
    _assert(m_equals(B, C, 1e-10));
    
    return 0;
}

int m_mult_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 0; B(1,2) = 6; B(1,3) = 24; B(1,4) = 0;
	B(2,1) = 3; B(2,2) = -3; B(2,3) = 0; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = 3; B(3,3) = 0; B(3,4) = 15;
	
	Matrix R = A * 3.0;
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int m_mult_02() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(4, 5);
	B(1,1) = 0; B(1,2) = 6; B(1,3) = 2; B(1,4) = 0; B(1,5) = 4;
	B(2,1) = 3; B(2,2) = -3; B(2,3) = 0; B(2,4) = 0; B(2,5) = 1;
	B(3,1) = 0; B(3,2) = 3; B(3,3) = 0; B(3,4) = 5; B(3,5) = 7;
	B(4,1) = 4; B(4,2) = 8; B(4,3) = 9; B(4,4) = 2; B(4,5) = 3;
	
	Matrix C(3, 5);
	C(1,1) = 6; C(1,2) = 18; C(1,3) = 0; C(1,4) = 40; C(1,5) = 58;
	C(2,1) = -3; C(2,2) = 9; C(2,3) = 2; C(2,4) = 0; C(2,5) = 3;
	C(3,1) = 23; C(3,2) = 37; C(3,3) = 45; C(3,4) = 10; C(3,5) = 16;
	
	Matrix R = A * B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_div_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 6; A(1,3) = 24; A(1,4) = 0;
	A(2,1) = 3; A(2,2) = -3; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 3; A(3,3) = 0; A(3,4) = 15;
	
	Matrix B(f, c);
	B(1,1) = 0; B(1,2) = 2; B(1,3) = 8; B(1,4) = 0;
	B(2,1) = 1; B(2,2) = -1; B(2,3) = 0; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = 1; B(3,3) = 0; B(3,4) = 5;
	
	Matrix R = A / 3.0;
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int m_div_02() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(4, 4);
	B(1,1) = 0; B(1,2) = 6; B(1,3) = 2; B(1,4) = 0;
	B(2,1) = 3; B(2,2) = -3; B(2,3) = 0; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = 3; B(3,3) = 0; B(3,4) = 5;
	B(4,1) = 4; B(4,2) = 8; B(4,3) = 9; B(4,4) = 2;
	
	Matrix C(3, 4);
	C(1,1) = -19.0/9.0; C(1,2) = -440.0/243.0; C(1,3) = -44.0/81.0; C(1,4) = 110.0/81.0;
	C(2,1) = 0; C(2,2) = 1.0/3.0; C(2,3) = 0; C(2,4) = 0;
	C(3,1) = -5.0/9.0; C(3,2) = -40.0/243.0; C(3,3) = 77.0/81.0; C(3,4) = 10.0/81.0;
	
	Matrix R = A / B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_equal_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix C(f, c);
	C(1,1) = 0; C(1,2) = 2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = 1; C(2,2) = -1; C(2,3) = 0; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = 1; C(3,3) = 0; C(3,4) = 5;
	
	Matrix B = A;
    
    _assert(m_equals(C, B, 1e-10));
    
    return 0;
}

int m_zeros_01() {
	int c = 4;
		
	Matrix A(1, c);
	A(1,1) = 0; A(1,2) = 0; A(1,3) = 0; A(1,4) = 0;
	Matrix B = zeros(4);
	
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_zeros_02() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 0; A(1,3) = 0; A(1,4) = 0;
	A(2,1) = 0; A(2,2) = 0; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 0; A(3,4) = 0;
	
	Matrix B = zeros(3, 4);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_eye_01() {
    int f = 3;
	
	Matrix A(f, f);
	A(1,1) = 1; A(1,2) = 0; A(1,3) = 0;
	A(2,1) = 0; A(2,2) = 1; A(2,3) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 1;
	
	Matrix B = eye(3);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_assign_row_01(){
	int f = 3;
	
	Matrix A(f, f);
	A(1,1) = 1; A(1,2) = 0; A(1,3) = 0;
	A(2,1) = 0; A(2,2) = 1; A(2,3) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 1;
	
	Matrix B(f);
	B(1) = 5;
	B(2) = 4;
	B(3) = 1;
	
	assign_row(A, B, 2);
	
	Matrix C(f, f);
	C(1,1) = 1; C(1,2) = 0; C(1,3) = 0;
	C(2,1) = 5; C(2,2) = 4; C(2,3) = 1;
	C(3,1) = 0; C(3,2) = 0; C(3,3) = 1;
	
    _assert(m_equals(A, C, 1e-10));
    
    return 0;
}

int m_assign_column_01(){
	int f = 3;
	
	Matrix A(f, f);
	A(1,1) = 1; A(1,2) = 0; A(1,3) = 0;
	A(2,1) = 0; A(2,2) = 1; A(2,3) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 1;
	
	Matrix B(f);
	B(1) = 5;
	B(2) = 4;
	B(3) = 1;
	
	assign_column(A, B, 2);
	
	Matrix C(f, f);
	C(1,1) = 1; C(1,2) = 5; C(1,3) = 0;
	C(2,1) = 0; C(2,2) = 4; C(2,3) = 0;
	C(3,1) = 0; C(3,2) = 1; C(3,3) = 1;
	
    _assert(m_equals(A, C, 1e-10));
    
    return 0;
}

int m_inv_01(){
   
   	Matrix A(3, 3);
	A(1,1) = 4; A(1,2) = 2; A(1,3) = 6;
	A(2,1) = 5; A(2,2) = 9; A(2,3) = 7;
	A(3,1) = 3; A(3,2) = 1; A(3,3) = 3;
	
	Matrix B = inv(A);
	
	
	Matrix C(3, 3);
	C(1,1) = -0.5; C(1,2) = 0; C(1,3) = 1;
	C(2,1) = -0.15; C(2,2) = 0.15; C(2,3) = -0.05;
	C(3,1) = 0.55; C(3,2) = -0.05; C(3,3) = -0.65;
	
    _assert(m_equals(B, C, 1e-10));
	
	return 0;
}

int m_transpose_01(){
	int f = 3;
	int c = 4;
	
	Matrix A(c, f);
	A(1,1) = 1; A(1,2) = 6; A(1,3) = 5;
	A(2,1) = 3; A(2,2) = 1; A(2,3) = 4;
	A(3,1) = 2; A(3,2) = 8; A(3,3) = 1;
	A(4,1) = 3; A(4,2) = 7; A(4,3) = 9;
	
	Matrix B = transpose(A);
	
	Matrix C(f, c);
	C(1,1) = 1; C(1,2) = 3; C(1,3) = 2; C(1,4) = 3;
	C(2,1) = 6; C(2,2) = 1; C(2,3) = 8; C(2,4) = 7;
	C(3,1) = 5; C(3,2) = 4; C(3,3) = 1; C(3,4) = 9;
	
	_assert(m_equals(B, C, 1e-10));
	
	return 0;
}

int m_dot_01() {
    int f = 4;
	
	Matrix A(f);
	A(1,1) = 3.6; A(1,2) =  2.0; A(1,3) = 8.0; A(1,4) = 2.0;
	
	Matrix B(f);
	B(1,1) = 2.0; B(1,2) =  6.0; B(1,3) = 1.0; B(1,4) = 9.0;
	
	double C = 45.2;
	
	double R = dot(A,B);
	
	int i;
	if(fabs(C-R) > 1e-10) {
		printf("%2.20lf %2.20lf\n",C,R);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
    
    return 0;
}

int m_norm_01() {
    int f = 4;
	
	Matrix A(f);
	A(1,1) = 1.0; A(1,2) =  -2.0; A(1,3) = 3.0;
	
	double C = 3.7417;
	
	double R = norm(A);
	
	int i;
	if(fabs(C-R) > 1e-4) {
		printf("%2.20lf %2.20lf\n",C,R);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
    
    return 0;
}

int m_extract_vector_01(){
	
	Matrix A(4);
	A(1) = 1; A(2) = 6; A(3) = 5; A(4) = 3;
	
	Matrix B = extract_vector(A, 2, 3);
	
	Matrix C(2);
	C(1) = 6; C(2) = 5;
	
	_assert(m_equals(B, C, 1e-10));
	
	return 0;
}

int m_union_vector_01(){
	
	Matrix A(3);
	A(1) = 1; A(2) = 6; A(3) = 5;
	
	Matrix B(3);
	B(1) = 3; B(2) = 1; B(3) = 4;
	
	Matrix C(6);
	C(1) = 1; C(2) = 6; C(3) = 5; C(4) = 3; C(5) = 1; C(6) = 4;
	
	Matrix D = union_vector(A, B);
	
	_assert(m_equals(C, D, 1e-10));
	
	return 0;
}

int all_tests()
{
    _verify(m_sum_01);
	_verify(m_sum_02);
    _verify(m_sub_01);
	_verify(m_sub_02);
	_verify(m_sub_03);
	_verify(m_mult_01);
	_verify(m_mult_02);
	_verify(m_div_01);
	_verify(m_div_02);
	_verify(m_equal_01);
	_verify(m_zeros_01);
    _verify(m_zeros_02);
	_verify(m_eye_01);
	_verify(m_assign_row_01);
	_verify(m_assign_column_01);
	_verify(m_inv_01);
	_verify(m_transpose_01);
	_verify(m_dot_01);
	_verify(m_norm_01);
	_verify(m_extract_vector_01);
	_verify(m_union_vector_01);

    return 0;
}


int main()
{
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
