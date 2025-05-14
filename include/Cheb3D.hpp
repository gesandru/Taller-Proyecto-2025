#ifndef _CHEB3D_
#define _CHEB3D_

#include "../include/Matrix.hpp"
#include <cmath>
#include <cstdio>

Matrix Cheb3D(const int t,const int N,const int Ta,const int Tb,Matrix &Cx,Matrix &Cy,Matrix &Cz);

#endif