#ifndef _CHEB3D_
#define _CHEB3D_

#include "../include/Matrix.hpp"
#include <cmath>
#include <cstdio>

Matrix& Cheb3D(const double t,const int N,const double Ta,const double Tb, Matrix &Cx, Matrix &Cy, Matrix &Cz);

#endif