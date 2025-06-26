#ifndef _DEINTEG_
#define _DEINTEG_

#include "../include/Matrix.hpp"
#include "../include/sign_.hpp"
#include <cmath>
#include <limits>

    Matrix& DEInteg(Matrix& (*func)(double, Matrix),double t,double tout,double relerr,
                          double abserr,int n_eqn,Matrix y);

#endif