#ifndef _ACCELHARMONIC_
#define _ACCELHARMONIC_

#include "../include/Matrix.hpp"
#include "../include/Legendre.hpp"
#include "../include/globals.hpp"
#include <cmath>
#include <tuple>

Matrix& AccelHarmonic(Matrix &r, Matrix &E, const int n_max, const int m_max);

#endif