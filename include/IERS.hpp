#ifndef _IERS_
#define _IERS_

#include "../include/Matrix.hpp"
#include "../include/const.hpp"
#include "../include/globals.hpp"
#include <cmath>
#include <tuple>

tuple< double, double, double, double, double, double, double, double, double> IERS(Matrix &eop,const double Mjd_UTC,const char interp = 'n');

#endif