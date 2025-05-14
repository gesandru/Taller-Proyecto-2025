#ifndef _TIMEUPDATE_
#define _TIMEUPDATE_

#include "../include/Matrix.hpp"
#include <cmath>

void TimeUpdate(Matrix &P, Matrix Phi, double Qdt = 0.0);

#endif