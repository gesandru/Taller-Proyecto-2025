#ifndef _MEASUPDATE_
#define _MEASUPDATE_

#include "../include/Matrix.hpp"
#include <cmath>
#include <tuple>

        tuple<Matrix&, Matrix&, Matrix&> MeasUpdate( Matrix &x, double z, double g, double s, Matrix G, Matrix &P, double n);

#endif