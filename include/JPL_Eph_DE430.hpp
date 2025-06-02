#ifndef _JPL_
#define _JPL_

#include "../include/Matrix.hpp"
#include "../include/globals.hpp"
#include "../include/Cheb3D.hpp"
#include <cmath>
#include <tuple>

    tuple<Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,
          Matrix&,Matrix&,Matrix&,Matrix&> JPL_Eph_DE430(double Mjd_TDB);

#endif