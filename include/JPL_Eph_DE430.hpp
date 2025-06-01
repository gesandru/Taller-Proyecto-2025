#ifndef _JPL_
#define _JPL_

#include "../include/Matrix.hpp"
#include <cmath>
#include <tuple>

    tuple<Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,
          Matrix&,Matrix&,Matrix&,Matrix&> JPL_Eph_DE430(Mjd_TDB);

#endif