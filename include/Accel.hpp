#ifndef _ACCEL_
#define _ACCEL_

#include "../include/Matrix.hpp"
#include "../include/const.hpp"
#include "../include/globals.hpp"
#include "../include/PrecMatrix.hpp"
#include "../include/NutMatrix.hpp"
#include "../include/PoleMatrix.hpp"
#include "../include/AccelHarmonic.hpp"
#include "../include/AccelPointMass.hpp"
#include "../include/IERS.hpp"
#include "../include/timediff.hpp"
#include "../include/GHAMatrix.hpp"
#include "../include/JPL_Eph_DE430.hpp"
#include "../include/Mjday_TDB.hpp"
#include <cmath>
#include <tuple>

        Matrix& Accel(double x, Matrix Y);

#endif