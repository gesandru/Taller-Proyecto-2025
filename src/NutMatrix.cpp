#include "../include/NutMatrix.hpp"

Matrix NutMatrix (double Mjd_TT){

//Mean obliquity of the ecliptic
double eps = MeanObliquity(Mjd_TT);

//Nutation in longitude and obliquity
tuple<double, double> result = NutAngles(Mjd_TT);

double dpsi = get<0>(result);
double deps = get<1>(result);

//Transformation from mean to true equator and equinox
return (R_x(-eps-deps)*R_z(-dpsi)*R_x(+eps));
}