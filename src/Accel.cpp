#include "../include/Accel.hpp"


Matrix& Accel(double x, Matrix Y){

tuple<double, double, double, double, double, double, double, double, double> result = 
IERS(eopdata,AuxParam.Mjd_UTC + x/86400.0,'l');

double x_pole = get<0>(result);
double y_pole = get<1>(result);
double UT1_UTC = get<2>(result);
double LOD = get<3>(result);
double dpsi = get<4>(result);
double deps = get<5>(result);
double dx_pole = get<6>(result);
double dy_pole = get<7>(result);
double TAI_UTC = get<8>(result);

tuple<double, double, double, double, double> result2 = timediff(UT1_UTC,TAI_UTC);

double UT1_TAI = get<0>(result2);
double UTC_GPS = get<1>(result2);
double UT1_GPS = get<2>(result2);
double TT_UTC = get<3>(result2);
double GPS_UTC = get<4>(result2);

double Mjd_UT1 = AuxParam.Mjd_UTC + x/86400.0 + UT1_UTC/86400.0;
double Mjd_TT = AuxParam.Mjd_UTC + x/86400.0 + TT_UTC/86400.0;

Matrix P = PrecMatrix(Const::MJD_J2000,Mjd_TT);
Matrix N = NutMatrix(Mjd_TT);
Matrix T = N * P;
Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

double MJD_TDB = Mjday_TDB(Mjd_TT);

 
tuple<Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,
          Matrix,Matrix,Matrix,Matrix> result3 = JPL_Eph_DE430(MJD_TDB);
 
Matrix r_Mercury = get<0>(result3);
Matrix r_Venus = get<1>(result3);
Matrix r_Earth = get<2>(result3);
Matrix r_Mars = get<3>(result3);
Matrix r_Jupiter = get<4>(result3);
Matrix r_Saturn = get<5>(result3);
Matrix r_Uranus = get<6>(result3);
Matrix r_Neptune = get<7>(result3);
Matrix r_Pluto = get<8>(result3);
Matrix r_Moon = get<9>(result3);
Matrix r_Sun = get<10>(result3);

//Acceleration due to harmonic gravity field
Matrix a = transpose(AccelHarmonic(transpose(extract_vector(Y,1,3)), E, AuxParam.n, AuxParam.m));

//Luni-solar perturbations
if (AuxParam.sun){
    a = a + AccelPointMass(extract_vector(Y,1,3),r_Sun,Const::GM_Sun);
}

if (AuxParam.moon){
    a = a + AccelPointMass(extract_vector(Y,1,3),r_Moon,Const::GM_Moon);
}

//Planetary perturbations
if (AuxParam.planets){
    a = a + AccelPointMass(extract_vector(Y,1,3),r_Mercury,Const::GM_Mercury);
    a = a + AccelPointMass(extract_vector(Y,1,3),r_Venus,Const::GM_Venus);
    a = a + AccelPointMass(extract_vector(Y,1,3),r_Mars,Const::GM_Mars);
    a = a + AccelPointMass(extract_vector(Y,1,3),r_Jupiter,Const::GM_Jupiter);
    a = a + AccelPointMass(extract_vector(Y,1,3),r_Saturn,Const::GM_Saturn);
    a = a + AccelPointMass(extract_vector(Y,1,3),r_Uranus,Const::GM_Uranus);    
    a = a + AccelPointMass(extract_vector(Y,1,3),r_Neptune,Const::GM_Neptune);
    a = a + AccelPointMass(extract_vector(Y,1,3),r_Pluto,Const::GM_Pluto);
}

return union_vector(extract_vector(Y,4,6),a);
}