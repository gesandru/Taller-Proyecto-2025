#include "../include/VarEqn.hpp"
/*
Matrix  VarEqn(double x, Matrix yPhi){

double UT1_UTC;
double TAI_UTC;
double UT1_TAI;
double UTC_GPS;
double UT1_GPS;
double TT_UTC;
double GPS_UTC;


[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,AuxParam.Mjd_UTC,'l');

Timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
double Mjd_UT1 = AuxParam.Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

// Transformation matrix
Matrix P = PrecMatrix(Const::MJD_J2000,AuxParam.Mjd_TT + x/86400.0);
Matrix N = NutMatrix(AuxParam.Mjd_TT + x/86400.0);
Matrix T = N * P;
Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

// State vector components
Matrix r = union_vector(yPhi,1,3);
Matrix v = union_vector(yPhi,4,6);
Matrix Phi = zeros(6,6);

// State transition matrix
for(int i = 1; i<=i;j++){
    for(int j = 1; j<=6;j++){
    Phi(i,j) = yPhi(6*i+j);
}

// Acceleration and gradient
Matrix a = AccelHarmonic( r, E, AuxParam.n, AuxParam.m );
Matrix G = G_AccelHarmonic( r, E, AuxParam.n, AuxParam.m );

// Time derivative of state transition matrix

//esto no sé si dejarlo como matriz o cómo vector
Matrix yPhip = zeros(42);
//Este seguro que es matriz
Matrix dfdy = zeros(6, 6);

for(int i = 1; i<=3;i++){
    for(int j = 1; j<=3;j++){
        dfdy(i,j) = 0.0;                 // dv/dr(i,j)
        dfdy(i+3,j) = G(i,j);            // da/dr(i,j)
        if(i==j){
            dfdy(i,j+3) = 1;
        }else{
            dfdy(i,j+3) = 0;             // dv/dv(i,j)
        }
        dfdy(i+3,j+3) = 0.0;             // da/dv(i,j)
    }
}

Matrix Phip = dfdy*Phi;

// Derivative of combined state vector and state transition matrix
for(int i = 1; i<=3;i++){
    yPhip(i)   = v(i);                 // dr/dt(i)
    yPhip(i+3) = a(i);                 // dv/dt(i)
}

for(int i = 1; i<=6;i++){
    for(int j = 1; j<=6;j++){
        yPhip(6*j+i) = Phip(i,j);     // dPhi/dt(i,j)
    }
}
  
return yPhip;
}
*/