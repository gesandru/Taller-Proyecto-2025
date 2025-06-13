#include "../include/VarEqn.hpp"

Matrix  VarEqn(double x, Matrix yPhi){
	
tuple< double, double, double, double, double, double, double, double, double> result = IERS(eopdata,AuxParam.Mjd_UTC,'l');

double x_pole = get<0>(result);
double y_pole = get<1>(result);
double UT1_UTC = get<2>(result);
double LOD = get<3>(result);
double dpsi = get<4>(result);
double deps = get<5>(result);
double dx_pole = get<6>(result);
double dy_pole = get<7>(result);
double TAI_UTC = get<8>(result);

tuple<double,double,double,double,double>result2 = timediff(UT1_UTC, TAI_UTC);

double UT1_TAI = get<0>(result2);
double UTC_GPS = get<1>(result2);
double UT1_GPS = get<2>(result2);
double TT_UTC = get<3>(result2);
double GPS_UTC = get<4>(result2);

double Mjd_UT1 = AuxParam.Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

// Transformation matrix
Matrix P = PrecMatrix(Const::MJD_J2000,AuxParam.Mjd_TT + x/86400.0);
Matrix N = NutMatrix(AuxParam.Mjd_TT + x/86400.0);
Matrix T = N * P;
Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

// State vector components
Matrix r = transpose(extract_vector(yPhi,1,3));
Matrix v = extract_vector(yPhi,4,6);
Matrix Phi = zeros(6,6);

// State transition matrix
for(int i = 1; i<=6;i++){
    for(int j = 1; j<=6;j++){
		Phi(j,i) = yPhi(6*i+j);
	}
}

// Acceleration and gradient
Matrix a = AccelHarmonic( r, E, AuxParam.n, AuxParam.m );

Matrix G = G_AccelHarmonic( r, E, AuxParam.n, AuxParam.m );

// Time derivative of state transition matrix
Matrix yPhip = zeros(42);
Matrix dfdy = zeros(6, 6);

for(int i = 1; i<=3;i++){
    for(int j = 1; j<=3;j++){
        dfdy(i,j) = 0.0;                 // dv/dr(i,j)
        dfdy(i+3,j) = G(i,j);            // da/dr(i,j)
        if(i==j){
            dfdy(i,j+3) = 1.0;
        }else{
            dfdy(i,j+3) = 0.0;             // dv/dv(i,j)
        }
        dfdy(i+3,j+3) = 0.0;             // da/dv(i,j)
    }
}

Matrix Phip = dfdy*Phi;

// Derivative of combined state vector and state transition matrix
for(int i = 1; i<=3;i++){
    yPhip(i)   = v(i);                 // dr/dt(i)
    yPhip(i+3) = a(i,1);                 // dv/dt(i)
}
for(int i = 1; i<=6;i++){
    for(int j = 1; j<=6;j++){
        yPhip(6*j+i) = Phip(i,j);     // dPhi/dt(i,j)
    }
}
return yPhip;
}
