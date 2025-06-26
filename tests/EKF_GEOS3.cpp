/**
 * @file EKF_GEOS3.cpp
 * @brief Archivo principal y de test completo de la aplicación
 *
 * @author George
 */


#include "../include/Matrix.hpp"
#include "../include/globals.hpp"
#include "../include/const.hpp"
#include "../include/DEInteg.hpp"
#include "../include/VarEqn.hpp"
#include "../include/Accel.hpp"
#include "../include/LTC.hpp"
#include "../include/Position.hpp"
#include "../include/IERS.hpp"
#include "../include/DEInteg.hpp"
#include "../include/timediff.hpp"
#include "../include/TimeUpdate.hpp"
#include "../include/AzElPa.hpp"
#include "../include/MeasUpdate.hpp"
#include <cmath>
#include <cstdio>
#include <tuple>

using namespace std;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)


/**
 * @brief Comprueba que los resultados sean correctos
 *
 * @param Y0 primera matríz resultado
 * @param Y_true Segunda matriz resultado
 */
int test(Matrix Y0, Matrix Y_true)
{
	int i;
	
    if(fabs(Y0(1)-Y_true(1)-39.0) > 1e-1) {
		printf("%2.20lf %2.20lf\n",Y0(1)-Y_true(1),39.0);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
	    if(fabs(Y0(2)-Y_true(2)-40.9) > 1e-1) {
		printf("%2.20lf %2.20lf\n",Y0(1)-Y_true(1),40.9);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
	    if(fabs(Y0(3)-Y_true(3)-(-14.6)) > 1e-1) {
		printf("%2.20lf %2.20lf\n",Y0(1)-Y_true(1),-14.6);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
	    if(fabs(Y0(4)-Y_true(4)-2.1) > 1e-1) {
		printf("%2.20lf %2.20lf\n",Y0(1)-Y_true(1),2.1);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
	    if(fabs(Y0(5)-Y_true(5)-(-2.4)) > 1e-1) {
		printf("%2.20lf %2.20lf\n",Y0(1)-Y_true(1),-2.4);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
	    if(fabs(Y0(6)-Y_true(6)-2.1) > 1e-1) {
		printf("%2.20lf %2.20lf\n",Y0(1)-Y_true(1),2.1);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
    return 0;
}

int main(){

DE430Coeff();

GGMO3S();

// Model parameters
AParam();

// read Earth orientation parameters
eop19620101(21413);

int nobs = 46;

GEOS3(nobs);

double sigma_range = 92.5;          // [m]
double sigma_az = 0.0224*Const::Rad; // [rad]
double sigma_el = 0.0139*Const::Rad; // [rad]

// Kaena Point station
double lat = Const::Rad*21.5748;     // [rad]
double lon = Const::Rad*(-158.2706); // [rad]
double alt = 300.20;                // [m]

Matrix Rs = transpose(Position(lon, lat, alt));

double Mjd1 = obs(1,1);
double Mjd2 = obs(9,1);
double Mjd3 = obs(18,1);

Matrix &r2 = zeros(3), &v2 = zeros(3);
r2(1) = 6221397.62857869; r2(2) = 2867713.77965738; r2(3) = 3006155.98509949;
v2(1) = 4645.04725161806; v2(2) = -2752.21591588204; v2(3) = -7507.99940987031;

//[r2,v2] = anglesg(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),
//                  Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);
// [r2,v2] = anglesdr(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),...
//                    Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);

Matrix Y0_apr = union_vector(r2,v2);

double Mjd0 = Mjday(1995,1,29,02,38,0);

double Mjd_UTC = obs(9,1);

AuxParam.Mjd_UTC = Mjd_UTC;
AuxParam.n      = 20;
AuxParam.m      = 20;
AuxParam.sun     = 1;
AuxParam.moon    = 1;
AuxParam.planets = 1;

int n_eqn  = 6;

Matrix Y = DEInteg(&Accel,0,-(obs(9,1)-Mjd0)*86400.0,1e-13,1e-6,6,Y0_apr);

Matrix P = zeros(6,6);
  
for(int i=1; i<=3; i++){
    P(i,i)=1e8;
}
for(int i=4; i<=6; i++){
    P(i,i)=1e3;
}

Matrix LT = LTC(lon,lat);

Matrix yPhi = zeros(42);
Matrix yPhi1 = zeros(42);
Matrix Phi  = zeros(6,6);

// Measurement loop
double t = 0.0, t_old;
Matrix Y_old, U, r, s, dAds, dEds, dAdY, K, dEdY, dDds, dDdY;
double Dist, UT1_UTC, TAI_UTC, TT_UTC, Mjd_TT, Mjd_UT1,
	theta, Azim, Elev;
tuple< double, double, double, double, double, double, double, double, double> result1;
tuple <double, double, double, double, double> result2;
tuple<double, double, Matrix, Matrix> result3;
tuple<Matrix, Matrix, Matrix> result4;

for(int i=1; i<=nobs; i++){
    // Previous step
    t_old = t;
    Y_old = Y;
    
    // Time increment and propagation
    Mjd_UTC = obs(i,1);                       // Modified Julian Date
    t       = (Mjd_UTC-Mjd0)*86400.0;         // Time since epoch [s]
	
	result1 = IERS(eopdata,Mjd_UTC,'l');
	UT1_UTC = get<2>(result1);
	TAI_UTC = get<8>(result1);
	
	result2 = timediff(UT1_UTC,TAI_UTC);
	TT_UTC = get<3>(result2);
	
    Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;
        
    for(int ii = 1; ii<=6; ii++){
        yPhi1(ii) = Y_old(ii);
        for (int j = 1; j<=6; j++){
            if (ii==j){ 
                yPhi1(6*j+ii) = 1.0;
            }else{
                yPhi1(6*j+ii) = 0.0;
            }
        }
    }
    
    yPhi = DEInteg (&VarEqn,0,t-t_old,1e-13,1e-6,42,yPhi1);
    
    // Extract state transition matrices
    for (int j = 1; j<=6; j++){
		for (int jj = 1; jj<=6; jj++){
			Phi(jj,j) = yPhi(6*j+jj);
		}
    }
    
    Y = DEInteg (&Accel,0,t-t_old,1e-13,1e-6,6,Y_old);
    
    // Topocentric coordinates
    theta = gmst(Mjd_UT1);                    // Earth rotation
    U = R_z(theta);
    r = transpose(extract_vector(Y,1,3));
    s = transpose(LT*(U*r-Rs));                          // Topocentric position [m]
    
    // Time update
    TimeUpdate(P, Phi);
        
    // Azimuth and partials  
	result3 = AzElPa(s);					// Azimuth, Elevation
	Azim = get<0>(result3);
	dAds = get<2>(result3);

	
    dAdY = union_vector(dAds*LT*U,zeros(3));
    
    // Measurement update
	result4 = MeasUpdate ( Y, obs(i,2), Azim, sigma_az, dAdY, P, 6 );
	Y = get<1>(result4);
	P = get<2>(result4);
    
    // Elevation and partials
    r = transpose(extract_vector(Y,1,3));
    s = transpose(LT*(U*r-Rs));                         // Topocentric position [m]
	
    result3 = AzElPa(s);					// Azimuth, Elevation
	Elev = get<1>(result3);
	dEds = get<3>(result3);
	
    dEdY = union_vector(dEds*LT*U,zeros(3));
    
    // Measurement update
    result4 = MeasUpdate ( Y, obs(i,3), Elev, sigma_el, dEdY, P, 6 );
	Y = get<1>(result4);
	P = get<2>(result4);
    
    // Range and partials
    r = transpose(extract_vector(Y,1,3));
    s = LT*(U*r-Rs);                       // Topocentric position [m]
    s = transpose(s);
    Dist = norm(s); dDds = s/Dist;         // Range
    dDdY = union_vector(dDds*LT*U,zeros(3));
    
    // Measurement update
    result4 = MeasUpdate ( Y, obs(i,4), Dist, sigma_range, dDdY, P, 6 );
	Y = get<1>(result4);
	P = get<2>(result4);
}

result1 = IERS(eopdata,obs(46,1),'l');
	UT1_UTC = get<2>(result1);
	TAI_UTC = get<8>(result1);
	
result2 = timediff(UT1_UTC,TAI_UTC);
	TT_UTC = get<3>(result2);
	
Mjd_TT = Mjd_UTC + TT_UTC/86400;
AuxParam.Mjd_UTC = Mjd_UTC;
AuxParam.Mjd_TT = Mjd_TT;

Matrix Y0 = DEInteg (&Accel,0,-(obs(46,1)-obs(1,1))*86400.0,1e-13,1e-6,6,Y);

Matrix Y_true(6);
Y_true(1) = 5753.173e3; Y_true(2) = 2673.361e3; Y_true(3) = 3440.304e3; Y_true(4) = 4.324207e3; Y_true(5) = -1.924299e3; Y_true(6) = -5.728216e3;

printf("\nError of Position Estimation\n");
printf("dX%10.1f [m]\n",Y0(1)-Y_true(1));
printf("dY%10.1f [m]\n",Y0(2)-Y_true(2));
printf("dZ%10.1f [m]\n",Y0(3)-Y_true(3));
printf("\nError of Velocity Estimation\n");
printf("dVx%8.1f [m/s]\n",Y0(4)-Y_true(4));
printf("dVy%8.1f [m/s]\n",Y0(5)-Y_true(5));
printf("dVz%8.1f [m/s]\n",Y0(6)-Y_true(6));

	int result = test(Y0, Y_true);

    if (result == 0)
        printf("PASSED\n");

    return result != 0;
}