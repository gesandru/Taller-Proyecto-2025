#include "..\include\Matrix.hpp"
#include "..\include\const.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include "../include/AzElPa.hpp"
#include "../include/AccelPointMass.hpp"
#include "../include/Mjday.hpp"
#include "../include/Mjday_TDB.hpp"
#include "../include/Cheb3D.hpp"
#include "../include/Frac.hpp"
#include "../include/MeanObliquity.hpp"
#include "../include/sign_.hpp"
#include "../include/Position.hpp"
#include "../include/TimeUpdate.hpp"
#include "../include/NutAngles.hpp"
#include "../include/timediff.hpp"
#include "../include/EccAnom.hpp"
#include "../include/Legendre.hpp"
#include "../include/IERS.hpp"
#include "../include/AccelHarmonic.hpp"

using namespace std;

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

int m_equals(Matrix A, Matrix B, double p) {
	if (A.n_row != B.n_row || A.n_column != B.n_column)
		return 0;
	else
		for(int i = 1; i <= A.n_row; i++)
			for(int j = 1; j <= A.n_column; j++)
				if(fabs(A(i,j)-B(i,j)) > p) {
					printf("%2.20lf %2.20lf\n",A(i,j),B(i,j));
					return 0;
				}
	
	return 1;
}

int Rx_01() {
    double f = 60;
	
	Matrix A(3,3);
	A(1,1) = 1.0; A(1,2) =  0; A(1,3) = 0;
	A(2,1) = 0; A(2,2) = -0.9524; A(2,3) = -0.3048;
	A(3,1) = 0; A(3,2) =  0.3048; A(3,3) = -0.9524;
	
	Matrix B = R_x(60);
    
    _assert(m_equals(A, B, 1e-4));
    
    return 0;
}

int Ry_01() {
    double f = 60;
	
	Matrix A(3,3);
	A(1,1) = -0.9524; A(1,2) =  0; A(1,3) = 0.3048;
	A(2,1) = 0; A(2,2) = 1.0; A(2,3) = 0;
	A(3,1) = -0.3048; A(3,2) =  0; A(3,3) = -0.9524;
	
	Matrix B = R_y(60);
    
    _assert(m_equals(A, B, 1e-4));
    
    return 0;
}

int Rz_01() {
    double f = 60;
	
	Matrix A(3,3);
	A(1,1) = -0.9524; A(1,2) =  -0.3048; A(1,3) = 0;
	A(2,1) = 0.3048; A(2,2) = -0.9524; A(2,3) = 0;
	A(3,1) = 0; A(3,2) =  0; A(3,3) = 1.0;
	
	Matrix B = R_z(60);
    
    _assert(m_equals(A, B, 1e-4));
    
    return 0;
}

int AzElPa_01() {
	
	Matrix A(3);
	A(1) = -147442.081321987; A(2) =  -2627423.95702519; A(3) = 340732.976228313;
	
	tuple<double, double, Matrix, Matrix> result = AzElPa(A);
	
	Matrix B(3);
	B(1) = -3.79406152108908e-07; B(2) =  2.1290980689938e-08; B(3) = 0;
	
	Matrix dAds = get<2>(result);

    _assert(m_equals(B, dAds, 1e-10));
	
	Matrix C(3);
	C(1) = 2.71129247023044e-09; C(2) =  4.83153434006981e-08; C(3) = 3.73737378628265e-07;
	
	Matrix dEds = get<3>(result);
    
    _assert(m_equals(C, dEds, 1e-4));
    
		int i;
		
		double D = 3.19765045275856;
		double Az = get<0>(result);
		
	if(fabs(D-Az) > 1e-10) {
		printf("%2.20lf %2.20lf\n",D,Az);
		i=0;
	}else{
		i=1;
	}
	
	_assert(i);
	
	double E = 0.128763196439668;
	double El = get<1>(result);
		
	if(fabs(E-El) > 1e-10) {
		printf("%2.20lf %2.20lf\n",E,El);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
    return 0;
}

int AccelpPointMass_01() {
	
Matrix A(3);
	A(1) = 5.5; A(2) =  5.5; A(3) = 5.5;
	
	Matrix B(3);
	B(1) = 1.0; B(2) =  1.0; B(3) = 1.0;
	
	double GM = 2.7;
    
	Matrix R = AccelPointMass(A,B,GM);
	
	Matrix C(3);
	C(1) = -0.5453; C(2) =  -0.5453; C(3) = -0.5453;
    
    _assert(m_equals(R, C, 1e-4));
    
    return 0;
}

int Mjday_01() {
	
	
    
	double E = 4.974610972222220e+04;
	double El = Mjday(1995,1,29,02,38,0);
	
	int i;
	if(fabs(E-El) > 1e-10) {
		printf("%2.20lf %2.20lf\n",E,El);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
    return 0;
}

int Mjday_02() {
	
	
    
	double E = 4.974610972222220e+04;
	double El = Mjday(1995,1,29,02,38);
	
	int i;
	if(fabs(E-El) > 1e-10) {
		printf("%2.20lf %2.20lf\n",E,El);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
    return 0;
}

int Mjday_TDB_01() {
	
	
    
	double E = 4.974611706232293e+04;
	double El = Mjday_TDB(4.974611706231468e+04);
	
	int i;
	if(fabs(E-El) > 1e-10) {
		printf("%2.20lf %2.20lf\n",E,El);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
    return 0;
}

int Cheb3D_01() {
	
Matrix A(3);
	A(1) = 5.0; A(2) =  3.0; A(3) = 8.0;
	
	Matrix B(3);
	B(1) = 13.0; B(2) =  74.0; B(3) = 29.0;
	
	Matrix C(3);
	C(1) = 45.0; C(2) =  63.0; C(3) = 17.0;
	
	int t = 34710;
	int N = 3;
	int Ta = 34704;
	int Tb = 34720;
    
	Matrix R = Cheb3D(t,N,Ta,Tb,A,B,C);
	
	Matrix D(3);
	D(1) = -2.75000000000000; D(2) =  -30.8750000000000; D(3) = 14.3750000000000;
    
    _assert(m_equals(R, D, 1e-4));
    
    return 0;
}

int Frac_01() {
	
	double E = 0.02;
	double El = Frac(1.02);
	
	int i;
	if(fabs(E-El) > 1e-10) {
		printf("%2.20lf %2.20lf\n",E,El);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
    return 0;
}

int MeanObliquity_01() {
	
	double E = 0.408791718328378;
	double El = MeanObliquity(100000.0);
	
	int i;
	if(fabs(E-El) > 1e-10) {
		printf("%2.20lf %2.20lf\n",E,El);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
    return 0;
}

int sign_01() {
	
	double E = -365.7;
	double El = sign_(-365.7, -1.0);
	
	int i;
	if(fabs(E-El) > 1e-10) {
		printf("%2.20lf %2.20lf\n",E,El);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
    return 0;
}

int sign_02() {
	
	double E = 365.7;
	double El = sign_(-365.7, 1.0);
	
	int i;
	if(fabs(E-El) > 1e-10) {
		printf("%2.20lf %2.20lf\n",E,El);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
    return 0;
}

int sign_03() {
	
	double E = -365.7;
	double El = sign_(365.7, -1.0);
	
	int i;
	if(fabs(E-El) > 1e-10) {
		printf("%2.20lf %2.20lf\n",E,El);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
    return 0;
}

int sign_04() {
	
	double E = 365.7;
	double El = sign_(365.7, 1.0);
	
	int i;
	if(fabs(E-El) > 1e-10) {
		printf("%2.20lf %2.20lf\n",E,El);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
    return 0;
}

int Position_01() {
	
Matrix A(3);
	A(1) = 1406307.88441300; A(2) =  965778.690137456; A(3) = 6125143.46191200;
	
	double lon = 76.0;
	double lat = 39.0;
	double h = 3.0;
	
	Matrix B = Position(lon,lat,h);
	
    _assert(m_equals(A, B, 1e-5));
    
    return 0;
}

int TimeUpdate_01() {
	
Matrix P = zeros(6,6);
	P(1,1) = 100000000.0; P(2,2) =  100000000.0; P(3,3) = 100000000.0;
	P(4,4) = 1000.0; P(5,5) =  1000.0; P(6,6) = 1000.0;
	
	Matrix Phi(6,6);
	Phi(1,1) = 1.00002525535511; Phi(1,2) =  7.08259834024473e-06; Phi(1,3) = 1.91609350460620e-07; Phi(1,4) = 5.00001498082565; Phi(1,5) = 1.17818628919961e-05; Phi(1,6) = 2.68390168073246e-07;
	Phi(2,1) = 7.08259815373561e-06; Phi(2,2) =  0.999988040046622; Phi(2,3) = 3.53016114843062e-08; Phi(2,4) = 1.17818628268260e-05; Phi(2,5) = 4.99995293819715; Phi(2,6) = 4.93631303180711e-08;
	Phi(3,1) = 1.91608861002907e-07; Phi(3,2) =  3.53015288644891e-08; Phi(3,3) = 0.999986704774626; Phi(3,4) = 2.68389762645616e-07; Phi(3,5) = 4.93630678596596e-08; Phi(3,6) = 4.99995072081276;
	Phi(4,1) = 1.01043851887223e-05; Phi(4,2) =  2.82768357826951e-06; Phi(4,3) = 6.44136325079115e-08; Phi(4,4) = 1.00002526606744; Phi(4,5) = 7.05571117883109e-06; Phi(4,6) = 1.30455621823661e-07;
	Phi(5,1) = 2.82768336557965e-06; Phi(5,2) =  -4.78603729288896e-06; Phi(5,3) = 1.18528331537947e-08; Phi(5,4) = 7.05571100144785e-06; Phi(5,5) = 0.999988029832331; Phi(5,6) = 2.39619698989173e-08;
	Phi(6,1) = 6.44131451075285e-08; Phi(6,2) =  1.18527461137171e-08; Phi(6,3) = -5.31820682446032e-06; Phi(6,4) = -5.31820682446032e-06; Phi(6,5) = 2.39618837211068e-08; Phi(6,6) = 0.999986704276552;
	
	Matrix B(6,6);
	B(1,1) = 100030051.289634; B(1,2) =  1416.64688492363; B(1,3) = 38.3247592060704; B(1,4) = 6010.60735344927; B(1,5) = 282.819148824426; B(1,6) = 6.41506095453518;
	B(2,1) = 1416.64688492363; B(2,2) =  100022607.558029; B(2,3) = 7.06085421838136; B(2,4) = 282.819193072307; B(2,5) = 4521.29708557101; B(2,6) = 1.18545645495969;
	B(3,1) = 38.3247592060704; B(3,2) =  7.06085421838136; B(3,3) = 100022340.479816; B(3,4) = 6.44240187111654; B(3,5) = 1.18547401528551; B(3,6) = 4468.07063235661;
	B(4,1) = 6010.60735344927; B(4,2) =  282.819193072307; B(4,3) = 6.44240187111654; B(4,4) = 1000.06154267741; B(4,5) = 0.0156154526297801; B(4,6) = -0.0051537065357492;
	B(5,1) = 282.819148824426; B(5,2) =  4521.29708557101; B(5,3) = 1.18547401528551; B(5,4) = 0.0156154526297801; B(5,5) = 999.979150066399; B(5,6) = 5.41618162396946e-05;
	B(6,1) = 6.41506095453518; B(6,2) =  1.18545645495969; B(6,3) = 4468.07063235661; B(6,4) = -0.0051537065357492; B(6,5) = 5.4123372055547e-05; B(6,6) = 999.976237491234;
	
	TimeUpdate(P, Phi, 0);
	
    _assert(m_equals(P, B, 1e-5));
    
    return 0;
}

int TimeUpdate_02() {
	
Matrix P = zeros(6,6);
	P(1,1) = 100000000.0; P(2,2) =  100000000.0; P(3,3) = 100000000.0;
	P(4,4) = 1000.0; P(5,5) =  1000.0; P(6,6) = 1000.0;
	
	Matrix Phi(6,6);
	Phi(1,1) = 1.00002525535511; Phi(1,2) =  7.08259834024473e-06; Phi(1,3) = 1.91609350460620e-07; Phi(1,4) = 5.00001498082565; Phi(1,5) = 1.17818628919961e-05; Phi(1,6) = 2.68390168073246e-07;
	Phi(2,1) = 7.08259815373561e-06; Phi(2,2) =  0.999988040046622; Phi(2,3) = 3.53016114843062e-08; Phi(2,4) = 1.17818628268260e-05; Phi(2,5) = 4.99995293819715; Phi(2,6) = 4.93631303180711e-08;
	Phi(3,1) = 1.91608861002907e-07; Phi(3,2) =  3.53015288644891e-08; Phi(3,3) = 0.999986704774626; Phi(3,4) = 2.68389762645616e-07; Phi(3,5) = 4.93630678596596e-08; Phi(3,6) = 4.99995072081276;
	Phi(4,1) = 1.01043851887223e-05; Phi(4,2) =  2.82768357826951e-06; Phi(4,3) = 6.44136325079115e-08; Phi(4,4) = 1.00002526606744; Phi(4,5) = 7.05571117883109e-06; Phi(4,6) = 1.30455621823661e-07;
	Phi(5,1) = 2.82768336557965e-06; Phi(5,2) =  -4.78603729288896e-06; Phi(5,3) = 1.18528331537947e-08; Phi(5,4) = 7.05571100144785e-06; Phi(5,5) = 0.999988029832331; Phi(5,6) = 2.39619698989173e-08;
	Phi(6,1) = 6.44131451075285e-08; Phi(6,2) =  1.18527461137171e-08; Phi(6,3) = -5.31820682446032e-06; Phi(6,4) = -5.31820682446032e-06; Phi(6,5) = 2.39618837211068e-08; Phi(6,6) = 0.999986704276552;
	
	Matrix B(6,6);
	B(1,1) = 100030051.289634; B(1,2) =  1416.64688492363; B(1,3) = 38.3247592060704; B(1,4) = 6010.60735344927; B(1,5) = 282.819148824426; B(1,6) = 6.41506095453518;
	B(2,1) = 1416.64688492363; B(2,2) =  100022607.558029; B(2,3) = 7.06085421838136; B(2,4) = 282.819193072307; B(2,5) = 4521.29708557101; B(2,6) = 1.18545645495969;
	B(3,1) = 38.3247592060704; B(3,2) =  7.06085421838136; B(3,3) = 100022340.479816; B(3,4) = 6.44240187111654; B(3,5) = 1.18547401528551; B(3,6) = 4468.07063235661;
	B(4,1) = 6010.60735344927; B(4,2) =  282.819193072307; B(4,3) = 6.44240187111654; B(4,4) = 1000.06154267741; B(4,5) = 0.0156154526297801; B(4,6) = -0.0051537065357492;
	B(5,1) = 282.819148824426; B(5,2) =  4521.29708557101; B(5,3) = 1.18547401528551; B(5,4) = 0.0156154526297801; B(5,5) = 999.979150066399; B(5,6) = 5.41618162396946e-05;
	B(6,1) = 6.41506095453518; B(6,2) =  1.18545645495969; B(6,3) = 4468.07063235661; B(6,4) = -0.0051537065357492; B(6,5) = 5.4123372055547e-05; B(6,6) = 999.976237491234;
	
	TimeUpdate(P, Phi);
	
    _assert(m_equals(P, B, 1e-5));
    
    return 0;
}

int NutAngles_01() {
	
	tuple <double, double> result = NutAngles(321.3);

	int i;
		
		double D = 5.31658395704401e-05;
		double Az = get<0>(result);
		
	if(fabs(D-Az) > 1e-10) {
		printf("%2.20lf %2.20lf\n",D,Az);
		i=0;
	}else{
		i=1;
	}
	
	_assert(i);
	
	double E = 3.49865427225745e-05;
	double El = get<1>(result);
		
	if(fabs(E-El) > 1e-10) {
		printf("%2.20lf %2.20lf\n",E,El);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
    return 0;
}

int timediff_01() {
	
	tuple <double, double, double, double, double> result = timediff(5216.2, 4653.1);

	int i;
		
		double A = 563.099999999999;
		double UT1_TAI = get<0>(result);
		
	if(fabs(A-UT1_TAI) > 1e-5) {
		printf("%2.20lf %2.20lf\n",A,UT1_TAI);
		i=0;
	}else{
		i=1;
	}
	
	_assert(i);
	
	double B = -4634.1;
	double UTC_GPS = get<1>(result);
		
	if(fabs(B-UTC_GPS) > 1e-5) {
		printf("%2.20lf %2.20lf\n",B,UTC_GPS);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
		double C = 582.099999999999;
	double UT1_GPS = get<2>(result);
		
	if(fabs(C-UT1_GPS) > 1e-5) {
		printf("%2.20lf %2.20lf\n",C,UT1_GPS);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
	double D = 4685.284;
	double TT_UTC = get<3>(result);
		
	if(fabs(D-TT_UTC) > 1e-5) {
		printf("%2.20lf %2.20lf\n",D,TT_UTC);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
		double E = 4634.1;
	double GPS_UTC = get<4>(result);
		
	if(fabs(E-GPS_UTC) > 1e-5) {
		printf("%2.20lf %2.20lf\n",E,UTC_GPS);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
    return 0;
}

int EccAnom_01() {

	int i;
		
		double D = 0.262936914041412;
		double Az = EccAnom(69.3,0.3);
		
	if(fabs(D-Az) > 1e-10) {
		printf("%2.20lf %2.20lf\n",D,Az);
		i=0;
	}else{
		i=1;
	}
	
	_assert(i);
	
    return 0;
}

int EccAnom_02() {

	int i;
	
	double E = 0.876771881596796;
	double El = EccAnom(69.3,0.9);
		
	if(fabs(E-El) > 1e-10) {
		printf("%2.20lf %2.20lf\n",E,El);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
    return 0;
}

int Legendre_01() {
	
	tuple<Matrix, Matrix> result = Legendre(3,3,5.1);
	
	Matrix A(4,4);
	A(1,1) = 1.0; A(1,2) =  0.0; A(1,3) = 0.0; A(1,4) = 0.0;
	A(2,1) = -1.603558068184872; A(2,2) =  0.654676654509079; A(2,3) = 0.0; A(2,4) = 0.0;
	A(3,1) = 1.756876908319470; A(3,2) =  -1.355301504740954; A(3,3) = 0.276661092784522; A(3,4) = 0.0;
	A(4,1) = -1.574601335784816; A(4,2) =  2.012120785940833; A(4,3) = -0.677674543560851; A(4,4) = 0.112950383527956;
	
		Matrix pnm = get<0>(result);

    _assert(m_equals(A, pnm, 1e-10));
	
	Matrix B(4,4);
	B(1,1) = 0.0; B(1,2) =  0.0; B(1,3) = 0.0; B(1,4) = 0.0;
	B(2,1) = 0.654676654509079; B(2,2) =  1.603558068184872; B(2,3) = 0.0; B(2,4) = 0.0;
	B(3,1) = -2.347451065785884; B(3,2) =  -2.766338975069329; B(3,3) = 1.355301504740954; B(3,4) = 0.0;
	B(4,1) = 4.928669226402895; B(4,2) =  2.785472285993957; B(4,3) = -3.043106902522991; B(4,4) = 0.829978421698789;
	
	Matrix dpnm = get<1>(result);

    _assert(m_equals(B, dpnm, 1e-10));
	
    return 0;
}

int IERS_01() {
	
	eop19620101(21413);
	
	tuple< double, double, double, double, double, double, double, double, double> result = IERS(eopdata, 4.974611635416653e+04);

	int i;
		
		double A = -5.595186212317044e-07;
		double x_pole = get<0>(result);
		
	if(fabs(A-x_pole) > 1e-15) {
		printf("%2.20lf %2.20lf\n",A,x_pole);
		i=0;
	}else{
		i=1;
	}
	
	_assert(i);
	
	double B = 2.334586344425293e-06;
	double y_pole = get<1>(result);
		
	if(fabs(B-y_pole) > 1e-15) {
		printf("%2.20lf %2.20lf\n",B,y_pole);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
		double C = 0.326067700000000;
	double UT1_UTC = get<2>(result);
		
	if(fabs(C-UT1_UTC) > 1e-10) {
		printf("%2.20lf %2.20lf\n",C,UT1_UTC);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
	double D = 0.002721300000000;
	double LOD = get<3>(result);
		
	if(fabs(D-LOD) > 1e-10) {
		printf("%2.20lf %2.20lf\n",D,LOD);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
	double E = -1.168643378314537e-07;
	double dpsi = get<4>(result);
		
	if(fabs(E-dpsi) > 1e-15) {
		printf("%2.20lf %2.20lf\n",E,dpsi);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
	double F = -2.487094184091920e-08;
	double deps = get<5>(result);
		
	if(fabs(F-deps) > 1e-20) {
		printf("%2.20lf %2.20lf\n",F,deps);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
	double G = -8.193351210751158e-10;
	double dx_pole = get<6>(result);
		
	if(fabs(G-dx_pole) > 1e-20) {
		printf("%2.20lf %2.20lf\n",G,dx_pole);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
	double H = -1.532011232306134e-09;
	double dy_pole = get<7>(result);
		
	if(fabs(H-dy_pole) > 1e-20) {
		printf("%2.20lf %2.20lf\n",H,dy_pole);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
	double I = 29.0;
	double TAI_UTC = get<8>(result);
		
	if(fabs(I-TAI_UTC) > 1e-2) {
		printf("%2.20lf %2.20lf\n",I,TAI_UTC);
		i=0;
	}else{
		i=1;
	}
    
    _assert(i);
	
    return 0;
}

int AccelHarmonic_01() {
	
	Matrix A(3,1);
	A(1,1) = 5.542555937228607e+06;
	A(2,1) = 3.213514867349196e+06;
	A(3,1) = 3.990892975876853e+06;
	
	Matrix B(3,3);
	B(1,1) = -0.978185453896254; B(1,2) =  0.207733066362260; B(1,3) = -0.000436950239569;
	B(2,1) = -0.207733028352522; B(2,2) =  -0.978185550768511; B(2,3) = -0.000131145697267;
	B(3,1) = -0.000454661708585; B(3,2) =  -0.000037515816903; B(3,3) = 0.999999895937642;
	
	Matrix R = AccelHarmonic(A, B, 20, 20);
	
	Matrix C(3,1);
	C(1,1) = -5.134836943433859;
	C(2,1) = -2.977177178054876;
	C(3,1) = -3.705918925034543;

    _assert(m_equals(C, R, 1e-10));
	
    return 0;
}

int all_tests()
{
    _verify(Rx_01);
	_verify(Ry_01);
	_verify(Rz_01);
	_verify(AzElPa_01);
	_verify(AccelpPointMass_01);
	_verify(Mjday_01);
	_verify(Mjday_02);
	_verify(Mjday_TDB_01);
	_verify(Cheb3D_01);
	_verify(Frac_01);
	_verify(MeanObliquity_01);
	_verify(sign_01);
	_verify(sign_02);
	_verify(sign_03);
	_verify(sign_04);
	_verify(Position_01);
	_verify(TimeUpdate_01);
	_verify(TimeUpdate_02);
	_verify(NutAngles_01);
	_verify(timediff_01);
	_verify(EccAnom_01);
	_verify(EccAnom_02);
	_verify(Legendre_01);
	_verify(IERS_01);
	_verify(AccelHarmonic_01);

    return 0;
}


int main()
{
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
