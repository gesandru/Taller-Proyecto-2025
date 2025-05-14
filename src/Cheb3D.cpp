#include "../include/Cheb3D.hpp"

Matrix& Cheb3D(const int t,const int N,const int Ta,const int Tb,Matrix &Cx,Matrix &Cy,Matrix &Cz){

//Check validity
if ( (t<Ta) || (Tb<t) ){
    printf("ERROR: Time out of range in Cheb3D::Value\n");
    exit(EXIT_FAILURE);
}

//Clenshaw algorithm
double tau = (2.0*t-Ta-Tb)/(Tb-Ta);  

Matrix f1 = zeros(1,3);
Matrix f2 = zeros(1,3);

Matrix &old_f1 = zeros(1);
Matrix &auxiliar = zeros(3);

for (int i = N; i>=2;i--){
    old_f1 = f1;
	
    auxiliar(1) = Cx(i); auxiliar(2) = Cy(i); auxiliar(3) = Cz(i);
    f1 = f1*2.0*tau-f2+auxiliar;
    f2 = old_f1;
}
    auxiliar(1) = Cx(1); auxiliar(2) = Cy(1); auxiliar(3) = Cz(1);

return (f1*tau-f2+auxiliar);

}