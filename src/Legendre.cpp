#include "../include/Legendre.hpp"

//fi [rad]
tuple<Matrix, Matrix> Legendre(const int n,const int m,const double fi){

Matrix pnm = zeros(n+1,m+1);
Matrix dpnm = zeros(n+1,m+1);

int i;

pnm(1,1)=1.0;
dpnm(1,1)=0.0;
pnm(2,2)=sqrt(3)*cos(fi);
dpnm(2,2)=-sqrt(3)*sin(fi);
// diagonal coefficients
for (i=2; i<=n; i++){
    pnm(i+1,i+1)= sqrt((2.0*i+1.0)/(2.0*i))*cos(fi)*pnm(i,i);
}
for (i=2; i<=n; i++){
    dpnm(i+1,i+1)= sqrt((2.0*i+1.0)/(2.0*i))*((cos(fi)*dpnm(i,i))- 
                  (sin(fi)*pnm(i,i)));
}
// horizontal first step coefficients
for (i=1; i<=n; i++){
    pnm(i+1,i)= sqrt(2.0*i+1.0)*sin(fi)*pnm(i,i);
}
for (i=1; i<=n; i++){
    dpnm(i+1,i)= sqrt(2.0*i+1.0)*((cos(fi)*pnm(i,i))+(sin(fi)*dpnm(i,i)));
}
// horizontal second step coefficients
int j=0;
int k=2;
while(true){
    for (i=k; i<=n; i++){
        pnm(i+1,j+1)=sqrt((2.0*i+1.0)/((i-j)*(i+j)))*((sqrt(2.0*i-1.0)*sin(fi)*pnm(i,j+1))
            -(sqrt(((i+j-1.0)*(i-j-1.0))/(2.0*i-3.0))*pnm(i-1,j+1)));
    }
    j = j+1;
    k = k+1;
    if (j>m){
        break;
    }
}
j = 0;
k = 2;
while(true){
    for (i=k; i<=n; i++){
        dpnm(i+1,j+1)=sqrt((2.0*i+1.0)/((i-j)*(i+j)))*((sqrt(2.0*i-1.0)*sin(fi)*dpnm(i,j+1))
             +(sqrt(2.0*i-1.0)*cos(fi)*pnm(i,j+1))-(sqrt(((i+j-1.0)*(i-j-1.0))/(2.0*i-3.0))*dpnm(i-1,j+1)));
    }
    j = j+1;
    k = k+1;
		if (j>m){
			break;
		}
	}
	return tie(pnm, dpnm);
}