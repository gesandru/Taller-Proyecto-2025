#include "../include/AzElPa.hpp"

using namespace std;

tuple<double, double, Matrix, Matrix> AzElPa(Matrix& s){

double pi2 = 2.0*numbers::pi;

double rho = sqrt(s(1)*s(1)+s(2)*s(2));

//Angles
double Az = atan2(s(1),s(2));

if (Az<0.0){
    Az = Az+pi2;
}

double El = atan ( s(3) / rho );

// Partials
Matrix dAds(3);
dAds(1)= s(2)/(rho*rho); dAds(2)= -s(1)/(rho*rho); dAds(3)= 0.0;
Matrix dEds(3);
dEds(1)= -s(1)*s(3)/rho / dot(s,s); dEds(2)= -s(2)*s(3)/rho / dot(s,s); dEds(3)= rho / dot(s,s);

return  tie(Az, El, dAds, dEds);
}