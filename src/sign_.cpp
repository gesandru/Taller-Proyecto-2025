#include "../include/sign_.hpp"

using namespace std;

// sign: returns absolute value of a with sign of b
double sign_(const double a,const double b){

double result = 0.0;

if (b>=0.0){
    result = abs(a);
}else{
    result = -abs(a);
}

return result;
}