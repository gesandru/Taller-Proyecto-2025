#include "../include/MeasUpdate.hpp"

tuple<Matrix&, Matrix&, Matrix&> MeasUpdate( Matrix &x, double z, double g, double s, Matrix G, Matrix &P, double n){

double Inv_W = s*s;

Matrix aux1 = (G*P*transpose(G)+Inv_W);

//Kalman gain
Matrix &K = P*transpose(G)*(1.0/aux1(1,1));

K = transpose(K);
//State update
x = x + K*(z-g);

K = transpose(K);

//Covariance update
P = (eye(n)-K*G)*P;

K = transpose(K);

return tie(K, x, P);
}