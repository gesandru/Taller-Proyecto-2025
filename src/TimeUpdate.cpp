#include "../include/TimeUpdate.hpp"

void TimeUpdate(Matrix &P, Matrix Phi, double Qdt){



P = Phi*P*transpose(Phi) + Qdt;
}

