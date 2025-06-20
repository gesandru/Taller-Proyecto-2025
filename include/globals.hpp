#ifndef _GLOBALS_
#define _GLOBALS_

#include "../include/Matrix.hpp"
#include <cmath>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <stdio.h>
#include <cstring>


//lo mismo haremos con GEOS3 y con GGMO3S

struct Param{
    double Mjd_UTC, Mjd_TT;
    int n, m, sun, moon, planets;
};

extern Param AuxParam;

void AParam();

/*
class Global{
    public:
    static int n_eqn;
    static Param *AuxParam;
    
    static Matrix *obs;
    static Matrix *GGMO3Sdata;

    static Matrix *PC;
	
    static void GEOS3();
    
    static Param* getParam();
};
*/
extern Matrix eopdata;
extern Matrix Cnm;
extern Matrix Snm;
extern Matrix PC;

void eop19620101(int c); //c = 21413

void GGMO3S();

void DE430Coeff();

#endif