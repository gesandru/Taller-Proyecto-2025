/**
 * @file globals.hpp
 * @brief Param y ficheros
 *
 * @author George
 *
 * @note Aquí se define la estructura para AuxParam y los métodos que leen de ficheros y lo guardan en memoria.
 */

#ifndef _GLOBALS_
#define _GLOBALS_

#include "../include/Matrix.hpp"
#include "../include/MjDay.hpp"
#include "../include/const.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <stdio.h>
#include <cstring>
#include <fstream>

struct Param{
    double Mjd_UTC, Mjd_TT;
    int n, m, sun, moon, planets;
};

extern Param AuxParam;

void AParam();

extern Matrix eopdata;
extern Matrix Cnm;
extern Matrix Snm;
extern Matrix PC;
extern Matrix obs;

/**
 * @brief Lee del fichero eop19620101.txt y lo guarda en la matríz eopdata
 *
 * @param c representa las lineas a leer del fichero, por defecto introducir 21413
 */
void eop19620101(int c);

/**
 * @brief Lee del fichero GGM03S.txt y lo guarda en las matrices Cnm y Snm
 *
 */
void GGMO3S();

/**
 * @brief Lee del fichero DE430Coeff.txt y lo guarda en la matríz PC
 *
 */
void DE430Coeff();

/**
 * @brief Lee del fichero GEOS3.txt y lo guarda en la matriz obs
 *
 * @param f representa las lineas a leer del fichero, por defecto introducir 46
 */
void GEOS3(int f);

#endif