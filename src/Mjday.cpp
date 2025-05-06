#include "../include/Mjday.h"

double Mjday(double yr, double mon, double day, double hr = 0.0, double min = 0.0, double sec = 0.0){

double jd = 367.0 * yr
    - floor( (7 * (yr + floor( (mon + 9) / 12.0) ) ) * 0.25 )
    + floor( 275 * mon / 9.0 )
    + day + 1721013.5
    + ( (sec/60.0 + min ) / 60.0 + hr ) / 24.0;

return jd-2400000.5;

}