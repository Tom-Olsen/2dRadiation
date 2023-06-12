#ifndef __INCLUDE_GUARD_Interpolation_h__
#define __INCLUDE_GUARD_Interpolation_h__
#include <math.h>
#include "DataTypes.hh"



// x€[0,1], f0=f(x=0), f1=f(x=1).
double LinearInterpolation(double x, double f0, double f1);



// x€[0,1], fm1=f(-1), fp0=f(0), fp1=f(1), fp2=f(2).
double CubicInterpolation(double x, double fm1, double fp0, double fp1, double fp2);



// x,y€[0,1]
// Calculates f(x,y) from f(0,0), f(0,1), f(1,0) and f(1,1):
// f(0,1) --- f(1,1)
//    |  (x,y)   |
// f(0,0) --- f(1,0)
// f10 <=> f(1,0)
double BilinearInterpolation
(double x  , double y  ,
 double f00, double f01,
 double f10, double f11);



#endif //__INCLUDE_GUARD_Interpolation_h__