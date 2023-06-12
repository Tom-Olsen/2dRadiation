#include "Interpolation.h"



// x€[0,1], f0=f(x=0), f1=f(x=1).
double LinearInterpolation(double x, double f0, double f1)
{
    return f0 * (1.0 - x) + f1 * x;
}



// x€[0,1], fm1=f(-1), fp0=f(0), fp1=f(1), fp2=f(2).
double CubicInterpolation(double x, double fm1, double fp0, double fp1, double fp2)
{
    double c0 = fp0;
    double c1 = 0.5 * (fp1 - fm1);
    double c2 = fm1 - 2.5 * fp0 + 2.0 * fp1 - 0.5 * fp2;
    double c3 = 0.5 * (fp2 - fm1) + 1.5 * (fp0 - fp1);
    return ((c3 * x + c2) * x + c1) * x + c0;
}



// x,y€[0,1]
// Calculates f(x,y) from f(0,0), f(0,1), f(1,0) and f(1,1):
// f(0,1) --- f(1,1)
//    |  (x,y)   |
// f(0,0) --- f(1,0)
// f10 <=> f(1,0)
double BilinearInterpolation
(double x  , double y  ,
 double f00, double f01,
 double f10, double f11)
{
    double f0 = f00 * (1.0 - x) + f10 * x;
    double f1 = f01 * (1.0 - x) + f11 * x;
    return f0 * (1.0 - y) + f1 * y;
}