#ifndef __INCLUDE_GUARD_Interpolation_hh__
#define __INCLUDE_GUARD_Interpolation_hh__
#include <math.h>
#include <iostream>
#include "Stencil.hh"
#include "DataTypes.hh"
#include "Profiler.hh"


/*
// x€[0,1], f0=f(x=0), f1=f(x=1).
INLINE double LinearInterpolation(double x, double f0, double f1)
{
    return f0*(1.0-x) + f1*x;
}
// x€[0,1], fm1=f(-1), fp0=f(0), fp1=f(1), fp2=f(2).
INLINE double CubicInterpolation(double x, double fm1, double fp0, double fp1, double fp2)
{
    double c0 = fp0;
    double c1 = 0.5 * (fp1 - fm1);
    double c2 = fm1 - 2.5 * fp0 + 2.0 * fp1 - 0.5 * fp2;
    double c3 = 0.5 * (fp2 - fm1) + 1.5 * (fp0 - fp1);

    return ((c3 * x + c2) * x + c1) * x + c0;
}



// i,j€[0,1]
// Calculates f(i,j) from f(0,0), f(0,1), f(1,0) and f(1,1):
// f(0,1) --- f(1,1)
//    |  (i,j)   |
// f(0,0) --- f(1,0)
// f10 <=> f(1,0)
INLINE double BilinearInterpolation
(double i  , double j  ,
 double f00, double f01,
 double f10, double f11)
{
    double f0 = f00*(1.0-i) + f10*i;
    double f1 = f01*(1.0-i) + f11*i;
    return f0*(1.0-j) + f1*j;
}
// Cubic Interpolation, x€[0,1],
// fm1m1=f(-1,-1), fm1p0=f(-1,0), fm1p1=f(-1,1), fm1p2=f(-1,2).
// fp0m1=f( 0,-1), fp0p0=f( 0,0), fp0p1=f( 0,1), fp0p2=f( 0,2).
// fp1m1=f( 1,-1), fp1p0=f( 1,0), fp1p1=f( 1,1), fp1p2=f( 1,2).
// fp2m1=f( 2,-1), fp2p0=f( 2,0), fp2p1=f( 2,1), fp2p2=f( 2,2).
INLINE double BicubicInterpolation
(double x, double y,
 double fm1m1, double fm1p0, double fm1p1, double fm1p2,
 double fp0m1, double fp0p0, double fp0p1, double fp0p2,
 double fp1m1, double fp1p0, double fp1p1, double fp1p2,
 double fp2m1, double fp2p0, double fp2p1, double fp2p2)
{
    double Fm1 = CubicInterpolation(x,fm1m1,fp0m1,fp1m1,fp2m1);
    double Fp0 = CubicInterpolation(x,fm1p0,fp0p0,fp1p0,fp2p0);
    double Fp1 = CubicInterpolation(x,fm1p1,fp0p1,fp1p1,fp2p1);
    double Fp2 = CubicInterpolation(x,fm1p2,fp0p2,fp1p2,fp2p2);
    return CubicInterpolation(y,Fm1,Fp0,Fp1,Fp2);
}*/
#endif //__INCLUDE_GUARD_Interpolation_hh__