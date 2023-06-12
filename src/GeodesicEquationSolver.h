#ifndef __INCLUDE_GUARD_GeodesicEquationSolver_h__
#define __INCLUDE_GUARD_GeodesicEquationSolver_h__
#include "Spacetimes.h"     // Metric data.



// Geodesic frequency Equation.
double Dnu(const double nu, const Coord& x, const Tensor2& v, Metric& metric);
// Geodesic position Equation.
Tensor2 Dx(const Coord& x, const Tensor2& v, Metric& metric);
// Geodesic velocity Equation.
Tensor2 Dv(const Coord& x, const Tensor2& v, Metric& metric);

// Euler solver in Lab Frame.
template<int timeDirection>
double Euler_GeodesicEquation(double dt, Coord& x, Tensor2& v, Metric& metric);
// Adaptive Runge-Kutta 45 Geodesic Equation solver in Lab Frame.
template<int timeDirection>
double RK45_GeodesicEquation(double dt, Coord& x, Tensor2& v, Metric& metric);


#endif //__INCLUDE_GUARD_GeodesicEquationSolver_h__