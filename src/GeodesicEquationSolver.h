#ifndef __INCLUDE_GUARD_GeodesicEquationSolver_h__
#define __INCLUDE_GUARD_GeodesicEquationSolver_h__
#include <fstream>                      // file input/output
#include "Utility.hh"                    // small useful functions.
#include "TensorTypes.hh"                // simple containers for rank 1-3 tensors
#include "Metric2D.h"                     // metric data
#include "AdvancedUtiltiy.h"            // utiliy depending on more complex classes

// Geodesic frequency Equation.
template<class Coord>
double Dnu(const double nu, const Coordinate2<Coord>& x, const Tensor2<Coord,LF>& v, Metric2D<Coord>& metric);
// Geodesic Space Equation.
template<class Coord>
Tensor2<Coord,LF> Dx(const Coordinate2<Coord>& x, const Tensor2<Coord,LF>& v, Metric2D<Coord>& metric);
// Geodesic Velocity Equation.
template<class Coord>
Tensor2<Coord,LF> Dv(const Coordinate2<Coord>& x, const Tensor2<Coord,LF>& v, Metric2D<Coord>& metric);

// Euler solver in Lab Frame.
template<int direction, class Coord>
double Euler_GeodesicEquation(double dt, Coordinate2<Coord>& x, Tensor2<Coord,LF>& v, Metric2D<Coord>& metric);
// Adaptive Runge-Kutta 45 Geodesic Equation solver in Lab Frame.
template<int direction, class Coord>
double RK45_GeodesicEquation(double dt, Coordinate2<Coord>& x, Tensor2<Coord,LF>& v, Metric2D<Coord>& metric);


#endif //__INCLUDE_GUARD_GeodesicEquationSolver_h__