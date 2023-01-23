#ifndef __INCLUDE_GUARD_Includes_hh__
#define __INCLUDE_GUARD_Includes_hh__

/* Files ending with .h have a corresponding .cpp file and can thus be compiled seperately.
 * Files ending with .hh use nested templates which cannot be precompiled. They must be included and compiled with the project.
 * The compilation steps are completely handeld by the Makefile. It is sufficient to include this 'Includes.hh' file in your cpp file to use the code.
 */

#include "ControlFlow.hh"                   // control flow classes and macros
#include "Utility.hh"                       // small useful functions.
#include "TensorTypes.hh"                   // simple containers for rank 1-3 tensors.
#include "Grid2D.h"                         // 2D Grid for numerical domain and maping to physical domain.
#include "Stencil.hh"                       // velocity stencils.
#include "Interpolation.hh"                 // several interpolation schemes which are needed throughout the code.
#include "Metric2D.h"                       // 2D Metric data. The numerical domain is defined by Grid3D.
#include "Spacetimes.h"                     // all spacetimes derive from the base class metric and define the actual components.
#include "FourierHarmonics.h"               // fourier interpolation thorugh discrete number of points in 2d plane.
#include "WriteData.h"                      // write stuff to vtr files.
#include "AdvancedUtiltiy.h"                // utiliy depending on more complex classes
#include "GeodesicEquationSolver.h"         // used to solve geodesic equation on spacelike hypersurface.
#include "Radiation.h"                      // Curved spacetime, rotating stencil.
#include "Profiler.hh"                      // Time measurements

// TODO:
// -interpolation.h/cpp
// -FourierHarmonics.h/cpp
// -SimulationData.h/cpp
// -WriteData.h/cpp <= delete this?
// -Introduce frame dependet tensors?
//  template<class CoordB>
//  Tensor2<CoordB,Frame> CoordinateTransformation(const Coordinate2<Coord>& x12);
//  template<class FrameB>
//  Tensor2<Coord,FrameB> FrameTransformation(const Tensor4x4<Coord,Frame>& tetrad);

/* Further notes on the Code:
 * In the test.cpp you can find functions which test every part of the code individually.
 * In the main.cpp you can find an implementation of a radiation transport simulation in a vacuum (no collision, just streaming).
 * Even though the grid and metric are formulated in 3D the current code only works in 2D.
 * For this a 3D grid+metric are initialized with a thickness in z-direction of 2 (1 is not possible due to dz=(zend-zstart)/(nz-1)=nan then).
 * In order to upgrade the code to a full 3D simulation only the following thins need to be updated:
 * -velocity stencils (directions+weights)
 * -streaming step (swap int k=0 with a loop over k and include k0,k1 Iat_000...111)
 * -AngularVectorInterpolation needs a 3D version
 */

#endif //__INCLUDE_GUARD_Includes_hh__