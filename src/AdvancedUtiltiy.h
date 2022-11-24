#ifndef __INCLUDE_GUARD_AdvancedUtility_h__
#define __INCLUDE_GUARD_AdvancedUtility_h__
#include "ControlFlow.hh"// used for template arguments
#include "Utility.hh"            // small useful functions (debug only).
#include "TensorTypes.hh"       // simple containers for rank 1-3 tensors
#include "Metric2D.h"           // metric data


template<class Coord, class FrameIn, class FrameOut>
Tensor2<Coord,FrameOut> Vec2ObservedByEulObs(const Tensor3<Coord,FrameIn>& u, const Coordinate2<Coord>& x, Metric2D<Coord>& metric);


#endif //__INCLUDE_GUARD_AdvancedUtility_h__