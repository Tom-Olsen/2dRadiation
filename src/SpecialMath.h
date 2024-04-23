#ifndef __INCLUDE_GUARD_AdvancedMath_h__
#define __INCLUDE_GUARD_AdvancedMath_h__
#include <vector>       // std::vector<T>
#include <algorithm>    // sort
#include "DataTypes.hh" // General relativity tensors.
#include "Spacetimes.h" // Metric data.

double Dot(const Tensor2 &vector0, const Tensor2 &vector1, const Tensor2x2 &gamma_ll);
double Dot(const Tensor3 &vector0, const Tensor3 &vector1, const Tensor3x3 &g_ll);
double Norm2(const Tensor2 &vector, const Tensor2x2 &gamma_ll);
double Norm2(const Tensor3 &vector, const Tensor3x3 &g_ll);
Tensor3 NullNormalize(const Tensor3 &vector, const Tensor3x3 &g_ll);

Tensor2 TransformIFtoLF(const Tensor2 &vector, const Tensor3x3 &tetrad);
Tensor2 TransformLFtoIF(const Tensor2 &vector, const Tensor3x3 &tetradInverse);

Tensor3 TransformIFtoLF(const Tensor3 &vector, const Tensor3x3 &tetrad);
Tensor3 TransformLFtoIF(const Tensor3 &vector, const Tensor3x3 &tetradInverse);

Tensor3x3 TransformIFtoLF(const Tensor3x3 &tensor, const Tensor3x3 &tetrad);
Tensor3x3 TransformLFtoIF(const Tensor3x3 &tensor, const Tensor3x3 &tetradInverse);

Tensor3x3 BoostMatrix(const Tensor2 & u);

template <class FrameIn, class FrameOut>
Tensor2 Vec2ObservedByEulObs(const Tensor3 &u, const Coord &xy, Metric &metric);
#endif //__INCLUDE_GUARD_AdvancedMath_h__