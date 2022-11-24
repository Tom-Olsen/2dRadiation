#ifndef __INCLUDE_GUARD_FourierHarmonics_h__
#define __INCLUDE_GUARD_FourierHarmonics_h__
#include <math.h>
#include <vector>
#include <iostream>
#include "Utility.hh"
#include "Stencil.hh"

// This namespace holds the Fourier Harmonics basis polynomials
namespace Fourier
{
    double Basis(const int N, const double phi);
    namespace Expansion
    {
        // data array length = fourierStencil.nDir
        std::vector<double> GetCoefficients(const Stencil& fourierStencil, const double* const data);
        double GetValue(double phi, const std::vector<double>& coefficients);
    };
};

#endif //__INCLUDE_GUARD_FourierHarmonics_h__