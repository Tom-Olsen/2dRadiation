#ifndef __INCLUDE_GUARD_FourierHarmonics_h__
#define __INCLUDE_GUARD_FourierHarmonics_h__
#include <math.h>
#include <vector>
#include <iostream>
#include "Utility.hh"
#include "Stencil.h"

// This namespace holds the Fourier Harmonics basis polynomials
namespace Fourier
{
    double Basis(int k, double phi);

    std::vector<double> GetCoefficients(const Stencil& stencil, const double* const data);
    void GetCoefficients(const Stencil& stencil, const double* const data, double* coefficients);

    double GetValue(double phi, const std::vector<double>& coefficients);
    double GetValue(double phi, double* coefficients, size_t nCoefficients);
};

#endif //__INCLUDE_GUARD_FourierHarmonics_h__