#include "FourierHarmonics.h"



double Fourier::Basis(int k, double phi)
{
    constexpr double sqrt2 = sqrt(2);
    if (k == 0)
        return 1.0;
    else if (k%2 != 0)
        return sqrt2 * cos(phi * (k+1.0)/2.0);
    else // (k%2 == 0)
        return sqrt2 * sin(phi * k/2.0);
}



std::vector<double> Fourier::GetCoefficients(const Stencil& stencil, const double* const data)
{
    std::vector<double> coefficients(stencil.nCoefficients);
    for(size_t i=0; i<stencil.nCoefficients; i++)
        coefficients[i] = 0;
    for(int d=0; d<stencil.nDir; d++)
    {
        double phi = stencil.Phi(d);
        double c = data[d] * stencil.W(d);
        
        for(int i=0; i<stencil.nCoefficients; i++)
            coefficients[i] += c * Basis(i,phi);
    }
    return coefficients;
}
void Fourier::GetCoefficients(const Stencil& stencil, const double* data, double* coefficients)
{
    for(size_t i=0; i<stencil.nCoefficients; i++)
        coefficients[i] = 0;
    for(size_t d=0; d<stencil.nDir; d++)
    {
        double phi = stencil.Phi(d);
        double c = data[d] * stencil.W(d);

        for(size_t i=0; i<stencil.nCoefficients; i++)
            coefficients[i] += c * Basis(i,phi);
    }
}



double Fourier::GetValue(double phi, const std::vector<double>& coefficients)
{
    double result = 0;
    for(size_t i=0; i<coefficients.size(); i++)
        result += coefficients[i] * Basis(i,phi);
    return result;
}
double Fourier::GetValue(double phi, double* coefficients, size_t nCoefficients)
{
    double result = 0;
    for(size_t i=0; i<nCoefficients; i++)
        result += coefficients[i] * Basis(i,phi);
    return result;
}