#include "FourierHarmonics.h"

double Fourier::Basis(const int k, const double phi)
{
    if (k == 0)
        return 1.0 / sqrt(2);
    else if (k%2 != 0)
        return cos(phi * (k+1.0)/2.0);
    else // (k%2 == 0)
        return sin(phi * k/2.0);
        // return 0;
}

// data array length = fourierStencil.nDir
std::vector<double> Fourier::Expansion::GetCoefficients(const Stencil& fourierStencil, const double* const data, double rotation)
{
    int N = fourierStencil.nDir;
    std::vector<double> coefficients(N);
    // Coefficients:
    for(int k=0; k<N; k++)
    {
        coefficients[k] = 0;
        // Integral:
        for(int i=0; i<N; i++)
            coefficients[k] += data[i] * Fourier::Basis(k,fourierStencil.Phi(i,rotation)) * fourierStencil.W(i);
    }
    return coefficients;
}
double Fourier::Expansion::GetValue(double phi, const std::vector<double>& coefficients)
{
    int N = coefficients.size();
    double result = 0;
    for(int j=0; j<N; j++)
        result += coefficients[j] * Fourier::Basis(j,phi);
    return result;
}