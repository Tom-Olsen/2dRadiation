#include <iostream>
#include <vector>
#include <array>
#include "../src/Utility.hh"
#include "../src/Stencil.h"
#include "../src/FourierHarmonics.h"
#include "../src/Interpolation.h"
using namespace std;

void UnitQuadrature(const Stencil &stencil)
{
    double E = 0;
    for (int d = 0; d < stencil.nDir; d++)
    {
        E += stencil.W(d) * 1.0;
    }
    PrintDouble(E, "E");
}

void Quadrature(const Stencil &stencil)
{
    cout << "nDir = " << stencil.nDir << endl;
    cout << "nOrder = " << stencil.nOrder << endl;
    cout << "nCoefficients = " << stencil.nCoefficients << endl;

    // Set original fourier coefficients:
    double *originalCoefficients = new double[stencil.nCoefficients]();
    originalCoefficients[0] = 5.0;
    for (size_t k = 1; k < stencil.nCoefficients; k++)
        originalCoefficients[k] = MySin(fmod(0.123 * M_PI * k, 2.0 * M_PI));

    // Calculate amplitude for each direction:
    double *data = new double[stencil.nDir]();
    for (size_t d = 0; d < stencil.nDir; d++)
    {
        double phi = stencil.Phi(d);
        data[d] = 0;
        for (int k = 0; k < stencil.nCoefficients; k++)
            data[d] += originalCoefficients[k] * Fourier::Basis(k, phi);
    }

    // Reconstruct coefficients with fourier transformation:
    double coefficients[stencil.nDir];
    Fourier::GetCoefficients(stencil, data, coefficients);

    // Print original coeffiients vs reconstructed coefficients:
    for (size_t k = 0; k < stencil.nCoefficients; k++)
        cout << "c" << k << ": " << Format(originalCoefficients[k]) << ", " << coefficients[k] << endl;
    cout << endl;

    // Use random test point to compare actual value to a reconstructed value:
    cout << "Phi, Real Value, Compressed Value:" << endl;
    double phi = RandomRange(0.0, 2.0 * M_PI);
    double value = 0;
    for (int k = 0; k < stencil.nCoefficients; k++)
        value += originalCoefficients[k] * Fourier::Basis(k, phi);
    cout << phi / M_PI << "pi, " << value << ", " << Fourier::GetValue(phi, coefficients, stencil.nCoefficients) << endl
         << endl;

    // Write original data points to file:
    std::ofstream file0(OUTPUTDIR + "sourcePoints.csv");
    file0 << "#x, y, z, value\n";
    for (int d = 0; d < stencil.nDir; d++)
        file0 << data[d] * cos(stencil.Phi(d)) << ", " << data[d] * sin(stencil.Phi(d)) << ", " << 0 << ", " << data[d] << "\n";
    file0.close();

    // Write many reconstructed points to file:
    std::ofstream file1(OUTPUTDIR + "fourierTransformPoints.csv");
    file1 << "#x, y, z, value\n";
    int m = 4 * stencil.nDir;
    for (int d = 0; d < m; d++)
    {
        double phi = 2.0 * M_PI * d / m;
        double value = Fourier::GetValue(phi, coefficients, stencil.nCoefficients);
        file1 << value * cos(phi) << ", " << value * sin(phi) << ", " << 0 << ", " << value << "\n";
    }
    file1.close();
}

void Interpolation(const Stencil &stencil)
{
    cout << "nDir = " << stencil.nDir << endl;
    cout << "nOrder = " << stencil.nOrder << endl;
    cout << "nCoefficients = " << stencil.nCoefficients << endl;

    // Set fourier coefficients:
    double *coefficients = new double[stencil.nCoefficients]();
    coefficients[0] = 10.0;
    for (size_t k = 1; k < stencil.nCoefficients; k++)
        coefficients[k] = 0.05 + MySin(fmod(0.123 * M_PI * sqrt(k), 2.0 * M_PI));

    // Calculate amplitude for each direction:
    double *data = new double[stencil.nDir]();
    for (size_t d = 0; d < stencil.nDir; d++)
    {
        double phi = stencil.Phi(d);
        data[d] = 0;
        for (int k = 0; k < stencil.nCoefficients; k++)
            data[d] += coefficients[k] * Fourier::Basis(k, phi);
    }

    // Calculate derivative of data:
    double *derivative = new double[stencil.nDir]();
    for (size_t d = 0; d < stencil.nDir; d++)
    {
        double delta = stencil.Phi((d + 1) % stencil.nDir) - stencil.Phi(d);
        double h = min(delta, 2.0 * M_PI - delta);
        derivative[d] = (data[(d + 1) % stencil.nDir] - data[d]) / h;
    }

    // Write original data points to file:
    std::ofstream file0(OUTPUTDIR + "sourcePoints.csv");
    file0 << "#x, y, z, value\n";
    for (int d = 0; d < stencil.nDir; d++)
        file0 << data[d] * cos(stencil.Phi(d)) << ", " << data[d] * sin(stencil.Phi(d)) << ", " << 0 << ", " << data[d] << "\n";
    file0.close();

    // Write derivatives to file:
    std::ofstream file1(OUTPUTDIR + "derivatives.csv");
    file1 << "#x, y, z, value\n";
    for (int d = 0; d < stencil.nDir; d++)
        file1 << derivative[d] * cos(stencil.Phi(d)) << ", " << derivative[d] * sin(stencil.Phi(d)) << ", " << 0 << ", " << derivative[d] << "\n";
    file1.close();

    // Write many reconstructed points to file:
    std::ofstream file2(OUTPUTDIR + "interpolatedPoints.csv");
    file2 << "#x, y, z, value\n";
    int m = stencil.interpolationGrid.nGrid;
    for (size_t i = 0; i < m; i++)
    {
        double phi = 2.0 * M_PI * (i + 0.5) / m;
        double d = stencil.interpolationGrid.d(phi);
        int d0 = std::floor(d);
        int d1 = (d0 + 1) % stencil.interpolationGrid.nGrid;
        std::array<size_t, 4> neighbourIndexesCubic0 = stencil.interpolationGrid.neighbourIndexesCubic[d0];
        std::array<size_t, 4> neighbourIndexesCubic1 = stencil.interpolationGrid.neighbourIndexesCubic[d1];
        std::array<double, 4> neighbourWeightsCubic0 = stencil.interpolationGrid.neighbourWeightsCubic[d0];
        std::array<double, 4> neighbourWeightsCubic1 = stencil.interpolationGrid.neighbourWeightsCubic[d1];

        double value0 = 0;
        double value1 = 0;
        for (size_t j = 0; j < 4; j++)
        {
            value0 += data[neighbourIndexesCubic0[j]] * neighbourWeightsCubic0[j];
            value1 += data[neighbourIndexesCubic1[j]] * neighbourWeightsCubic1[j];
        }
        double value = std::max(0.0, LinearInterpolation(d - d0, value0, value1));
        file2 << value * cos(phi) << ", " << value * sin(phi) << ", " << 0 << ", " << value << "\n";
    }
    file2.close();
}

void StencilAnalysis()
{
    std::ofstream fileOut("../output/fluxMax.csv");
    int nOrders[5] = {50, 100, 150, 200, 250};
    int nGhosts[5] = {0, 5, 10, 15, 20};

    for (int j = 0; j < 5; j++)
        for (int i = 0; i < 5; i++)
        {
            int nGhost = nGhosts[i];
            int nOrder = nOrders[j];
            Stencil stencil(nOrder, nGhost);

            std::string directory = "../output/" + std::to_string(nOrder) + "/";
            CreateDirectory(directory);

            stencil.fluxToSigmaTable.WriteToCsv(directory + "sigmaTable " + stencil.name, -1);
            fileOut << stencil.name << " " << stencil.relativeFluxMax << "\n";
        }
}

int main()
{
    UnitQuadrature(Stencil(20, 5));
    // Quadrature(Stencil(5,0));
    // Interpolation(Stencil(9,5));
    // Interpolation(Stencil(15,5));
    // Interpolation(Stencil(30,10));

    // Stencil(20, 5).interpolationGrid.Print();
    // Stencil(100, 0).Print();
    // Stencil(50,10).WriteToCsv();
    // Stencil(9,5).interpolationGrid.Print();

    // StencilAnalysis();

    // Stencil s(50, 0);
    // double sum = 0;
    // for (int d = 0; d < s.nDir; d++)
    //     sum += s.W(d);
    // cout << sum << endl;
    // s.fluxToSigmaTable.WriteToCsv("../output/sigmaTable", -1);
    // s.fluxToNormalizationTable.WriteToCsv("../output/normalizationTable", -1);
}