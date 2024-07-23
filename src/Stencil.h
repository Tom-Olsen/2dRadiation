#ifndef __INCLUDE_GUARD_Stencil_h__
#define __INCLUDE_GUARD_Stencil_h__
#include <vector>
#include <fstream>
#include "Utility.hh"
#include "DataTypes.hh"

struct Stencil;

struct InterpolationGrid
{
public:
    // Public Properties:
    size_t nGrid;
    Index2Buffer neighbourIndexesLinear;
    Double2Buffer neighbourWeightsLinear;
    Index4Buffer neighbourIndexesCubic;
    Double4Buffer neighbourWeightsCubic;

    // Constructor:
    InterpolationGrid() = default;
    InterpolationGrid(size_t nGrid, const Stencil &stencil);

    // Grid Index:
    double d(double phi) const;
    double Phi(size_t d) const;
    double Cx(size_t d) const;
    double Cy(size_t d) const;
    Tensor2 C(size_t d) const;

private:
    // Initialization:
    std::array<size_t, 2> NeighbourIndexesLinear(double phi, const Stencil &stencil);
    std::array<double, 2> LagrangeWeightsLinear(double phi, std::array<size_t, 2> neighbours, const Stencil &stencil);
    std::array<size_t, 4> NeighbourIndexesCubic(double phi, const Stencil &stencil);
    std::array<double, 4> LagrangeWeightsCubic(double phi, std::array<size_t, 4> neighbours, const Stencil &stencil);

public:
    // Debugging:
    void Print() const;
};

// Polar Coordinates convention:
// phi € [0,2pi]
// phi =     0 => right
// phi =  pi/2 => up
// phi =    pi => left
// phi = 3pi/2 => down
struct Stencil
{
public:
    // Public Properties:
    std::string name;     // Name of the stencel, e.g. "Stencil9.13"
    size_t nDir;          // Number of directions in stencil.
    size_t nGhost;        // Number of ghost directions in stencil.
    size_t nOrder;        // Quadrature integration.
    size_t nCoefficients; // Number of exact Fourier Harmonics that cna be integrated.
    InterpolationGrid interpolationGrid;
    LookUpTable<double> fluxToSigmaTable;
    LookUpTable<double> fluxToNormalizationTable;
    double sigmaMax;
    double relativeFluxMax;
    static constexpr double deltaPhi = M_PI / 8.0 + 1e-8; // +- angle range in which ghost directions are arranged.
    static constexpr double maxInterpolationError = 0.01; // 1%

private:
    // Internal Buffers:
    DoubleBuffer w;
    DoubleBuffer phi;
    DoubleBuffer cx;
    DoubleBuffer cy;
    // Initialization:
    void SortDirections();
    void PopulateLookUpTable();

public:
    // Constructor:
    Stencil(size_t nOrder, int nGhost = 0, bool isIntensityStencil = true);
    // Getters:
    double W(size_t d) const;
    double Phi(size_t d) const;
    double Cx(size_t d) const;
    double Cy(size_t d) const;
    Tensor2 C(size_t d) const;

public:
    // Debugging:
    void Print() const;
    void PrintAll() const;
    void WriteToCsv() const;
};
#endif //__INCLUDE_GUARD_Stencil_h__