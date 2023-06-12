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
    Index4Buffer  neighbourIndexes;
    Double4Buffer neighbourWeights;

    // Constructor:
    InterpolationGrid() = default;
    InterpolationGrid(size_t nGrid, const Stencil& stencil);

    // Grid Index:
    double d(double phi) const;
    
private:
    // Initialization:
    std::array<size_t,4> SurroundingNeighbourIndexes(double phi, const Stencil& stencil);
    std::array<double,4> LagrangeWeights(double phi, std::array<size_t,4> neighbours, const Stencil& stencil);

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
    std::string name;       // Name of the stencel, e.g. "Stencil9.13"
    size_t nDir;            // Number of directions in stencil.
    size_t nGhost;          // Number of ghost directions in stencil.
    size_t nOrder;          // Quadrature integration.
    size_t nCoefficients;   // Number of exact Fourier Harmonics that cna be integrated.
    InterpolationGrid interpolationGrid;

private:
    // Internal Buffers:
    RealBuffer w;
    RealBuffer phi;
    RealBuffer cx;
    RealBuffer cy;
    // Initialization:
    void SortDirections();

public:
    // Constructor:
    Stencil(size_t nOrder, int nGhost = 0);
    // Getters:
    double W(size_t d) const;
    double Phi(size_t d) const;
    double Cx(size_t d) const;
    double Cy(size_t d) const;
    Tensor2 C(size_t d) const;

public:
    // Debugging:
    void Print() const;
    void WriteToCsv() const;
};
#endif //__INCLUDE_GUARD_Stencil_h__