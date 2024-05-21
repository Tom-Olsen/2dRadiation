#ifndef __INCLUDE_GUARD_Grid_h__
#define __INCLUDE_GUARD_Grid_h__
#include <fstream>        // File input/output.
#include "ControlFlow.hh" // Template arguments and profiling macros.
#include "Utility.hh"     // Utility functions.
#include "DataTypes.hh"   // General relativity tensors.
#include "Profiler.hh"    // Time measurement profiler.

class Grid
{
private:
    double m_cfl = 1;

public:
    // Members:
    size_t nx, ny, nxy;
    double dx, dy, dt;
    double startx, starty;
    double endx, endy;

    // Constructors:
    Grid() = delete;
    Grid(size_t nx_, size_t ny_, Coord start_, Coord end_);
    Grid(const Grid &grid);

    // Setters/Getters:
    void SetCFL(double cfl);
    double GetCFL();

    // Grid Access Tools:
    size_t Index(size_t i, size_t j);
    double x(size_t i);
    double y(size_t j);
    double x(double i);
    double y(double j);
    Coord xy(size_t i, size_t j);
    Coord xy(double i, double j);
    Coord xy(size_t ij);
    double i(double x);
    double j(double y);
    size_t i(size_t ij);
    size_t j(size_t ij);
    Coord ij(const Coord &xy);

    // Domain Checks:
    bool OutsideDomain(const Coord &xy);
    bool OutsideDomain(double i, double j);

    // Write Data to file:
    void WriteFrametoCsv(float time, const RealBuffer &r, const RealBuffer &g, const RealBuffer &b, const RealBuffer &a, std::string directory, std::string name = "");

    // Debugging:
    void Print();
};
#endif //__INCLUDE_GUARD_Grid_h__