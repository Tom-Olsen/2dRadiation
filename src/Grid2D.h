#ifndef __INCLUDE_GUARD_Grid2D_h__
#define __INCLUDE_GUARD_Grid2D_h__
#include <iomanip>              // std::setprecision
#include <fstream>              // file input/output
#include <jsoncpp/json/json.h>  // everything about json files
#include "ControlFlow.hh"// used for template arguments
#include "Utility.hh"            // exit_on_error(char* msg)
#include "TensorTypes.hh"       // simple containers for rank 1-3 tensors
#include "Profiler.hh"

template<class Coord>
class Grid2D
{
public:
    // Members:
    int n1;
    int n2;
    int n12;
    double cfl=1;
    double d1;
    double d2;
    double dt;
    double start1;
    double start2;
    double end1;
    double end2;

    // Constructors:
    Grid2D() = delete;
    Grid2D(const Int2& n, const Coordinate2<Coord>& start, const Coordinate2<Coord>& end);
    Grid2D(const Grid2D& grid);

    // Setters:
    void SetCFL(double cfl_);

    // Grid2D Access Tools:
    int Index(int i, int j);
    int Index(int i, int j, int d);
    double xCoord(double i, double j);
    double yCoord(double i, double j);
    double rCoord(double i, double j);
    double rCoord(int ij);
    double phCoord(double i, double j);
    double x1Coord(double i, double j);
    double x2Coord(double i, double j);
    Coordinate2<xy> xyCoord(double i, double j);
    Coordinate2<rph> rphCoord(double i, double j);
    Coordinate2<Coord> x12Coord(double i, double j);
    double i(double x1);
    double j(double x2);
    Double2 ij(double x1, double x2);
    Double2 ij(const Coordinate2<Coord>& x12);

    double Integrate(double* values);

    // Domain Checks:
    bool OutsideDomain(const Coordinate2<Coord>& x12);
    bool OutsideDomain(const int i, const int j);

    // Write Data to json file:
    void WriteFrametoJson
    (float time, double* r, double* g, double* b, double* a,
     const int frameNumber, std::string directory, std::string name="");
};
#endif //__INCLUDE_GUARD_Grid2D_h__