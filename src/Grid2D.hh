#ifndef __INCLUDE_GUARD_Grid2D_hh__
#define __INCLUDE_GUARD_Grid2D_hh__
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
    Grid2D(const Int2& n, const Coordinate2<Coord>& start, const Coordinate2<Coord>& end) :
    n1(n[0]), n2(n[1]), start1(start[1]), start2(start[2]), end1(end[1]), end2(end[2])
    {
       if constexpr(std::is_same<Coord,xy>::value)
       {
           if(end[1] < start[1])
               exit_on_error("Grid_xy x€[a,b] has b<a!");
           if(end[2] < start[2])
               exit_on_error("Grid_xy y€[a,b] has b<a!");
       }
       if constexpr(std::is_same<Coord,rph>::value)
       {
           if(start[1] < 0 or end[1] < 0)
               exit_on_error("Grid_rph radius range includes negative numbers!");
           if(end[2] - start[2] > 2.0*M_PI)
               exit_on_error("Grid_rph angle range is over 2pi!");
       }
       d1  = (end1-start1) / (n1-1);
       d2  = (end2-start2) / (n2-1);
       dt  = cfl*std::min(d1,d2);
       n12 = n1*n2;
       if(n1==1 or n2==1)
           exit_on_error("Grid2D must have at least 2 LP in each Dimension.");
    }
    Grid2D(const Grid2D& grid) :
    n1(grid.n1), n2(grid.n2), start1(grid.start1), start2(grid.start2), end1(grid.end1), end2(grid.end2)
    {
        d1 = (end1-start1) / (n1-1);
        d2 = (end2-start2) / (n2-1);
        dt = cfl*std::min(d1,d2);
        n12 =n1*n2;
    }

    // Setters:
    void SetCFL(double cfl_)
    {
        cfl = cfl_;
        dt = cfl*std::min(d1,d2);
    }

    // Indexing:
    INLINE int Index(int i, int j)
    { return i + j*n1; }
    INLINE int Index(int i, int j, int d)
    { return i + j*n1 + d*n12; }
    INLINE double i(double x1)
    { return (x1-start1)/d1; }
    INLINE double j(double x2)
    { return (x2-start2)/d2; }
    INLINE Double2 ij(double x1, double x2)
    {
        double i = (x1-start1)/d1;
        double j = (x2-start2)/d2;
        return Double2(i,j);
    }
    INLINE Double2 ij(const Coordinate2<Coord>& x12)
    {
        double i = (x12[1]-start1)/d1;
        double j = (x12[2]-start2)/d2;
        return Double2(i,j);
    }

    // Coordinates:
    INLINE double xCoord(double i, double j)
    {
        if constexpr(std::is_same<Coord,xy>::value)
            return start1 + i*d1;
        if constexpr(std::is_same<Coord,rph>::value)
        {
            double r  = start1 + i*d1;
            double ph = start2 + j*d2;
            return r*cos(ph);
        }
    }
    INLINE double yCoord(double i, double j)
    {
        if constexpr(std::is_same<Coord,xy>::value)
            return start2 + j*d2;
        if constexpr(std::is_same<Coord,rph>::value)
        {
            double r  = start1 + i*d1;
            double ph = start2 + j*d2;
            return r*sin(ph);
        }
    }
    INLINE double rCoord(double i, double j)
    {
        if constexpr(std::is_same<Coord,xy>::value)
        {
            double x = start1 + i*d1;
            double y = start2 + j*d2;
            return sqrt(x*x + y*y);
        }
        if constexpr(std::is_same<Coord,rph>::value)
            return start1 + i*d1;
    }
    INLINE double rCoord(int ij)
    {
        int i = ij % n1;
        int j = ij / n1;
        if constexpr(std::is_same<Coord,xy>::value)
        {
            double x = start1 + i*d1;
            double y = start2 + j*d2;
            return sqrt(x*x + y*y);
        }
        if constexpr(std::is_same<Coord,rph>::value)
            return start1 + i*d1;
    }
    INLINE double phCoord(double i, double j)
    {
        if constexpr(std::is_same<Coord,xy>::value)
        {
            double x = start1 + i*d1;
            double y = start2 + j*d2;
            return fmod(atan2(y,x) + 2.0*M_PI, 2.0*M_PI);
        }
        if constexpr(std::is_same<Coord,rph>::value)
            return start2 + j*d2;
    }
    INLINE double x1Coord(double i, double j)
    { return start1 + i*d1; }
    INLINE double x2Coord(double i, double j)
    { return start2 + j*d2; }
    INLINE Coordinate2<xy> xyCoord(double i, double j)
    {
        if constexpr(std::is_same<Coord,xy>::value)
            return Coordinate2<xy>(start1 + i*d1, start2 + j*d2);
        if constexpr(std::is_same<Coord,rph>::value)
        {
            double r  = start1 + i*d1;
            double ph = start2 + j*d2;
            return Coordinate2<xy>(r*cos(ph),r*sin(ph));
        }
    }
    INLINE Coordinate2<rph> rphCoord(double i, double j)
    {
        if constexpr(std::is_same<Coord,xy>::value)
        {
            double x = start1 + i*d1;
            double y = start2 + j*d2;
            return Coordinate2<rph>(sqrt(x*x + y*y),fmod(atan2(y,x) + 2.0*M_PI, 2.0*M_PI));
        }
        if constexpr(std::is_same<Coord,rph>::value)
            return Coordinate2<rph>(start1 + i*d1, start2 + j*d2);
    }
    INLINE Coordinate2<Coord> x12Coord(double i, double j)
    { return Coordinate2<Coord>(start1 + i*d1, start2 + j*d2); }

    // Domain Checks:
    INLINE bool OutsideDomain(const Coordinate2<Coord>& x12)
    {
        return
        x12[1]>end1 or x12[1]<start1 or
        x12[2]>end2 or x12[2]<start2;
    }
    INLINE bool OutsideDomain(const int i, const int j)
    {
        return
        i>n1-1 or i<0 or
        j>n2-1 or j<0;
    }

    // Write Data to json file:
    void WriteFrametoJson
    (float time, double* r, double* g, double* b, double* a,
     const int frameNumber, std::string directory, std::string name="")
    {
        PROFILE_FUNCTION();
        // main body:
        Json::Value jsonData;

        // primitive types:
        if constexpr(std::is_same<Coord,xy>::value)
            jsonData["meshType"] = 0;
        if constexpr(std::is_same<Coord,rph>::value)
            jsonData["meshType"] = 1;
        jsonData["time"] = time;

        // arrays:
        Json::Value colors(Json::arrayValue);
        for(int ij=0; ij<n12; ij++)
        {
            Json::Value Color;
            Color["r"] = r[ij];
            Color["g"] = g[ij];
            Color["b"] = b[ij];
            Color["a"] = a[ij];
            colors.append(Color);
        }
        jsonData["colors"] = colors;

        // structs:
        Json::Value start;
        start["x"] = start1;
        start["y"] = start2;
        jsonData["start"] = start;

        Json::Value end;
        end["x"] = end1;
        end["y"] = end2;
        jsonData["end"] = end;

        Json::Value resolution;
        resolution["x"] = n1;
        resolution["y"] = n2;
        jsonData["resolution"] = resolution;

        // write json to file:
        // TODO: check if directory exists and dont create it if not necessary
        CreateDirectory(directory);
        if(name == "")
        {
            std::ofstream fileOut(directory + "/data" + FrameNumber(frameNumber) + ".json");
            fileOut << jsonData;
            fileOut.close();
        }
        else
        {
            std::ofstream fileOut(directory + "/" + name + ".json");
            fileOut << jsonData;
            fileOut.close();
        }
    }
};
#endif //__INCLUDE_GUARD_Grid2D_hh__