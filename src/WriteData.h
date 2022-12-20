#ifndef __INCLUDE_GUARD_WriteData_h__
#define __INCLUDE_GUARD_WriteData_h__
#include <fstream>	    // std::file
#include "Utility.hh"    // basic utility functions
#include "Grid2D.hh"	    // structure of numerical Grid
#include "Metric2D.h"     // metric data

// Write nFrames frames of the zeroth and first Momenta to directory/data_0...nFrames as .vtr file.
template<class Coord>
void Write_Radiation_vtr(double* E, double* Fx, double* Fy, Grid2D<Coord>& grid, int nFrames, std::string directory);

// Write nFrames frames of the zeroth and first Momenta to name.json file.
template<class Coord>
void Write_Radiation_json(double* E, double* Fx, double* Fy, Grid2D<Coord>& grid, int nFrames, std::string name);

#endif //__INCLUDE_GUARD_WriteData_h__