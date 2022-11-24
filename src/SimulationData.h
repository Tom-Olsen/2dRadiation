#ifndef __INCLUDE_GUARD_SimulationData_h__
#define __INCLUDE_GUARD_SimulationData_h__
#include <fstream>                  // std::file
#include <vector>                   // std::vector datastructure (analog to list in c#)
#include "ControlFlow.hh"	        // used for template arguments
#include "Stencil.hh"               // velocity stencil
#include "Grid2D.h"                 // structure of numerical Grid
#include "Metric2D.h"	            // Metric data
#include "Utility.hh"	            // basic utility


template<class Coord>
struct SimulationData
{
    // General simulation parameters:
    float simTime;
    int nDir;
    int nMom;
    int nFourier;
    float sigma;
    Distribution& distribution;
    Metric2D<Coord>& metric;

    // Derived from simulation parameters:
    int timeSteps;
	UniformStencil  uniStencil;
	DirectedStencil dirStencil;
	UniformStencil  momentStencil;
	UniformStencil  fourierStencil;

    // Initial data:
    double* initialE;
    double* initialAngle;
    double* initialKappa0;
    double* initialEta;

    // Data management:
    std::string name;
    std::string date;
    std::string directoryPath;
    std::vector<std::string> timeNames;
    std::vector<double> timeMeasurements;

    // Constructor/Destructor:
    SimulationData() = delete;
    SimulationData(double simTime_, int nDir_, int nMom_, int nFourier_, double sigma_,
    Distribution& distribution_, Metric2D<Coord>& metric_, std::string name_);
    //(double sigma_, int nDir_, int nMom_, int nFourier_, Metric2D<Coord>& metric_, double simTime_, std::string name_);
    ~SimulationData();

    // Methods:
    void AddTimeMeasurement(std::string timeName, double timeMeasurement);
    void LogSimulationParameters();
};


#endif //__INCLUDE_GUARD_SimulationData_h__