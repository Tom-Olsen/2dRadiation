#ifndef __INCLUDE_GUARD_Log_hh__
#define __INCLUDE_GUARD_Log_hh__
#include <math.h>
#include "Grid2D.h"
#include "Metric2D.h"
#include "Stencil.hh"


template<class Coord>
class Log
{
public:
    // General simulation parameters:
    int timeSteps;
    float simTime;
    Stencil& stencil;
    Stencil& fourierStencil;
    Metric2D<Coord>& metric;

    // Data management:
    std::string name;
    std::string date;
    std::string directoryPath;
    std::vector<std::string> timeNames;
    std::vector<double> timeMeasurements;

    Log(std::string name_, double simTime_, Stencil& stencil_, Stencil& fourierStencil_, Metric2D<Coord>& metric_) :
    name(name_), stencil(stencil_), fourierStencil(fourierStencil_), metric(metric_)
    {
        // Derived from simulation parameters:
	    timeSteps = ceil(simTime_/metric_.grid.dt);
	    simTime = timeSteps * metric_.grid.dt;

	    // Creation time:
        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        std::ostringstream oss;
        oss << std::put_time(&tm, "%H.%M.%S - %d.%m.%Y");
	    date = oss.str();

	    // File system overhead:
	    directoryPath = "output/" + name;
	    if(!CreateDirectory(directoryPath))
	    	double a = 0;
	    	// exit_on_error("Failed to create output directory.");
    }
    ~Log()
    {
        std::ofstream file(directoryPath+"/Log.txt");

	    file << "Creation Date: " << date << std::endl << std::endl;

	    file << "Spacetime: " << metric.Name() << std::endl;
	    file << "Black Hole  Mass,   M = " << metric.m << std::endl;
	    file << "Black Hole  Spin,   a = " << metric.a << std::endl << std::endl;

	    file << "Grid Structure:" << std::endl;
	    file << "Simulation Time = " << simTime << std::endl;
	    file << "(start0,start1) = (" << metric.grid.start1 << "," << metric.grid.start2 << ");" << std::endl;
	    file << "(end0  ,end1  ) = (" << metric.grid.end1   << "," << metric.grid.end2   << ");" << std::endl;
	    file << "(N0    ,N1    ) = (" << metric.grid.n1     << "," << metric.grid.n2     << ");" << std::endl;
	    file << "N12 = " << metric.grid.n12 << std::endl;
	    file << "Nt  = " << timeSteps << std::endl;
	    file << "dx  = " << metric.grid.d1  << std::endl;
	    file << "dy  = " << metric.grid.d2  << std::endl;
	    file << "dt  = " << metric.grid.dt  << std::endl;
	    file << "cfl = " << metric.grid.cfl << std::endl << std::endl;
    
	    file << "Stencil Properties:" << std::endl;
	    file << "nDir     = " << stencil.nDir << ";" << std::endl;
	    file << "nFourier = " << fourierStencil.nDir << ";" << std::endl;
	    file << "Normal distribution sigma = " << stencil.sigma << std::endl << std::endl;

	    file << "Time measurements:" << std::endl;
	    for(int i=0; i<timeMeasurements.size(); i++)
            file << timeNames[i] << ": " << timeMeasurements[i] << "s" << std::endl;
    }
    


    void AddTimeMeasurement(std::string timeName, double timeMeasurement)
    {
        timeNames.push_back(timeName);
        timeMeasurements.push_back(timeMeasurement);
    }
};

#endif //__INCLUDE_GUARD_Log_hh__