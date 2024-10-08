#ifndef __INCLUDE_GUARD_Logger_hh__
#define __INCLUDE_GUARD_Logger_hh__
#include <math.h>       // Basic math.
#include "Spacetimes.h" // Metric data.
#include "Stencil.h"    // Velocity stencil.

class Logger
{
public:
    // General simulation parameters:
    Stencil &intensityStencil;
    Stencil &streamingStencil;
    Metric &metric;
    int timeSteps;
    double simTime;
    int maxItterationCount;
    double averageItterationCount;

    // Data management:
    std::string name;
    std::string date;
    std::string directoryPath;
    std::vector<std::string> timeNames;
    std::vector<double> timeMeasurements;

    Logger(Stencil &intensityStencil, Stencil &streamingStencil, Metric &metric)
        : intensityStencil(intensityStencil), streamingStencil(streamingStencil), metric(metric){};

    void SetValues(std::string name, double simTime)
    {
        this->name = name;

        // Derived from simulation parameters:
        timeSteps = ceil(simTime / metric.grid.dt);
        this->simTime = simTime = timeSteps * metric.grid.dt;

        // Creation time:
        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        std::ostringstream oss;
        oss << std::put_time(&tm, "%H.%M.%S - %d.%m.%Y");
        date = oss.str();

        // File system overhead:
        directoryPath = OUTPUTDIR + name;
    }

    void Write()
    {
        CreateDirectory(directoryPath);
        std::ofstream file(directoryPath + "/Log.txt");

        file << "Creation Date: " << date << std::endl << std::endl;

        // Spacetime:
        file << "Spacetime: " << metric.Name() << std::endl;
        file << "Black Hole  Mass,   M = " << metric.m << std::endl;
        file << "Black Hole  Spin,   a = " << metric.a << std::endl << std::endl;

        // Grid:
        file << "Grid Structure:" << std::endl;
        file << "Simulation Time = " << simTime << std::endl;
        file << "(startx,starty) = (" << metric.grid.startx << "," << metric.grid.starty << ")" << std::endl;
        file << "(endx  ,  endy) = (" << metric.grid.endx << "," << metric.grid.endy << ")" << std::endl;
        file << "(Nx    ,    Ny) = (" << metric.grid.nx << "," << metric.grid.ny << ")" << std::endl;
        file << "Nt         = " << timeSteps << std::endl;
        file << "dx         = " << metric.grid.dx << std::endl;
        file << "dy         = " << metric.grid.dy << std::endl;
        file << "dt         = " << metric.grid.dt << std::endl;
        file << "cfl        = " << metric.grid.GetCFL() << std::endl << std::endl;

        // Stencil:
        file << "Stencil Properties:" << std::endl;
        file << "Intensity Stencil: nDir          = " << intensityStencil.nDir << std::endl;
        file << "Intensity Stencil: nOrder        = " << intensityStencil.nOrder << std::endl;
        file << "Intensity Stencil: nCoefficients = " << intensityStencil.nCoefficients << std::endl;
        file << "Intensity Stencil: max Interpolation Error = " << 100 * intensityStencil.maxInterpolationError << "%" << std::endl;
        file << "Intensity Stencil: sigma max               = " << intensityStencil.sigmaMax << std::endl;
        file << "Intensity Stencil: relative flux max       = " << intensityStencil.relativeFluxMax << std::endl;
        file << "Streaming Stencil: nDir          = " << streamingStencil.nDir << std::endl;
        file << "Streaming Stencil: nOrder        = " << streamingStencil.nOrder << std::endl;
        file << "Streaming Stencil: nCoefficients = " << streamingStencil.nCoefficients << std::endl << std::endl;

        // Lambda Itteration:
        file << "Lambda Itteration:" << std::endl;
        file << "Max Lambda Itteration Count     = " << maxItterationCount << std::endl;
        file << "Average Lambda Itteration Count = " << averageItterationCount << std::endl << std::endl;

        // Times:
        file << "Time measurements:" << std::endl;
        for (int i = 0; i < timeMeasurements.size(); i++)
            file << timeNames[i] << ": " << timeMeasurements[i] << "s" << std::endl;
        file << "Computation Time: " << ComputationTime() << "s" << std::endl;
        
        // Benchmark:
        file << "Benchmark: " << metric.grid.nxy * timeSteps / (1e6 * ComputationTime()) << " MLUPS" << std::endl;
    }

    void AddTimeMeasurement(std::string timeName, double timeMeasurement)
    {
        timeNames.push_back(timeName);
        timeMeasurements.push_back(timeMeasurement);
    }

    double TotalTime()
    {
        for (int i = 0; i < timeMeasurements.size(); i++)
        {
            if (timeNames[i] == "Total Time")
                return timeMeasurements[i];
        }
        return 0;
    }
    double WritingTime()
    {
        for (int i = 0; i < timeMeasurements.size(); i++)
        {
            if (timeNames[i] == "void Grid::WriteFrametoCsv(float, const DoubleBuffer&, const DoubleBuffer&, const DoubleBuffer&, const DoubleBuffer&, std::string, std::string)")
                return timeMeasurements[i];
        }
        return 0;
    }
    double ComputationTime()
    {
        return TotalTime() - WritingTime();
    }
};

#endif //__INCLUDE_GUARD_Logger_hh__