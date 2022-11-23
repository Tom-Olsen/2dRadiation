#include "SimulationData.h"



template<class Coord>
SimulationData<Coord>::SimulationData
(double simTime_, int nDir_, int nMom_, int nFourier_, double sigma_,
 Distribution& distribution_, Metric2D<Coord>& metric_, std::string name_):
 nDir(nDir_),
 nMom(nMom_),
 nFourier(nFourier_),
 sigma(sigma_),
 distribution(distribution_),
 metric(metric_),
 name(name_),
 uniStencil(nDir_),
 dirStencil(nDir_,distribution_),
 momentStencil(nMom_),
 fourierStencil(nFourier_)
{
    // Derived from simulation parameters:
	timeSteps = ceil(simTime_/metric.grid.dt);
	simTime = timeSteps * metric.grid.dt;

    // Initial data:
	initialE = new double[metric.grid.n12]();
	initialAngle = new double[metric.grid.n12]();
	initialKappa0 = new double[metric.grid.n12]();
	initialEta = new double[metric.grid.n12]();

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



template<class Coord>
SimulationData<Coord>::~SimulationData()
{
	delete[] initialE;
	delete[] initialAngle;
	delete[] initialKappa0;
	delete[] initialEta;
}



// formats a int number to 008 etc.
std::string formatNumber(int number)
{
	std::string formatedNumber;
	int maxDigits = 3;
	int digitsOfNumber = numDigits(number);

	for(int i=0; i<(maxDigits-digitsOfNumber); i++)
		formatedNumber += "0";
    if(number!=0)
    	formatedNumber += std::to_string(number);

	return formatedNumber;
}



template<class Coord>
void SimulationData<Coord>::AddTimeMeasurement(std::string timeName, double timeMeasurement)
{
    timeNames.push_back(timeName);
    timeMeasurements.push_back(timeMeasurement);
}



template<class Coord>
void SimulationData<Coord>::LogSimulationParameters()
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
	file << "nDir     = " << nDir << ";" << std::endl;
	file << "nMom     = " << nMom << ";" << std::endl;
	file << "nFourier = " << nFourier << ";" << std::endl;
	file << "Distribution Type : " << distribution.Name() << std::endl;
	file << "Normal distribution sigma = " << sigma << std::endl << std::endl;

	file << "Time measurements:" << std::endl;
	for(int i=0; i<timeMeasurements.size(); i++)
        file << timeNames[i] << ": " << timeMeasurements[i] << "s" << std::endl;
}



template struct SimulationData<xy>;
template struct SimulationData<rph>;