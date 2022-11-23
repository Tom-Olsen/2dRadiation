#include <iostream> 		// cin/cout
#include <math.h>   		// basic maths
#include <fstream>  		// file input/output
#include <stdlib.h> 		// srand, rand
#include <chrono>   		// time
#include <filesystem>		// folder/file management
#include <iomanip>			// std::put:time
#include "src/Includes.hh"	// all src files for this project

template<class Coord, int FourierOrder>
void SetInitialData1(SimulationData<Coord,FourierOrder>& simData)
{
	Metric2D<Coord>& metric = simData.metric;
	Grid2D<Coord>& grid = metric.grid;
	Stencil& stencil = simData.stencil;
	for(int j=0; j<grid.n2; j++)
	{
		int i = 1;
		if(3.0 < grid.yCoord(i,j) and grid.yCoord(i,j) < 3.5)
		{
			// Set Initial Data in LF, single directions:
			Coordinate2<Coord> x = grid.xyCoord(i,j);
			Tensor3<Coord,LF> u(1,1,0);
			u.NullNormalize(metric.GetMetric_ll(x));

			Tensor2<Coord,IF> v = Vec2ObservedByEulObs<Coord,LF,IF>(u,x,simData.metric);
			Tensor2<xy,IF> vxy = v.template Transform<xy>(x);
			int d = GetNearestVector(vxy,stencil);


			double I = 1.0;
			simData.Iinitial[grid.Index(i,j,d)] = I;
		}
	}
}



template<class Coord, int FourierOrder>
void SetInitialData2(SimulationData<Coord,FourierOrder>& simData)
{
	Metric2D<Coord>& metric = simData.metric;
	Grid2D<Coord>& grid = metric.grid;
	Stencil& stencil = simData.stencil;
	for(int i=0; i<grid.n1; i++)
	{
		int j = grid.n2-2;
		if(3.0 < grid.rCoord(i,j) and grid.rCoord(i,j) < 3.5)
		{
			// Set Initial Data in LF, single directions:
			Coordinate2<Coord> x = grid.x12Coord(i,j);
			Tensor3<Coord,LF> u(1,0,-1);
			u.NullNormalize(metric.GetMetric_ll(x));

			Tensor2<Coord,IF> v = Vec2ObservedByEulObs<Coord,LF,IF>(u,x,simData.metric);
			Tensor2<xy,IF> vxy = v.template Transform<xy>(x);
			int d = GetNearestVector(vxy,stencil);

			double I = 1.0;
			simData.Iinitial[grid.Index(i,j,d)] = I;
		}
	}
}



// TODOTOM:
// -Transformation of moments from IF to LF 
// -Adaptive TE, use Kretschmann scalar
// -Dispersion Plot
// -Momenta integral corrct? 1/2pi?
int main(int argc, char *argv[])
{
	if(argc < 4)
		exit_on_error("Too few input arguments!\nInput: (int)setup (int)nGrid1 (int)nGrid2 (int)simTime");
	int setup = atoi(argv[1]);
	int n1 = atoi(argv[2]);
	int n2 = atoi(argv[3]);
	float simTime = atof(argv[4]);

	Stencil& stencil = DirectedStencil<50>::GetInstance();
	Stencil& uniStencil = UniformStencil<50>::GetInstance();
	UniformStencil<5>& fourierStencil = UniformStencil<5>::GetInstance();

	if (setup == 1)
	{
		Int2 N(n1,n2);
    	Coordinate2<xy> Start(0,0);
    	Coordinate2<xy> End(4,4);
    	Grid2D<xy> grid(N,Start,End);
		SchwarzSchild<xy> metric(grid,1.0,0.0);
		SimulationData simData(metric,stencil,uniStencil,fourierStencil,simTime,"run1");
		SetInitialData1(simData);

		Radiation radiation(simData);
		radiation.RunSimulation();
	}
	if (setup == 2)
	{
		Int2 N(n1,n2);
    	Coordinate2<rph> Start(2,0);
    	Coordinate2<rph> End(4,M_PI_2);
    	Grid2D<rph> grid(N,Start,End);
    	SchwarzSchild<rph> metric(grid,1.0,0.0);
		SimulationData simData(metric,stencil,uniStencil,fourierStencil,simTime,"run2");
		SetInitialData2(simData);

		Radiation radiation(simData);
		radiation.RunSimulation();
	}
}