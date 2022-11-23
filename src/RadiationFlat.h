#ifndef __INCLUDE_GUARD_RadiationFlat_h__
#define __INCLUDE_GUARD_RadiationFlat_h__
#include <math.h>
#include <omp.h>
#include <fstream>
#include "ControlFlowClasses.hh"	// used for template arguments
#include "Utility.h"				// basic utility functions
#include "Stencil.hh"				// velocity stencils.
#include "TensorTypes.hh" 			// simple containers for rank 1-3 tensors.
#include "Grid2D.h"					// 2D Grid for numerical domain and maping to physical domain.
#include "Metric2D.h"				// 2D Metric data. The numerical domain is defined by Grid3D.
#include "FourierHarmonics.h"		// fourier interpolation thorugh discrete number of points in 2d plane.
#include "GeodesicEquationSolver.h"	// used to solve geodesic equation on spacelike hypersurface.
#include "SimulationData.h"			// Summary of simulation parameters.



template<class Coord>
class RadiationFlat
{
public:
    SimulationData<Coord>& simData;
	Grid2D<Coord>& grid;
	int nDir;
	Stencil& uniStencil;
	Stencil& dirStencil;
	Stencil& momentStencil;
	double* rotation;
	double* E;
	double* Fx;
	double* Fy;
	double* Pxx;
	double* Pxy;
	double* Pyy;
	double* Fnorm;
    double* kappa0;
    double* kappa1;
    double* kappaA;
    double* eta;
	double* I;
	double* temp;

private:
	inline double InterpolateIntensity(int ij, double phi);

public:
	RadiationFlat() = delete;
	RadiationFlat(SimulationData<Coord>& simData_);
	~RadiationFlat();

	void LoadInitialData();
	void NormalizeInitialData();
	void InitVakuumMicrophysics();
	void ComputeMoments();
	void RotateStencil();
	void Stream();
	void Collide();

	void RunSimulation(int writeFrequency = 1, bool updateFourierHarmonics = false, bool keepSourceNodesActive = true, bool writeData = true, bool printToTerminal = true);
};
#endif //__INCLUDE_GUARD_RadiationFlat_h__