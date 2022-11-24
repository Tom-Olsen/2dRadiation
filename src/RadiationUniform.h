#ifndef __INCLUDE_GUARD_RadiationUniform_h__
#define __INCLUDE_GUARD_RadiationUniform_h__
#include <math.h>
#include <omp.h>
#include <fstream>
#include "ControlFlow.hh"			// used for template arguments
#include "Utility.hh"				// basic utility functions
#include "Stencil.hh"				// velocity stencils.
#include "TensorTypes.hh" 			// simple containers for rank 1-3 tensors.
#include "Grid2D.h"					// 2D Grid for numerical domain and maping to physical domain.
#include "Metric2D.h"				// 2D Metric data. The numerical domain is defined by Grid3D.
#include "FourierHarmonics.h"		// fourier interpolation thorugh discrete number of points in 2d plane.
#include "GeodesicEquationSolver.h"	// used to solve geodesic equation on spacelike hypersurface.
#include "SimulationData.h"			// Summary of simulation parameters.



template<class Coord>
class RadiationUniform
{
public:
    SimulationData<Coord>& simData;
	Grid2D<Coord>& grid;
	Metric2D<Coord>& metric;
	int nDir;
	Stencil& stencil;
	Stencil& fourierStencil;
	double* E;
	double* Fx;
	double* Fy;
	double* Pxx;
	double* Pxy;
	double* Pyy;
	double* E_LF;
	double* Fx_LF;
	double* Fy_LF;
	double* Pxx_LF;
	double* Pxy_LF;
	double* Pyy_LF;
	double* Fnorm_LF;
    double* kappa0;
    double* kappa1;
    double* kappaA;
    double* eta;
	double* I;
	double* temp;
	double* coefficientsS;
	double* coefficientsX1;
	double* coefficientsX2;
	double* coefficientsCx1;
	double* coefficientsCx2;

private:
	Coordinate2<Coord> GetTempCoordinate(int i, int j, double phi);
	Tensor2<Coord,LF> GetTemp2Velocity(int i, int j, double phi);
	double GetFrequencyShift(int i, int j, double phi);
	double InterpolateIntensity(int ij, double phi);
	void SaveStreamingFourierExtrapolation();

public:
	RadiationUniform() = delete;
	RadiationUniform(SimulationData<Coord>& simData_);
	~RadiationUniform();

	void LoadInitialData();
	void NormalizeInitialData();
	void UpdateFourierCoefficients();
	void InitVakuumMicrophysics();
	void ComputeMomentsIF();
	void ComputeMomentsLF();
	void Stream();
	void GeodesicStream();
	void FlatStream();
	void Collide();

	void RunSimulation(int writeFrequency = 1, bool updateFourierHarmonics = false, bool keepSourceNodesActive = true, bool writeData = true, bool printToTerminal = true);
};
#endif //__INCLUDE_GUARD_RadiationUniform_h__