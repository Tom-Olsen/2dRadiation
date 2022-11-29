#ifndef __INCLUDE_GUARD_Radiation_h__
#define __INCLUDE_GUARD_Radiation_h__
#include <math.h>
#include <omp.h>
#include <fstream>
#include <vector>
#include <string>
#include "ControlFlow.hh"			// used for template arguments
#include "Utility.hh"				// basic utility functions
#include "Stencil.hh"				// velocity stencils.
#include "TensorTypes.hh" 			// simple containers for rank 1-3 tensors.
#include "Grid2D.h"					// 2D Grid for numerical domain and maping to physical domain.
#include "Metric2D.h"				// 2D Metric data. The numerical domain is defined by Grid3D.
#include "FourierHarmonics.h"		// fourier interpolation thorugh discrete number of points in 2d plane.
#include "GeodesicEquationSolver.h"	// used to solve geodesic equation on spacelike hypersurface.
#include "Log.hh"
#include "Utility.hh"



// Input system:
enum StreamingType { StreamCurvedRot, StreamCurvedNoRot, StreamFlatRot, StreamFlatNoRot };
std::string StreamingName(int n);


struct RunParameters
{
    std::string name;
    double simTime;
    int writeFrequency;
    bool updateFourierHarmonics;
    bool keepSourceNodesActive;
    bool writeData;
    bool printToTerminal;
};



template<class Coord>
class Radiation
{
public:
	StreamingType streamingType;
	Grid2D<Coord>& grid;
	Metric2D<Coord>& metric;
	Stencil& stencil;
	Stencil& fourierStencil;
	int nDir;

	bool* isInitialGridPoint;
	double* initialE;
	double* initialRotation;
	double* initialKappa0;
	double* initialKappa1;
	double* initialKappaA;
	double* initialEta;

	double* rotation;
	double* rotationNew;
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
	inline __attribute__((always_inline)) Coordinate2<Coord> GetTempCoordinate(int i, int j, double phi);
	inline __attribute__((always_inline)) Tensor2<Coord,LF> GetTemp2Velocity(int i, int j, double phi);
	inline __attribute__((always_inline)) double GetFrequencyShift(int i, int j, double phi);
	inline __attribute__((always_inline)) double IntensityAt(int ij, double phi);
	inline __attribute__((always_inline)) Tensor2<Coord,IF> AverageF(int i, int j);

public:
	Radiation() = delete;
	Radiation(Grid2D<Coord>& grid_, Metric2D<Coord>& metric_, Stencil& stencil_, Stencil& fourierStencil_, StreamingType streamingType_);
	~Radiation();

	void LoadInitialData();
	void NormalizeInitialData();
	void UpdateFourierCoefficients();
	void ComputeMomentsIF();
	void ComputeMomentsLF();
	void StreamCurvedRot();
	void StreamCurvedNoRot();
	void StreamFlatRot();
	void StreamFlatNoRot();
	void Collide();

	void RunSimulation(RunParameters params);
};
#endif //__INCLUDE_GUARD_Radiation_h__