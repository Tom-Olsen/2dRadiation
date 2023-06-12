#ifndef __INCLUDE_GUARD_Radiation_h__
#define __INCLUDE_GUARD_Radiation_h__
#include "GeodesicEquationSolver.h" // Solves geodesic equation, given xyz coordinates, LF 3 velocity, and metric.
#include "SpecialMath.h"            // More specific tensor operations, nullNormalize etc.
#include "FourierHarmonics.h"       // Fourier harmonic functions and fourier transforms.
#include "Log.hh"                   // log final results.
#include "Config.hh"                // Config for simulation parameters.



class Radiation
{
public:
	StreamingType streamingType;
	Grid& grid;
	Metric& metric;
	Stencil& stencil;
	Stencil& streamingStencil;

	double sigma = 1.0;

	bool* isInitialGridPoint;
	RealBuffer initialE;
	RealBuffer initialNx;
	RealBuffer initialNy;
	RealBuffer initialKappa0;
	RealBuffer initialKappa1;
	RealBuffer initialKappaA;
	RealBuffer initialEta;

	RealBuffer rotationAngle;
	RealBuffer rotationAngleNew;
	RealBuffer E;
	RealBuffer Fx;
	RealBuffer Fy;
	RealBuffer Pxx;
	RealBuffer Pxy;
	RealBuffer Pyy;
	RealBuffer E_LF;
	RealBuffer Fx_LF;
	RealBuffer Fy_LF;
	RealBuffer Pxx_LF;
	RealBuffer Pxy_LF;
	RealBuffer Pyy_LF;
    RealBuffer kappa0;
    RealBuffer kappa1;
    RealBuffer kappaA;
    RealBuffer eta;
	RealBuffer I;
	RealBuffer Inew;
	RealBuffer coefficientsS;
	RealBuffer coefficientsX;
	RealBuffer coefficientsY;
	RealBuffer coefficientsCx;
	RealBuffer coefficientsCy;

	Radiation() = delete;
	Radiation(Metric& metric, Stencil& stencil, Stencil& streamingStencil, StreamingType streamingType);
	~Radiation();

	size_t Index(size_t ij, size_t d);
	size_t Index(size_t i, size_t j, size_t d);

	size_t HarmonicIndex(size_t f, size_t ij);

	Coord GetTempCoordinate(size_t ij, double angle);
	Tensor2 GetTemp2VelocityIF(size_t ij, double angle);
	double GetFrequencyShift(size_t ij, double angle);
	double IntensityAt(size_t ij, Tensor2 vTempIF);
	Tensor2 AverageF(size_t i, size_t j);

	void NormalizeInitialDirections();
	void LoadInitialData();
	void NormalizeInitialIntensities();
	void UpdateFourierCoefficients();
	void ComputeMomentsIF();
	void ComputeMomentsLF();
	void UpdateRotationMatrizes();

	void StreamFlatFixed();
	void StreamFlatAdaptive();
	void StreamCurvedFixed();
	void StreamCurvedAdaptive();
	void Collide();

	void WriteIntensitiesToCsv(double time, const int frameNumber, std::string directory, std::string name);

	void RunSimulation(Config config);
};
#endif //__INCLUDE_GUARD_Radiation_h__