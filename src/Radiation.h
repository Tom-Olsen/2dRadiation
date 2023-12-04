#ifndef __INCLUDE_GUARD_Radiation_h__
#define __INCLUDE_GUARD_Radiation_h__
#include "GeodesicEquationSolver.h" // Solves geodesic equation, given xyz coordinates, LF 3 velocity, and metric.
#include "SpecialMath.h"            // More specific tensor operations, nullNormalize etc.
#include "FourierHarmonics.h"       // Fourier harmonic functions and fourier transforms.
#include "Logger.hh"                // logs final results.
#include "Config.hh"                // Config for simulation parameters.

class Radiation
{
private:
    // Constants:
    static constexpr int HALO = 1;
    static constexpr double MIN_FLUX_NORM = 1e-16;
    static constexpr double MIN_ENERGY_DENSITY = 1e-16;
    static constexpr double LAMBDA_ITTERATION_TOLERENCE = 1e-6;
    static constexpr int MAX_LAMBDA_ITERATIONS = 100;
    static constexpr double MAX_INTERPOLATION_ERROR = 0.01; // 1%

public:
    Grid &grid;
    Metric &metric;
    Stencil &stencil;
    Stencil &streamingStencil;
    Config config;
    Logger logger;

    // Set from the outside:
    bool *isInitialGridPoint;
    RealBuffer initialE_LF;
    RealBuffer initialFx_LF;
    RealBuffer initialFy_LF;
    RealBuffer initialPxx_LF;
    RealBuffer initialPxy_LF;
    RealBuffer initialPyy_LF;
    RealBuffer initialKappa0;
    RealBuffer initialKappa1;
    RealBuffer initialKappaA;
    RealBuffer initialEta;
    RealBuffer initialI;
    RealBuffer initialFluxAngle_IF;
    // kappa* and eta must be givne in CGS units
    // and will be converted to code units in the LoadInitialData() method.

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
    RealBuffer F_LF;
    RealBuffer Pxx_LF;
    RealBuffer Pxy_LF;
    RealBuffer Pyy_LF;
    RealBuffer kappa0;
    RealBuffer kappa1;
    RealBuffer kappaA;
    RealBuffer eta;
    RealBuffer ux;
    RealBuffer uy;
    RealBuffer I;
    RealBuffer Inew;
    RealBuffer coefficientsS;
    RealBuffer coefficientsX;
    RealBuffer coefficientsY;
    RealBuffer coefficientsCx;
    RealBuffer coefficientsCy;

    Radiation() = delete;
    Radiation(Metric &metric, Stencil &stencil, Stencil &streamingStencil, Config config);
    ~Radiation();

    size_t Index(size_t ij, size_t d);
    size_t Index(size_t i, size_t j, size_t d);

    size_t HarmonicIndex(size_t f, size_t ij);

    Coord GetTempCoordinate(size_t ij, double angle);
    Tensor2 GetTemp2VelocityIF(size_t ij, double angle);
    double GetFrequencyShift(size_t ij, double angle);
    double IntensityAt(size_t ij, Tensor2 vTempIF);
    Tensor2 AverageF(size_t i, size_t j);

    Tensor3 InitialDataLFtoIF(size_t ij);
    double SigmaMax();
    double FluxMax();
    void LoadInitialData();
    void UpdateFourierCoefficients();
    void ComputeMomentsIF();
    void ComputeMomentsLF();
    void UpdateRotationMatrizes();

    void StreamFlatFixed();
    void StreamFlatAdaptive();
    void StreamCurvedFixed();
    void StreamCurvedAdaptive();

    void CollideStaticFluidForwardEuler();
    void CollideStaticFluidBackwardEuler();
    void CollideForwardEuler();
    void CollideBackwardEuler();
    void CollideBackwardEulerTest();

    void RunSimulation();
};
#endif //__INCLUDE_GUARD_Radiation_h__