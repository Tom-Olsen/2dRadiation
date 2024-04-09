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

public:
    Grid &grid;
    Metric &metric;
    Stencil &stencil;
    Stencil &streamingStencil;
    Config config;
    Logger logger;

    // Initial data is set from the outside (in code units):
    bool *isInitialGridPoint;
    // config.initialDataType = InitialDataType::Intensities
    RealBuffer initialI;
    RealBuffer initialFluxAngle_IF;
    // config.initialDataType = InitialDataType::Moments
    RealBuffer initialE_LF;
    RealBuffer initialFx_LF;
    RealBuffer initialFy_LF;
    RealBuffer initialPxx_LF;
    RealBuffer initialPxy_LF;
    RealBuffer initialPyy_LF;

    // Rotation angle of stencils:
    RealBuffer rotationAngle;
    RealBuffer rotationAngleNew;

    // Inertial frame moments:
    RealBuffer E;
    RealBuffer Fx;
    RealBuffer Fy;
    RealBuffer Pxx;
    RealBuffer Pxy;
    RealBuffer Pyy;

    // Lab frame moments:
    RealBuffer E_LF;
    RealBuffer Fx_LF;
    RealBuffer Fy_LF;
    RealBuffer Pxx_LF;
    RealBuffer Pxy_LF;
    RealBuffer Pyy_LF;
    RealBuffer F_LF; // only for saving

    // Fluid properties:
    RealBuffer kappa0;
    RealBuffer kappa1;
    RealBuffer kappaA;
    RealBuffer eta;
    RealBuffer ux;
    RealBuffer uy;

    // Population intensities:
    RealBuffer I;
    RealBuffer Inew;

    // Fourier coefficients for geodesic streaming:
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
    void LoadInitialData();
    void UpdateFourierCoefficients();
    void ComputeMomentsIF();
    void ComputeMomentsLF();
    void UpdateRotationMatrizes();

    void StreamFlatFixed();
    void StreamFlatAdaptive();
    void StreamCurvedFixed();
    void StreamCurvedAdaptive();
    void StreamGeodesicFixed();

    void Collide();

    void RunSimulation();
};
#endif //__INCLUDE_GUARD_Radiation_h__