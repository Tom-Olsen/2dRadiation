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
    static constexpr double LAMBDA_ITTERATION_TOLERENCE = 1e-4; // 0.01% error
    static constexpr int MAX_LAMBDA_ITERATIONS = 100;

public:
    Grid &grid;
    Metric &metric;
    Stencil &stencil;
    Stencil &streamingStencil;
    Config config;
    Logger logger;
    double sigmaOverwrite = -1;

    // Initial data is set from the outside (in code units):
    bool *isInitialGridPoint;
    // config.initialDataType = InitialDataType::Intensities
    DoubleBuffer initialI;
    DoubleBuffer initialFluxAngle_IF;
    // config.initialDataType = InitialDataType::Moments
    DoubleBuffer initialE_LF;
    DoubleBuffer initialFx_LF;
    DoubleBuffer initialFy_LF;
    DoubleBuffer initialPxx_LF;
    DoubleBuffer initialPxy_LF;
    DoubleBuffer initialPyy_LF;

    // Rotation angle of stencils:
    DoubleBuffer rotationAngle;
    DoubleBuffer rotationAngleNew;

    // Inertial frame moments:
    DoubleBuffer E;
    DoubleBuffer Fx;
    DoubleBuffer Fy;
    DoubleBuffer Pxx;
    DoubleBuffer Pxy;
    DoubleBuffer Pyy;

    // Lab frame moments:
    DoubleBuffer E_LF;
    DoubleBuffer Fx_LF;
    DoubleBuffer Fy_LF;
    DoubleBuffer Pxx_LF;
    DoubleBuffer Pxy_LF;
    DoubleBuffer Pyy_LF;
    DoubleBuffer F_LF; // only for saving

    // Fluid properties:
    DoubleBuffer kappa0;
    DoubleBuffer kappa1;
    DoubleBuffer kappaA;
    DoubleBuffer eta;
    DoubleBuffer ux;
    DoubleBuffer uy;

    // Population intensities:
    DoubleBuffer I;
    DoubleBuffer Inew;

    // Fourier coefficients for geodesic streaming:
    DoubleBuffer coefficientsS;
    DoubleBuffer coefficientsX;
    DoubleBuffer coefficientsY;
    DoubleBuffer coefficientsCx;
    DoubleBuffer coefficientsCy;

    // Testing:
    IntBuffer itterationCount;
    int maxItterationCount = 0;
    double averageItterationCount = 0;

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
    double NormalizationOverwrite();
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