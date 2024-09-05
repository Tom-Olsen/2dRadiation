#include <iostream>
#include "../src/Radiation.h"
using namespace std;

// Macros:
#define WRITE_DATA true
#define PRINT_SETUP true
#define PRINT_PROGRESS true
#define PRINT_RESULTS true
#define SAVE_ID false

Logger SphereWave(Stencil stencil, StreamingType streamingType, double cfl)
{
    // Create Radiation object:
    size_t nx = 301;
    size_t ny = 301;
    Coord start(-1, -1);
    Coord end(1, 1);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil streamingStencil(5, 0, false);

    // Config:
    Config config =
        {
            .name = "Sphere Wave 2d/" + StreamingName(streamingType) + " " + stencil.name + Format(cfl, 2) + "cfl",
            .t0 = 0,
            .simTime = 0.7,
            .writePeriod = 1.0,
            .updateFourierHarmonics = false,
            .keepSourceNodesActive = false,
            .writeData = WRITE_DATA,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .saveInitialData = SAVE_ID,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Intensities,
        };

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, config);

    // Initial Data:
    double emissionRadius = 0.05;
    for (size_t j = 0; j < grid.ny; j++)
        for (size_t i = 0; i < grid.nx; i++)
        {
            size_t ij = grid.Index(i, j);
            Coord xy = grid.xy(i, j);
            double r = xy.EuklNorm();
            // Fluid:
            radiation.kappa0[ij] = 0;
            radiation.kappa1[ij] = 0;
            radiation.kappaA[ij] = 0;
            radiation.eta[ij] = 0;
            radiation.ux[ij] = 0;
            radiation.uy[ij] = 0;
            if (r < emissionRadius)
            {
                radiation.isInitialGridPoint[ij] = true;
                // Initial data given by intensities:
                for(int d=0; d<stencil.nDir; d++)
                    radiation.initialI[radiation.Index(ij, d)] = 1;
                radiation.initialFluxAngle_IF[ij] = (r > 1e-6) ? fmod(MyAtan2(xy[2], xy[1]) + 2.0 * M_PI, 2.0 * M_PI) : 0;
            }
        }
    radiation.RunSimulation();
    return radiation.logger;
}
void SphereWaveAnalysis(int n)
{
    double cfl = 0.9;
    if(n==0) SphereWave(Stencil(20, 0), StreamingType::FlatFixed   , cfl);
    if(n==0) SphereWave(Stencil(16, 4), StreamingType::FlatAdaptive, cfl);
     
    if(n==1) SphereWave(Stencil(50, 0), StreamingType::FlatFixed   , cfl);
    if(n==1) SphereWave(Stencil(40,10), StreamingType::FlatAdaptive, cfl);
     
    if(n==1) SphereWave(Stencil(100,  0), StreamingType::FlatFixed   , cfl);
    if(n==1) SphereWave(Stencil( 80, 20), StreamingType::FlatAdaptive, cfl);
     
    if(n==0) SphereWave(Stencil(200, 0), StreamingType::FlatFixed   , cfl);
    if(n==0) SphereWave(Stencil(160,40), StreamingType::FlatAdaptive, cfl);
}

Logger Shadow(Stencil stencil, StreamingType streamingType, double cfl)
{
    // Create Radiation object:
    size_t nx = 191;
    size_t ny = 191;
    Coord start(-0.2, -0.2);
    Coord end(1.7, 1.7);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil streamingStencil(5, 0, false);

    // Config:
    Config config =
        {
            .name = "Shadow 2d/" + StreamingName(streamingType) + " " + stencil.name + Format(cfl, 2) + "cfl",
            .t0 = 0,
            .simTime = 1.5,
            .writePeriod = 2,
            .updateFourierHarmonics = false,
            .keepSourceNodesActive = true,
            .writeData = WRITE_DATA,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .saveInitialData = SAVE_ID,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, config);

    // Initial Data:
    Coord planetPos(0.75, 0.75);
    double sunRadius = 0.10;
    double planetRadius = 0.25;
    for (size_t j = 0; j < grid.ny; j++)
        for (size_t i = 0; i < grid.nx; i++)
        {
            size_t ij = grid.Index(i, j);
            Coord xy = grid.xy(i, j);
            double r = xy.EuklNorm();
            radiation.kappa0[ij] = 0;
            radiation.kappa1[ij] = 0;
            radiation.kappaA[ij] = 0;
            radiation.eta[ij] = 0;
            radiation.ux[ij] = 0;
            radiation.uy[ij] = 0;
            if (r < sunRadius)
            {
                radiation.isInitialGridPoint[ij] = true;
                radiation.initialE_LF[ij] = 1;
                radiation.initialFx_LF[ij] = 0;
                radiation.initialFy_LF[ij] = 0;
            }
            double dist = (planetPos - xy).EuklNorm();
            if (dist <= planetRadius)
                radiation.kappaA[ij] = 1e10;
        }
    radiation.RunSimulation();
    return radiation.logger;
}
void ShadowAnalysis(int n)
{
    double cfl = 0.9;
    if (n == 0) Shadow(Stencil(50,  0), StreamingType::FlatFixed, cfl);
    if (n == 1) Shadow(Stencil(100,  0), StreamingType::FlatFixed, cfl);
    if (n == 2) Shadow(Stencil(200,  0), StreamingType::FlatFixed, cfl);
    if (n == 3) Shadow(Stencil( 40, 10), StreamingType::FlatAdaptive, cfl);
    if (n == 4) Shadow(Stencil( 80, 20), StreamingType::FlatAdaptive, cfl);
    if (n == 5) Shadow(Stencil(160, 40), StreamingType::FlatAdaptive, cfl);
}

Logger Star(Stencil stencil, StreamingType streamingType, double cfl, double kappaA)
{
    // Create Radiation object:
    size_t nx = 401;
    size_t ny = 401;
    Coord start(-4, -4);
    Coord end(4, 4);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil streamingStencil(5, 0, false);
    
    // Config:
    Config config =
        {
            .name = "Star 2d/" + StreamingName(streamingType) + " " + stencil.name + Format(cfl, 2) + "cfl" + Format(kappaA, 0) + "kappaA",
            .t0 = 0,
            .simTime = 10,
            .writePeriod = 11,
            .updateFourierHarmonics = false,
            .keepSourceNodesActive = false,
            .writeData = WRITE_DATA,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .saveInitialData = SAVE_ID,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, config);

    // Initial Data:
    double starRadius = 1;
    double E0 = 1.0 - exp(-kappaA);
    double F1 = 0.325;
    for (size_t j = 0; j < grid.ny; j++)
        for (size_t i = 0; i < grid.nx; i++)
        {
            size_t ij = grid.Index(i, j);
            Coord xy = grid.xy(i, j);
            double r = xy.EuklNorm();
            radiation.kappa0[ij] = 0;
            radiation.kappa1[ij] = 0;
            radiation.ux[ij] = 0;
            radiation.uy[ij] = 0;
            if (r <= starRadius + grid.dx / 2.0)
            {
                radiation.isInitialGridPoint[ij] = true;
                radiation.kappaA[ij] = radiation.eta[ij] = kappaA / starRadius;
                radiation.initialE_LF[ij] = E0;
                radiation.initialFx_LF[ij] = 0;
                radiation.initialFy_LF[ij] = 0;
            }
            else
            {
                radiation.isInitialGridPoint[ij] = true;
                radiation.kappaA[ij] = radiation.eta[ij] = 0;
                radiation.initialE_LF[ij]  = E0 / (r * r);
                radiation.initialFx_LF[ij] = F1 / (r * r) * xy[1] / r;
                radiation.initialFy_LF[ij] = F1 / (r * r) * xy[2] / r;
            }
        }
    radiation.RunSimulation();
    return radiation.logger;
}
void StarAnalysis(int n)
{
    double cfl = 0.9;
    if (n ==  0) Star(Stencil(100,  0), StreamingType::FlatFixed   , cfl, 1);
    if (n ==  1) Star(Stencil(100,  0), StreamingType::FlatFixed   , cfl, 10);
    if (n ==  2) Star(Stencil(100,  0), StreamingType::FlatFixed   , cfl, 1e10);
    if (n ==  3) Star(Stencil(200,  0), StreamingType::FlatFixed   , cfl, 1);
    if (n ==  4) Star(Stencil(200,  0), StreamingType::FlatFixed   , cfl, 10);
    if (n ==  5) Star(Stencil(200,  0), StreamingType::FlatFixed   , cfl, 1e10);
    if (n ==  6) Star(Stencil( 80, 20), StreamingType::FlatAdaptive, cfl, 1);
    if (n ==  7) Star(Stencil( 80, 20), StreamingType::FlatAdaptive, cfl, 10);
    if (n ==  8) Star(Stencil( 80, 20), StreamingType::FlatAdaptive, cfl, 1e10);
    if (n ==  9) Star(Stencil(160, 40), StreamingType::FlatAdaptive, cfl, 1);
    if (n == 10) Star(Stencil(160, 40), StreamingType::FlatAdaptive, cfl, 10);
    if (n == 11) Star(Stencil(160, 40), StreamingType::FlatAdaptive, cfl, 1e10);
}

Logger BeamCrossing(Stencil stencil, StreamingType streamingType, double cfl)
{
    // Create Radiation object:
    size_t nx = 201;
    size_t ny = 101;
    Coord start(-0.5, -0.25);
    Coord end(0.5, 0.25);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil streamingStencil(5, 0, false);

    InitialDataType initialDataType = InitialDataType::Moments;
    if (streamingType == StreamingType::FlatFixed)
        initialDataType = InitialDataType::Intensities;

    // Config:
    Config config =
        {
            .name = "Beam Crossing 2d/" + StreamingName(streamingType) + "_" + stencil.name + Format(cfl, 2) + "cfl",
            .t0 = 0,
            .simTime = 0.85,
            .writePeriod = 1,
            .updateFourierHarmonics = false,
            .keepSourceNodesActive = true,
            .writeData = WRITE_DATA,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .saveInitialData = SAVE_ID,
            .streamingType = streamingType,
            .initialDataType = initialDataType,
        };

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, config);

    // Initial Data:
    Tensor2 dir0 = Tensor2(0.3,  0.1).EuklNormalized();
    Tensor2 dir1 = Tensor2(0.3, -0.1).EuklNormalized();

    // Find nearest directions in stencil:
    float dist0 = 1e10;
    float dist1 = 1e10;
    float d0 = -1;
    float d1 = -1;
    for(int d=0; d<stencil.nDir; d++)
    {
        float dist = (dir0 - stencil.C(d)).EuklNorm();
        if(dist < dist0)
        {
            dist0 = dist;
            d0 = d;
        }
        dist = (dir1 - stencil.C(d)).EuklNorm();
        if(dist < dist1)
        {
            dist1 = dist;
            d1 = d;
        }
    }

    for (size_t j = 0; j < grid.ny; j++)
        for (size_t i = 0; i < grid.nx; i++)
        {
            size_t ij = grid.Index(i, j);
            Coord xy = grid.xy(i, j);
            double x = xy[1];
            double y = xy[2];
            radiation.kappa0[ij] = 0;
            radiation.kappa1[ij] = 0;
            radiation.kappaA[ij] = 0;
            radiation.eta[ij] = 0;

            // Beam 0, from bottom to top:
            if ( -0.45 - grid.dx < x && x <= -0.45 && -0.20 < y && y < -0.15)
            {
                radiation.isInitialGridPoint[ij] = true;
                radiation.initialE_LF[ij] = 1;
                radiation.initialFx_LF[ij] = dir0[1];
                radiation.initialFy_LF[ij] = dir0[2];
                // Only used for FlatFixed streaming:
                radiation.initialI[radiation.Index(ij, d0)] = 1.0 / stencil.W(d0);
            }
            // Beam 1, from bottom to top:
            if ( -0.45 - grid.dx < x && x <= -0.45 && 0.15 < y && y < 0.20)
            {
                radiation.isInitialGridPoint[ij] = true;
                radiation.initialE_LF[ij] = 1;
                radiation.initialFx_LF[ij] = dir1[1];
                radiation.initialFy_LF[ij] = dir1[2];
                // Only used for FlatFixed streaming:
                radiation.initialI[radiation.Index(ij, d1)] = 1.0 / stencil.W(d1);
            }
        }
    radiation.RunSimulation();
    return radiation.logger;
}
void BeamCrossingAnalysis(int n)
{
    double cfl = 0.9;
    if(n == 0) BeamCrossing(Stencil(200,  0), StreamingType::FlatFixed   , cfl);
    if(n == 1) BeamCrossing(Stencil( 80, 20), StreamingType::FlatAdaptive, cfl);
    if(n == 2) BeamCrossing(Stencil(160, 40), StreamingType::FlatAdaptive, cfl);
}

Logger Diffusion(Stencil stencil, StreamingType streamingType, double kappaS, double lambda, double cfl, double correctionFactor)
{
    // Create Radiation object:
    size_t nx = 201;
    size_t ny = 201;
    Coord start(-0.5, -0.5);
    Coord end(0.5, 0.5);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil streamingStencil(5, 0, false);

    // Initial Data:
    // double lambda = 0.0;  // = 3kappa1 / kappa0
    double kappa0 = kappaS / (1.0 - lambda / 9.0);
    double kappa1 = kappa0 * lambda / 3.0;
    double PE = kappa0 * grid.dx;
    double D = 1.0 / (2.0 * kappa0) * (1.0 + correctionFactor * PE);

    // Config:
    double t0 = 1;
    std::string name = "Diffusion 2d/" + StreamingName(streamingType) + " " + stencil.name + " " + Format(cfl, 2, true) + "cfl " + std::to_string(nx) + "nx " + std::to_string(ny) + "ny " + std::to_string((int)kappa0) + "kappa0 " + std::to_string((int)kappa1) + "kappa1 " + Format(PE, 1, true) + "PE " + Format(correctionFactor, 3);
    Config config =
        {
            .name = name,
            .t0 = t0,
            .simTime = 2.0,
            .writePeriod = 1.0,
            .updateFourierHarmonics = false,
            .keepSourceNodesActive = false,
            .writeData = WRITE_DATA,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .saveInitialData = true,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, config);

    for (size_t j = 0; j < grid.ny; j++)
        for (size_t i = 0; i < grid.nx; i++)
        {
            size_t ij = grid.Index(i, j);
            Coord xy = grid.xy(i, j);
            double x = xy[1];
            double y = xy[2];
            double r = xy.EuklNorm();
            radiation.kappa0[ij] = kappa0;
            radiation.kappa1[ij] = kappa1;
            radiation.kappaA[ij] = 0;
            radiation.eta[ij] = 0;
            radiation.isInitialGridPoint[ij] = true;
            radiation.ux[ij] = 0.0;
            radiation.uy[ij] = 0.0;

            double E = 1.0 / t0 * exp(-r * r / (4.0 * D * t0));
            radiation.initialE_LF[ij] = E;
            radiation.initialFx_LF[ij] = (x * E) / (2.0 * t0 * (1.0 + correctionFactor * PE));
            radiation.initialFy_LF[ij] = (y * E) / (2.0 * t0 * (1.0 + correctionFactor * PE));
        }

    radiation.RunSimulation();
    return radiation.logger;
}
void DiffusionAnalysis(int n)
{
    double lambda = 0;
    double correctionFactor = 0.64;
    double cfl = 0.9;
    if (n == 0) Diffusion(Stencil(40, 0), StreamingType::FlatFixed   ,    100.0, lambda, cfl, correctionFactor);
    if (n == 1) Diffusion(Stencil(40, 0), StreamingType::FlatFixed   , 100000.0, lambda, cfl, correctionFactor);
    if (n == 2) Diffusion(Stencil(32, 8), StreamingType::FlatAdaptive,    100.0, lambda, cfl, correctionFactor);
    if (n == 3) Diffusion(Stencil(32, 8), StreamingType::FlatAdaptive, 100000.0, lambda, cfl, correctionFactor);
}

Logger MovingDiffusion(Stencil stencil, StreamingType streamingType, double kappaS, double lambda, double cfl, double correctionFactor, double ux)
{
    // Create Radiation object:
    size_t nx = 301;
    size_t ny = 101;
    Coord start(-1.0, -0.5);
    Coord end(2.0, 0.5);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil streamingStencil(5, 0, false);

    // Initial Data:
    // double lambda = 0.0;  // = 3kappa1 / kappa0
    double kappa0 = kappaS / (1.0 - lambda / 9.0);
    double kappa1 = kappa0 * lambda / 3.0;
    double PE = kappa0 * grid.dx;
    double D = 1.0 / (2.0 * kappa0) * (1.0 + correctionFactor * PE);

    // Account for time dilation:
    // We use t=1,2,3 when looking at the system from the lab frame (fluid is moving light gets dragged by the fluid)
    // We use t=(1,2,3)/gamma(ux) when looking at the system from the fluid frame (fluid is at rest)
    double simTime = 2.0;
    if (ux == 0)
    {
        double gamma = 1.0 / sqrt(1.0 - 0.5 * 0.5);
        simTime = simTime / gamma;
    }

    // Config:
    double t0 = 1;
    std::string name = "Moving Diffusion 2d/" + StreamingName(streamingType) + " " + stencil.name + " " + Format(cfl, 2, true) + "cfl " + std::to_string(nx) + "nx " + std::to_string(ny) + "ny " + std::to_string((int)kappa0) + "kappa0 " + std::to_string((int)kappa1) + "kappa1 " + Format(PE, 1, true) + "PE " + Format(ux, 3) + "ux";
    Config config =
        {
            .name = name,
            .t0 = t0,
            .simTime = simTime,
            .writePeriod = simTime / 2.0,
            .updateFourierHarmonics = false,
            .keepSourceNodesActive = false,
            .writeData = WRITE_DATA,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .saveInitialData = true,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, config);

    for (size_t j = 0; j < grid.ny; j++)
        for (size_t i = 0; i < grid.nx; i++)
        {
            size_t ij = grid.Index(i, j);
            Coord xy = grid.xy(i, j);
            double r = xy.EuklNorm();
            double x = xy[1];
            double y = xy[2];
            radiation.kappa0[ij] = kappa0;
            radiation.kappa1[ij] = kappa1;
            radiation.kappaA[ij] = 0;
            radiation.eta[ij] = 0;
            radiation.isInitialGridPoint[ij] = true;
            radiation.ux[ij] = ux;
            radiation.uy[ij] = 0.0;

            // Fluid Frame Moments:
            double E = 1.0 / t0 * exp(-r * r / (4.0 * D * t0));
            double Fx = (x * E) / (2.0 * t0 * (1.0 + correctionFactor * PE));
            double Fy = (y * E) / (2.0 * t0 * (1.0 + correctionFactor * PE));
            double Pxx = E / (2.0 * (1.0 + correctionFactor * PE));
            double Pxy = 0.0;
            double Pyx = 0.0;
            double Pyy = E / (2.0 * (1.0 + correctionFactor * PE));

            // Lab Frame Moments:
            Tensor3x3 boost = BoostMatrix(Tensor2(-ux, 0));
            Tensor3x3 Tff = Tensor3x3(E, Fx, Fy, Fx, Pxx, Pxy, Fy, Pyx, Pyy);
            Tensor3x3 Tlf(0);
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int I = 0; I < 3; I++)
                        for (int J = 0; J < 3; J++)
                            Tlf[{i, j}] += Tff[{I, J}] * boost[{i, I}] * boost[{j, J}];

            // Use FF E for both so that the scaling is the same
            radiation.initialE_LF[ij] = E;//Tlf[{0, 0}];
            radiation.initialFx_LF[ij] = Tlf[{1, 0}];
            radiation.initialFy_LF[ij] = Tlf[{2, 0}];
        }
        
    radiation.RunSimulation();
    return radiation.logger;
}
void MovingDiffusionAnalysis(int n)
{
    double lambda = 0;
    double cfl = 0.9;
    double correctionFactor = 0.64;

    if (n == 0) MovingDiffusion(Stencil(40, 0), StreamingType::FlatFixed   , 1000.0, lambda, cfl, correctionFactor, 0.0);
    if (n == 1) MovingDiffusion(Stencil(40, 0), StreamingType::FlatFixed   , 1000.0, lambda, cfl, correctionFactor, 0.5);
    if (n == 2) MovingDiffusion(Stencil(32, 8), StreamingType::FlatAdaptive, 1000.0, lambda, cfl, correctionFactor, 0.0);
    if (n == 3) MovingDiffusion(Stencil(32, 8), StreamingType::FlatAdaptive, 1000.0, lambda, cfl, correctionFactor, 0.5);
}

Logger CurvedBeam(Stencil stencil, StreamingType streamingType, double cfl, int nF = 5, bool updateFourierHarmonics = false)
{
    size_t nx = 251;
    size_t ny = 201;
    Coord start(0, 0);
    Coord end(5, 4);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(cfl);
    SchwarzSchild metric(grid, 1.0, 0.0); // needs at least stencil with order 5
    Stencil streamingStencil(nF, 0, false);
    
    // Config:
    string getsUpdated = (updateFourierHarmonics) ? " 1" : "";
    Config config =
        {
            .name = "Curved Beam 2d/" + StreamingName(streamingType) + " " + stencil.name + Format(cfl, 2) + "cfl " + std::to_string(nx) + "nx " + std::to_string(ny) + "ny " + std::to_string(nF) + "nF" + getsUpdated,
            .t0 = 0,
            .simTime = 10,
            .writePeriod = 11,
            .updateFourierHarmonics = updateFourierHarmonics,
            .keepSourceNodesActive = true,
            .writeData = WRITE_DATA,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .saveInitialData = SAVE_ID,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, config);
    // Turn this on to overwrite sigma (needed to create 'wrong' simultion for Appendix).
    // radiation.sigmaOverwrite = std::min(1.5 * stencil.sigmaMax, 500.0);

    // Initial Data:
    PARALLEL_FOR(2)
    for (size_t j = 0; j < grid.ny; j++)
        for (size_t i = 0; i < grid.nx; i++)
        {
            size_t ij = grid.Index(i, j);
            Coord xy = grid.xy(i, j);
            double y = xy[2];
            radiation.kappa0[ij] = 0;
            radiation.kappa1[ij] = 0;
            radiation.kappaA[ij] = 0;
            radiation.eta[ij] = 0;
            if (3.0 < y && y < 3.5 && i <= 1)
            {
                Tensor3 uLF(1, 1, 0);
                uLF = NullNormalize(uLF, metric.GetMetric_ll(ij));
                Tensor2 vLF = Vec2ObservedByEulObs<LF, LF>(uLF, xy, metric);

                radiation.isInitialGridPoint[ij] = true;
                radiation.initialE_LF[ij] = 1;
                radiation.initialFx_LF[ij] = 10 * vLF[1];
                radiation.initialFy_LF[ij] = 10 * vLF[2];
            }
        }
    radiation.RunSimulation();
    return radiation.logger;
}
void CurvedBeamAnalysis(int n)
{
    double cfl = 0.9;
    if (n == 0) CurvedBeam(Stencil( 50,  0), StreamingType::CurvedFixed  , cfl, 5);
    if (n == 1) CurvedBeam(Stencil(100,  0), StreamingType::CurvedFixed  , cfl, 5);
    if (n == 2) CurvedBeam(Stencil(200,  0), StreamingType::CurvedFixed  , cfl, 5);

    if (n == 3) CurvedBeam(Stencil( 40, 10), StreamingType::CurvedAdaptive, cfl, 5);
    if (n == 4) CurvedBeam(Stencil( 80, 20), StreamingType::CurvedAdaptive, cfl, 5);
    if (n == 5) CurvedBeam(Stencil(160, 40), StreamingType::CurvedAdaptive, cfl, 5);
}

Logger CurvedDiffusion(Stencil stencil, StreamingType streamingType, double kappaS, double lambda, double cfl)
{
    size_t nx = 401;
    size_t ny = 401;
    Coord start(-4, -4);
    Coord end(4, 4);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(cfl);
    SchwarzSchild metric(grid, 1.0, 0.0); // needs at least stencil with order 5
    Stencil streamingStencil(5, 0, false);
    
    // Initial Data:
    // double lambda = 0.0;  // = 3kappa1 / kappa0
    double kappa0 = kappaS / (1.0 - lambda / 9.0);
    double kappa1 = kappa0 * lambda / 3.0;
    double PE = kappa0 * grid.dx;
    double D = 1.0 / (2.0 * kappa0);

    // Config:
    Config config =
        {
            .name = "Curved Diffusion 2d/" + StreamingName(streamingType) + " " + stencil.name + Format(cfl, 2) + "cfl " + std::to_string(nx) + "nx " + std::to_string(ny) + "ny",
            .t0 = 0,
            .simTime = 3,
            .writePeriod = 4,
            .updateFourierHarmonics = false,
            .keepSourceNodesActive = false,
            .writeData = WRITE_DATA,
            .printSetup = PRINT_SETUP,
            .printProgress = PRINT_PROGRESS,
            .printResults = PRINT_RESULTS,
            .saveInitialData = SAVE_ID,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, config);

    // Initial Data:
    PARALLEL_FOR(2)
    for (size_t j = 0; j < grid.ny; j++)
        for (size_t i = 0; i < grid.nx; i++)
        {
            size_t ij = grid.Index(i, j);
            Coord xy = grid.xy(i, j);
            double r = xy.EuklNorm();
            double x = xy[1];
            double y = xy[2];
            radiation.kappa0[ij] = kappa0;
            radiation.kappa1[ij] = kappa1;
            radiation.kappaA[ij] = 0;
            radiation.eta[ij] = 0;
            radiation.ux[ij] = 0;//-0.5 * x / r;
            radiation.uy[ij] = 0;//-0.5 * y / r;
            if (2.2 < r && r < 2.5)
            {
                Tensor3 uLF(1, x, y);
                uLF = NullNormalize(uLF, metric.GetMetric_ll(ij));
                Tensor2 vLF = Vec2ObservedByEulObs<LF, LF>(uLF, xy, metric);

                radiation.isInitialGridPoint[ij] = true;
                radiation.initialE_LF[ij] = 1;
                radiation.initialFx_LF[ij] = 10 * vLF[1];
                radiation.initialFy_LF[ij] = 10 * vLF[2];
            }
        }
    radiation.RunSimulation();
    return radiation.logger;
}

void StencilAnalysis()
{
    Stencil( 20, 0).Print();
    Stencil( 50, 0).Print();
    Stencil(100, 0).Print();
    Stencil(200, 0).Print();

    Stencil( 16,  4).Print();
    Stencil( 40, 10).Print();
    Stencil( 80, 20).Print();
    Stencil(160, 40).Print();
}

void StreamingTypePerformanceAnalysis(int n)
{
    if (n == 0)
        CurvedBeam(Stencil(200, 0), StreamingType::FlatFixed, 0.75, 5, false);
    if (n == 1)
        CurvedBeam(Stencil(200, 0), StreamingType::CurvedFixed, 0.75, 3, false);
    if (n == 2)
        CurvedBeam(Stencil(200, 0), StreamingType::CurvedFixed, 0.75, 5, false);
    if (n == 3)
        CurvedBeam(Stencil(200, 0), StreamingType::CurvedFixed, 0.75, 7, false);
    if (n == 4)
        CurvedBeam(Stencil(200, 0), StreamingType::CurvedFixed, 0.75, 9, false);
    if (n == 5)
        CurvedBeam(Stencil(200, 0), StreamingType::CurvedFixed, 0.75, 3, true);
    if (n == 6)
        CurvedBeam(Stencil(200, 0), StreamingType::CurvedFixed, 0.75, 5, true);
    if (n == 7)
        CurvedBeam(Stencil(200, 0), StreamingType::CurvedFixed, 0.75, 7, true);
    if (n == 8)
        CurvedBeam(Stencil(200, 0), StreamingType::CurvedFixed, 0.75, 9, true);
    if (n == 9)
        CurvedBeam(Stencil(200, 0), StreamingType::GeodesicFixed, 0.75, 5, false);
}
void MetricDataForLukasCurvedBeam()
{
    // Needed for plots of Lukas M1 curved beam data:
    // Fx = alphas * Fx_data - Bx * E_data
    // Fy = alphas * Fy_data - By * E_data
    size_t nx = 200;
    size_t ny = 200;
    Coord start(-0.1, -1.0);
    Coord end(5.9695, 4.97);
    Grid grid(nx, ny, start, end);
    KerrSchild metric(grid, 1.0, 0.0);
    
    std::ofstream alphaOut("../output/alpha.txt");
    std::ofstream betaxOut("../output/betaX.txt");
    std::ofstream betayOut("../output/betaY.txt");

    for(size_t j = 0; j < grid.ny; j++)
    {
        size_t ij = grid.Index(0, j);
        alphaOut << metric.alpha[ij];
        betaxOut << metric.beta1_u[ij];
        betayOut << metric.beta2_u[ij];

        for(size_t i = 1; i < grid.nx; i++)
        {
            ij = grid.Index(i, j);
            alphaOut << " " << metric.alpha[ij];
            betaxOut << " " << metric.beta1_u[ij];
            betayOut << " " << metric.beta2_u[ij];
        }
        alphaOut << std::endl;
        betaxOut << std::endl;
        betayOut << std::endl;
    }

    alphaOut.close();
    betaxOut.close();
    betayOut.close();
}


int main(int argc, char *argv[])
{
    int n = 1;
    if (argc > 1)
        n = atoi(argv[1]);
        
    // Paper:
    // SphereWaveAnalysis(n);               // Done
    ShadowAnalysis(n);                   // Done
    // StarAnalysis(n);                     // Done
    // BeamCrossingAnalysis(n);             // Done
    // DiffusionAnalysis(n);                // Done
    // MovingDiffusionAnalysis(n);          // Done
    // CurvedBeamAnalysis(n);               // Done
    // StencilAnalysis();

    // StreamingTypePerformanceAnalysis(n); // Done
    // MetricDataForLukasCurvedBeam();
    
    CurvedBeam(Stencil( 80, 20), StreamingType::CurvedAdaptive, 0.9, 5);
}