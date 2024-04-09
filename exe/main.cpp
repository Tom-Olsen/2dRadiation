#include <iostream>
#include "../src/Radiation.h"
using namespace std;

void StraightBeam(Stencil stencil, StreamingType streamingType, double cfl)
{
    // Create Radiation object:
    size_t nx = 900 + 1 + 2;
    size_t ny = 300 + 1 + 2;
    double dx = 3.0 / (nx - 1.0 - 2.0);
    double dy = 1.0 / (ny - 1.0 - 2.0);
    Coord start(0 - dx, -0.5 - dy);
    Coord end(3 + dx, 0.5 + dy);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil streamingStencil(5, 0, false);

    // Config:
    Config config =
        {
            .name = "Straight Beam 2d/" + StreamingName(streamingType) + " " + stencil.name + Format(cfl, 2) + "cfl",
            .t0 = 0,
            .simTime = 2.8,
            .writePeriod = 3, // write first and last frame only.
            .updateFourierHarmonics = false,
            .keepSourceNodesActive = true,
            .writeData = true,
            .printToTerminal = true,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, config);

    // Initial Data:
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
            if (-0.25 < y && y < 0.25 && i <= 1)
            {
                radiation.isInitialGridPoint[ij] = true;
                radiation.initialE_LF[ij] = 1;
                radiation.initialFx_LF[ij] = 1;
                radiation.initialFy_LF[ij] = 0;
            }
        }
    radiation.RunSimulation();
}

void BeamCrossing(Stencil stencil, StreamingType streamingType, double cfl)
{
    // Create Radiation object:
    size_t nx = 200 + 1 + 2;
    size_t ny = 100 + 1 + 2;
    double dx = 1.0 / (nx - 1.0 - 2.0);
    double dy = 0.5 / (ny - 1.0 - 2.0);
    Coord start(-0.5 - dx, -0.25 - dy);
    Coord end(0.5 + dx, 0.25 + dy);
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
            .writePeriod = 1, // write first and last frame only.
            .updateFourierHarmonics = false,
            .keepSourceNodesActive = true,
            .writeData = true,
            .printToTerminal = true,
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
            if ( -0.45 - dx < x && x <= -0.45 && -0.20 < y && y < -0.15)
            {
                radiation.isInitialGridPoint[ij] = true;
                radiation.initialE_LF[ij] = 1;
                radiation.initialFx_LF[ij] = dir0[1];
                radiation.initialFy_LF[ij] = dir0[2];
                // Only used for FlatFixed streaming:
                radiation.initialI[radiation.Index(ij, d0)] = 1.0 / stencil.W(d0);
            }
            // Beam 1, from bottom to top:
            if ( -0.45 - dx < x && x <= -0.45 && 0.15 < y && y < 0.20)
            {
                size_t ij = grid.Index(i, j);
                radiation.isInitialGridPoint[ij] = true;
                radiation.initialE_LF[ij] = 1;
                radiation.initialFx_LF[ij] = dir1[1];
                radiation.initialFy_LF[ij] = dir1[2];
                // Only used for FlatFixed streaming:
                radiation.initialI[radiation.Index(ij, d1)] = 1.0 / stencil.W(d1);
            }
        }
    radiation.RunSimulation();
}

void StraightBeamShadow(Stencil stencil, StreamingType streamingType, double cfl)
{
    // Create Radiation object:
    size_t nx = 900 + 1 + 2;
    size_t ny = 300 + 1 + 2;
    double dx = 3.0 / (nx - 1.0 - 2.0);
    double dy = 1.0 / (ny - 1.0 - 2.0);
    Coord start(0 - dx, -0.5 - dy);
    Coord end(3 + dx, 0.5 + dy);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil streamingStencil(5, 0, false);

    // Config:
    Config config =
        {
            .name = "Straight Beam Shadow 2d/" + StreamingName(streamingType) + " " + stencil.name + Format(cfl, 2) + "cfl",
            .t0 = 0,
            .simTime = 2.8,
            .writePeriod = 3, // write first and last frame only.
            .updateFourierHarmonics = false,
            .keepSourceNodesActive = true,
            .writeData = true,
            .printToTerminal = true,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, config);

    // Initial Data:
    Coord center(1.0 / 3.0, 0.0);
    double radius = 0.125;
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
            if (-0.25 < y && y < 0.25 && i <= 1)
            {
                radiation.isInitialGridPoint[ij] = true;
                radiation.initialE_LF[ij] = 1;
                radiation.initialFx_LF[ij] = 1;
                radiation.initialFy_LF[ij] = 0;
            }
            double dist = (center - xy).EuklNorm();
            if (dist <= radius)
                radiation.kappaA[ij] = 1e10;
        }
    radiation.RunSimulation();
}

void SphereWave(Stencil stencil, StreamingType streamingType, double cfl)
{
    // Create Radiation object:
    size_t nx = 300 + 1 + 2;
    size_t ny = 300 + 1 + 2;
    double dx = 2.0 / (nx - 1.0 - 2.0);
    double dy = 2.0 / (ny - 1.0 - 2.0);
    Coord start(-1 - dx, -1 - dy);
    Coord end(1 + dx, 1 + dy);
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
            .writePeriod = 1.0, // write first and last frame only.
            .updateFourierHarmonics = false,
            .keepSourceNodesActive = false,
            .writeData = true,
            .printToTerminal = true,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Intensities,
        };

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, config);

    // Initial Data:
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
            if (r < 0.1)
            {
                radiation.isInitialGridPoint[ij] = true;
                // Initial data given by moments:
                radiation.initialE_LF[ij] = 1;
                radiation.initialFx_LF[ij] = (r > 1e-6) ? xy[1] / r : 0;
                radiation.initialFy_LF[ij] = (r > 1e-6) ? xy[2] / r : 0;
                // Initial data given by intensities:
                for(int d=0; d<stencil.nDir; d++)
                    radiation.initialI[radiation.Index(ij, d)] = 1;
                radiation.initialFluxAngle_IF[ij] = 0;
            }
        }
    radiation.RunSimulation();
}

void SphereWaveShadow(Stencil stencil, StreamingType streamingType, double cfl)
{
    // Create Radiation object:
    size_t nx, ny;
    nx = ny = 300;
    Coord start(-0.5, -0.5);
    Coord end(2, 2);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil streamingStencil(5, 0, false);

    // Config:
    Config config =
        {
            .name = "Sphere Wave Shadow 2d/" + StreamingName(streamingType) + " " + stencil.name + Format(cfl, 2) + "cfl",
            .t0 = 0,
            .simTime = 1.55,
            .writePeriod = 2, // write first and last frame only.
            .updateFourierHarmonics = false,
            .keepSourceNodesActive = true,
            .writeData = true,
            .printToTerminal = true,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, config);

    // Initial Data:
    Coord center(0.75, 0.75);
    double sourceRadius = 0.25;
    double shadowRadius = 0.25;
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
            if (r < sourceRadius)
            {
                radiation.isInitialGridPoint[ij] = true;
                radiation.initialE_LF[ij] = 1;
                radiation.initialFx_LF[ij] = xy[1] / r;
                radiation.initialFy_LF[ij] = xy[2] / r;
            }
            double dist = (center - xy).EuklNorm();
            if (dist <= shadowRadius)
                radiation.kappaA[ij] = 1e10;
        }
    radiation.RunSimulation();
}

void CurvedBeam(Stencil stencil, StreamingType streamingType, double cfl, int nF = 5, bool updateFourierHarmonics = false)
{
    size_t nx = 250 + 1 + 2;
    size_t ny = 200 + 1 + 2;
    double dx = 5.0 / (nx - 1.0 - 2.0);
    double dy = 4.0 / (ny - 1.0 - 2.0);
    Coord start(0 - dx, 0 - dy);
    Coord end(5 + dx, 4 + dy);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(cfl);
    SchwarzSchild metric(grid, 1.0, 0.0); // needs at least stencil with order 5
    // KerrSchild metric(grid, 1.0, 0.0);   // initial direction is somehow wrong
    Stencil streamingStencil(nF, 0, false);
    
    // Config:
    string getsUpdated = (updateFourierHarmonics) ? " 1" : "";
    Config config =
        {
            .name = "Curved Beam 2d/" + StreamingName(streamingType) + " " + stencil.name + Format(cfl, 2) + "cfl " + std::to_string(nF) + "nF" + getsUpdated,
            .t0 = 0,
            .simTime = 10.1,
            .writePeriod = 10, // write first and last frame only.
            .updateFourierHarmonics = updateFourierHarmonics,
            .keepSourceNodesActive = true,
            .writeData = true,
            .printToTerminal = true,
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
                // High F, will be normalized to |F| = E anyway:
                radiation.initialFx_LF[ij] = 10 * vLF[1];
                radiation.initialFy_LF[ij] = 10 * vLF[2];
            }
        }
    radiation.RunSimulation();
}

void Diffusion(Stencil stencil, StreamingType streamingType, double kappaS, double lambda, double cfl)
{
    // Create Radiation object:
    size_t nx, ny;
    nx = ny = 200 + 1 + 2;
    double dx = 1.0 / (nx - 1.0 - 2.0);
    double dy = 1.0 / (ny - 1.0 - 2.0);
    Coord start(-0.5 - dx, -0.5 - dy);
    Coord end(0.5 + dx, 0.5 + dy);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil streamingStencil(5, 0, false);

    // Initial Data:
    // double lambda = 0.0;  // = 3kappa1 / kappa0
    double kappa0 = kappaS / (1.0 - lambda / 9.0);
    double kappa1 = kappa0 * lambda / 3.0;
    double correctionFactor = 0.5;
    double PE = kappa0 * grid.dx;
    double D = 1.0 / (2.0 * kappa0) * (1.0 + correctionFactor * PE);

    // Config:
    double t0 = 1;
    std::string name = "Diffusion 2d/" + StreamingName(streamingType) + " " + stencil.name + " " + Format(cfl, 2, true) + "cfl " + std::to_string(nx) + "nx " + std::to_string(ny) + "ny " + std::to_string((int)kappa0) + "kappa0 " + std::to_string((int)kappa1) + "kappa1 " + Format(PE, 1, true) + "PE";
    Config config =
        {
            .name = name,
            .t0 = t0,
            .simTime = 0.5,
            .writePeriod = 0.25,
            .updateFourierHarmonics = false,
            .keepSourceNodesActive = false,
            .writeData = true,
            .printToTerminal = true,
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
            // radiation.initialFx_LF[ij] = (x * E) / (2.0 * t0);
            // radiation.initialFy_LF[ij] = (y * E) / (2.0 * t0);
            radiation.initialFx_LF[ij] = (x * E) / (2.0 * t0 * (1.0 + correctionFactor * PE));
            radiation.initialFy_LF[ij] = (y * E) / (2.0 * t0 * (1.0 + correctionFactor * PE));
        }

    radiation.RunSimulation();
}

void MovingDiffusion(Stencil stencil, StreamingType streamingType, double kappaS, double lambda, double cfl)
{
    // Create Radiation object:
    size_t nx = 601;
    size_t ny = 51;
    Coord start(-3.0, -0.5);
    Coord end(3.0, 0.5);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil streamingStencil(5, 0, false);

    // Initial Data:
    // double lambda = 0.0;  // = 3kappa1 / kappa0
    double kappa0 = kappaS / (1.0 - lambda / 9.0);
    double kappa1 = kappa0 * lambda / 3.0;
    double PE = kappa0 * grid.dx;

    // Config:
    double t0 = 0;
    std::string name = "MovingDiffusion 2d/" + StreamingName(streamingType) + " " + stencil.name + " " + Format(cfl, 2, true) + "cfl " + std::to_string(nx) + "nx " + std::to_string(ny) + "ny " + std::to_string((int)kappa0) + "kappa0 " + std::to_string((int)kappa1) + "kappa1 " + Format(PE, 1, true) + "PE";
    Config config =
        {
            .name = name,
            .t0 = t0,
            .simTime = 2.0,
            .writePeriod = 1.0,
            .updateFourierHarmonics = false,
            .keepSourceNodesActive = false,
            .writeData = true,
            .printToTerminal = true,
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
            radiation.ux[ij] = 0.5;
            radiation.uy[ij] = 0.0;

            double E = exp(-9 * x * x);
            radiation.initialE_LF[ij] = E;
            radiation.initialFx_LF[ij] = 0.655 * E;
            radiation.initialFy_LF[ij] = 0;
        }

    radiation.RunSimulation();
}

void TestBeam(Stencil stencil, StreamingType streamingType, double cfl, Tensor2 F)
{
    // Create Radiation object:
    size_t nx = 200 + 1 + 2;
    size_t ny = 200 + 1 + 2;
    double dx = 1.0 / (nx - 1.0 - 2.0);
    double dy = 1.0 / (ny - 1.0 - 2.0);
    Coord start(-1 - dx, -1 - dy);
    Coord end(1 + dx, 1 + dy);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil streamingStencil(5, 0, false);

    // Config:
    Config config =
        {
            .name = "Test Beam 2d/" + StreamingName(streamingType) + " " + stencil.name + "(" + Format(F[1], 1) + "," + Format(F[2], 1) + ")F" + Format(cfl, 2) + "cfl",
            .t0 = 0,
            .simTime = 2.8,
            .writePeriod = 3, // write first and last frame only.
            .updateFourierHarmonics = false,
            .keepSourceNodesActive = true,
            .writeData = true,
            .printToTerminal = true,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, config);

    // Initial Data:
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
            if (r < 0.2)
            {
                radiation.isInitialGridPoint[ij] = true;
                radiation.initialE_LF[ij] = 1;
                radiation.initialFx_LF[ij] = F[1];
                radiation.initialFy_LF[ij] = F[2];
            }
        }
    radiation.RunSimulation();
}


void CflTestCurvedBeam(int n)
{
    if (n == 0)
        CurvedBeam(Stencil(200, 0), StreamingType::CurvedFixed, 1.00, 5, false);
    if (n == 1)
        CurvedBeam(Stencil(200, 0), StreamingType::CurvedFixed, 0.75, 5, false);
    if (n == 2)
        CurvedBeam(Stencil(200, 0), StreamingType::CurvedFixed, 0.50, 5, false);
    if (n == 3)
        CurvedBeam(Stencil(200, 0), StreamingType::CurvedAdaptive, 1.00, 5, false);
    if (n == 4)
        CurvedBeam(Stencil(200, 0), StreamingType::CurvedAdaptive, 0.75, 5, false);
    if (n == 5)
        CurvedBeam(Stencil(200, 0), StreamingType::CurvedAdaptive, 0.50, 5, false);
}
void CflTestDiffusion(int n)
{
    double lambda = 0;
    double kappaS = 10000;
    if (n == 0)
        MovingDiffusion(Stencil(20, 0), StreamingType::FlatFixed, kappaS, lambda, 1.00);
    if (n == 1)
        MovingDiffusion(Stencil(20, 0), StreamingType::FlatFixed, kappaS, lambda, 0.75);
    if (n == 2)
        MovingDiffusion(Stencil(20, 0), StreamingType::FlatFixed, kappaS, lambda, 0.50);
    if (n == 3)
        MovingDiffusion(Stencil(20, 0), StreamingType::FlatAdaptive, kappaS, lambda, 1.00);
    if (n == 4)
        MovingDiffusion(Stencil(20, 0), StreamingType::FlatAdaptive, kappaS, lambda, 0.75);
    if (n == 5)
        MovingDiffusion(Stencil(20, 0), StreamingType::FlatAdaptive, kappaS, lambda, 0.50);
}



// Paper Analysises:
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
void SphereWaveAnalysis(int n)
{
    if(n==0) SphereWave(Stencil(20,0), StreamingType::FlatFixed   , 0.75);
    if(n==1) SphereWave(Stencil(20,0), StreamingType::FlatAdaptive, 0.75);
    if(n==2) SphereWave(Stencil(16,4), StreamingType::FlatAdaptive, 0.75);
    
    if(n==0) SphereWave(Stencil(30,0), StreamingType::FlatFixed   , 0.75);
    if(n==1) SphereWave(Stencil(30,0), StreamingType::FlatAdaptive, 0.75);
    if(n==2) SphereWave(Stencil(24,6), StreamingType::FlatAdaptive, 0.75);
    
    if(n==0) SphereWave(Stencil(40,0), StreamingType::FlatFixed   , 0.75);
    if(n==1) SphereWave(Stencil(40,0), StreamingType::FlatAdaptive, 0.75);
    if(n==2) SphereWave(Stencil(32,8), StreamingType::FlatAdaptive, 0.75);
    
    if(n==0) SphereWave(Stencil(50,0), StreamingType::FlatFixed   , 0.75);
    if(n==1) SphereWave(Stencil(50,0), StreamingType::FlatAdaptive, 0.75);
    if(n==2) SphereWave(Stencil(40,10), StreamingType::FlatAdaptive, 0.75);
}
void BeamCrossingAnalysis(int n)
{
    if(n == 0) BeamCrossing(Stencil(200 ,0), StreamingType::FlatFixed   , 0.95);
    if(n == 1) BeamCrossing(Stencil(140, 0), StreamingType::FlatAdaptive, 0.95);
    if(n == 2) BeamCrossing(Stencil(100,40), StreamingType::FlatAdaptive, 0.95);
}
void DiffusionAnalysis(int n)
{
    // Diffusion:
    double lambda = 0;
    double cfl = 0.1;   // low cfl needed to reproduce correct flux density. With higher cfl (e.g 0.25) only energy density is reproduced correctly.
    if (n == 0) Diffusion(Stencil(20, 0), StreamingType::FlatFixed, 100.0, lambda, cfl);
    if (n == 1) Diffusion(Stencil(20, 0), StreamingType::FlatFixed, 1000.0, lambda, cfl);
    if (n == 2) Diffusion(Stencil(20, 0), StreamingType::FlatFixed, 10000.0, lambda, cfl);
    if (n == 3) Diffusion(Stencil(20, 0), StreamingType::FlatFixed, 100000.0, lambda, cfl);
    if (n == 4) Diffusion(Stencil(20, 0), StreamingType::FlatFixed, 1000000.0, lambda, cfl);
    
    // if (n == 0) MovingDiffusion(Stencil(20, 0), StreamingType::FlatFixed, 10000.0, lambda, cfl);
}


int main(int argc, char *argv[])
{
    int n = 1;
    if (argc > 1)
        n = atoi(argv[1]);

    // Straight Beam:
    // if (n == 0)
    // StraightBeam(Stencil(200, 0), StreamingType::FlatFixed, 0.75, Tensor2(1, 0));
    // if (n == 1)
    // StraightBeam(Stencil(110, 0), StreamingType::FlatAdaptive, 0.75, Tensor2(1, 0));
    // if (n == 2)
    // StraightBeam(Stencil(100, 10), StreamingType::FlatAdaptive, 0.75, Tensor2(1, 0));

    // Straight Beam Shadow:
    // if (n == 0)
    //     StraightBeamShadow(Stencil(200, 0), StreamingType::FlatFixed, 0.75);
    // if (n == 1)
    //     StraightBeamShadow(Stencil(100, 10), StreamingType::FlatAdaptive, 0.75);
    // if (n == 2)
    //     StraightBeamShadow(Stencil(110, 0), StreamingType::FlatAdaptive, 0.75);

    // Sphere Wave
    // if(n==1) SphereWave(Stencil(60,0), StreamingType::FlatFixed   , 0, 0.75);
    // if(n==2) SphereWave(Stencil(20,5), StreamingType::FlatAdaptive, 0, 0.75);
    // if(n==3) SphereWave(Stencil(25,0), StreamingType::FlatAdaptive, 0, 0.75);

    // Sphere Wave Shadow
    // if(n==2) SphereWaveShadow(Stencil(200, 0), StreamingType::FlatFixed   ,  10, 0.75);
    // if(n==2) SphereWaveShadow(Stencil(110, 0), StreamingType::FlatAdaptive,  10, 0.75);
    // if(n==2) SphereWaveShadow(Stencil(100,10), StreamingType::FlatAdaptive,  10, 0.75);

    // Curved Beam:
    // if (n == 3)
    //     CurvedBeam(Stencil(200, 0), StreamingType::CurvedFixed, 0.75);
    // if (n == 4)
    //     CurvedBeam(Stencil(110, 0), StreamingType::CurvedAdaptive, 0.75);
    // if (n == 5)
    //     CurvedBeam(Stencil(100, 10), StreamingType::CurvedAdaptive, 0.75);
    // CurvedBeam(Stencil(200, 0), StreamingType::CurvedAdaptive, 0.75);

    // Diffusion:
    // double lambda = 0;
    // double cfl = 0.25;
    // if (n == 0)
        // Diffusion(Stencil(20, 0), StreamingType::FlatFixed, 100.0, lambda, 0.25);
        // MovingDiffusion(Stencil(20, 0), StreamingType::FlatFixed, 10000.0, lambda, 0.25);
    // if (n == 0)
    // Diffusion(Stencil(20, 0), StreamingType::FlatFixed, 100.0, lambda, cfl);
    // if (n == 1)
    // Diffusion(Stencil(20, 0), StreamingType::FlatFixed, 500.0, lambda, cfl);
    // if (n == 2)
    // Diffusion(Stencil(20, 0), StreamingType::FlatFixed, 1000.0, lambda, cfl);
    // if (n == 3)
    // Diffusion(Stencil(20, 0), StreamingType::FlatFixed, 10000.0, lambda, cfl);
    // if (n == 4)
    // Diffusion(Stencil(20, 0), StreamingType::FlatFixed, 100000.0, lambda, cfl);
    // if (n == 5)
    // Diffusion(Stencil(20, 0), StreamingType::FlatFixed, 1000000.0, lambda, cfl);

    // if (n == 0)
    // Diffusion1D(Stencil(20, 0), StreamingType::FlatFixed, 100.0, lambda, cfl);

    // Testing:
    // CflTestCurvedBeam(n);
    // CflTestDiffusion(n);
    // if (n == 0)
    //    TestBeam(Stencil(200), StreamingType::FlatFixed, 0.75, Tensor2(1, 1));
    // if (n == 1)
    //    TestBeam(Stencil(110), StreamingType::FlatAdaptive, 0.75, Tensor2(1, 1));
    // if (n == 2)
    //    TestBeam(Stencil(100, 10), StreamingType::FlatAdaptive, 0.75, Tensor2(1, 1));


    // Paper:
    // StreamingTypePerformanceAnalysis(n);
    // SphereWaveAnalysis(n); // not used in paper, I use 3D version instead.
    // BeamCrossingAnalysis(n);
    DiffusionAnalysis(n);
}