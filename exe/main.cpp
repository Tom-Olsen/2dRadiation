#include <iostream>
#include "../src/Radiation.h"
using namespace std;

void StraightBeam(Stencil stencil, StreamingType streamingType, double cfl, Tensor2 F)
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
    Stencil streamingStencil(5);

    // Config:
    Config config =
        {
            .name = "Straight Beam 2d/" + StreamingName(streamingType) + " " + stencil.name + "(" + Format(F[1], 1) + "," + Format(F[2], 1) + ")F" + Format(cfl, 2) + "cfl",
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
            radiation.initialKappa0[ij] = 0;
            radiation.initialKappa1[ij] = 0;
            radiation.initialKappaA[ij] = 0;
            radiation.initialEta[ij] = 0;
            if (-0.25 < y && y < 0.25 && i <= 1)
            {
                radiation.isInitialGridPoint[ij] = true;
                radiation.initialE_LF[ij] = 1;
                radiation.initialFx_LF[ij] = F[1];
                radiation.initialFy_LF[ij] = F[2];
            }
        }
    radiation.RunSimulation();
}

void BeamCrossing(Stencil stencil, StreamingType streamingType, double cfl)
{
    // Create Radiation object:
    size_t nx = 300 + 1 + 2;
    size_t ny = 300 + 1 + 2;
    double dx = 3.0 / (nx - 1.0 - 2.0);
    double dy = 3.0 / (ny - 1.0 - 2.0);
    Coord start(0 - dx, 0 - dy);
    Coord end(3 + dx, 3 + dy);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil streamingStencil(5);

    // Config:
    Config config =
        {
            .name = "Beam Crossing 2d/" + StreamingName(streamingType) + "_" + stencil.name + Format(cfl, 2) + "cfl",
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
            double x = xy[1];
            double y = xy[2];
            radiation.initialKappa0[ij] = 0;
            radiation.initialKappa1[ij] = 0;
            radiation.initialKappaA[ij] = 0;
            radiation.initialEta[ij] = 0;
            if (1 < y && y < 2 && i <= 1)
            {
                radiation.isInitialGridPoint[ij] = true;
                radiation.initialE_LF[ij] = 1;
                radiation.initialFx_LF[ij] = 1;
                radiation.initialFy_LF[ij] = 0;
            }
            if (1 < x && x < 2 && j <= 1)
            {
                size_t ij = grid.Index(i, j);
                radiation.isInitialGridPoint[ij] = true;
                radiation.initialE_LF[ij] = 1;
                radiation.initialFx_LF[ij] = 0;
                radiation.initialFy_LF[ij] = 1;
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
    Stencil streamingStencil(5);

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
            radiation.initialKappa0[ij] = 0;
            radiation.initialKappa1[ij] = 0;
            radiation.initialKappaA[ij] = 0;
            radiation.initialEta[ij] = 0;
            if (-0.25 < y && y < 0.25 && i <= 1)
            {
                radiation.isInitialGridPoint[ij] = true;
                radiation.initialE_LF[ij] = 1;
                radiation.initialFx_LF[ij] = 1;
                radiation.initialFy_LF[ij] = 0;
            }
            double dist = (center - xy).EuklNorm();
            if (dist <= radius)
                radiation.initialKappaA[ij] = 1e10;
        }
    radiation.RunSimulation();
}

void SphereWave(Stencil stencil, StreamingType streamingType, double F, double cfl)
{
    // Create Radiation object:
    size_t nx, ny;
    nx = ny = 300;
    Coord start(-1, -1);
    Coord end(1, 1);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil streamingStencil(5);

    // Config:
    Config config =
        {
            .name = "Sphere Wave 2d/" + StreamingName(streamingType) + " " + stencil.name + Format(F, 1) + "F" + Format(cfl, 2) + "cfl",
            .t0 = 0,
            .simTime = 0.7,
            .writePeriod = 1, // write first and last frame only.
            .updateFourierHarmonics = false,
            .keepSourceNodesActive = false,
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
            radiation.initialKappa0[ij] = 0;
            radiation.initialKappa1[ij] = 0;
            radiation.initialKappaA[ij] = 0;
            radiation.initialEta[ij] = 0;
            if (r < 0.1)
            {
                radiation.isInitialGridPoint[ij] = true;
                radiation.initialE_LF[ij] = 1;
                radiation.initialFx_LF[ij] = xy[1] * F;
                radiation.initialFy_LF[ij] = xy[2] * F;
            }
        }
    radiation.RunSimulation();
}

void SphereWaveShadow(Stencil stencil, StreamingType streamingType, double F, double cfl)
{
    // Create Radiation object:
    size_t nx, ny;
    nx = ny = 300;
    Coord start(-0.5, -0.5);
    Coord end(2, 2);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil streamingStencil(5);

    // Config:
    Config config =
        {
            .name = "Sphere Wave Shadow 2d/" + StreamingName(streamingType) + " " + stencil.name + Format(F, 1) + "F" + Format(cfl, 2) + "cfl",
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
            radiation.initialKappa0[ij] = 0;
            radiation.initialKappa1[ij] = 0;
            radiation.initialKappaA[ij] = 0;
            radiation.initialEta[ij] = 0;
            if (r < sourceRadius)
            {
                radiation.isInitialGridPoint[ij] = true;
                radiation.initialE_LF[ij] = 1;
                radiation.initialFx_LF[ij] = F * xy[1];
                radiation.initialFy_LF[ij] = F * xy[2];
            }
            double dist = (center - xy).EuklNorm();
            if (dist <= shadowRadius)
                radiation.initialKappaA[ij] = 1e10;
        }
    radiation.RunSimulation();
}

void CurvedBeam(Stencil stencil, StreamingType streamingType, double cfl)
{
    size_t nx = 500 + 1 + 2;
    size_t ny = 400 + 1 + 2;
    double dx = 5.0 / (nx - 1.0 - 2.0);
    double dy = 4.0 / (ny - 1.0 - 2.0);
    Coord start(0 - dx, 0 - dy);
    Coord end(5 + dx, 4 + dy);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(cfl);
    SchwarzSchild metric(grid, 1.0, 0.0); // needs at least stencil with order 5
    // KerrSchild metric(grid, 1.0, 0.0);   // initial direction is somehow wrong
    Stencil streamingStencil(5);

    // Config:
    Config config =
        {
            .name = "Curved Beam 2d/" + StreamingName(streamingType) + " " + stencil.name + Format(cfl, 2) + "cfl",
            .t0 = 0,
            .simTime = 10,
            .writePeriod = 11, // write first and last frame only.
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
    PARALLEL_FOR(2)
    for (size_t j = 0; j < grid.ny; j++)
        for (size_t i = 0; i < grid.nx; i++)
        {
            size_t ij = grid.Index(i, j);
            Coord xy = grid.xy(i, j);
            double y = xy[2];
            radiation.initialKappa0[ij] = 0;
            radiation.initialKappa1[ij] = 0;
            radiation.initialKappaA[ij] = 0;
            radiation.initialEta[ij] = 0;
            if (3.0 < y && y < 3.5 && i <= 1)
            {
                Tensor3 uLF(1, 1, 0);
                uLF = NullNormalize(uLF, metric.GetMetric_ll(ij));
                Tensor2 vLF = Vec2ObservedByEulObs<LF, LF>(uLF, xy, metric);
                radiation.isInitialGridPoint[ij] = true;
                radiation.initialE_LF[ij] = 1;
                radiation.initialFx_LF[ij] = 1 * vLF[1];
                radiation.initialFy_LF[ij] = 1 * vLF[2];

                // Tensor3 uLF(1,1,0);
                // uLF = NullNormalize(uLF,metric.GetMetric_ll(ij));
                // Tensor2 vLF = Vec2ObservedByEulObs<LF,LF>(uLF, xy, metric);
                // double E = 2;
                // Tensor3 vec(E, vLF[1], vLF[2]);
                // vec = NullNormalize(vec,metric.GetMetric_ll(ij));
                // radiation.isInitialGridPoint[ij] = true;
                // radiation.initialE_LF[ij]  = vec[0];
                // radiation.initialFx_LF[ij] = vec[1];
                // radiation.initialFy_LF[ij] = vec[2];
                //  vec.Print("vec");
                //  PrintDouble(Norm2(vec,metric.GetMetric_ll(ij)), "|vec|");
                //  vLF.Print("vLF");
                //  PrintDouble(Norm2(vLF,metric.GetGamma_ll(ij)), "|vLF|");
            }
        }
    radiation.RunSimulation();
}

void Diffusion(Stencil stencil, StreamingType streamingType, double cfl)
{
    // Create Radiation object:
    size_t nx, ny;
    nx = ny = 201;
    Coord start(-0.5, -0.5);
    Coord end(0.5, 0.5);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(cfl);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil streamingStencil(5);

    // Config:
    int d = 2;
    double t = 0;
    // double t = 1;
    Config config =
        {
            .name = "Diffusion 2d/" + StreamingName(streamingType) + " " + stencil.name + Format(cfl, 2) + "cfl " + std::to_string(nx) + "nx " + std::to_string(ny) + "ny",
            .t0 = t,
            .simTime = 3,
            .writePeriod = 1,
            .updateFourierHarmonics = false,
            .keepSourceNodesActive = false,
            .writeData = true,
            .printToTerminal = true,
            .streamingType = streamingType,
            .initialDataType = InitialDataType::Moments,
        };

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, config);

    // Initial Data Carlo:
    // double kappaS = 100; // = kappa0 - kappa1/3
    // double lambda = 0.0; // = 3kappa1 / kappa0
    // double kappa0 = kappaS / (1.0 - lambda / 9.0);
    // double kappa1 = kappa0 * lambda / 3.0;

    // Initial Data Lucas:
    double kappa0 = 100;
    double lambda = 0.0;
    double kappa1 = kappa0 * lambda / 3.0;
    double A = 1;
    double sigma0 = 0.1;
    double D = (1.0 - grid.dt * kappa0 / 2.0) / (d * kappa0);
    double sigmaD = sqrt(2 * D * t);

    for (size_t j = 0; j < grid.ny; j++)
        for (size_t i = 0; i < grid.nx; i++)
        {
            size_t ij = grid.Index(i, j);
            Coord xy = grid.xy(i, j);
            double x = xy[1];
            double y = xy[2];
            double r = xy.EuklNorm();
            radiation.initialKappa0[ij] = kappa0;
            radiation.initialKappa1[ij] = kappa1;
            radiation.initialKappaA[ij] = 0;
            radiation.initialEta[ij] = 0;

            radiation.isInitialGridPoint[ij] = true;
            // Carlo:
            // double E = pow(kappaS / t, d / 2.0) * exp(-3.0 * kappaS * r * r / (4.0 * t));
            // radiation.initialE_LF[ij] = E;
            // radiation.initialFx_LF[ij] = x / (2.0 * t) * E;
            // radiation.initialFy_LF[ij] = y / (2.0 * t) * E;
            // Lucas:
            double E = A * (sigma0 * sigma0) / (sigma0 * sigma0 + sigmaD * sigmaD) * exp(-r * r / (2 * sigma0 * sigma0 + 2 * sigmaD * sigmaD));
            radiation.initialE_LF[ij] = E;
            radiation.initialFx_LF[ij] = (x * E) / (3.0 * kappa0 * (sigma0 * sigma0 + sigmaD * sigmaD));
            radiation.initialFy_LF[ij] = (y * E) / (3.0 * kappa0 * (sigma0 * sigma0 + sigmaD * sigmaD));
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
    Stencil streamingStencil(5);

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
            radiation.initialKappa0[ij] = 0;
            radiation.initialKappa1[ij] = 0;
            radiation.initialKappaA[ij] = 0;
            radiation.initialEta[ij] = 0;
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
// Note:
// -TE changed from 1e-6 to 1e-4

// TODO:
// -fix kerr metric initial direction?
// -try higher cfl condition
// -test Halo=1 vs Halo=2
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

    // Beam Crossing:
    // if(n==0) BeamCrossing(Stencil(200 ,0), StreamingType::FlatFixed   , 2000, 0.75);
    // if(n==1) BeamCrossing(Stencil(110, 0), StreamingType::FlatAdaptive,  600, 0.75);
    // if(n==2) BeamCrossing(Stencil(100,10), StreamingType::FlatAdaptive,  600, 0.75);

    // Straight Beam Shadow:
    // if (n == 0)
    // StraightBeamShadow(Stencil(200, 0), StreamingType::FlatFixed, 0.75);
    // if (n == 1)
    // StraightBeamShadow(Stencil(100, 10), StreamingType::FlatAdaptive, 0.75);
    // if (n == 2)
    // StraightBeamShadow(Stencil(110, 0), StreamingType::FlatAdaptive, 0.75);

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

    // Diffusion:
    // if (n == 0)
    //     Diffusion(Stencil(200, 0), StreamingType::FlatFixed, 0.2);
    // if (n == 1)
    //     Diffusion(Stencil(110, 0), StreamingType::FlatAdaptive, 0.2);
    // if (n == 2)
    //     Diffusion(Stencil(100, 10), StreamingType::FlatAdaptive, 0.2);

    // Testing:
    if (n == 0)
        TestBeam(Stencil(200), StreamingType::FlatFixed, 0.75, Tensor2(1, 1));
    if (n == 1)
        TestBeam(Stencil(110), StreamingType::FlatAdaptive, 0.75, Tensor2(1, 1));
    if (n == 2)
        TestBeam(Stencil(100, 10), StreamingType::FlatAdaptive, 0.75, Tensor2(1, 1));
}