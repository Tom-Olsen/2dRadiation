#include <iostream>
#include "../src/Radiation.h"
using namespace std;



void SphereWave(Stencil stencil, StreamingType streamingType)
{
    // Create Radiation object:
    size_t nx, ny;
    nx = ny = 200;
    Coord start(-1,-1);
    Coord end(1,1);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(0.5);
    Minkowski metric(grid, 1.0, 0.0);
    Stencil streamingStencil(5);

    Radiation radiation(metric, stencil, streamingStencil, streamingType);
    radiation.sigma = 1.0;

    // Initial Data:
    for(size_t j=0; j<grid.ny; j++)
    for(size_t i=0; i<grid.nx; i++)
    {
        size_t ij = grid.Index(i,j);
        Coord xy = grid.xy(i,j);
        double r = xy.EuklNorm();
        if (r < 0.1)
        {
            radiation.isInitialGridPoint[ij] = true;
            radiation.initialE[ij] = 1;
        }
    }

    // Start simulation:
    Config config =
    {
        .name = "SphereWave_" + StreamingName(streamingType) + "_" + stencil.name,
        .simTime = 0.7,
        .writeFrequency = 5,
        .updateFourierHarmonics = false,
        .keepSourceNodesActive = false,
        .writeData = true,
        .printToTerminal = true,
    };
    radiation.RunSimulation(config);
}
void CurvedBeam(size_t nx, size_t ny, Stencil stencil, int sigma, int simTime, StreamingType streamingType, std::string comment = "")
{
    // Grid, Metric, Stencil:
    Coord start(-0.1, -1.1);
    Coord end  (6.0, 5.0);
    Grid grid(nx, ny, start, end);
    grid.SetCFL(0.5);
    SchwarzSchild metric(grid, 1.0, 0.0);   // needs at least stencil with order 5
    // KerrSchild metric(grid, 1.0, 0.0);   // initial direction is somehow wrong
    Stencil streamingStencil(5);

    std::cout << "ij max: " << grid.Index(grid.nx-1, grid.ny-1) << std::endl;

    // Radiation:
    Radiation radiation(metric, stencil, streamingStencil, streamingType);
    radiation.sigma = sigma;

    // Initial Data:
    PARALLEL_FOR(2)
    for(size_t j=0; j<grid.ny; j++)
    for(size_t i=0; i<grid.nx; i++)
    {
        size_t ij = grid.Index(i,j);
        Coord xy = grid.xy(i,j);
        double x = xy[1];
        double y = xy[2];
        if(3.0 < y && y < 3.5 && x < 0)
        {
            Tensor3 uLF(1,1,0);
            uLF = NullNormalize(uLF,metric.GetMetric_ll(ij));
            Tensor2 vIF = Vec2ObservedByEulObs<LF,IF>(uLF, xy, metric);
            
            radiation.isInitialGridPoint[ij] = true;
            radiation.initialE[ij] = 1;
            radiation.initialNx[ij] = vIF[1];
            radiation.initialNy[ij] = vIF[2];
            radiation.initialKappa0[ij] = 0;
            radiation.initialKappa1[ij] = 0;
            radiation.initialKappaA[ij] = 0;
            radiation.initialEta[ij] = 0;
        }
    }

    // Start simulation:
    Config config =
    {
        .name = "CurvedBeam_" + metric.Name() + "_" + stencil.name + "_"
              + std::to_string(nx) + "x" + std::to_string(ny) + "y_"
              + std::to_string(sigma) + "s_" + std::to_string(simTime) + "t_" + StreamingName(streamingType) + "_" + comment,
        .simTime = (double)simTime,
        .writeFrequency = 20,
        .updateFourierHarmonics = false,
        .keepSourceNodesActive = true,
        .writeData = true,
        .printToTerminal = true,
    };
    radiation.RunSimulation(config);
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
    if(argc > 1)
        n = atoi(argv[1]);

    // Runs for paper: Sphere Wave
    SphereWave(Stencil(10,0), StreamingType::FlatFixed);
    SphereWave(Stencil(20,0), StreamingType::FlatFixed);
    SphereWave(Stencil(30,0), StreamingType::FlatFixed);
    SphereWave(Stencil(40,0), StreamingType::FlatFixed);
    SphereWave(Stencil(50,0), StreamingType::FlatFixed);
    SphereWave(Stencil(60,0), StreamingType::FlatFixed);
    SphereWave(Stencil(10,0), StreamingType::FlatAdaptive);
    SphereWave(Stencil(10,5), StreamingType::FlatAdaptive);
    SphereWave(Stencil(20,0), StreamingType::FlatAdaptive);
    SphereWave(Stencil(20,5), StreamingType::FlatAdaptive);

    // Runs for paper: Curved Beam
    // CurvedBeam(200,200, Stencil( 50, 0), 100, 10, StreamingType::CurvedFixed);
    // CurvedBeam(200,200, Stencil(100, 0), 250, 10, StreamingType::CurvedFixed);
    // CurvedBeam(200,200, Stencil(200, 0), 700, 10, StreamingType::CurvedFixed);
    // CurvedBeam(200,200, Stencil( 50, 0), 100, 10, StreamingType::CurvedAdaptive);
    // CurvedBeam(200,200, Stencil( 50,10), 100, 10, StreamingType::CurvedAdaptive);
    // CurvedBeam(200,200, Stencil( 60, 0), 100, 10, StreamingType::CurvedAdaptive);
    // CurvedBeam(200,200, Stencil(100, 0), 250, 10, StreamingType::CurvedAdaptive);
    // CurvedBeam(200,200, Stencil(100,10), 250, 10, StreamingType::CurvedAdaptive);
    // CurvedBeam(200,200, Stencil(110, 0), 250, 10, StreamingType::CurvedAdaptive);
    // CurvedBeam(200,200, Stencil(200, 0), 700, 10, StreamingType::CurvedAdaptive);
    // CurvedBeam(200,200, Stencil(200,20), 700, 10, StreamingType::CurvedAdaptive);
    // CurvedBeam(200,200, Stencil(220, 0), 700, 10, StreamingType::CurvedAdaptive);
}