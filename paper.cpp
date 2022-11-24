#include <iostream>     // cin/cout
#include <math.h>       // basic maths
#include <fstream>      // file input/output
#include <stdlib.h>     // srand, rand
#include <time.h>       // time
#include <filesystem>   // folder/file management
#include <iomanip>      // std::put:time
#include "src/Includes.hh"

using namespace std;



void RingWaveFlatRot(double sigma_, int nDir_, int nGrid_, double cfl_)
{
    double simTime = 1;
    int nDir = nDir_;
    int nGrid = nGrid_;
    double sigma = sigma_;
    string fileName = "RingWaveFlatRot NN" + std::to_string(nGrid) +
        " N" + std::to_string(nDir) +
        " cfl"  + std::to_string(cfl_);

    Int2 N(nGrid,nGrid);
    Coordinate2<xy> Start(-1.1, -1.1);
    Coordinate2<xy> End(1.1, 1.1);
    Grid2D<xy> grid(N,Start,End);
    grid.SetCFL(cfl_);
    Minkowski<xy> metric(grid,1.0,0.0);
    RotStencil stencil(nDir_,sigma_);
    UniformStencil fourierStencil(5);
    
    // Initialize simulation:
    RadiationNew radiation(grid,metric,stencil,fourierStencil);

    double r0 = 0.1;
    double rSigma = 0.04;
    for(int j=0; j<grid.n2; j++)
    {
        for(int i=0; i<grid.n1; i++)
        {
            double r = grid.rCoord(i,j);
            double ph = grid.phCoord(i,j);
            int ij = grid.Index(i,j);
            if(r < 4*r0)
            {
                radiation.isInitialGridPoint[ij] = true;
                radiation.initialE[ij] = exp(-0.5 * pow((r - r0)/rSigma,2));
                radiation.initialRotation[ij] = ph;
                radiation.initialKappa0[ij] = 0;
                radiation.initialKappa1[ij] = 0;
                radiation.initialKappaA[ij] = 0;
                radiation.initialEta[ij] = 0;
            }
        }
    }
    
    // Run simulation:
    int writeFrequency = 2;
    bool updateFourierHarmonics = false;
    bool keepSourceNodesActive = false;
    bool writeData = true;
    // bool writeData = false;
    bool printToTerminal = true;
    radiation.RunSimulation(fileName,simTime,writeFrequency,updateFourierHarmonics,keepSourceNodesActive,writeData,printToTerminal);
}
void CurvedBeamCloseRot(double sigma_, int nDir_, int nGrid_, int nFourier_, double cfl_)
{
    double simTime = 10;
    int nDir = nDir_;
    int nFourier = nFourier_;
    int nGrid = nGrid_;
    double sigma = sigma_;
    string fileName = "CurvedBeamCloseRot NN" + std::to_string(nGrid) +
        " N" + std::to_string(nDir) + 
        " F" + std::to_string(nFourier) + 
        " s" + std::to_string(sigma) +
        " cfl" + std::to_string(cfl_);

    Int2 N(nGrid,nGrid);
    Coordinate2<xy> Start(-0.1, -1.0);
    Coordinate2<xy> End(6.0, 5.0);
    Grid2D<xy> grid(N,Start,End);
    grid.SetCFL(cfl_);
    // KerrSchild<xy> metric(grid,1.0,0.5);
    SchwarzSchild<xy> metric(grid,1.0,0.0);
    // Minkowski<xy> metric(grid,1.0,0.0);
    RotStencil stencil(nDir_,sigma_);
    UniformStencil fourierStencil(5);
    
    // Initialize simulation:
    RadiationNew radiation(grid,metric,stencil,fourierStencil);

    for(int j=0; j<grid.n2; j++)
    {
        for(int i=0; i<grid.n1; i++)
        {
            if(3.0 < grid.yCoord(i,j) && grid.yCoord(i,j) < 3.5 && grid.xCoord(i,j) < 0.0)
            {
                int ij = grid.Index(i,j);
                Coordinate2<xy> x = grid.xyCoord(i,j);
                Tensor3<xy,LF> u(1,1,0);
                u.NullNormalize(metric.GetMetric_ll(x));
                Tensor2<xy,IF> v = Vec2ObservedByEulObs<xy,LF,IF>(u,x,metric);

                radiation.isInitialGridPoint[ij] = true;
                radiation.initialE[ij] = 1;
                radiation.initialRotation[ij] = v.Angle();
                radiation.initialKappa0[ij] = 0;
                radiation.initialKappa1[ij] = 0;
                radiation.initialKappaA[ij] = 0;
                radiation.initialEta[ij] = 0;
            }
        }
    }

    // Run simulation:
    int writeFrequency = 5;
    bool updateFourierHarmonics = false;
    bool keepSourceNodesActive = true;
    bool writeData = true;
    bool printToTerminal = true;
    radiation.RunSimulation(fileName,simTime,writeFrequency,updateFourierHarmonics,keepSourceNodesActive,writeData,printToTerminal);
}



void RingWaveFlatDirected(double sigma_, int nDir_, int nGrid_, double cfl_)
{
    double simTime = 1;
    int nDir = nDir_;
    int nMom = 4 * nDir;
    int nGrid = nGrid_;
    // X2 distribution;
    // PieceWiseLinear0 distribution;
    PieceWiseLinear1 distribution;
    double sigma = sigma_;
    string fileName = "RingWaveFlatDirected cubic x^2 NN" + std::to_string(nGrid) +
        " N" + std::to_string(nDir) +
        " Mom" + std::to_string(nMom) +
        " cfl"  + std::to_string(cfl_);

    Int2 N(nGrid,nGrid);
    Coordinate2<xy> Start(-1.1, -1.1);
    Coordinate2<xy> End(1.1, 1.1);
    Grid2D<xy> grid(N,Start,End);
    grid.SetCFL(cfl_);
    Minkowski<xy> metric(grid,1.0,0.0);
    SimulationData simData(simTime, nDir, nMom, nDir, sigma, distribution, metric, fileName);

    double r0 = 0.1;
    double rSigma = 0.04;
    for(int j=0; j<grid.n2; j++)
    {
        for(int i=0; i<grid.n1; i++)
        {
            double r = grid.rCoord(i,j);
            double ph = grid.phCoord(i,j);
            int ij = grid.Index(i,j);
            simData.initialE[ij] = exp(-0.5 * pow((r - r0)/rSigma,2));
            simData.initialAngle[ij] = ph;
        }
    }

    // Initialize simulation:
    RadiationFlat radiation(simData);
    
    // Run simulation:
    int writeFrequency = 2;
    bool updateFourierHarmonics = false;
    bool keepSourceNodesActive = false;
    // bool writeData = true;
    bool writeData = false;
    bool printToTerminal = true;
    radiation.RunSimulation(writeFrequency,updateFourierHarmonics,keepSourceNodesActive,writeData,printToTerminal);
}



void RingWaveFlatUniform(int nDir_, int nGrid_, double cfl_)
{
    double simTime = 1;
    int nDir = nDir_;
    int nGrid = nGrid_;
    X2 distribution;
    string fileName = "RingWaveFlatUniform NN" + std::to_string(nGrid) +
        " N" + std::to_string(nDir) +
        " cfl"  + std::to_string(cfl_);

    Int2 N(nGrid,nGrid);
    Coordinate2<xy> Start(-1.1, -1.1);
    Coordinate2<xy> End(1.1, 1.1);
    Grid2D<xy> grid(N,Start,End);
    grid.SetCFL(cfl_);
    Minkowski<xy> metric(grid,1.0,0.0);
    SimulationData simData(simTime, nDir, nDir, nDir, 0, distribution, metric, fileName);
    
    for(int j=0; j<grid.n2; j++)
    {
        for(int i=0; i<grid.n1; i++)
        {
            double r = grid.rCoord(i,j);
            double phi = grid.phCoord(i,j);
            if(r < 5.0 * grid.d1)
            {
                int ij = grid.Index(i,j);
                simData.initialE[ij] = 1;
                simData.initialAngle[ij] = phi;
            }
        }
    }

    // Initialize simulation:
    RadiationFlatUniform radiation(simData);
    // Run simulation:
    int writeFrequency = 2;
    bool updateFourierHarmonics = false;
    bool keepSourceNodesActive = false;
    // bool writeData = true;
    bool writeData = false;
    bool printToTerminal = true;
    radiation.RunSimulation(writeFrequency,updateFourierHarmonics,keepSourceNodesActive,writeData,printToTerminal);
}



void BeamCrossing(double sigma_, int nDir_, int nGrid_, double cfl_)
{
    double simTime = 2.5;
    int nDir = nDir_;
    int nMom = 4 * nDir;
    int nFourier = 5;
    int nGrid = nGrid_;
    X2 distribution;
    double sigma = sigma_;
    string fileName = "BeamCrossing cubic x^2 NN" + std::to_string(nGrid) +
        " N" + std::to_string(nDir) +
        " s" + std::to_string(sigma) +
        " Mom" + std::to_string(nMom) +
        " cfl" + std::to_string(cfl_);

    Int2 N(nGrid,nGrid);
    Coordinate2<xy> Start(-1, -1);
    Coordinate2<xy> End(1, 1);
    Grid2D<xy> grid(N,Start,End);
    grid.SetCFL(cfl_);
    Minkowski<xy> metric(grid,1.0,0.0);

    SimulationData simData(simTime, nDir, nMom, nFourier, sigma, distribution, metric, fileName);

    //for(int j=0; j<grid.n2; j++)
    //{
    //    for(int i=0; i<grid.n1; i++)
    //    {
    //        if(-0.2 < grid.xCoord(i,j) && grid.xCoord(i,j) < 0.2 && grid.yCoord(i,j) < -0.8)
    //        {// x€[-0.2,0.2] y<-0.8
    //            int ij = grid.Index(i,j);
    //            simData.initialE[ij] = 1;
    //            simData.initialAngle[ij] = M_PI / 2.0;
    //        }
    //        if(-0.2 < grid.yCoord(i,j) && grid.yCoord(i,j) < 0.2 && grid.xCoord(i,j) < -0.8)
    //        {// y€[-0.2,0.2] x<-0.8
    //            int ij = grid.Index(i,j);
    //            simData.initialE[ij] = 1;
    //            simData.initialAngle[ij] = 2.0 * M_PI;
    //        }
    //    }
    //}
    
    for(int j=0; j<grid.n2; j++)
    {
        for(int i=0; i<grid.n1; i++)
        {
            if(-0.5 < grid.yCoord(i,j) && grid.yCoord(i,j) < -0.3 && grid.xCoord(i,j) < -0.8)
            {// x€[-0.2,0.2] y<-0.8
                int ij = grid.Index(i,j);
                simData.initialE[ij] = 1;
                simData.initialAngle[ij] = M_PI / 8.0;
            }
            if(0.3 < grid.yCoord(i,j) && grid.yCoord(i,j) < 0.5 && grid.xCoord(i,j) < -0.8)
            {// y€[-0.2,0.2] x<-0.8
                int ij = grid.Index(i,j);
                simData.initialE[ij] = 1;
                simData.initialAngle[ij] = 2.0 * M_PI - M_PI / 8.0;
            }
        }
    }

    // Initialize simulation:
    Radiation radiation(simData);
    // Run simulation:
    int writeFrequency = 1;
    bool updateFourierHarmonics = false;
    bool keepSourceNodesActive = true;
    bool writeData = true;
    bool printToTerminal = true;
    radiation.RunSimulation(writeFrequency,updateFourierHarmonics,keepSourceNodesActive,writeData,printToTerminal);
}



void CurvedBeamClose(double sigma_, int nDir_, int nGrid_, double cfl_)
{
    double simTime = 10;
    int nDir = nDir_;
    int nMom = 4 * nDir;
    int nFourier = 5;
    int nGrid = nGrid_;
    X2 distribution;
    double sigma = sigma_;
    string fileName = "CurvedBeamClose cubic x^2 NN" + std::to_string(nGrid) +
        " N" + std::to_string(nDir) + 
        " s" + std::to_string(sigma) +
        " Mom" + std::to_string(nMom) +
        " cfl" + std::to_string(cfl_);

    Int2 N(nGrid,nGrid);
    Coordinate2<xy> Start(-0.1, -1.0);
    Coordinate2<xy> End(6.0, 5.0);
    Grid2D<xy> grid(N,Start,End);
    grid.SetCFL(cfl_);
    // KerrSchild<xy> metric(grid,1.0,0.5);
    SchwarzSchild<xy> metric(grid,1.0,0.0);
    // Minkowski<xy> metric(grid,1.0,0.0);
    SimulationData simData(simTime, nDir, nMom, nFourier, sigma, distribution, metric, fileName);

    for(int j=0; j<grid.n2; j++)
    {
        for(int i=0; i<grid.n1; i++)
        {
            if(3.0 < grid.yCoord(i,j) && grid.yCoord(i,j) < 3.5 && grid.xCoord(i,j) < 0.0)
            {
                Coordinate2<xy> x = grid.xyCoord(i,j);
                Tensor3<xy,LF> u(1,1,0);
                u.NullNormalize(metric.GetMetric_ll(x));
                Tensor2<xy,IF> v = Vec2ObservedByEulObs<xy,LF,IF>(u,x,simData.metric);

                int ij = grid.Index(i,j);
                simData.initialE[ij] = 1;
                simData.initialAngle[ij] = v.Angle();
            }
        }
    }
    // Initialize simulation:
    Radiation radiation(simData);

    // Run simulation:
    int writeFrequency = 5;
    bool updateFourierHarmonics = false;
    bool keepSourceNodesActive = true;
    bool writeData = true;
    bool printToTerminal = true;
    radiation.RunSimulation(writeFrequency,updateFourierHarmonics,keepSourceNodesActive,writeData,printToTerminal);
}



void CurvedBeamFar(double sigma_, int nDir_, int nGrid_, double cfl_)
{
    double simTime = 15;
    int nDir = nDir_;
    int nMom = 4 * nDir;
    int nFourier = 5;
    int nGrid = nGrid_;
    X2 distribution;
    double sigma = sigma_;
    string fileName = "CurvedBeamFar cubic x^2 NN" + std::to_string(nGrid) +
        " N" + std::to_string(nDir) + 
        " s" + std::to_string(sigma) +
        " Mom" + std::to_string(nMom) +
        " cfl" + std::to_string(cfl_);

	Int2 N(nGrid,nGrid);
    Coordinate2<xy> Start(-0.1, 3.9);
    Coordinate2<xy> End(8.1, 9.1);
    Grid2D<xy> grid(N,Start,End);
    grid.SetCFL(cfl_);
    // KerrSchild<xy> metric(grid,1.0,0.0);
    SchwarzSchild<xy> metric(grid,1.0,0.0);
    // Minkowski<xy> metric(grid,1.0,0.0);
	SimulationData simData(simTime, nDir, nMom, nFourier, sigma, distribution, metric, fileName);
	
	for(int j=0; j<grid.n2; j++)
	{
	    for(int i=0; i<grid.n1; i++)
        {
		    if(7.0 < grid.yCoord(i,j) && grid.yCoord(i,j) < 8.0 && grid.xCoord(i,j) < 0.0)
		    {
		    	Coordinate2<xy> x = grid.xyCoord(i,j);
		    	Tensor3<xy,LF> u(1,1,0);
		    	u.NullNormalize(metric.GetMetric_ll(x));
		    	Tensor2<xy,IF> v = Vec2ObservedByEulObs<xy,LF,IF>(u,x,simData.metric);

                int ij = grid.Index(i,j);
                simData.initialE[ij] = 1;
                simData.initialAngle[ij] = v.Angle();
		    }
        }
	}
    // Initialize simulation:
    Radiation radiation(simData);

    // Run simulation:
    int writeFrequency = 5;
    bool updateFourierHarmonics = false;
    bool keepSourceNodesActive = true;
    bool writeData = true;
    bool printToTerminal = true;
    radiation.RunSimulation(writeFrequency,updateFourierHarmonics,keepSourceNodesActive,writeData,printToTerminal);
}



void RecreateMetricDataForLukasM1Data()
{
    {
	    Int2 N(200,200);
        Coordinate2<xy> Start(-0.1000000000000000056,-1.000000000000000000);
        Coordinate2<xy> End(5.969500000000000028,4.970000000000001528);
        Grid2D<xy> grid(N,Start,End);
        KerrSchild<xy> metric(grid,1.0,0.0);

        ofstream fileAlpha("output/alphas_close.txt");
        ofstream fileBetaX("output/BXs_close.txt");
        ofstream fileBetaY("output/BYs_close.txt");
        fileAlpha << fixed << setprecision(16);
        fileBetaX << fixed << setprecision(16);
        fileBetaY << fixed << setprecision(16);
        for(int j=0; j<grid.n2; j++)
        {
            for(int i=0; i<grid.n1; i++)
            {
                int ij = grid.Index(i,j);
                Tensor2<xy,LF> beta = metric.GetBeta_l(ij);
                fileAlpha << metric.GetAlpha(ij);
                fileBetaX << beta[1];
                fileBetaY << beta[2];
                if(i != grid.n1-1)
                {
                    fileAlpha << " ";
                    fileBetaX << " ";
                    fileBetaY << " ";
                }
            }
            fileAlpha << "\n";
            fileBetaX << "\n";
            fileBetaY << "\n";
        }
        fileAlpha.close();
        fileBetaX.close();
        fileBetaY.close();
    }

    {
	    Int2 N(200,200);
        Coordinate2<xy> Start(0.0,3.999999999999999556);
        Coordinate2<xy> End(7.960000000000001741,8.974999999999999645);
        Grid2D<xy> grid(N,Start,End);
        KerrSchild<xy> metric(grid,1.0,0.0);

        ofstream fileAlpha("output/alphas_far.txt");
        ofstream fileBetaX("output/BXs_far.txt");
        ofstream fileBetaY("output/BYs_far.txt");
        fileAlpha << fixed << setprecision(16);
        fileBetaX << fixed << setprecision(16);
        fileBetaY << fixed << setprecision(16);
        for(int j=0; j<grid.n2; j++)
        {
            for(int i=0; i<grid.n1; i++)
            {
                int ij = grid.Index(i,j);
                Tensor2<xy,LF> beta = metric.GetBeta_l(ij);
                fileAlpha << metric.GetAlpha(ij);
                fileBetaX << beta[1];
                fileBetaY << beta[2];
                if(i != grid.n1-1)
                {
                    fileAlpha << " ";
                    fileBetaX << " ";
                    fileBetaY << " ";
                }
            }
            fileAlpha << "\n";
            fileBetaX << "\n";
            fileBetaY << "\n";
        }
        fileAlpha.close();
        fileBetaX.close();
        fileBetaY.close();
    }
}



int main()
{
    RingWaveFlatRot(1,12,200,0.5);
    RingWaveFlatRot(1,24,200,0.5);
    RingWaveFlatRot(1,36,200,0.5);
    RingWaveFlatRot(1,48,200,0.5);
    RingWaveFlatRot(1,60,200,0.5);
    // CurvedBeamCloseRot(0.01,60,200,5,0.5);
    // RingWaveFlatDirected(1,60,200,0.5);
    // RingWaveFlatUniform(60,200,0.5);
    return 0;

    //// Ring Wave:
    //{
    //    double sigma = 1;
    //    int nDirs[3] =
    //    { 10, 20, 30 };
    //    int nGrids[2] =
    //    { 100, 200 };
    //    for(int nD=0; nD<3; nD++)
    //        for(int nG=0; nG<2; nG++)
    //            RingWave(1,nDirs[nD],nGrids[nG]);
    //}
    //return 0;

    //// Beam Crossing:
    //{
    //    int nDirs[6] =
    //    { 10, 20, 30, 40, 50, 60 };
    //    double sigmas[21] =
    //    { 0.010, 0.012, 0.014, 0.016, 0.018, 0.020, 0.022, 0.024, 0.026, 0.028, 0.030, 0.032, 0.034, 0.036, 0.038, 0.040, 0.042, 0.044, 0.046, 0.048, 0.050 };
    //    double cfls[3] =
    //    { 1.00, 0.50, 0.25 };
    //    for(int c=0; c<3; c++)
    //    for(int s=0; s<21; s++)
    //        BeamCrossing(sigmas[s],60,100,cfls[c]);
    //}
    //return 0;

    //// Curved Beam Close:
    //{
    //    int nDirs[4] =
    //    { 20, 30, 40, 50 };
    //    double sigmas[3] =
    //    { 0.02, 0.03, 0.04 };
    //    double nGrids[2] =
    //    { 100, 200 };
    //    double cfls[3] =
    //    { 1.0, 0.5, 0.25 };
    //
    //    // void CurvedBeamClose(double sigma_, int nDir_, int nGrid_, double cfl_)
    //    for(int nG=0; nG<2; nG++)
    //        for(int nCfl=0; nCfl<3; nCfl++)
    //            for(int nD=0; nD<4; nD++)
    //                for(int nS=0; nS<3; nS++)
    //                    CurvedBeamClose(sigmas[nS],nDirs[nD],nGrids[nG],cfls[nCfl]);
    //}
    //return 0;


}