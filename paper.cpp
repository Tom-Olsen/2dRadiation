#include <iostream>     // cin/cout
#include <math.h>       // basic maths
#include <fstream>      // file input/output
#include <stdlib.h>     // srand, rand
#include <time.h>       // time
#include <filesystem>   // folder/file management
#include <iomanip>      // std::put:time
#include "src/Includes.hh"

using namespace std;



void RingWaveFlat(double simTime_, int nDir_, int nGrid_, double cfl_, double sigma_, bool writeData_, StreamingType streamingType_)
{
    string fileName = "RingWaveFlat NN" + std::to_string(nGrid_) +
        " N" + std::to_string(nDir_) +
        " s" + std::to_string(sigma_) +
        " cfl"  + std::to_string(cfl_) + " " + StreamingName(streamingType_);

    Int2 N(nGrid_,nGrid_);
    double halfWidth = simTime_ + 0.1;
    Coordinate2<xy> Start(-halfWidth, -halfWidth);
    Coordinate2<xy> End(halfWidth, halfWidth);
    Grid2D<xy> grid(N,Start,End);
    grid.SetCFL(cfl_);
    Minkowski<xy> metric(grid,1.0,0.0);
    DynamicStencil stencil(nDir_,sigma_);
    StaticStencil fourierStencil(5);
    
    // Initialize simulation:
    Radiation radiation(grid,metric,stencil,fourierStencil, streamingType_);

    // Set initial data:
    {
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
    }
    
    // Run simulation:
    RunParameters params =
    {
        .name = fileName,
        .simTime = simTime_,
        .writeFrequency = 2,
        .updateFourierHarmonics = false,
        .keepSourceNodesActive = false,
        .writeData = writeData_,
        .printToTerminal = true
    };
    radiation.RunSimulation(params);
}



void BeamCrossing(double simTime_, int nDir_, int nGrid_, double cfl_, double sigma_, StreamingType streamingType_)
{
    string fileName = "BeamCrossing NN" + std::to_string(nGrid_) +
        " N" + std::to_string(nDir_) +
        " s" + std::to_string(sigma_) +
        " cfl"  + std::to_string(cfl_) + " " + StreamingName(streamingType_);

    Int2 N(nGrid_,nGrid_);
    Coordinate2<xy> Start(-1.1, -1.1);
    Coordinate2<xy> End(1.1, 1.1);
    Grid2D<xy> grid(N,Start,End);
    grid.SetCFL(cfl_);
    Minkowski<xy> metric(grid,1.0,0.0);
    DynamicStencil stencil(nDir_,sigma_);
    StaticStencil fourierStencil(5);
 
    // Initialize simulation:
    Radiation radiation(grid,metric,stencil,fourierStencil, streamingType_);

    for(int j=0; j<grid.n2; j++)
    for(int i=0; i<grid.n1; i++)
    {
        if(-0.2 < grid.xCoord(i,j) && grid.xCoord(i,j) < 0.2 && grid.yCoord(i,j) < -0.8)
        {// x€[-0.2,0.2] y<-0.8
            int ij = grid.Index(i,j);
            radiation.isInitialGridPoint[ij] = true;
            radiation.initialE[ij] = 1;
            radiation.initialRotation[ij] = M_PI / 2.0;
            radiation.initialKappa0[ij] = 0;
            radiation.initialKappa1[ij] = 0;
            radiation.initialKappaA[ij] = 0;
            radiation.initialEta[ij] = 0;
        }
        if(-0.2 < grid.yCoord(i,j) && grid.yCoord(i,j) < 0.2 && grid.xCoord(i,j) < -0.8)
        {// y€[-0.2,0.2] x<-0.8
            int ij = grid.Index(i,j);
            radiation.isInitialGridPoint[ij] = true;
            radiation.initialE[ij] = 1;
            radiation.initialRotation[ij] = 2.0 * M_PI;
            radiation.initialKappa0[ij] = 0;
            radiation.initialKappa1[ij] = 0;
            radiation.initialKappaA[ij] = 0;
            radiation.initialEta[ij] = 0;
        }
    }
    
    //for(int j=0; j<grid.n2; j++)
    //{
    //    for(int i=0; i<grid.n1; i++)
    //    {
    //        if(-0.5 < grid.yCoord(i,j) && grid.yCoord(i,j) < -0.3 && grid.xCoord(i,j) < -0.8)
    //        {// x€[-0.2,0.2] y<-0.8
    //            int ij = grid.Index(i,j);
    //            radiation.isInitialGridPoint[ij] = true;
    //            radiation.initialE[ij] = 1;
    //            radiation.initialRotation[ij] = M_PI / 8.0;
    //            radiation.initialKappa0[ij] = 0;
    //            radiation.initialKappa1[ij] = 0;
    //            radiation.initialKappaA[ij] = 0;
    //            radiation.initialEta[ij] = 0;
    //        }
    //        if(0.3 < grid.yCoord(i,j) && grid.yCoord(i,j) < 0.5 && grid.xCoord(i,j) < -0.8)
    //        {// y€[-0.2,0.2] x<-0.8
    //            int ij = grid.Index(i,j);
    //            radiation.isInitialGridPoint[ij] = true;
    //            radiation.initialE[ij] = 1;
    //            radiation.initialRotation[ij] = 2.0 * M_PI - M_PI / 8.0;
    //            radiation.initialKappa0[ij] = 0;
    //            radiation.initialKappa1[ij] = 0;
    //            radiation.initialKappaA[ij] = 0;
    //            radiation.initialEta[ij] = 0;
    //        }
    //    }
    //}

    // Run simulation:
    RunParameters params =
    {
        .name = fileName,
        .simTime = simTime_,
        .writeFrequency = 5,
        .updateFourierHarmonics = false,
        .keepSourceNodesActive = true,
        .writeData = true,
        // .writeData = false,
        .printToTerminal = true
    };
    radiation.RunSimulation(params);
}


    
void CurvedBeamClose(double simTime_, int nDir_, int nGrid_, int nFourier_, double cfl_, double sigma_, StreamingType streamingType_)
{
    string fileName = "CurvedBeamClose NN" + std::to_string(nGrid_) +
        " N" + std::to_string(nDir_) + 
        " F" + std::to_string(nFourier_) + 
        " s" + std::to_string(sigma_) +
        " cfl" + std::to_string(cfl_) + " " + StreamingName(streamingType_);

    Int2 N(nGrid_,nGrid_);
    Coordinate2<xy> Start(-0.1, -1.0);
    Coordinate2<xy> End(6.0, 5.0);
    Grid2D<xy> grid(N,Start,End);
    grid.SetCFL(cfl_);
    // KerrSchild<xy> metric(grid,1.0,0.0);
    SchwarzSchild<xy> metric(grid,1.0,0.0);
    // Minkowski<xy> metric(grid,1.0,0.0);
    DynamicStencil stencil(nDir_,sigma_);
    StaticStencil fourierStencil(nFourier_);
    
    // Initialize simulation:
    Radiation radiation(grid,metric,stencil,fourierStencil, streamingType_);

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
    RunParameters params =
    {
        .name = fileName,
        .simTime = simTime_,
        .writeFrequency = 20,
        // .updateFourierHarmonics = false,
        .updateFourierHarmonics = true,
        .keepSourceNodesActive = true,
        .writeData = true,
        // .writeData = false,
        .printToTerminal = true
    };
    radiation.RunSimulation(params);
}



void CurvedBeamFar(double simTime_, int nDir_, int nGrid_, int nFourier_, double cfl_, double sigma_, StreamingType streamingType_)
{
    string fileName = "CurvedBeamFar NN" + std::to_string(nGrid_) +
        " N" + std::to_string(nDir_) + 
        " F" + std::to_string(nFourier_) + 
        " s" + std::to_string(sigma_) +
        " cfl" + std::to_string(cfl_) + " " + StreamingName(streamingType_);

	Int2 N(nGrid_,nGrid_);
    Coordinate2<xy> Start(-0.1, 3.9);
    Coordinate2<xy> End(8.1, 9.1);
    Grid2D<xy> grid(N,Start,End);
    grid.SetCFL(cfl_);
    // KerrSchild<xy> metric(grid,1.0,0.0);
    SchwarzSchild<xy> metric(grid,1.0,0.0);
    // Minkowski<xy> metric(grid,1.0,0.0);
    DynamicStencil stencil(nDir_,sigma_);
    StaticStencil fourierStencil(nFourier_);

    // Initialize simulation:
    Radiation radiation(grid,metric,stencil,fourierStencil, streamingType_);

	for(int j=0; j<grid.n2; j++)
	{
	    for(int i=0; i<grid.n1; i++)
        {
		    if(7.0 < grid.yCoord(i,j) && grid.yCoord(i,j) < 8.0 && grid.xCoord(i,j) < 0.0)
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
    RunParameters params =
    {
        .name = fileName,
        .simTime = simTime_,
        .writeFrequency = 20,
        .updateFourierHarmonics = false,
        .keepSourceNodesActive = true,
        .writeData = true,
        // .writeData = false,
        .printToTerminal = true
    };
    radiation.RunSimulation(params);
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



void FourierHarmonicsError(int fourierOrder)
{
    string fileName = "output/FourierHarmonicsError" + std::to_string(fourierOrder);
    std::ofstream fileOut(fileName + ".txt");
    fileOut << " s      \t x1      \t x2      \t v1      \t v2\n";


    Coordinate2<xy> x0(sqrt(2)+0.1, sqrt(2)+0.1);
    // Coordinate2<xy> x0(1,3);
    double t = 0.1;
    int steps = 100;
    int nDir = 50;
    
    Int2 N(200,200);
    Coordinate2<xy> Start(x0[1]-0.2, x0[2]-0.2);
    Coordinate2<xy> End(x0[1]+0.2, x0[2]+0.2);
    Grid2D<xy> grid(N,Start,End);
    KerrSchild<xy> metric(grid,1.0,0.0);
    // SchwarzSchild<xy> metric(grid,1.0,0.0);
    // Minkowski<xy> metric(grid,1.0,0.0);
    StaticStencil fourierStencil(fourierOrder);
    double alpha = metric.GetAlpha(x0);

	double dataS [fourierStencil.nDir];
	double dataX1[fourierStencil.nDir];
	double dataX2[fourierStencil.nDir];
    double dataV1[fourierStencil.nDir];
    double dataV2[fourierStencil.nDir];

    // Get geodesic data and save to file:
    for(int d=0; d<fourierStencil.nDir; d++)
    {
        double s = 1;
        Coordinate2<xy> x = x0;
        Tensor2<xy,IF> vIF = fourierStencil.Cxy(d);
        Tensor3<xy,IF> uIF(alpha, vIF[1]*alpha, vIF[2]*alpha);
        Tensor2<xy,LF> vLF = Vec2ObservedByEulObs<xy,IF,LF>(uIF,x0,metric);
        
        for(int k=0; k<steps; k++)
            s *= RK45_GeodesicEquation<-1>(t/steps,x,vLF,metric);
            
		dataS[d]  = s;
		dataX1[d] = x[1];
		dataX2[d] = x[2];
		dataV1[d] = vLF[1];
		dataV2[d] = vLF[2];
    }
    
    // Create fourier coefficients:
	std::vector<double> cS  = Fourier::Expansion::GetCoefficients(fourierStencil,dataS );
	std::vector<double> cX1 = Fourier::Expansion::GetCoefficients(fourierStencil,dataX1);
	std::vector<double> cX2 = Fourier::Expansion::GetCoefficients(fourierStencil,dataX2);
	std::vector<double> cV1 = Fourier::Expansion::GetCoefficients(fourierStencil,dataV1);
	std::vector<double> cV2 = Fourier::Expansion::GetCoefficients(fourierStencil,dataV2);
    
    // Write fourier coefficients to file:
    fileOut << "Coefficients:\n";
    for(int d=0; d<fourierStencil.nDir; d++)
    {
        fileOut << Format(cS[d] ) << "\t";
        fileOut << Format(cX1[d]) << "\t";
        fileOut << Format(cX2[d]) << "\t";
        fileOut << Format(cV1[d]) << "\t";
        fileOut << Format(cV2[d]) << "\n";
    }fileOut << "\n";

    // Calculate nDir Photons via Geodesic and via Fourier and compare:
    double absErrorS  = 0;
    double absErrorX1 = 0;
    double absErrorX2 = 0;
    double absErrorV1 = 0;
    double absErrorV2 = 0;
    double relErrorS  = 0;
    double relErrorX1 = 0;
    double relErrorX2 = 0;
    double relErrorV1 = 0;
    double relErrorV2 = 0;

    DynamicStencil stencil(nDir);
    for(int d=0; d<stencil.nDir; d++)
    {
        double phi = stencil.Phi(d);
        double s = 1;
        double s0 = s;
        Coordinate2<xy> x = x0;
        Tensor2<xy,IF> vIF = stencil.Cxy(d);
        Tensor3<xy,IF> uIF(alpha, vIF[1]*alpha, vIF[2]*alpha);
        Tensor2<xy,LF> vLF = Vec2ObservedByEulObs<xy,IF,LF>(uIF,x0,metric);
        Tensor2<xy,LF> vLF0 = vLF;
        
        for(int k=0; k<steps; k++)
            s *= RK45_GeodesicEquation<-1>(t/steps,x,vLF,metric);
            
        //errorS  = max(errorS , 100 * abs(Fourier::Expansion::GetValue(phi,cS ) - s     ) / abs(s      - s0     ));
        //errorX1 = max(errorX1, 100 * abs(Fourier::Expansion::GetValue(phi,cX1) - x[1]  ) / abs(x[1]   - x0[1]  ));
        //errorX2 = max(errorX2, 100 * abs(Fourier::Expansion::GetValue(phi,cX2) - x[2]  ) / abs(x[2]   - x0[2]  ));
        //errorV1 = max(errorV1, 100 * abs(Fourier::Expansion::GetValue(phi,cV1) - vLF[1]) / abs(vLF[1] - vLF0[1]));
        //errorV2 = max(errorV2, 100 * abs(Fourier::Expansion::GetValue(phi,cV2) - vLF[2]) / abs(vLF[2] - vLF0[2]));

        double errorS  = abs(Fourier::Expansion::GetValue(phi,cS ) - s     );
        double errorX1 = abs(Fourier::Expansion::GetValue(phi,cX1) - x[1]  );
        double errorX2 = abs(Fourier::Expansion::GetValue(phi,cX2) - x[2]  );
        double errorV1 = abs(Fourier::Expansion::GetValue(phi,cV1) - vLF[1]);
        double errorV2 = abs(Fourier::Expansion::GetValue(phi,cV2) - vLF[2]);
        if (absErrorS < errorS)
        {
            absErrorS = errorS;
            relErrorS = 100 * absErrorS / abs(s);
        }
        if (absErrorX1 < errorX1)
        {
            absErrorX1 = errorX1;
            relErrorX1 = 100 * absErrorX1 / abs(x[1]);
        }
        if (absErrorX2 < errorX2)
        {
            absErrorX2 = errorX2;
            relErrorX2 = 100 * absErrorX2 / abs(x[2]);
        }
        if (absErrorV1 < errorV1)
        {
            absErrorV1 = errorV1;
            relErrorV1 = 100 * absErrorS / abs(vLF[1]);
        }
        if (absErrorV2 < errorV2)
        {
            absErrorV2 = errorV2;
            relErrorV2 = 100 * absErrorV2 / abs(vLF[2]);
        }
    }

    fileOut << "Errors:\n";
    fileOut << Format(relErrorS ) << "\t";
    fileOut << Format(relErrorX1) << "\t";
    fileOut << Format(relErrorX2) << "\t";
    fileOut << Format(relErrorV1) << "\t";
    fileOut << Format(relErrorV2);

    fileOut.close();
}

// TODO:
// -test curved beam with the new idea (see Radiation::StreamCurvedDynamic())
int main()
{
    CurvedBeamClose(10, 100, 200, 5, 0.5, 0.03, StreamingType::CurvedDynamic);

    // Geodesic streaming does not give same results as fourier streaming atm!
    // CurvedBeamClose(10, 100, 200, 5, 0.5, 0.03, StreamingType::GeodesicDynamic);
    // CurvedBeamClose(10, 100, 200, 3, 0.5, 0.03, StreamingType::CurvedDynamic);
    // CurvedBeamFar  (12, 100, 200, 5, 0.5, 0.02, StreamingType::GeodesicDynamic);
    // CurvedBeamFar  (12, 100, 200, 5, 0.5, 0.02, StreamingType::CurvedDynamic);

    // These give high contrast between static and dynamic streaming with same initial data:
    // RingWaveFlat(1.0, 20, 200, 0.5, 1.25, true, StreamingType::FlatDynamic);
    // RingWaveFlat(1.0, 30, 200, 0.5, 1.25, true, StreamingType::FlatDynamic);
    // RingWaveFlat(1.0, 40, 200, 0.5, 1.25, true, StreamingType::FlatDynamic);
    // RingWaveFlat(1.0, 50, 200, 0.5, 1.25, true, StreamingType::FlatDynamic);
    // RingWaveFlat(1.0, 20, 200, 0.5, 1.25, true, StreamingType::FlatStatic);
    // RingWaveFlat(1.0, 30, 200, 0.5, 1.25, true, StreamingType::FlatStatic);
    // RingWaveFlat(1.0, 40, 200, 0.5, 1.25, true, StreamingType::FlatStatic);
    // RingWaveFlat(1.0, 50, 200, 0.5, 1.25, true, StreamingType::FlatStatic);
    
    // Long Ring Wave Simulations:
    // RingWaveFlat(5.0, 50, 200, 0.5, 1.25, true, StreamingType::FlatStatic);
    // RingWaveFlat(5.0, 50, 200, 0.5, 1.25, true, StreamingType::FlatDynamic);

    // These do not oversaturate:
    // BeamCrossing(2.0, 20, 100, 0.5, 0.001, StreamingType::FlatStatic); // comparison to dynamic stencil
    // BeamCrossing(2.0, 30, 100, 0.5, 0.150, StreamingType::FlatDynamic);
    // BeamCrossing(2.0, 40, 100, 0.5, 0.130, StreamingType::FlatDynamic);
    // BeamCrossing(2.0, 50, 100, 0.5, 0.120, StreamingType::FlatDynamic);
    // BeamCrossing(2.0, 60, 100, 0.5, 0.110, StreamingType::FlatDynamic);
    // BeamCrossing(2.0, 70, 100, 0.5, 0.105, StreamingType::FlatDynamic);
    // BeamCrossing(2.0, 80, 100, 0.5, 0.100, StreamingType::FlatDynamic);
    // BeamCrossing(2.0, 90, 100, 0.5, 0.090, StreamingType::FlatDynamic);
    // BeamCrossing(2.0,100, 100, 0.5, 0.090, StreamingType::FlatDynamic);
    
    // These give good resuls:
    // CurvedBeamClose(10,  50, 200, 5, 0.5, 0.05, StreamingType::CurvedDynamic);
    // CurvedBeamClose(10, 100, 200, 5, 0.5, 0.03, StreamingType::CurvedDynamic);
    // CurvedBeamClose(10, 200, 200, 5, 0.5, 0.02, StreamingType::CurvedDynamic);
    // CurvedBeamClose(10,  50, 200, 5, 0.5, 0.12, StreamingType::CurvedStatic);
    // CurvedBeamClose(10, 100, 200, 5, 0.5, 0.08, StreamingType::CurvedStatic);
    // CurvedBeamClose(10, 200, 200, 5, 0.5, 0.05, StreamingType::CurvedStatic);

    // These give good resuls:
    // CurvedBeamFar(12,  50, 200, 5, 0.5, 0.04, StreamingType::CurvedDynamic);
    // CurvedBeamFar(12, 100, 200, 5, 0.5, 0.02, StreamingType::CurvedDynamic);
    // CurvedBeamFar(12, 200, 200, 5, 0.5, 0.01, StreamingType::CurvedDynamic);
    // CurvedBeamFar(12,  50, 200, 5, 0.5, 0.08, StreamingType::CurvedStatic);
    // CurvedBeamFar(12, 100, 200, 5, 0.5, 0.06, StreamingType::CurvedStatic);
    // CurvedBeamFar(12, 200, 200, 5, 0.5, 0.04, StreamingType::CurvedStatic);

    // Performance comparison, static vs dynamic stancil in flat spacetime:
    // for(int i=0; i<5; i++)
    // {
        // double sigma = 1 + 0.1*i;
        // RingWaveFlat(10.0, 60, 400, 0.5, sigma, false, StreamingType::FlatStatic);
        // RingWaveFlat(10.0, 60, 400, 0.5, sigma, false, StreamingType::FlatDynamic);
        // RingWaveFlat(10.0, 30, 400, 0.5, sigma, false, StreamingType::FlatDynamic);
    // }

    // Benchmark Fourier Harmonics:
    // FourierHarmonicsError(3);
    // FourierHarmonicsError(5);
    // FourierHarmonicsError(7);
    // FourierHarmonicsError(9);

    // Performance comparison, Flat streaming vs Fourier (static & dynamic) curved streaming:
    // CurvedBeamClose(10, 200, 200, 3, 0.5, 0.02, StreamingType::CurvedDynamic);
    // CurvedBeamClose(10, 200, 200, 5, 0.5, 0.02, StreamingType::CurvedDynamic);
    // CurvedBeamClose(10, 200, 200, 7, 0.5, 0.02, StreamingType::CurvedDynamic);
    // CurvedBeamClose(10, 200, 200, 9, 0.5, 0.02, StreamingType::CurvedDynamic);
    // CurvedBeamClose(10, 200, 200, 3, 0.5, 0.02, StreamingType::FlatStatic);
}