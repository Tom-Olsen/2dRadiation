#include <iostream>     // cin/cout
#include <math.h>       // basic maths
#include <fstream>      // file input/output
#include <stdlib.h>     // srand, rand
#include <time.h>       // time
#include <filesystem>   // folder/file management
#include <iomanip>      // std::put:time
#include "src/Includes.hh"

using namespace std;



//void FourierExpansionCoefficients(int nDir, int Solver)
//{
//    // Grid & Metric:
//    using Coord = xy;
//    Int2 N  = {41,41};
//    Coordinate2<Coord> Start = {0,0};
//    Coordinate2<Coord> End   = {4,4};
//    Grid2D grid(N,Start,End);
//    KerrSchild<Coord> metric(grid);
//    // using Coord = rph;
//    // using SpaceTime = SchwarzSchild;
//    // Int2 N  = {100,100};
//    // Coordinate2<Coord> Start = {1,0};
//    // Coordinate2<Coord> End   = {4,M_PI_2};
//    // Grid2D grid(N,Start,End);
//    // Metric2D<Coord,SpaceTime> metric(grid);
//
//    // Stencil Streaming:
//    Stencil stencil = UniformStencil(nDir);
//    // Coordinate2<Coord> x0 = Coordinate2<xy>{sqrt(2)+0.1,sqrt(2)+0.1}.Transform<Coord>();
//    Coordinate2<Coord> x0 = Coordinate2<xy>{3,3}.Transform<Coord>();
//    double s0 = 1;
//    double xData[nDir];
//    double yData[nDir];
//    double sData[nDir];
//    for(int d=0; d<stencil.nDir; d++)
//    {
//        Coordinate2<Coord> x = x0;
//        double s = s0;
//
//        double alpha = metric.GetAlpha(x);
//        Tensor2<xy,IF> cxy = stencil.Cxy(d);
//        Tensor2<Coord,IF> c = cxy.Transform<Coord>(x);
//        Tensor3<Coord,IF> u(alpha, c[1]*alpha, c[2]*alpha);
//        Tensor2<Coord,LF> v = Vec2ObservedByEulObs<Coord,IF,LF>(u,x,metric);
//
//    	if(metric.InsideBH(x) or grid.OutsideDomain(x))
//            continue;
//        if(Solver == 0)
//            s *= RK45_GeodesicEquation<-1>(0.5*grid.dt,x,v,metric);
//        if(Solver == 1)
//            s *= Euler_GeodesicEquation<-1>(0.5*grid.dt,x,v,metric);
//        xData[d] = x[1];
//        yData[d] = x[2];
//        sData[d] = s;
//    }
//
//    FourierHarmonicsExpansionN<nDir> xFourExp; xFourExp.SetCoefficients(stencil,xData);
//    FourierHarmonicsExpansionN<nDir> yFourExp; yFourExp.SetCoefficients(stencil,yData);
//    FourierHarmonicsExpansionN<nDir> sFourExp; sFourExp.SetCoefficients(stencil,sData);
//
//    if(Solver == 0)
//        cout << "RK45 solver:" << endl;
//    if(Solver == 1)
//        cout << "Euler solver:" << endl;
//    xFourExp.Print("x Expansion",true,8);
//    yFourExp.Print("y Expansion",true,8);
//    sFourExp.Print("s Expansion",true,8);
//}



//void FourierExpansionRelativeError(int nDir)
//{
//    // Grid & Metric:
//    using Coord = xy;
//    Int2 N  = {41,41};
//    Coordinate2<Coord> Start = {0,0};
//    Coordinate2<Coord> End   = {4,4};
//    Grid2D grid(N,Start,End);
//    SchwarzSchild<Coord> metric(grid);
//    // KerrSchild<Coord> metric(grid);
//    // using Coord = rph;
//    // using SpaceTime = SchwarzSchild;
//    // Int2 N  = {100,100};
//    // Coordinate2<Coord> Start = {1,0};
//    // Coordinate2<Coord> End   = {4,M_PI_2};
//    // Grid2D grid(N,Start,End);
//    // Metric2D<Coord,SpaceTime> metric(grid);
//
//    // Stencil Streaming:
//    Stencil stencil = UniformStencil(nDir);
//    Coordinate2<Coord> x0 = Coordinate2<xy>{sqrt(2)+0.1,sqrt(2)+0.1}.Transform<Coord>();
//    double s0 = 1;
//    double xData[nDir];
//    double yData[nDir];
//    double sData[nDir];
//    for(int d=0; d<stencil.nDir; d++)
//    {
//        Coordinate2<Coord> x = x0;
//        double s = s0;
//
//        double alpha = metric.GetAlpha(x);
//        Tensor2<xy,IF> cxy = stencil.Cxy(d);
//        Tensor2<Coord,IF> c = cxy.Transform<Coord>(x);
//        Tensor3<Coord,IF> u(alpha, c[1]*alpha, c[2]*alpha);
//        Tensor2<Coord,LF> v = Vec2ObservedByEulObs<Coord,IF,LF>(u,x,metric);
//
//    	if(metric.InsideBH(x) or grid.OutsideDomain(x))
//            continue;
//        s *= RK45_GeodesicEquation<-1>(grid.dt,x,v,metric);
//        // s *= Euler_GeodesicEquation<-1>(grid.dt,x,v,metric);
//        xData[d] = x[1];
//        yData[d] = x[2];
//        sData[d] = s;
//    }
//
//    FourierHarmonicsExpansionN<nDir> xFourExp; xFourExp.SetCoefficients(stencil,xData);
//    FourierHarmonicsExpansionN<nDir> yFourExp; yFourExp.SetCoefficients(stencil,yData);
//    FourierHarmonicsExpansionN<nDir> sFourExp; sFourExp.SetCoefficients(stencil,sData);
//
//    cout << fixed;
//    cout.precision(8);
//
//    // Reference Stencil Streaming:
//    constexpr int nDirRef = 100;
//    Stencil stencilRef = UniformStencil(nDirRef);
//    double xDataRef[nDirRef];
//    double yDataRef[nDirRef];
//    double sDataRef[nDirRef];
//    for(int d=0; d<nDirRef; d++)
//    {
//        Coordinate2<Coord> x = x0;
//        double s = s0;
//
//        double alpha = metric.GetAlpha(x);
//        Tensor2<xy,IF> cxy = stencilRef.Cxy(d);
//        Tensor2<Coord,IF> c = cxy.Transform<Coord>(x);
//        Tensor3<Coord,IF> u(alpha, c[1]*alpha, c[2]*alpha);
//        Tensor2<Coord,LF> v = Vec2ObservedByEulObs<Coord,IF,LF>(u,x,metric);
//
//    	if(metric.InsideBH(x) or grid.OutsideDomain(x))
//            continue;
//        s *= RK45_GeodesicEquation<-1>(grid.dt,x,v,metric);
//        xDataRef[d] = x[1];
//        yDataRef[d] = x[2];
//        sDataRef[d] = s;
//    }
//
//    // max rel error:
//    double maxRelErrorX = 0;
//    double maxRelErrorY = 0;
//    double maxRelErrorS = 0;
//    ofstream fileRef("output/stencil streaming reference.txt");
//    ofstream fileExp("output/stencil streaming " + to_string(nDir) + ".txt");
//    fileRef << "x,y,z,s" << endl;
//    fileExp << "x,y,z,s" << endl;
//    for(int d=0; d<nDirRef; d++)
//    {
//        double phi = stencilRef.Phi(d);
//        double xRef = xDataRef[d];
//        double xExp = xFourExp.GetValue(phi);
//        double yRef = yDataRef[d];
//        double yExp = yFourExp.GetValue(phi);
//        double sRef = sDataRef[d];
//        double sExp = sFourExp.GetValue(phi);
//
//        Coordinate2<Coord> x12Ref{xRef,yRef};
//        Coordinate2<xy> xyRef = x12Ref.Transform<xy>();
//        Coordinate2<Coord> x12Exp{xExp,yExp};
//        Coordinate2<xy> xyExp = x12Exp.Transform<xy>();
//
//        fileRef << xyRef[1] << "," << xyRef[2] << ","  << 0 << ","<< sRef << endl;
//        fileExp << xyExp[1] << "," << xyExp[2] << ","  << 0 << ","<< sExp << endl;
//
//        double relErrorX = abs((xRef - xExp) / xRef);
//        double relErrorY = abs((yRef - yExp) / yRef);
//        double relErrorS = abs((sRef - sExp) / sRef);
//        if (maxRelErrorX < relErrorX)
//            maxRelErrorX = relErrorX;
//        if (maxRelErrorY < relErrorY)
//            maxRelErrorY = relErrorY;
//        if (maxRelErrorS < relErrorS)
//            maxRelErrorS = relErrorS;
//    }
//    fileRef.close();
//    fileExp.close();
//
//    cout << 100*maxRelErrorX << "\t" << 100*maxRelErrorY << "\t" << 100*maxRelErrorS << endl;
//}
//void FourierExpansionRelativeError()
//{
//    cout << "x,y,s errors:" << endl;
//    for(int i=3; i<=35; i+=2)
//    {
//        FourierExpansionRelativeError(i);
//    }
//}



//template<class Coord>
//void StartSimulation(int nGrid1, int nGrid2, float simTime,
//bool updateFourierHarmonics, bool keepSourceNodesActive, bool writeData, bool printToTerminal)
//{
//	// Stencil:
//	Stencil uniStencil = UniformStencil<nDir>::GetInstance();
//	Stencil dirStencil = DirectedStencil<nDir>::GetInstance();
//	Stencil momentStencil = DirectedStencil<2*nDir>::GetInstance();
//	UniformStencil<FourierOrder>& fourierStencil = UniformStencil<FourierOrder>::GetInstance();
//
//	Int2 N(nGrid1,nGrid2);
//    Coordinate2<Coord> Start;
//    Coordinate2<Coord> End;
//    if constexpr(std::is_same<Coord,xy>::value)
//    {
//        Start = Coordinate2<Coord>(0,0);
//        End = Coordinate2<Coord>(4,4);
//    }
//    if constexpr(std::is_same<Coord,rph>::value)
//    {
//        Start = Coordinate2<Coord>(1,0);
//        End = Coordinate2<Coord>(4,M_PI_2);
//    }
//    Grid2D<Coord> grid(N,Start,End);
//    SchwarzSchild<Coord> metric(grid,1.0,0.0);
//	SimulationData simData(metric,uniStencil,dirStencil,momentStencil,fourierStencil,simTime,"a");
//
//	int d;
//	bool init = false;
//	for(int j=0; j<grid.n2; j++)
//	{
//		int i = 1;
//		if(3.0 < grid.yCoord(i,j) and grid.yCoord(i,j) < 3.5)
//		{
//			// Set Initial Data in LF, single directions:
//			Coordinate2<Coord> x = grid.xyCoord(i,j);
//			Tensor3<Coord,LF> u(1,1,0);
//			u.NullNormalize(metric.GetMetric_ll(x));
//			Tensor2<Coord,IF> v = Vec2ObservedByEulObs<Coord,LF,IF>(u,x,simData.metric);
//			Tensor2<xy,IF> vxy = v.template Transform<xy>(x);
//
//			simData.initialE[grid.Index(i,j)] = 1.0;
//			simData.initialAngle[grid.Index(i,j)] = vxy.Angle();
//		}
//	}
//	Radiation radiation(simData);
//    radiation.RunSimulation(1,updateFourierHarmonics,keepSourceNodesActive,writeData,printToTerminal);
//}
//void UpdateFourierExpansionTimeMeasurements()
//{
//    int N = 10;
//    constexpr int nX = 100;
//    constexpr int nY = 100;
//    constexpr double duration = 2;
//
//
//    bool updateFourierHarmonics = false;
//    bool keepSourceNodesActive = true;
//    bool writeData = false;
//    bool printToTerminal = false;
//
//    // In: void Radiation<Coord,FourierOrder>::RunSimulation(), choose Flat or Geodesic streaming.
//    //cout << "Reference" << endl;
//    //for(int i=0; i<N; i++)
//    //{
//    //    StartSimulation<xy,5,200>(nX,nY,duration,updateFourierHarmonics,keepSourceNodesActive,writeData,printToTerminal);
//    //}
//
//    cout << "Fourier Order = 5:" << endl;
//    for(int i=0; i<N; i++)
//    {
//        StartSimulation<xy,5,200>(nX,nY,duration,updateFourierHarmonics,keepSourceNodesActive,writeData,printToTerminal);
//    }
//    cout << endl << endl;
//    
//    cout << "Fourier Order = 7:" << endl;
//    for(int i=0; i<N; i++)
//    {
//        StartSimulation<xy,7,200>(nX,nY,duration,updateFourierHarmonics,keepSourceNodesActive,writeData,printToTerminal);
//    }
//    cout << endl << endl;
//
//    cout << "Fourier Order = 9:" << endl;
//    for(int i=0; i<N; i++)
//    {
//        StartSimulation<xy,9,200>(nX,nY,duration,updateFourierHarmonics,keepSourceNodesActive,writeData,printToTerminal);
//    }
//    cout << endl << endl;
//
//    cout << "Fourier Order = 11:" << endl;
//    for(int i=0; i<N; i++)
//    {
//        StartSimulation<xy,11,200>(nX,nY,duration,updateFourierHarmonics,keepSourceNodesActive,writeData,printToTerminal); 
//    }
//    cout << endl << endl;
//
//    cout << "Fourier Order = 13:" << endl;
//    for(int i=0; i<N; i++)
//    {
//        StartSimulation<xy,13,200>(nX,nY,duration,updateFourierHarmonics,keepSourceNodesActive,writeData,printToTerminal); 
//    }
//    cout << endl << endl;
//
//    cout << "Fourier Order = 15:" << endl;
//    for(int i=0; i<N; i++)
//    {
//        StartSimulation<xy,15,200>(nX,nY,duration,updateFourierHarmonics,keepSourceNodesActive,writeData,printToTerminal);
//    }
//    cout << endl << endl;
//
//    cout << "Fourier Order = 17:" << endl;
//    for(int i=0; i<N; i++)
//    {
//        StartSimulation<xy,17,200>(nX,nY,duration,updateFourierHarmonics,keepSourceNodesActive,writeData,printToTerminal);
//    }
//    cout << endl << endl;
//
//    cout << "Fourier Order = 19:" << endl;
//    for(int i=0; i<N; i++)
//    {
//        StartSimulation<xy,19,200>(nX,nY,duration,updateFourierHarmonics,keepSourceNodesActive,writeData,printToTerminal); 
//    }
//    cout << endl << endl;
//
//    cout << "Fourier Order = 21:" << endl;
//    for(int i=0; i<N; i++)
//    {
//        StartSimulation<xy,21,200>(nX,nY,duration,updateFourierHarmonics,keepSourceNodesActive,writeData,printToTerminal);
//    }
//    cout << endl << endl;
//
//    cout << "Fourier Order = 23:" << endl;
//    for(int i=0; i<N; i++)
//    {
//        StartSimulation<xy,23,200>(nX,nY,duration,updateFourierHarmonics,keepSourceNodesActive,writeData,printToTerminal);
//    }
//    cout << endl << endl;
//
//    cout << "Fourier Order = 25:" << endl;
//    for(int i=0; i<N; i++)
//    {
//        StartSimulation<xy,25,200>(nX,nY,duration,updateFourierHarmonics,keepSourceNodesActive,writeData,printToTerminal);
//    }
//}



void RingWaveFlatDirected(double sigma_, int nDir_, int nGrid_, float cfl_)
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



void RingWaveFlatUniform(int nDir_, int nGrid_, float cfl_)
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
    // FourierExpansionCoefficients<9,0>();
    // FourierExpansionCoefficients<9,1>();
    // FourierExpansionRelativeError();
    // UpdateFourierExpansionTimeMeasurements();

    // SourceSinkTest();
    // CurvedBeamClose();
    // CurvedBeamFar();
    // PolarEnergyConservation();
    // GeodesicPolarEnergyConservation();
    // CartesianEnergyConservation();
    // GeodesicCartEnergyConservation();
    // RecreateMetricDataForLukasM1Data();
    
    // BeamCrossing(0.024,60,200);
    //cout << "Directed 12" << endl;
    //for(int i=0; i<10; i++)
    //    RingWaveFlatDirected(1,12,200,0.5);
    //cout << "Directed 24" << endl;
    //for(int i=0; i<10; i++)
    //    RingWaveFlatDirected(1,24,200,0.5);
    //cout << "Directed 36" << endl;
    //for(int i=0; i<10; i++)
    //    RingWaveFlatDirected(1,36,200,0.5);
    //cout << "Directed 48" << endl;
    //for(int i=0; i<10; i++)
    //    RingWaveFlatDirected(1,48,200,0.5);
    //cout << "Directed 60" << endl;
    //for(int i=0; i<10; i++)
    //    RingWaveFlatDirected(1,60,200,0.5);
    //cout << "Uniform 12" << endl;
    //for(int i=0; i<10; i++)
    //    RingWaveFlatUniform(12,200,0.5);
    //cout << "Uniform 24" << endl;
    //for(int i=0; i<10; i++)
    //    RingWaveFlatUniform(24,200,0.5);
    //cout << "Uniform 36" << endl;
    //for(int i=0; i<10; i++)
    //    RingWaveFlatUniform(36,200,0.5);
    //cout << "Uniform 48" << endl;
    //for(int i=0; i<10; i++)
    //    RingWaveFlatUniform(48,200,0.5);
    //cout << "Uniform 60" << endl;
    //for(int i=0; i<10; i++)
    //    RingWaveFlatUniform(60,200,0.5);
    // CurvedBeamClose(0.04,50,100,0.5);
    // CurvedBeamFar(0.03,70,200,0.5);

    RingWaveFlatDirected(1,60,200,0.5);
    RingWaveFlatUniform(60,200,0.5);
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

    // Curved Beam Close:
    {
        int nDirs[4] =
        { 20, 30, 40, 50 };
        double sigmas[3] =
        { 0.02, 0.03, 0.04 };
        double nGrids[2] =
        { 100, 200 };
        double cfls[3] =
        { 1.0, 0.5, 0.25 };

        // void CurvedBeamClose(double sigma_, int nDir_, int nGrid_, double cfl_)
        for(int nG=0; nG<2; nG++)
            for(int nCfl=0; nCfl<3; nCfl++)
                for(int nD=0; nD<4; nD++)
                    for(int nS=0; nS<3; nS++)
                        CurvedBeamClose(sigmas[nS],nDirs[nD],nGrids[nG],cfls[nCfl]);
    }
    return 0;


}