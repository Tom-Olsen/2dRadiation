#include <iostream>     // cin/cout
#include <math.h>       // basic maths
#include <fstream>      // file input/output
#include <stdlib.h>     // srand, rand
#include <time.h>       // time
#include <filesystem>   // folder/file management
#include <iomanip>      // std::put:time
#include <jsoncpp/json/json.h>
#include "src/Includes.hh"
#include "src/eigen/Eigen/Dense"
// #include "src/Stencil.hh"

using namespace std;
//using namespace tensors;

void Test_Utiliy()
{
    cout << "------------------------------------------------------------" << endl;
    cout << "Utility Test:" << endl;
    cout << Format(1.23456789) << endl;
    cout << Format(-15) << endl;
    cout << Format(0) << endl;
    cout << Format(-1.0*0.0) << endl;
    // exit_on_error("Test!");

    cout << "------------------------------------------------------------" << endl << endl;
}



void Test_TensorTypes()
{
    cout << "------------------------------------------------------------" << endl;
    cout << "TensorTypes Test:" << endl;

    cout << "Coordinate2:" << endl;
    {
        Coordinate2<xy>  xyCoord(1.0,2.0);
        Coordinate2<rph> rphCoord = xyCoord.Transform<rph>();
        Coordinate2<xy>  xyCoord_ = rphCoord.Transform<xy>();
        xyCoord.Print("xy");
        rphCoord.Print("rph");
        xyCoord_.Print("xy_", true);
    }
    cout << "Tensor2:" << endl;
    {
        Coordinate2<xy>  xyCoord(1.1, 0.9);
        Coordinate2<rph> rphCoord = xyCoord.Transform<rph>();
        Tensor2<xy,IF> vxy(1.0,2.0);
        Tensor2<rph,IF> vrph = vxy.Transform<rph>(xyCoord);
        Tensor2<xy,IF> vxy_ = vrph.Transform<xy>(xyCoord.Transform<rph>());
        vxy.Print("vxy");
        vrph.Print("vrph");
        vxy_.Print("vxy_");
        PrintDouble(vxy.Norm(Tensor2x2<xy,IF>(1,0, 0,1)),"|vxy| ");
        PrintDouble(vrph.Norm(Tensor2x2<rph,IF>(1,0, 0,rphCoord[1]*rphCoord[1])),"|vrph|", true);
    }
    cout << "Tensor3:" << endl;
    {
        Coordinate2<xy>  xyCoord(1.5, 2.3);
        Coordinate2<rph> rphCoord = xyCoord.Transform<rph>();
        Tensor3<xy,IF> vxy(-1.0,0.4,sqrt(1-0.4*0.4));
        Tensor3<rph,IF> vrph = vxy.Transform<rph>(xyCoord);
        Tensor3<xy,IF> vxy_ = vrph.Transform<xy>(xyCoord.Transform<rph>());
        vxy.Print("vxy");
        vrph.Print("vrph");
        vxy_.Print("vxy_");
        PrintDouble(vxy.Norm(Tensor3x3<xy,IF>(-1,0,0, 0,1,0, 0,0,1)),"|vxy| ");
        PrintDouble(vrph.Norm(Tensor3x3<rph,IF>(-1,0,0, 0,1,0, 0,0,rphCoord[1]*rphCoord[1])),"|vrph|", true);
    }
    //cout << "Double<N>:" << endl;
    //{
    //    Double<5> a(0.0,1.0,2.0,3.0,4.0,5.0);
    //    a.Print("a");
    //    for(int i=0; i<5; i++)
    //        cout << a[i] << ",";
    //    cout << endl;
    //}
    cout << "------------------------------------------------------------" << endl << endl;
}



void Test_Grid()
{
    cout << "------------------------------------------------------------" << endl;
    cout << "Grid2D Test:" << endl;
    cout << "The file 'grid.json' has been created in the folder 'output'." << endl;
    // using Coord = xy;
    using Coord = rph;
    Int2 N = {100,100};
    // Coordinate2<Coord> Start = {1,1};
    // Coordinate2<Coord> End = {4,4};
    Coordinate2<Coord> Start = {1,0};
    Coordinate2<Coord> End = {4,M_PI_2};
    Grid2D<Coord> grid(N,Start,End);
    double one[grid.n12];
    double E[grid.n12];
    double Fx[grid.n12];
    double Fy[grid.n12];
    
    for(int j=0; j<grid.n2; j++)
    for(int i=0; i<grid.n1; i++)
    {
        int ij = grid.Index(i,j);
        one[ij] = 1;
        E[ij]  = grid.rCoord(i,j);
        Fx[ij] = grid.xCoord(i,j);
        Fy[ij] = grid.yCoord(i,j);
    }
    
    cout << "Integral: " << grid.Integrate(one) << endl;

    std::filesystem::path currentPath = std::filesystem::current_path();
	std::string directoryPath = currentPath/"output";
    grid.WriteFrametoJson(0,E,Fx,Fy,one,0,directoryPath,"grid");
    cout << "------------------------------------------------------------" << endl << endl;
}



void Test_Eigen()
{
    cout << "------------------------------------------------------------" << endl;
    cout << "LGS_Solver Test:" << endl;
    using namespace Eigen;

    // Data arrays:
    double A[] = {1,2,3,  4,5,6,  7,8,10};
    double b[] = {3,3,4};
    double x[3];

    // Expression templates:
    Map<Matrix<double,3,3,RowMajor>> mA(A);
    Map<Vector3d> vb(b);
    Map<Vector3d> vx(x);

    // Solve linear system of equations:
    vx = mA.partialPivLu().solve(vb);

    // Output results:
    cout << "Here is the matrix A:\n" << mA << endl;
    cout << "Here is the vector b:\n" << vb << endl;
    cout << "The solution to Ax = b is:\n" << vx << endl;

    // Invert Matrix:
    Matrix<double,3,3,RowMajor> mAinv = mA.inverse();
    cout << "The inverse Matrix A^-1 is:" << endl;
    cout << mAinv << endl;

    // 4x4 matrix det and inv test:
    double B[] = {1,1,1,1,  1,2,3,4,  1,1,2,2, 1,2,2,1};
    Map<Matrix<double,4,4,RowMajor>> mB(B);
    cout << "     B:\n" << mB << endl;
    cout << "det(B):\n" << mB.determinant() << endl;
    cout << "B^(-1):\n" << mB.inverse() << endl;

    cout << "------------------------------------------------------------" << endl << endl;
}



void Test_Stencil()
{
    cout << "------------------------------------------------------------" << endl;
    cout << "Stencil Test:" << endl;
    constexpr int N = 6;
    X2 distribution;
    UniformStencil uniStencil = UniformStencil(N);
    DirectedStencil dirStencil = DirectedStencil(N,distribution);

    for(int k=0; k<uniStencil.nDir; k++)
        uniStencil.Cxy(k).Print("c");
    for(int k=0; k<uniStencil.nDir; k++)
        Double2(uniStencil.Phi(k),uniStencil.W(k)).Print("phi, w");

    for(int k=0; k<dirStencil.nDir; k++)
        dirStencil.Cxy(k).Print("c");
    for(int k=0; k<dirStencil.nDir; k++)
        Double2(dirStencil.Phi(k),dirStencil.W(k)).Print("phi, w");
        
    // Stencil rotation test:
    {
        constexpr int N = 6;
        double alpha = -M_PI / 8.0;
        DirectedStencil stencil = DirectedStencil(N,distribution);

        for(int i=0; i<N; i++)
            PrintDouble(stencil.Phi(i),"phi");
        for(int i=0; i<N; i++)
            stencil.Cxy(i).Print("cxy");
        cout << stencil.Index(0.087266) << endl << endl; // should be 0

            
        for(int i=0; i<N; i++)
            PrintDouble(stencil.Phi(i,alpha),"phi");
        for(int i=0; i<N; i++)
            stencil.Cxy(i,alpha).Print("cxy");
        cout << stencil.Index(0.087266,alpha) << endl;  // should be between 0,1 or 4,5 depending of sign of alpha
    }

    // Grid of rotated stencils:
//    {
//        using Coord = xy;
//        Int2 N = {10,10};
//        Coordinate2<Coord> Start = {1,1};
//        Coordinate2<Coord> End = {4,4};
//        Grid2D<Coord> grid(N,Start,End);
//
//        double alpha[grid.n12];
//        for(int i=0; i<grid.n1; i++)
//        for(int j=0; j<grid.n2; j++)
//        {
//            int ij = grid.Index(i,j);
//            alpha[ij] = grid.phCoord(i,j) - M_PI / 2.0;
//        }
//
//        constexpr int nDir = 10;
//        DirectedStencil<nDir>& stencil = DirectedStencil<nDir>::GetInstance();
//
//
//        // main body:
//        Json::Value jsonData;
//
//        // arrays:
//        Json::Value positions(Json::arrayValue);
//        Json::Value directions(Json::arrayValue);
//        Json::Value colors(Json::arrayValue);
//
//        // set data:
//        for(int i=0; i<grid.n1; i++)
//        for(int j=0; j<grid.n2; j++)
//        {
//            int ij = grid.Index(i,j);
//            Coordinate2<xy> x = grid.xyCoord(i,j);
//            for(int d=0; d<nDir; d++)
//            {
//                Tensor2<xy,IF> v = stencil.Cxy(d,alpha[ij]);
//                
//                Json::Value position;
//                position["x"] = x[1];
//                position["y"] = x[2];
//                positions.append(position);
//
//                Json::Value direction;
//                double length = min(grid.d1,grid.d2);
//                direction["x"] = length * v[1];
//                direction["y"] = length * v[2];
//                directions.append(direction);
//
//                Json::Value color;
//                color["r"] = 0.8;
//                color["g"] = 0;
//                color["b"] = 0;
//                color["a"] = 1;
//                colors.append(color);
//            }
//        }
//        jsonData["positions"] = positions;
//        jsonData["directions"] = directions;
//        jsonData["colors"] = colors;
//
//        // write json to file:
//        ofstream file("output/RotatedStencilOnGrid.json");
//        file << jsonData;
//        file.close();
//    }

    //integral with both stencils:
    {
        constexpr int N = 50;
        X2 distribution;
        UniformStencil uniStencil = UniformStencil(N);
        DirectedStencil dirStencil = DirectedStencil(N,distribution);

        double uniI[N];
        for(int d=0; d<N; d++)
            uniI[d] = ((double) rand() / (RAND_MAX));
            
        
		double dirI[N];
		for(int k=0; k<N; k++)
		{
			double phi = dirStencil.Phi(k);	// angle of new stencil
			double d = uniStencil.Index(phi);      // d value of above angle in old stencil

			int dFloor = floor(d);
			int m1 = (dFloor - 1 + N) % N;
			int p0 = (dFloor + 0 + N) % N;
			int p1 = (dFloor + 1 + N) % N;
			int p2 = (dFloor + 2 + N) % N;

			dirI[k] = CubicInterpolation(d - dFloor, uniI[m1], uniI[p0], uniI[p1], uniI[p2]);
        }

        double uniIntegral = 0;
        double dirIntegral = 0;
        for(int d=0; d<N; d++)
        {
            uniIntegral += uniI[d] * uniStencil.W(d);
            dirIntegral += dirI[d] * dirStencil.W(d);
        }
        Double2(uniIntegral, dirIntegral).Print("Integrals");
    }

    cout << "------------------------------------------------------------" << endl << endl;
}


// TODOTOM:
// -Bilinear Interpolation not fully adepted to polar grid, yet.
void Test_Interpolation1()
{
    cout << "Bilinear Interpolation:" << endl;
    cout << "Use 'Table to Points' filter in Paraview on the file 'output/BilInt.csv'." << endl;
    // using Coord = xy;
    using Coord = rph;

    Int2 N = {2,2};
    Coordinate2<Coord> Start = {1,1};
    Coordinate2<Coord> End = {3,3};
    Grid2D<Coord> grid(N,Start,End);
    double f00 = 0;
    double f01 = 1;
    double f10 = 2;
    double f11 = 3;
    ofstream file("output/BilInt.csv");
    file << "x,y,z,color" << endl;
    int res = 20;
    for(int j=0; j<res; j++)
    for(int i=0; i<res; i++)
    {
        double color = BilinearInterpolation(i/(res-1.0),j/(res-1.0),f00,f01,f10,f11);
        Coordinate2 xy = grid.xyCoord(i/(res-1.0),j/(res-1.0));
        file << xy[1] << "," << xy[2] << "," << 0 << "," << color << endl;
    }
    file.close();
    cout << endl;
}
void Test_Interpolation2()
{
    cout << "Bicubic Interpolation:" << endl;
    cout << "Use 'Table to Points' filter in Paraview on the file 'output/BicuInt.csv'." << endl;
    using Coord = xy;
    // using Coord = rph;

    Int2 N = {4,4};
    Coordinate2<Coord> Start = {1,1};
    Coordinate2<Coord> End = {3,3};
    Grid2D<Coord> grid(N,Start,End);
    double fm1m1 = 1.0; double fm1p0 = 0.0; double fm1p1 = 0.0; double fm1p2 = 1.0;
    double fp0m1 = 0.0; double fp0p0 = 1.0; double fp0p1 = 1.0; double fp0p2 = 0.0;
    double fp1m1 = 0.0; double fp1p0 = 1.0; double fp1p1 = 1.0; double fp1p2 = 0.0;
    double fp2m1 = 1.0; double fp2p0 = 0.0; double fp2p1 = 0.0; double fp2p2 = 1.0;
    ofstream file("output/BicuInt.csv");
    file << "x,y,z,color" << endl;
    int res = 20;
    for(int j=0; j<res; j++)
    for(int i=0; i<res; i++)
    {
        double color = BicubicInterpolation(i/(res-1.0),j/(res-1.0),
        fm1m1, fm1p0, fm1p1, fm1p2,
        fp0m1, fp0p0, fp0p1, fp0p2,
        fp1m1, fp1p0, fp1p1, fp1p2,
        fp2m1, fp2p0, fp2p1, fp2p2);
        Coordinate2 xy = grid.xyCoord(i/(res-1.0),j/(res-1.0));
        file << xy[1] << "," << xy[2] << "," << 0 << "," << color << endl;
    }
    file.close();
    cout << endl;
}
void Test_Interpolation5()
{
    cout << "Vector Interpolation Test: open 'output/VectorInterpolation.json' with unity." << endl;

    // Give the LBM velocity set random intensities/lengths:
    constexpr int N = 6;
    UniformStencil stencil = UniformStencil(N);
    // DirectedStencil stencil = DirectedStencil(N);
    double I[N];
    srand(time(NULL));
    for(int d=0; d<N; d++)
        I[d] = 4.0 * ((double) rand() / (RAND_MAX)) + 1.0;


    // New stencil set for interpolation:
    constexpr int M = 2*N;
    // UniformStencil stencilNew = UniformStencil(M);
    X2 distribution;
    DirectedStencil stencilNew = DirectedStencil(M,distribution);
    double INew[M];
    for(int d=0; d<M; d++)
    {
        double phi = stencilNew.Phi(d);
	    double k = stencil.Index(phi);
        cout << k << ", " << phi << endl;

	    int kFloor = floor(k);
	    int m1 = (kFloor - 1 + N) % N;
	    int p0 = (kFloor + 0 + N) % N;
	    int p1 = (kFloor + 1 + N) % N;
	    int p2 = (kFloor + 2 + N) % N;
        // cout << k << ", " << k-kFloor << ": " << m1 << ", " << p0 << ", " << p1 << ", " << p2 << endl;
        INew[d] = CubicInterpolation(k-kFloor,I[m1],I[p0],I[p1],I[p2]);
    }

    // main body:
    Json::Value jsonData;

    // arrays:
    Json::Value positions(Json::arrayValue);
    Json::Value directions(Json::arrayValue);
    Json::Value colors(Json::arrayValue);
    for(int i=0; i<M; i++)
    {
        Json::Value position;
        position["x"] = 0;
        position["y"] = 0;
        positions.append(position);

        Json::Value direction;
        direction["x"] = INew[i] * stencilNew.Cx(i);
        direction["y"] = INew[i] * stencilNew.Cy(i);
        directions.append(direction);
        
        Json::Value color;
        color["r"] = 1;
        color["g"] = 0;
        color["b"] = 0;
        color["a"] = 1;
        colors.append(color);
    }
    for(int i=0; i<N; i++)
    {
        Json::Value position;
        position["x"] = 0;
        position["y"] = 0;
        positions.append(position);

        Json::Value direction;
        direction["x"] = I[i] * stencil.Cx(i);
        direction["y"] = I[i] * stencil.Cy(i);
        directions.append(direction);
        
        Json::Value color;
        color["r"] = 0;
        color["g"] = 1;
        color["b"] = 0;
        color["a"] = 1;
        colors.append(color);
    }
    jsonData["positions"] = positions;
    jsonData["directions"] = directions;
    jsonData["colors"] = colors;

    // write json to file:
    ofstream file("output/VectorInterpolation.json");
    file << jsonData;
    file.close();
}
//void Test_Interpolation6()
//{
//    cout << "Cubic Angular Vector Interpolation Test: open 'output/VectorInterpolation.vtk' with paraview." << endl;
//    // Create direction vectors that don't match the velocity set:
//    UniformStencil stencil = UniformStencil(8);
//    int dir = 4*stencil.nDir;
//    Tensor2<xy,IF> v[dir];
//    double vIntensity[dir];
//    for(int d=0; d<dir; d++)
//    {
//        double phi = (2*M_PI/dir)*d;
//        v[d][1]  = cos(phi);
//        v[d][2]  = sin(phi);
//        v[d].Print("v");
//    }
//
//    // Give the LBM velocity set random intensities/lengths:
//    srand(time(NULL));
//    double I[stencil.nDir];
//    for(int d=0; d<stencil.nDir; d++)
//        I[d] = ((double) rand() / (RAND_MAX));
//        
//    for(int d=0; d<dir; d++)
//    {
//        Int2 index = GetNearestVectors(v[d],stencil);
//        double phiLower  = stencil.Phi(index[0]);
//        double phiHigher = stencil.Phi(index[1]);
//        double phi = (2*M_PI/dir)*d;
//        if(phiLower>phiHigher)
//        {
//            phi       = std::fmod(phi       + M_PI, 2*M_PI); 
//            phiLower  = std::fmod(phiLower  + M_PI, 2*M_PI); 
//            phiHigher = std::fmod(phiHigher + M_PI, 2*M_PI); 
//        }
//        double x = Map01(phi,phiLower,phiHigher);
//        vIntensity[d] = CubicInterpolation(x,I[index[0]-1],I[index[0]],I[index[0]+1],I[index[0]+2]);
//        
//        cout.precision(6);
//        cout << std::fixed;
//        cout << "index0:" << index[0] << ", index1:" << index[1] << ", I[i0]:" << I[index[0]] << ", I[d]:" << vIntensity[d] << ", I[i1]:" << I[index[1]] << endl;
//    }
//
//    // Output results to file:
//    ofstream fileVTK("output/VectorInterpolationCubic.vtk");
//
//    // Overhead:
//    fileVTK << "# vtk DataFile Version 3.0" << endl;
//    fileVTK << "Velocity vectors based on table" << endl;
//    fileVTK << "ASCII" << endl << endl;
//
//    // Vector Startpunkte:
//    fileVTK << "DATASET POLYDATA" << endl;
//    fileVTK << "POINTS " << dir+stencil.nDir << " float" << endl;
//    for(int i=0; i<dir+stencil.nDir; i++)
//        fileVTK << "0.0 0.0 0.0" << endl;
//    fileVTK << endl;
//
//    // Vector Endpunkte
//    fileVTK << "POINT_DATA " << dir+stencil.nDir << endl;
//    fileVTK << "VECTORS point_vectors float" << endl;
//    for(int d=0; d<stencil.nDir; d++)
//        fileVTK << I[d]*stencil.Cx(d) << " " << I[d]*stencil.Cy(d) << " " << 0 << endl;
//    for(int d=0; d<dir; d++)
//        fileVTK << vIntensity[d]*v[d][1] << " " << vIntensity[d]*v[d][2] << " " << 0 << endl;
//    fileVTK << endl;
//
//    // Farben:
//    fileVTK << "SCALARS colour float 1" << endl;
//    fileVTK << "LOOKUP_TABLE default" << endl;
//    for(int d=0; d<stencil.nDir; d++)
//        fileVTK << 1 << endl;
//    for(int d=0; d<dir; d++)
//        fileVTK << 0 << endl;
//    fileVTK.close();
//
//    cout << "The file 'VectorInterpolationCubic.vtk' has been created in the folder 'output'." << endl;
//    cout << "Load it with Paraview, add Glyph filter, set options to:" << endl;
//    cout << "Scalar Array: point_vectors" << endl;
//    cout << "Glyph Mode: All Points" << endl;
//    cout << endl;
//}
void Test_Interpolation()
{
    cout << "------------------------------------------------------------" << endl;
    cout << "Interpolation Test:" << endl;
    // Test_Interpolation1();
    // Test_Interpolation2();
    // Test_Interpolation3();
    // Test_Interpolation4();
    Test_Interpolation5();
    // Test_Interpolation6();
    cout << "------------------------------------------------------------" << endl << endl;
}



void Test_Metric()
{
    cout << "------------------------------------------------------------" << endl;
    cout << "Metric Test:" << endl;
    cout << "All following quantities are evaluated at the following coordinates:" << endl;
    Coordinate2<xy> coord_xy(2.5,2.5);
    Coordinate2<rph> coord_rph = coord_xy.Transform<rph>();
    coord_xy.Print("xy");
    coord_rph.Print("rph");

    cout << "---------- Cartesian Tests ----------" << endl;
    {
        Int2 N(201,201);
        Coordinate2<xy> Start(1,1);
        Coordinate2<xy> End  (4,4);
        Grid2D<xy> grid(N,Start,End);
        KerrSchild<xy> metric = KerrSchild(grid,1.0,0.0);
        // SchwarzSchild<xy> metric = SchwarzSchild(grid,1.0,0.0);

        // Test Point:
        grid.ij(coord_xy).Print("ij");

        // Some Tensors:
        double alpha = metric.GetAlpha(coord_xy);
        Tensor3x3 g_ll = metric.GetMetric_ll(coord_xy);
        Tensor3x3 g_uu = metric.GetMetric_uu(coord_xy);
        Tensor3x3 tetrad = metric.GetTetrad(coord_xy);
        Tensor2x2 dl_Beta_u = metric.GetDerivBeta_lu(coord_xy);
        Tensor3x3x3 dl_g_ll = metric.GetDerivMetric_lll(coord_xy);
        PrintDouble(alpha,"alpha");
        g_ll.Print("g_ll");
        g_uu.Print("g_uu");
        tetrad.Print("tetrad",true);
        dl_g_ll.Print("dl_g_ll");
        dl_Beta_u.Print("dl_Beta_u");

        // Vectors:
        Tensor3<xy,IF> vIFxy(alpha,0.6*alpha,sqrt(1-0.6*0.6)*alpha);
        Tensor3<rph,IF> vIFrph = vIFxy.Transform<rph>(coord_xy);
        vIFxy.Print("vxy ");
        vIFrph.Print("vrph");
        PrintDouble(vIFxy.Norm(metric.GetMinkowskiMetric_ll(coord_xy)),"|vxy|(IF)");
        PrintDouble(vIFrph.Norm(Tensor3x3<rph,IF>(-1,0,0, 0,1,0, 0,0,coord_rph[1]*coord_rph[1])),"|vrph|(IF)",true);

        Tensor3<xy,LF> vLFxy = vIFxy.Transform<LF>(tetrad);
        vLFxy.Print("v");
        PrintDouble(vLFxy.Norm(g_ll),"|vxy|(LF)");

        // Json:
        std::filesystem::path currentPath = std::filesystem::current_path();
	    std::string directoryPath = currentPath/"output";
        grid.WriteFrametoJson(0, metric.g00_ll, metric.g11_ll, metric.g22_ll, metric.g00_ll, 0, directoryPath, "metric_Cart0");
        grid.WriteFrametoJson(0, metric.g01_ll, metric.g02_ll, metric.g12_ll, metric.g01_ll, 0, directoryPath, "metric_Cart1");
        grid.WriteFrametoJson(0, metric.tetrad00_ul, metric.tetrad11_ul, metric.tetrad22_ul, metric.tetrad00_ul, 0, directoryPath, "tetrad_Cart0");
        grid.WriteFrametoJson(0, metric.tetrad10_ul, metric.tetrad20_ul, metric.tetrad21_ul, metric.tetrad10_ul, 0, directoryPath, "tetrad_Cart1");
    }
    cout << "-------------------------------------" << endl;

    cout << "------------ Polar Tests ------------" << endl;
    {
        Int2 N  = {201,201};
        Coordinate2<rph> Start(1,0);
        Coordinate2<rph> End  (4,M_PI_2);
        Grid2D<rph> grid(N,Start,End);
        SchwarzSchild<rph> metric = SchwarzSchild(grid,1.0,0.0);

        // Test Point:
        grid.ij(coord_rph).Print("ij");
        
        // Some Tensors:
        double alpha = metric.GetAlpha(coord_rph);
        Tensor3x3 g_ll = metric.GetMetric_ll(coord_rph);
        Tensor3x3 g_uu = metric.GetMetric_uu(coord_rph);
        Tensor3x3 tetrad = metric.GetTetrad(coord_rph);
        Tensor2x2 dl_Beta_u = metric.GetDerivBeta_lu(coord_rph);
        Tensor3x3x3 dl_g_ll = metric.GetDerivMetric_lll(coord_rph);
        PrintDouble(alpha,"alpha");
        g_ll.Print("g_ll");
        g_uu.Print("g_uu");
        tetrad.Print("tetrad",true);
        dl_g_ll.Print("dl_g_ll");
        dl_Beta_u.Print("dl_Beta_u");

        // Vectors:
        Tensor3<xy,IF> vIFxy(alpha,0.6*alpha,sqrt(1-0.6*0.6)*alpha);
        Tensor3<rph,IF> vIFrph = vIFxy.Transform<rph>(coord_xy);
        vIFxy.Print("v");
        vIFrph.Print("v");
        PrintDouble(vIFxy.Norm(Tensor3x3<xy,IF>(-1,0,0, 0,1,0, 0,0,1)),"|vxy|(IF)");
        PrintDouble(vIFrph.Norm(metric.GetMinkowskiMetric_ll(coord_rph)),"|vrph|(IF)",true);

        Tensor3<rph,LF> vLFrph = vIFrph.Transform<LF>(metric.GetTetrad(coord_rph));
        vLFrph.Print("v");
        PrintDouble(vLFrph.Norm(g_ll),"|vrph|(LF)");
        
        // Json:
        std::filesystem::path currentPath = std::filesystem::current_path();
	    std::string directoryPath = currentPath/"output";
        grid.WriteFrametoJson(0,metric.g00_ll,metric.g11_ll, metric.g22_ll,metric.g12_ll,0,directoryPath,"metric_Polar");
    }
    cout << "-------------------------------------" << endl;
    cout << "------------------------------------------------------------" << endl << endl;
}



void Test_AdvancedUtility()
{
    cout << "------------------------------------------------------------" << endl;
    cout << "Mostly Tests for transformations between different frames:" << endl << endl;

    // Cartesian Grid and Metric:
    Int2 Nxy  = {100,100};
    Coordinate2<xy> Startxy = {2,2};
    Coordinate2<xy> Endxy   = {3,3};
    Grid2D<xy> gridxy(Nxy,Startxy,Endxy);
    SchwarzSchild<xy> metricxy(gridxy,1.0,0.0);
    
    // Polar Grid and Metric:
    Int2 Nrph  = {100,100};
    Coordinate2<rph> Startrph = {2,0};
    Coordinate2<rph> Endrph   = {3,M_PI_2};
    Grid2D<rph> gridrph(Nrph,Startrph,Endrph);
    SchwarzSchild<rph> metricrph(gridrph,1.0,0.0);

    cout << "Cartesian Tests:" << endl;
    {
        Coordinate2<xy> x(2.2,2.6);
        double alpha = metricxy.GetAlpha(x);
        Tensor3<xy,IF> u(alpha, 0.6*alpha, 0.8*alpha);
        Tensor2<xy,IF> v = Vec2ObservedByEulObs<xy,IF,IF>(u,x,metricxy);
        u.Print("u");
        v.Print("v");
        PrintDouble(u.Norm(metricxy.GetMinkowskiMetric_ll(x)),"|u|");
        PrintDouble(v.Norm(metricxy.GetMinkowskiGamma_uu(x)),"|v|",true);
    }
    {
        Coordinate2<xy> x(2.2,2.6);
        double alpha = metricxy.GetAlpha(x);
        Tensor3<xy,IF> u(alpha, 0.6*alpha, 0.8*alpha);
        Tensor2<xy,LF> v = Vec2ObservedByEulObs<xy,IF,LF>(u,x,metricxy);
        u.Print("u");
        v.Print("v");
        PrintDouble(u.Norm(metricxy.GetMinkowskiMetric_ll(x)),"|u|");
        PrintDouble(v.Norm(metricxy.GetGamma_ll(x)),"|v|",true);
    }
    {
        Coordinate2<xy> x(2.2,2.6);
        double alpha = metricxy.GetAlpha(x);
        Tensor3<xy,IF> uIF(alpha, 0.6*alpha, 0.8*alpha);
        Tensor3<xy,LF> u = uIF.Transform<LF>(metricxy.GetTetrad(x));
        Tensor2<xy,IF> v = Vec2ObservedByEulObs<xy,LF,IF>(u,x,metricxy);
        u.Print("u");
        v.Print("v");
        PrintDouble(u.Norm(metricxy.GetMetric_ll(x)),"|u|");
        PrintDouble(v.Norm(metricxy.GetMinkowskiGamma_uu(x)),"|v|",true);
    }
    {
        Coordinate2<xy> x(2.2,2.6);
        double alpha = metricxy.GetAlpha(x);
        Tensor3<xy,IF> uIF(alpha, 0.6*alpha, 0.8*alpha);
        Tensor3<xy,LF> u = uIF.Transform<LF>(metricxy.GetTetrad(x));
        Tensor2<xy,LF> v = Vec2ObservedByEulObs<xy,LF,LF>(u,x,metricxy);
        u.Print("u");
        v.Print("v");
        PrintDouble(u.Norm(metricxy.GetMetric_ll(x)),"|u|");
        PrintDouble(v.Norm(metricxy.GetGamma_ll(x)),"|v|",true);
    }

    cout << "Polar Tests:" << endl;
    {
        Coordinate2<rph> x(2.3,1.1);
        double alpha = metricrph.GetAlpha(x);
        Tensor3<rph,IF> u(alpha, 0.6*alpha, 0.8*alpha);
        u.NullNormalize(metricrph.GetMinkowskiMetric_ll(x));
        Tensor2<rph,IF> v = Vec2ObservedByEulObs<rph,IF,IF>(u,x,metricrph);
        u.Print("u");
        v.Print("v");
        PrintDouble(u.Norm(metricrph.GetMinkowskiMetric_ll(x)),"|u|");
        PrintDouble(v.Norm(metricrph.GetMinkowskiGamma_uu(x)),"|v|",true);
    }
    {
        Coordinate2<rph> x(2.3,1.4);
        double alpha = metricrph.GetAlpha(x);
        Tensor3<rph,IF> u(alpha, 0.6*alpha, 0.8*alpha);
        u.NullNormalize(metricrph.GetMinkowskiMetric_ll(x));
        Tensor2<rph,LF> v = Vec2ObservedByEulObs<rph,IF,LF>(u,x,metricrph);
        u.Print("u");
        v.Print("v");
        PrintDouble(u.Norm(metricrph.GetMinkowskiMetric_ll(x)),"|u|");
        PrintDouble(v.Norm(metricrph.GetGamma_ll(x)),"|v|",true);
    }
    {
        Coordinate2<rph> x(2.3,1.4);
        double alpha = metricrph.GetAlpha(x);
        Tensor3<rph,IF> uIF(alpha, 0.6*alpha, 0.8*alpha);
        uIF.NullNormalize(metricrph.GetMinkowskiMetric_ll(x));
        Tensor3<rph,LF> u = uIF.Transform<LF>(metricrph.GetTetrad(x));
        Tensor2<rph,IF> v = Vec2ObservedByEulObs<rph,LF,IF>(u,x,metricrph);
        u.Print("u");
        v.Print("v");
        PrintDouble(u.Norm(metricrph.GetMetric_ll(x)),"|u|");
        PrintDouble(v.Norm(metricrph.GetMinkowskiGamma_uu(x)),"|v|",true);
    }
    {
        Coordinate2<rph> x(2.3,1.4);
        double alpha = metricrph.GetAlpha(x);
        Tensor3<rph,IF> uIF(alpha, 0.6*alpha, 0.8*alpha);
        uIF.NullNormalize(metricrph.GetMinkowskiMetric_ll(x));
        Tensor3<rph,LF> u = uIF.Transform<LF>(metricrph.GetTetrad(x));
        Tensor2<rph,LF> v = Vec2ObservedByEulObs<rph,LF,LF>(u,x,metricrph);
        u.Print("u");
        v.Print("v");
        PrintDouble(u.Norm(metricrph.GetMetric_ll(x)),"|u|");
        PrintDouble(v.Norm(metricrph.GetGamma_ll(x)),"|v|",true);
    }
    cout << "------------------------------------------------------------" << endl << endl;
}



void Test_GeodesicRay()
{
    cout << "------------------------------------------------------------" << endl;
    cout << "Create Geodesic photon trajectories." << endl;
    cout << "Open 'output/GeodesicRay.csv' with Paraview and add Table to Points filter." << endl;
    
    ofstream filexy("output/GeodesicRayXY.csv");
    ofstream filerph("output/GeodesicRayRPH.csv");
    filexy << "x, y, z, color" << endl;
    filerph << "x, y, z, color" << endl;

    // Grid & Metrics:
    Int2 Nxy  = {100,100};
    Coordinate2<xy> Startxy = {-10,-10};
    Coordinate2<xy> Endxy   = {10,10};
    Grid2D gridxy(Nxy,Startxy,Endxy);
    SchwarzSchild<xy> metricxy(gridxy,1.0,0.0);
    // KerrSchild<xy> metricxy(gridxy,1.0,0.0);

    Int2 Nrph  = {100,100};
    Coordinate2<rph> Startrph = {0,0};
    Coordinate2<rph> Endrph   = {10*sqrt(2),2*M_PI};
    Grid2D gridrph(Nrph,Startrph,Endrph);
    SchwarzSchild<rph> metricrph(gridrph,1.0,0.0);

    // Horizon:
    int NdirPhi = 1024;
    for(int i=0; i<NdirPhi; i++)
    {
        double phi     = (2*M_PI/NdirPhi)*(i+0.5);
        double xSphere = 2*metricxy.m*cos(phi);
        double ySphere = 2*metricxy.m*sin(phi);
        filexy << xSphere << "," << ySphere << "," << 0 << "," << 1 << endl;
        filerph << xSphere << "," << ySphere << "," << 0 << "," << 1 << endl;
    }
    // Geodesic Rays:
    int NRays = 80;
    for(int i=0; i<NRays; i++)
    {
        Coordinate2<xy> x;
        x[1] = Startxy[1] + (Endxy[1]-Startxy[1])/(NRays+1)*(i+1);
        x[2] = Startxy[2];

        Tensor3<xy,LF> u(1,0,1);
        u.NullNormalize(metricxy.GetMetric_ll(x));
        Tensor2<xy,LF> v = Vec2ObservedByEulObs<xy,LF,LF>(u,x,metricxy);
        
        double s = 1;
        for(int k=0; k<1000; k++)
        {
            // s *= Euler_GeodesicEquation<1>(200.0/1000.0,x,v,metric);
            s *= RK45_GeodesicEquation<1>(200.0/1000.0,x,v,metricxy);
            if(metricxy.InsideBH(x) or gridxy.OutsideDomain(x))
                break;
            filexy << x[1] << ", " << x[2] << ", " << 0 << "," << s << endl;
        }
    }
    for(int i=0; i<NRays; i++)
    {
        Coordinate2<xy> xCart;
        xCart[1] = Startxy[1] + (Endxy[1]-Startxy[1])/(NRays+1)*(i+1);
        xCart[2] = Startxy[2];
        Coordinate2<rph> x = xCart.Transform<rph>();

        Tensor3<xy,LF> uxy(1,0,1);
        Tensor3<rph,LF> u = uxy.Transform<rph>(xCart);
        u.NullNormalize(metricrph.GetMetric_ll(x));
        Tensor2<rph,LF> v = Vec2ObservedByEulObs<rph,LF,LF>(u,x,metricrph);
        
        double s = 1;
        for(int k=0; k<1000; k++)
        {
            // s *= Euler_GeodesicEquation<1>(200.0/1000.0,x,v,metric);
            s *= RK45_GeodesicEquation<1>(200.0/1000.0,x,v,metricrph);
            if(metricrph.InsideBH(x) or gridrph.OutsideDomain(x))
                break;
            xCart = x.Transform<xy>();
            filerph << xCart[1] << ", " << xCart[2] << ", " << 0 << "," << s << endl;
        }
    }
    filexy.close();
    filerph.close();
    cout << "------------------------------------------------------------" << endl << endl;
}



void Test_PhotonSphere()
{
    cout << "------------------------------------------------------------" << endl;
    cout << "Create Geodesic photon trajectory on photon sphere." << endl;
    cout << "Open 'output/PhotonSphere.csv' with Paraview and add Table to Points filter." << endl;
    
    ofstream filexy("output/PhotonSphereXY.csv");
    ofstream filerph("output/PhotonSphereRPH.csv");
    filexy  << "x, y, z, color" << endl;
    filerph << "x, y, z, color" << endl;


    // Grid & Metrics:
    Int2 Nxy  = {100,100};
    Coordinate2<xy> Startxy = {-10,-10};
    Coordinate2<xy> Endxy   = {10,10};
    Grid2D gridxy(Nxy,Startxy,Endxy);
    SchwarzSchild<xy> metricxy(gridxy,1.0,0.0);

    Int2 Nrph  = {100,100};
    Coordinate2<rph> Startrph = {0,0};
    Coordinate2<rph> Endrph   = {10*sqrt(2),2*M_PI};
    Grid2D gridrph(Nrph,Startrph,Endrph);
    SchwarzSchild<rph> metricrph(gridrph,1.0,1.0);

    // Horizon:
    int NdirPhi = 1024;
    for(int i=0; i<NdirPhi; i++)
    {
        double phi     = (2*M_PI/NdirPhi)*(i+0.5);
        double xSphere = 2*metricxy.m*cos(phi);
        double ySphere = 2*metricxy.m*sin(phi);
        filexy  << xSphere << "," << ySphere << "," << 0 << "," << 1 << endl;
        filerph << xSphere << "," << ySphere << "," << 0 << "," << 1 << endl;
    }
    
    // Geodesic Rays:
    {
        Coordinate2<xy> x(0.0,3.0);
        Tensor3<xy,LF> u(1,1,0);
        u.NullNormalize(metricxy.GetMetric_ll(x));
        Tensor2<xy,LF> v = Vec2ObservedByEulObs<xy,LF,LF>(u,x,metricxy);

        double s = 1;
        for(int k=0; k<1000; k++)
        {
            //Euler_GeodesicEquation<1>(200.0/1000.0,x,vLF,metricxy);
            s *= RK45_GeodesicEquation<1>(200.0/1000.0,x,v,metricxy);
            if(metricxy.InsideBH(x) or gridxy.OutsideDomain(x))
                break;
            filexy << x[1] << ", " << x[2] << ", " << 0 << "," << s << endl;
        }
    }
    {
        Coordinate2<rph> x(3.0,M_PI_2);
        Tensor3<rph,LF> u(1,0,-1);
        u.NullNormalize(metricrph.GetMetric_ll(x));
        Tensor2<rph,LF> v = Vec2ObservedByEulObs<rph,LF,LF>(u,x,metricrph);

        double s = 1;
        for(int k=0; k<1000; k++)
        {
            //Euler_GeodesicEquation<1>(200.0/1000.0,x,vLF,metricrph);
            s *= RK45_GeodesicEquation<1>(200.0/1000.0,x,v,metricrph);
            if(metricrph.InsideBH(x) or gridrph.OutsideDomain(x))
                break;
            Coordinate2<xy> xCart = x.Transform<xy>();
            filerph << xCart[1] << ", " << xCart[2] << ", " << 0 << "," << s << endl;
        }
    }
    filexy.close();
    cout << "------------------------------------------------------------" << endl << endl;
}



void Test_FourierQuadrature()
{
    //{
    //    constexpr int N = 16;
    //    UniformStencil stencil = UniformStencil(N);
    //
    //    double data[N];
    //    for(int i=0; i<N; i++)
    //    {
    //        double phi = stencil.Phi(i);
    //        data[i] = 0;
    //        for(int k=0; k<N; k++)
    //            data[i] += (1.0+0.1*k) * Fourier::Basis(k,phi);
    //    }
    //    vector<double> c = Fourier::Expansion::GetCoefficients(stencil,data);
    //    for(int i=0; i<c.size(); i++)
    //        cout << "c" << to_string(i) << ": " << c[i] << endl;
    //
    //    cout << "Test reconstruction of function values:" << endl;
    //    for(int i=0; i<N; i++)
    //    {
    //        double phi = stencil.Phi(i);
    //        Double2(data[i],Fourier::Expansion::GetValue(phi,c)).Print(to_string(i));
    //    }
    //    stencil.Print();
    //}

    {
        constexpr int N = 16;
        RotStencil stencil = RotStencil(N);
        double rotation = 0.2;//stencil.Phi(0);
    
        double data[N];
        for(int i=0; i<N; i++)
        {
            double phi = stencil.Phi(i,rotation);
            data[i] = 0;
            for(int k=0; k<N-1; k++)
                data[i] += (1.0+0.1*k) * Fourier::Basis(k,phi);
        }
        vector<double> c = Fourier::Expansion::GetCoefficients(stencil,data,rotation);
        for(int i=0; i<c.size(); i++)
            cout << "c" << to_string(i) << ": " << c[i] << endl;
    
        cout << "Test reconstruction of function values:" << endl;
        for(int i=0; i<N; i++)
        {
            double phi = stencil.Phi(i,rotation);
            Double2(data[i],Fourier::Expansion::GetValue(phi,c)).Print(to_string(i));
        }
        for(int i=0; i<N; i++)
        {
            Double2(i,stencil.Index(stencil.Phi(i,rotation),rotation)).Print("index test");
        }
        stencil.Print(rotation);
    }
}


// this tests uses old fourier expansion.
void Test_FourierExpansion()
{
    cout << "------------------------------------------------------------" << endl;
    cout << "Test fourier harmonic expansion with stencil streaming:" << endl;
    
    ofstream fileRef("output/FourierReference.csv");
    ofstream fileExp("output/FourierExpansion.csv");
    fileRef << "x, y, z, s" << endl;
    fileExp << "x, y, z, s" << endl;

    // Grid & Metric:
    using Coord = xy;
    Int2 N  = {20,20};
    Coordinate2<Coord> Start = {0,0};
    Coordinate2<Coord> End   = {4,4};
    Grid2D grid(N,Start,End);
    KerrSchild<Coord> metric(grid);
    // using Coord = rph;
    // Int2 N  = {20,20};
    // Coordinate2<Coord> Start = {1,0};
    // Coordinate2<Coord> End   = {4,M_PI_2};
    // Grid2D grid(N,Start,End);
    // SchwarzSchild<Coord> metric(grid);

    // Horizon:
    int NdirPhi = 1024;
    for(int i=0; i<NdirPhi; i++)
    {
        double phi     = (2*M_PI/NdirPhi)*(i+0.5)/4;
        double xSphere = 2*metric.m*cos(phi);
        double ySphere = 2*metric.m*sin(phi);
        fileRef << xSphere << ", " << ySphere << ", " << 0 << ", " << 1 << endl;
        fileExp << xSphere << ", " << ySphere << ", " << 0 << ", " << 1 << endl;
    }

    // Reference Stencil Streaming:
    UniformStencil stencilRef = UniformStencil(32);
    for(int i=1; i<grid.n1-1; i++)
    for(int j=1; j<grid.n2-1; j++)
    for(int d=0; d<stencilRef.nDir; d++)
    {
        double s = 1;
        Coordinate2<Coord> x = grid.x12Coord(i,j);
        Coordinate2<xy> xyCoord = grid.xyCoord(i,j);

        double alpha = metric.GetAlpha(x);
        Tensor2<xy,IF> cxy = stencilRef.Cxy(d);
        Tensor2<Coord,IF> c = cxy.Transform<Coord>(xyCoord);
        Tensor3<Coord,IF> u(alpha, c[1]*alpha, c[2]*alpha);
        Tensor2<Coord,LF> v = Vec2ObservedByEulObs<Coord,IF,LF>(u,x,metric);

    	if(metric.InsideBH(x))
            continue;
        fileRef << xyCoord[1] << ", " << xyCoord[2] << ", " << 0 << ", " << s << std::endl;
        for(int k=0; k<5; k++)
        {
            s *= RK45_GeodesicEquation<-1>(0.5*grid.dt/5.0,x,v,metric);
            if(metric.InsideBH(x) or grid.OutsideDomain(x))
                break;
            xyCoord = x.Transform<xy>();
            fileRef << xyCoord[1] << ", " << xyCoord[2] << ", " << 0 << ", " << s << std::endl;
        }
    }
    fileRef.close();
    
    // Fourier Expansion Stencil Streaming:
    // n=5 needed for redshift due to 4 roots instead of 2.
    constexpr int n = 5;
    UniformStencil stencilExp = UniformStencil(n);
    for(int i=1; i<grid.n1-1; i++)
    for(int j=1; j<grid.n2-1; j++)
    {
        double sData[stencilExp.nDir];
        double x1Data[stencilExp.nDir];
        double x2Data[stencilExp.nDir];
        for(int d=0; d<stencilExp.nDir; d++)
        {
            double s = 1;
            Coordinate2<Coord> x = grid.x12Coord(i,j);
            Coordinate2<xy> xyCoord = grid.xyCoord(i,j);

            double alpha = metric.GetAlpha(x);
            Tensor2<xy,IF> cxy = stencilExp.Cxy(d);
            Tensor2<Coord,IF> c = cxy.Transform<Coord>(xyCoord);
            Tensor3<Coord,IF> u(alpha, c[1]*alpha, c[2]*alpha);
            Tensor2<Coord,LF> v = Vec2ObservedByEulObs<Coord,IF,LF>(u,x,metric);

        	if(metric.InsideBH(x))
            {
                sData[d] = 1;
                x1Data[d] = 0;
                x2Data[d] = 0;
                continue;
            }
            fileExp << xyCoord[1] << ", " << xyCoord[2] << ", " << 0 << ", " << s << std::endl;
            for(int k=0; k<5; k++)
                s *= RK45_GeodesicEquation<-1>(0.5*grid.dt/5.0,x,v,metric);
            sData[d] = s;
            x1Data[d] = x[1];
            x2Data[d] = x[2];
        }
        vector<double> cS  = Fourier::Expansion::GetCoefficients(stencilExp,sData);
        vector<double> cX1 = Fourier::Expansion::GetCoefficients(stencilExp,x1Data);
        vector<double> cX2 = Fourier::Expansion::GetCoefficients(stencilExp,x2Data);

        for(int d=0; d<stencilRef.nDir; d++)
        {
            double phi = stencilRef.Phi(d);
            double s = Fourier::Expansion::GetValue(phi,cS);
            Coordinate2<Coord> x(Fourier::Expansion::GetValue(phi,cX1), Fourier::Expansion::GetValue(phi,cX2));
            Coordinate2<xy> xyCoord = x.Transform<xy>();
            fileExp << xyCoord[1] << ", " << xyCoord[2] << ", " << 0 << ", " << s << std::endl;

            if(std::isnan(xyCoord[1]) || std::isnan(xyCoord[2]) || std::isnan(s))
            {
                Int2(i,j).Print("ij");
                grid.xyCoord(i,j).Print("xy");
                PrintDouble(d,"d");
                //sFourExp.Print("s",true);
                //x1FourExp.Print("x1",true);
                //x2FourExp.Print("x2",true);
                for(int k=0; k<n; k++)
                    Double3(sData[k],x1Data[k],x2Data[k]).Print("s,x1,x2 data");
                exit_on_error();
            }
        }
    }
    fileExp.close();
    cout << "------------------------------------------------------------" << endl << endl;
}




void Test_JsonWriter()
{
    Json::Value jsonData;   
    Json::Value E(Json::arrayValue);
    Json::Value Fx(Json::arrayValue);
    Json::Value Fy(Json::arrayValue);
    for(int i=0; i<10; i++)
    {
        E.append(Json::Value(i));
        Fx.append(Json::Value(2*i));
        Fy.append(Json::Value(i/2.0));
    }

    jsonData["dimensions"] = 2;
    jsonData["numberOfFrames"] = 150;
    jsonData["nx"] = 101;
    jsonData["ny"] = 101;
    jsonData["dx"] = 0.05;
    jsonData["dy"] = 0.05;
    jsonData["xStart"] = 0;
    jsonData["yStart"] = 0;
    jsonData["xEnd"] = 5;
    jsonData["yEnd"] = 5;
    jsonData["frames"]["E"] = E;
    jsonData["frames"]["Fx"] = Fx;
    jsonData["frames"]["Fy"] = Fy;

    ofstream fileOut("output/example.json");
    fileOut << jsonData << endl;

    ifstream fileIn("output/example.json");
    Json::Value value;
    Json::Reader reader;

    reader.parse(fileIn,value);
    cout << value << endl;
}



void Test_FrequencyTransform()
{
    cout << "------------------------------------------------------------" << endl;
    
    Int2 N  = {100,100};
    Coordinate2<xy> Start = {1,1};
    Coordinate2<xy> End   = {3,3};
    Grid2D grid(N,Start,End);
    KerrSchild<xy> metric(grid,1.0,0.0);

    constexpr int Ndir = 8;
    UniformStencil stencil = UniformStencil(Ndir);

    Coordinate2<xy> x(2.0,2.0);
    Tensor3x3<xy,LF> g_ll = metric.GetMetric_ll(x);
    Tensor2x2<xy,LF> gamma_ll = metric.GetGamma_ll(x);

    Json::Value jsonData;
    jsonData["numberOfVectors"] = Ndir;
    Json::Value xValue(Json::arrayValue);
    Json::Value yValue(Json::arrayValue);
    Json::Value vxValue(Json::arrayValue);
    Json::Value vyValue(Json::arrayValue);
    for(int i=0; i<Ndir; i++)
    {
        Tensor2<xy,IF> c = stencil.Cxy(i);

        Tensor3<xy,LF> u(1, c[1], c[2]);
        u.NullNormalize(g_ll);
        Tensor2<xy,LF> v = Vec2ObservedByEulObs<xy,LF,LF>(u,x,metric);

        xValue.append(Json::Value(x[1]));
        yValue.append(Json::Value(x[2]));
        vxValue.append(Json::Value(v[1]));
        vyValue.append(Json::Value(v[2]));
    }
    jsonData["x"] = xValue;
    jsonData["y"] = yValue;
    jsonData["vx"] = vxValue;
    jsonData["vy"] = vyValue;

    ofstream file("output/stencilLF.json");
    file << jsonData;
    file.close();

    // Test: Eulerian 4-vel in comoving IF
    Tensor3x3<xy,Tetrad> tetrad = metric.GetTetrad(x);
    double alpha = metric.GetAlpha(x);
    Tensor2<xy,LF> beta = metric.GetBeta_u(x);
    Tensor3<xy,LF> nLF(1.0/alpha, -beta[1]/alpha, -beta[2]/alpha);
    Tensor3<xy,IF> nIF = nLF.Transform<IF>(tetrad);

    nLF.Print("n");
    nIF.Print("n");

    // Test: Frequency in LF and IF:
    double nu = 3.1415;
    Tensor3<xy,LF> pLF(nu, 1, 1);
    pLF.NullNormalize(g_ll);
    pLF.Print("p");
    PrintDouble(pLF.Norm(g_ll), "|p|");

    Tensor3<xy,IF> pIF = pLF.Transform<IF>(tetrad);
    pIF.Print("p");
    PrintDouble(pIF.Norm(metric.GetMinkowskiMetric_ll(x)),"|p|");

    double nuLF = 0;
    for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
        nuLF += g_ll[{i,j}] * pLF[i] * nLF[j];
        
    double nuIF = 0;
    for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
        nuIF += metric.GetMinkowskiMetric_ll(x)[{i,j}] * pIF[i] * nIF[j];

    PrintDouble(nuLF,"nuLF");
    PrintDouble(nuIF,"nuIF");

    cout << "------------------------------------------------------------" << endl << endl;
}

template<class Coord>
void CheckInsideBH(Metric2D<Coord>& metric, Coordinate2<Coord> point)
{
    cout << metric.InsideBH(point) << endl;
}
void Test_Virtual_Override()
{
    Int2 N(10,10);
    Coordinate2<xy> Start(-1, -1);
    Coordinate2<xy> End(1, 1);
    Grid2D<xy> grid(N,Start,End);
    Minkowski<xy> metric(grid,1.0,0.0);

    Coordinate2<xy> point(0,0);
    cout << metric.InsideBH(point) << endl;
    CheckInsideBH(metric, point);
}



int main()
{
    // Test_Utiliy();
    // Test_TensorTypes();
    // Test_Grid();
    // Test_Eigen();
    // Test_Stencil();
    // Test_Interpolation();
    // Test_Metric();
    // Test_AdvancedUtility();
    // Test_GeodesicRay();
    // Test_PhotonSphere();
    Test_FourierQuadrature();
    // Test_FourierExpansion();
    // Test_JsonWriter();
    // Test_FrequencyTransform();
    // Test_Virtual_Override();

    // UniformStencil b = UniformStencil(10);
    // Stencil& a = b;
    // a.Print();
    // b.Print();
}