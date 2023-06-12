#include <iostream>
#include "../src/Spacetimes.h"
#include "../src/SpecialMath.h"
using namespace std;



void Tetrad()
{
    size_t nx, ny;
    nx = ny = 200;
    Coord start(1,1);
    Coord end(3,3);
    Grid grid(nx, ny, start, end);
    KerrSchild metric(grid, 1.0, 0.0);
    // SchwarzSchild metric(grid, 1.0, 0.0);

    Coord xy(1.5, 2);
    Tensor3x3 g_ll = metric.GetMetric_ll(xy);
    Tensor3x3 tetrad = metric.GetTetrad(xy);
    Tensor3x3 eta(0);

    for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    for(int I=0; I<3; I++)
    for(int J=0; J<3; J++)
        eta[{i,j}] += tetrad[{I,i}] * tetrad[{J,j}] * g_ll[{I,J}];

    g_ll.Print("  g_ll",1);
    tetrad.Print("tetrad",1);
    cout << "eta recreated from tetrad^T g tetrad:" << endl;
    eta.Print("   eta",1);
}



void PhotonVelocity()
{
    size_t nx, ny;
    nx = ny = 200;
    Coord start(1,1);
    Coord end(3,3);
    Grid grid(nx, ny, start, end);
    // SchwarzSchild metric(grid, 1.0, 0.0);
    KerrSchild metric(grid, 1.0, 0.0);
    
    Coord xy(1.5, 2);
    double alpha = metric.GetAlpha(xy);
    Tensor3x3 g_ll = metric.GetMetric_ll(xy);
    Tensor2x2 eta2x2_ll = metric.GetMinkowskiGamma_ll(xy);
    Tensor3x3 eta3x3_ll = metric.GetMinkowskiMetric_ll(xy);
    Tensor2x2 gamma_ll = metric.GetGamma_ll(xy);

    {
        cout << "Four velocity defined in LF and transformed to three velocity in IF" << endl;
        Tensor3 uLF(1,0.2,0.4);
        uLF = NullNormalize(uLF,g_ll);
        Tensor2 vIF = Vec2ObservedByEulObs<LF,IF>(uLF,xy,metric);

        uLF.Print(" uLF ");
        PrintDouble(Norm2(uLF,g_ll),"|uLF|");
        vIF.Print(" vIF ");
        PrintDouble(Norm2(vIF,eta2x2_ll),"|vIF|");
    }cout << endl;
    {
        cout << "Four velocity defined in LF and transformed to three velocity in LF" << endl;
        Tensor3 uLF(1,0.2,0.4);
        uLF = NullNormalize(uLF,g_ll);
        Tensor2 vLF = Vec2ObservedByEulObs<LF,LF>(uLF, xy, metric);

        uLF.Print(" uLF ");
        PrintDouble(Norm2(uLF,g_ll),"|uLF|");
        vLF.Print(" vLF ");
        PrintDouble(Norm2(vLF,gamma_ll),"|vLF|");
    }cout << endl;
    {
        cout << "Four velocity defined in IF and transformed to three velocity in LF:" << endl;
        Tensor3 uIF(1*alpha,0.2*alpha,0.4*alpha);
        uIF = NullNormalize(uIF,eta3x3_ll);
        Tensor2 vLF = Vec2ObservedByEulObs<IF,LF>(uIF, xy, metric);

        uIF.Print(" uIF ");
        PrintDouble(Norm2(uIF,eta3x3_ll),"|uIF|");
        vLF.Print(" vLF ");
        PrintDouble(Norm2(vLF,gamma_ll),"|vLF|");
    }cout << endl;
    {
        cout << "Four velocity defined in IF and transformed to three velocity in IF:" << endl;
        Tensor3 uIF(1*alpha,0.2*alpha,0.4*alpha);
        uIF = NullNormalize(uIF,eta3x3_ll);
        Tensor2 vIF = Vec2ObservedByEulObs<IF,IF>(uIF, xy, metric);

        uIF.Print(" uIF ");
        PrintDouble(Norm2(uIF,eta3x3_ll),"|uIF|");
        vIF.Print(" vIF ");
        PrintDouble(Norm2(vIF,eta2x2_ll),"|vIF|");
    }cout << endl;
}



int main()
{
    // Tetrad();
    PhotonVelocity();
}