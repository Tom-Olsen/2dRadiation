#include <iostream>
#include "../src/GeodesicEquationSolver.h"
#include "../src/SpecialMath.h"
using namespace std;



int main()
{
    std::ofstream fileOut((std::string)OUTPUTDIR + (std::string)"Test_GeodesicEquationSolver.txt");
    fileOut << "#x, y, z, s \n";

    size_t nx, ny;
    nx = ny = 200;
    Coord start(-6,-6);
    Coord end(6,6);
    Grid grid(nx, ny, start, end);
    // Minkowski metric(grid, 1.0, 0.0);
    // KerrSchild metric(grid, 1.0, 0.0);
    SchwarzSchild metric(grid, 1.0, 0.0);
    cout << "Initialization complete." << endl;

    int n = 30;
    for(int i=0; i<n; i++)
    {
        double s = 1.0;
        Coord x(start[1] + (i+0.5) * (end[1] - start[1]) / n, start[2]);
        Tensor3x3 g_ll = metric.GetMetric_ll(x);
        Tensor2x2 gamma_ll = metric.GetGamma_ll(x);
        Tensor2x2 delta_ll = metric.GetMinkowskiGamma_ll(x);

        Tensor3 uLF(1,0,1);
        uLF = NullNormalize(uLF, g_ll);
        Tensor2 vLF = Vec2ObservedByEulObs<LF,LF>(uLF,x,metric);

        fileOut << x[1] << ", " << x[2] << ", " << 0 << ", " << s << "\n";
        while(true)
        {
            if (grid.OutsideDomain(x) || metric.InsideBH(x))
                break;
            s *= RK45_GeodesicEquation<1>(5 * grid.dt, x, vLF, metric);
            fileOut << x[1] << ", " << x[2] << ", " << 0 << ", " << s << "\n";
        }
        cout << "Photon(" << i << ") complete." << endl;
    }

    fileOut.close();
}