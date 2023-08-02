#include <iostream>
#include "../src/Grid.h"
using namespace std;

int main()
{
    size_t nx = 20;
    size_t ny = 20;
    Coord start(0, 0);
    Coord end(1, 1);
    Grid grid(nx, ny, start, end);

    auto dist = [](Coord a, Coord b)
    {
        double x = (a[1] - b[1]);
        double y = (a[2] - b[2]);
        return sqrt(x * x + y * y);
    };

    Coord rOrb(0.3, 0.3);
    Coord gOrb(0.3, 0.7);
    Coord bOrb(0.7, 0.3);
    Coord aOrb(0.5, 0.5);
    RealBuffer data0(grid.nxy);
    RealBuffer data1(grid.nxy);
    RealBuffer data2(grid.nxy);
    RealBuffer data3(grid.nxy);

    for (size_t j = 0; j < grid.ny; j++)
        for (size_t i = 0; i < grid.nx; i++)
        {
            size_t ij = grid.Index(i, j);
            Coord xy = grid.xy(i, j);
            data0[ij] = dist(rOrb, xy);
            data1[ij] = dist(gOrb, xy);
            data2[ij] = dist(bOrb, xy);
            data3[ij] = dist(aOrb, xy);
            ;
        }
    grid.WriteFrametoCsv(0, data0, data1, data2, data3, OUTPUTDIR, "GridTest");

    cout << "output/GridTest... has been created. Plot it with ParaView (Filter:Table to Points)." << endl;
}