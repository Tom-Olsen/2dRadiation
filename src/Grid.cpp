#include "Grid.h"

// Constructors:
Grid::Grid(size_t nx, size_t ny, Coord start, Coord end) : nx(nx), ny(ny), startx(start[1]), starty(start[2]), endx(end[1]), endy(end[2])
{
    if (end[1] < start[1])
        ExitOnError("Grid x€[a,b] has b<a!");
    if (end[2] < start[2])
        ExitOnError("Grid y€[a,b] has b<a!");
    if (nx == 1 or ny == 1)
        ExitOnError("Grid must have at least 2 LP in each Dimension.");

    dx = (endx - startx) / (nx - 1.0);
    dy = (endy - starty) / (ny - 1.0);
    dt = m_cfl * std::min(dx, dy);
    nxy = nx * ny;
}
Grid::Grid(const Grid &grid) : nx(grid.nx), ny(grid.ny), startx(grid.startx), starty(grid.starty), endx(grid.endx), endy(grid.endy)
{
    dx = (endx - startx) / (nx - 1);
    dy = (endy - starty) / (ny - 1);
    dt = m_cfl * std::min(dx, dy);
    nxy = nx * ny;
}

// Setters/Getters:
void Grid::SetCFL(double cfl)
{
    m_cfl = cfl;
    dt = m_cfl * std::min(dx, dy);
}
double Grid::GetCFL()
{
    return m_cfl;
}

// Grid Access Tools:
size_t Grid::Index(size_t i, size_t j)
{
    return i + j * nx;
}
Coord Grid::xy(size_t i, size_t j)
{
    return Coord(startx + i * dx, starty + j * dy);
}
Coord Grid::xy(double i, double j)
{
    return Coord(startx + i * dx, starty + j * dy);
}
Coord Grid::xy(size_t ij)
{
    return xy(i(ij), j(ij));
}
double Grid::i(double x)
{
    return (x - startx) / dx;
}
double Grid::j(double y)
{
    return (y - starty) / dy;
}
size_t Grid::i(size_t ij)
{
    return ij % nx;
}
size_t Grid::j(size_t ij)
{
    return ij / nx;
}
Coord Grid::ij(const Coord &xy)
{
    return Coord(i(xy[1]), j(xy[2]));
}

// Domain Checks:
bool Grid::OutsideDomain(const Coord &xyz)
{
    return xyz[1] > endx or xyz[1] < startx or
           xyz[2] > endy or xyz[2] < starty;
}
bool Grid::OutsideDomain(double i, double j)
{
    return i > nx - 1 or i < 0 or
           j > ny - 1 or j < 0;
}

// Write Data to file:
void Grid::WriteFrametoCsv(float time, const RealBuffer &r, const RealBuffer &g, const RealBuffer &b, const RealBuffer &a, std::string directory, std::string name)
{
    PROFILE_FUNCTION();
    CreateDirectory(directory);

    name = (name == "") ? "data" : name;
    std::ofstream fileOut(directory + name + Format(time, 3) + "t " + std::to_string(nx) + "x " + std::to_string(ny) + "y" + ".csv");

    fileOut << "#nx=" << nx << "\n";
    fileOut << "#ny=" << ny << "\n";
    fileOut << "#startx=" << startx << "\n";
    fileOut << "#starty=" << starty << "\n";
    fileOut << "#endx=" << endx << "\n";
    fileOut << "#endy=" << endy << "\n";
    fileOut << "#t=" << time << "\n";
    fileOut << "#x,y,r,g,b,a\n";

    for (size_t j = 0; j < ny; j++)
        for (size_t i = 0; i < nx; i++)
        {
            size_t ij = Index(i, j);
            Coord x = xy(i, j);
            fileOut << x[1] << "," << x[2] << ",";
            fileOut << r[ij] << "," << g[ij] << "," << b[ij] << "," << a[ij] << "\n";
        }

    fileOut.close();
}