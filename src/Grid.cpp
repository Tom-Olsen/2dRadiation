#include "Grid.h"

// Constructors:
Grid::Grid(size_t nx_, size_t ny_, Coord start_, Coord end_)
{
    if (end_[1] < start_[1])
        ExitOnError("Grid x€[a,b] has b<a!");
    if (end_[2] < start_[2])
        ExitOnError("Grid y€[a,b] has b<a!");
    if (nx_ <= 1 or ny_ <= 1)
        ExitOnError("Grid must have at least 2 LP in each Dimension.");

    // Add 2 ghost cells and 1 extra cell due to off by one quirk.
    nx = nx_ + 1 + 2;
    ny = ny_ + 1 + 2;
    nxy = nx * ny;
    dx = (end_[1] - start_[1]) / (nx - 1.0 - 2.0);
    dy = (end_[2] - start_[2]) / (ny - 1.0 - 2.0);
    dt = m_cfl * std::min(dx, dy);
    startx = start_[1] - dx;
    starty = start_[2] - dy;
    endx = end_[1] + dx;
    endy = end_[2] + dy;
}
Grid::Grid(const Grid &grid) :
nx(grid.nx), ny(grid.ny), nxy(grid.nxy),
dx(grid.dx), dy(grid.dy), dt(grid.dt),
startx(grid.startx), starty(grid.starty),
endx(grid.endx), endy(grid.endy),
m_cfl(grid.m_cfl) {}

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
double Grid::x(size_t i)
{
    return startx + i * dx;
}
double Grid::y(size_t j)
{
    return starty + j * dy;
}
double Grid::x(double i)
{
    return startx + i * dx;
}
double Grid::y(double j)
{
    return starty + j * dy;
}
Coord Grid::xy(size_t i, size_t j)
{
    return Coord(x(i), y(j));
}
Coord Grid::xy(double i, double j)
{
    return Coord(x(i), y(j));
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

// Debugging:
void Grid::Print()
{
    std::cout << "Grid: nx=" << nx << "\n";
    std::cout << "Grid: ny=" << ny << "\n";
    std::cout << "Grid: nxy=" << nxy << "\n";
    std::cout << "Grid: dx=" << dx << "\n";
    std::cout << "Grid: dy=" << dy << "\n";
    std::cout << "Grid: dt=" << dt << "\n";
    std::cout << "Grid: startx=" << startx << "\n";
    std::cout << "Grid: starty=" << starty << "\n";
    std::cout << "Grid: endx=" << endx << "\n";
    std::cout << "Grid: endy=" << endy << "\n";
    std::cout << "Grid: cfl=" << m_cfl << "\n";
}