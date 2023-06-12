#include "Spacetimes.h"

// ------------------------------ Minkowski ------------------------------
Minkowski::Minkowski(Grid& grid_, double m_, double a_) : Metric(grid_, m_, a_)
{
    this->InitializeMetricOnGrid();
    this->InitializeMetricDerivativesOnGrid();
    this->InitializeAdmComponentsOnGrid();
    this->InitializeBoostedTetradOnGrid();
}

bool Minkowski::InsideBH(const Coord& xy)
{
    return false;
}
std::string Minkowski::Name()
{
    return "Minkowski";
}

Tensor3x3 Minkowski::MetricFunction(const Coord& xy)
{
    return Tensor3x3(-1,0,0, 0,1,0, 0,0,1);
}
// -----------------------------------------------------------------------



// ------------------------------ SchwarzSchild ------------------------------
SchwarzSchild::SchwarzSchild(Grid& grid_, double m_, double a_) : Metric(grid_, m_, a_)
{
    this->InitializeMetricOnGrid();
    this->InitializeMetricDerivativesOnGrid();
    this->InitializeAdmComponentsOnGrid();
    this->InitializeBoostedTetradOnGrid();
}

bool SchwarzSchild::InsideBH(const Coord& xy)
{
    // Buffer zone must be bigger here or geodesic equations dont converge.
    return xy.EuklNormSquared() <= 2.1 * 2.1 * this->m * this->m;
}
std::string SchwarzSchild::Name()
{
    return "SchwarzSchild";
}

Tensor3x3 SchwarzSchild::MetricFunction(const Coord& xy)
{
    double rs = 2.0 * this->m;
    double r2 = xy.EuklNormSquared();
    double r = sqrt(r2);
    if(r > 2.0*this->m)
    {
        Tensor3x3 g_ll(0.0);
        g_ll[{0,0}] = -1.0 + rs / r;
        g_ll[{1,1}] = xy[1]*xy[1]/(r2-rs*r) + xy[2]*xy[2]/r2;
        g_ll[{2,2}] = xy[1]*xy[1]/r2 + xy[2]*xy[2]/(r2-rs*r);
        g_ll[{2,1}] = g_ll[{1,2}] = rs*xy[1]*xy[2] / (r2*(r-rs));
        return g_ll;
    }
    else
        return Tensor3x3(-1,0,0, 0,1,0, 0,0,1);
}
// ---------------------------------------------------------------------------



// ------------------------------ KerrSchild ------------------------------
KerrSchild::KerrSchild(Grid& grid_, double m_, double a_) : Metric(grid_, m_, a_)
{
    this->InitializeMetricOnGrid();
    this->InitializeMetricDerivativesOnGrid();
    this->InitializeAdmComponentsOnGrid();
    this->InitializeBoostedTetradOnGrid();
}

bool KerrSchild::InsideBH(const Coord& xy)
{
    return xy.EuklNormSquared() <= 2.1 * 2.1 * this->m * this->m;
}
std::string KerrSchild::Name()
{
    return "KerrSchild";
}

Tensor3x3 KerrSchild::MetricFunction(const Coord& xy)
{
    double a2 = this->a * this->a;
    double R2 = xy.EuklNormSquared();
    if(R2 > a2)
    {
        double r2 = R2 - a2;
        double r  = sqrt(r2);
        double rho2 = r2 + a2;
        double H = this->m / r;
        Tensor3 l_l( 1, (r * xy[1] + this->a * xy[2]) / rho2, (r * xy[2] - this->a * xy[1]) / rho2);
        Tensor3 l_u(-1, (r * xy[1] + this->a * xy[2]) / rho2, (r * xy[2] - this->a * xy[1]) / rho2);
        Tensor3x3 g_ll;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                if (i == 0 && j == 0)
                    g_ll[{i,j}] = -1.0 + 2.0 * H * l_l[i] * l_l[j];
                else
                    g_ll[{i,j}] = (i==j) + 2.0 * H * l_l[i] * l_l[j];
            }
        return g_ll;
    }
    else
        return Tensor3x3(-1,0,0, 0,1,0, 0,0,1);
}
// ------------------------------------------------------------------------