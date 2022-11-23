#include "Spacetimes.h"

// ------------------------------ Minkowski ------------------------------
template<class Coord>
Minkowski<Coord>::Minkowski(Grid2D<Coord>& grid_, double m_, double a_) : Metric2D<Coord>(grid_, m_, a_)
{
    this->InitializeMetricOnGrid();
    this->InitializeMetricDerivativesOnGrid();
    this->InitializeAdmComponentsOnGrid();
    this->InitializeBoostedTetradOnGrid();
}

template<class Coord>
bool Minkowski<Coord>::InsideBH(const int i, const int j)
{
    return false;
}
template<class Coord>
bool Minkowski<Coord>::InsideBH(const Coordinate2<Coord>& x12)
{
    return false;
}
template<class Coord>
std::string Minkowski<Coord>::Name()
{
    return "Minkowski";
}

template<class Coord>
Tensor3x3<Coord,LF> Minkowski<Coord>::MetricFunction(Coordinate2<Coord> x)
{
    if constexpr(std::is_same<Coord,xy>::value)
        return Tensor3x3<Coord,LF>(-1,0,0 ,0,1,0, 0,0,1);

    else //if  constexpr(std::is_same<Coord,xy>::value)
        return Tensor3x3<Coord,LF>(-1,0,0 ,0,1,0, 0,0,x[1]*x[1]);
}
// -----------------------------------------------------------------------



// ------------------------------ SchwarzSchild ------------------------------
template<class Coord>
SchwarzSchild<Coord>::SchwarzSchild(Grid2D<Coord>& grid_, double m_, double a_) : Metric2D<Coord>(grid_, m_, a_)
{
    this->InitializeMetricOnGrid();
    this->InitializeMetricDerivativesOnGrid();
    this->InitializeAdmComponentsOnGrid();
    this->InitializeBoostedTetradOnGrid();
}

template<class Coord>
bool SchwarzSchild<Coord>::InsideBH(const int i, const int j)
{
    return this->grid.rCoord(i,j) <= 2.1*this->m;
}
template<class Coord>
bool SchwarzSchild<Coord>::InsideBH(const Coordinate2<Coord>& x12)
{
    if constexpr(std::is_same<Coord,xy>::value)
        return x12[1]*x12[1] + x12[2]*x12[2] <= 2.1*2.1*this->m*this->m;
    if constexpr(std::is_same<Coord,rph>::value)
        return x12[1] <= 2.1*this->m;
}
template<class Coord>
std::string SchwarzSchild<Coord>::Name()
{
    return "SchwarzSchild";
}

template<class Coord>
Tensor3x3<Coord,LF> SchwarzSchild<Coord>::MetricFunction(Coordinate2<Coord> x)
{
    if constexpr(std::is_same<Coord,xy>::value)
    {
        double rs = 2.0 * this->m;
        double r2 = x[1]*x[1] + x[2]*x[2];
        double r = sqrt(r2);
        if(r > 2.0*this->m)
        {
            Tensor3x3<Coord,LF> g_ll(0.0);
            g_ll[{0,0}] = -1.0 + rs / r;
            g_ll[{1,1}] = x[1]*x[1]/(r2-rs*r) + x[2]*x[2]/r2;
            g_ll[{2,1}] = g_ll[{1,2}] = rs*x[1]*x[2] / (r2*(r-rs));
            g_ll[{2,2}] = x[1]*x[1]/r2 + x[2]*x[2]/(r2-rs*r);
            return g_ll;
        }
        else
            return Tensor3x3<Coord,LF>(-1,0,0 ,0,1,0, 0,0,1);
    }
    
    else //if constexpr(std::is_same<Coord,xy>::value)
    {
        if(x[1] > 2.0*this->m)
        {
            Tensor3x3<Coord,LF> g_ll(0.0);
            g_ll[{0,0}] = -(1.0 - 2.0*this->m/x[1]);
            g_ll[{1,1}] = 1.0/(1.0 - 2.0*this->m/x[1]);
            g_ll[{2,2}] = x[1]*x[1];
            return g_ll;
        }
        else
            return Tensor3x3<Coord,LF>(-1,0,0 ,0,1,0, 0,0,x[1]*x[1]);
    }
}
// ---------------------------------------------------------------------------



// ------------------------------ KerrSchild ------------------------------
template<class Coord>
KerrSchild<Coord>::KerrSchild(Grid2D<Coord>& grid_, double m_, double a_) : Metric2D<Coord>(grid_, m_, a_)
{
    this->InitializeMetricOnGrid();
    this->InitializeMetricDerivativesOnGrid();
    this->InitializeAdmComponentsOnGrid();
    this->InitializeBoostedTetradOnGrid();
}

template<class Coord>
bool KerrSchild<Coord>::InsideBH(const int i, const int j)
{
    return this->grid.rCoord(i,j) <= 2.1*this->m;
}
template<class Coord>
bool KerrSchild<Coord>::InsideBH(const Coordinate2<Coord>& x12)
{
    if constexpr(std::is_same<Coord,xy>::value)
        return x12[1]*x12[1] + x12[2]*x12[2] <= 2.1*2.1*this->m*this->m;
    if constexpr(std::is_same<Coord,rph>::value)
        return x12[1] <= 2.1*this->m;
}
template<class Coord>
std::string KerrSchild<Coord>::Name()
{
    return "KerrSchild";
}

template<class Coord>
Tensor3x3<Coord,LF> KerrSchild<Coord>::MetricFunction(Coordinate2<Coord> x)
{
    if constexpr(std::is_same<Coord,xy>::value)
    {
        double a2 = this->a * this->a;
        double R2 = x[1]*x[1] + x[2]*x[2];
        if(R2 > a2)
        {
            double r2 = R2 - a2;
            double r  = sqrt(r2);
            double rho2 = r2 + a2;
            double H = this->m / r;
            Tensor3<Coord,LF> l_l( 1, (r * x[1] + this->a * x[2]) / rho2, (r * x[2] - this->a * x[1]) / rho2);
            Tensor3<Coord,LF> l_u(-1, (r * x[1] + this->a * x[2]) / rho2, (r * x[2] - this->a * x[1]) / rho2);
            Tensor3x3<Coord,LF> g_ll;
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
            return Tensor3x3<Coord,LF>(-1,0,0 ,0,1,0, 0,0,1);
    }
    
    else //if constexpr(std::is_same<Coord,xy>::value)
    {
        exit_on_error("KerrSchild_rph not implemented yet!");
            return Tensor3x3<Coord,LF>(-1,0,0 ,0,1,0, 0,0,1);
    }
}
// ------------------------------------------------------------------------



template class Minkowski<xy>;
template class Minkowski<rph>;
template class SchwarzSchild<xy>;
template class SchwarzSchild<rph>;
template class KerrSchild<xy>;
template class KerrSchild<rph>;