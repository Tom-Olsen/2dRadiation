#include "Metric2D.h"



// Constructor/Desctructor
template<class Coord>
Metric2D<Coord>::Metric2D(Grid2D<Coord>& grid_, double m_, double a_) : grid(grid_), m(m_), a(a_)
{
    g00_ll = new double[grid.n12]();    g00_uu = new double[grid.n12]();
    g01_ll = new double[grid.n12]();    g01_uu = new double[grid.n12]();
    g02_ll = new double[grid.n12]();    g02_uu = new double[grid.n12]();
    g11_ll = new double[grid.n12]();    g11_uu = new double[grid.n12]();
    g12_ll = new double[grid.n12]();    g12_uu = new double[grid.n12]();
    g22_ll = new double[grid.n12]();    g22_uu = new double[grid.n12]();
    d0_g00_lll = new double[grid.n12]();    d1_g00_lll = new double[grid.n12]();    d2_g00_lll = new double[grid.n12]();
    d0_g01_lll = new double[grid.n12]();    d1_g01_lll = new double[grid.n12]();    d2_g01_lll = new double[grid.n12]();
    d0_g02_lll = new double[grid.n12]();    d1_g02_lll = new double[grid.n12]();    d2_g02_lll = new double[grid.n12]();
    d0_g11_lll = new double[grid.n12]();    d1_g11_lll = new double[grid.n12]();    d2_g11_lll = new double[grid.n12]();
    d0_g12_lll = new double[grid.n12]();    d1_g12_lll = new double[grid.n12]();    d2_g12_lll = new double[grid.n12]();
    d0_g22_lll = new double[grid.n12]();    d1_g22_lll = new double[grid.n12]();    d2_g22_lll = new double[grid.n12]();
    d0_g00_luu = new double[grid.n12]();    d1_g00_luu = new double[grid.n12]();    d2_g00_luu = new double[grid.n12]();
    d0_g01_luu = new double[grid.n12]();    d1_g01_luu = new double[grid.n12]();    d2_g01_luu = new double[grid.n12]();
    d0_g02_luu = new double[grid.n12]();    d1_g02_luu = new double[grid.n12]();    d2_g02_luu = new double[grid.n12]();
    d0_g11_luu = new double[grid.n12]();    d1_g11_luu = new double[grid.n12]();    d2_g11_luu = new double[grid.n12]();
    d0_g12_luu = new double[grid.n12]();    d1_g12_luu = new double[grid.n12]();    d2_g12_luu = new double[grid.n12]();
    d0_g22_luu = new double[grid.n12]();    d1_g22_luu = new double[grid.n12]();    d2_g22_luu = new double[grid.n12]();
    alpha = new double[grid.n12]();
    beta1_u = new double[grid.n12]();    beta1_l = new double[grid.n12]();
    beta2_u = new double[grid.n12]();    beta2_l = new double[grid.n12]();
    gamma11_ll = new double[grid.n12]();    gamma11_uu = new double[grid.n12]();
    gamma12_ll = new double[grid.n12]();    gamma12_uu = new double[grid.n12]();
    gamma22_ll = new double[grid.n12]();    gamma22_uu = new double[grid.n12]();
    d1_alpha_l = new double[grid.n12]();
    d2_alpha_l = new double[grid.n12]();
    d1_beta1_lu = new double[grid.n12]();    d1_beta2_lu = new double[grid.n12]();
    d2_beta1_lu = new double[grid.n12]();    d2_beta2_lu = new double[grid.n12]();
    d1_beta1_ll = new double[grid.n12]();    d1_beta2_ll = new double[grid.n12]();
    d2_beta1_ll = new double[grid.n12]();    d2_beta2_ll = new double[grid.n12]();
    d1_gamma11_lll = new double[grid.n12]();    d2_gamma11_lll = new double[grid.n12]();
    d1_gamma12_lll = new double[grid.n12]();    d2_gamma12_lll = new double[grid.n12]();
    d1_gamma22_lll = new double[grid.n12]();    d2_gamma22_lll = new double[grid.n12]();
    d1_gamma11_luu = new double[grid.n12]();    d2_gamma11_luu = new double[grid.n12]();
    d1_gamma12_luu = new double[grid.n12]();    d2_gamma12_luu = new double[grid.n12]();
    d1_gamma22_luu = new double[grid.n12]();    d2_gamma22_luu = new double[grid.n12]();
    K11_ll = new double[grid.n12];
    K12_ll = new double[grid.n12];
    K22_ll = new double[grid.n12];
    tetrad00_ul = new double[grid.n12]();    tetrad01_ul = new double[grid.n12]();    tetrad02_ul = new double[grid.n12]();
    tetrad10_ul = new double[grid.n12]();    tetrad11_ul = new double[grid.n12]();    tetrad12_ul = new double[grid.n12]();
    tetrad20_ul = new double[grid.n12]();    tetrad21_ul = new double[grid.n12]();    tetrad22_ul = new double[grid.n12]();
}
//template<class Coord>
//Metric2D<Coord>::Metric2D(const Metric2D& metric) : grid(metric.grid), m(metric.a), a(metric.a)
//{
//    g00_ll = new double[grid.n12]();    g00_uu = new double[grid.n12]();
//    g01_ll = new double[grid.n12]();    g01_uu = new double[grid.n12]();
//    g02_ll = new double[grid.n12]();    g02_uu = new double[grid.n12]();
//    g11_ll = new double[grid.n12]();    g11_uu = new double[grid.n12]();
//    g12_ll = new double[grid.n12]();    g12_uu = new double[grid.n12]();
//    g22_ll = new double[grid.n12]();    g22_uu = new double[grid.n12]();
//    d0_g00_lll = new double[grid.n12]();    d1_g00_lll = new double[grid.n12]();    d2_g00_lll = new double[grid.n12]();
//    d0_g01_lll = new double[grid.n12]();    d1_g01_lll = new double[grid.n12]();    d2_g01_lll = new double[grid.n12]();
//    d0_g02_lll = new double[grid.n12]();    d1_g02_lll = new double[grid.n12]();    d2_g02_lll = new double[grid.n12]();
//    d0_g11_lll = new double[grid.n12]();    d1_g11_lll = new double[grid.n12]();    d2_g11_lll = new double[grid.n12]();
//    d0_g12_lll = new double[grid.n12]();    d1_g12_lll = new double[grid.n12]();    d2_g12_lll = new double[grid.n12]();
//    d0_g22_lll = new double[grid.n12]();    d1_g22_lll = new double[grid.n12]();    d2_g22_lll = new double[grid.n12]();
//    d0_g00_luu = new double[grid.n12]();    d1_g00_luu = new double[grid.n12]();    d2_g00_luu = new double[grid.n12]();
//    d0_g01_luu = new double[grid.n12]();    d1_g01_luu = new double[grid.n12]();    d2_g01_luu = new double[grid.n12]();
//    d0_g02_luu = new double[grid.n12]();    d1_g02_luu = new double[grid.n12]();    d2_g02_luu = new double[grid.n12]();
//    d0_g11_luu = new double[grid.n12]();    d1_g11_luu = new double[grid.n12]();    d2_g11_luu = new double[grid.n12]();
//    d0_g12_luu = new double[grid.n12]();    d1_g12_luu = new double[grid.n12]();    d2_g12_luu = new double[grid.n12]();
//    d0_g22_luu = new double[grid.n12]();    d1_g22_luu = new double[grid.n12]();    d2_g22_luu = new double[grid.n12]();
//    alpha = new double[grid.n12]();
//    beta1_u = new double[grid.n12]();    beta1_l = new double[grid.n12]();
//    beta2_u = new double[grid.n12]();    beta2_l = new double[grid.n12]();
//    gamma11_ll = new double[grid.n12]();    gamma11_uu = new double[grid.n12]();
//    gamma12_ll = new double[grid.n12]();    gamma12_uu = new double[grid.n12]();
//    gamma22_ll = new double[grid.n12]();    gamma22_uu = new double[grid.n12]();
//    d1_alpha_l = new double[grid.n12]();
//    d2_alpha_l = new double[grid.n12]();
//    d1_beta1_lu = new double[grid.n12]();    d1_beta2_lu = new double[grid.n12]();
//    d2_beta1_lu = new double[grid.n12]();    d2_beta2_lu = new double[grid.n12]();
//    d1_beta1_ll = new double[grid.n12]();    d1_beta2_ll = new double[grid.n12]();
//    d2_beta1_ll = new double[grid.n12]();    d2_beta2_ll = new double[grid.n12]();
//    d1_gamma11_lll = new double[grid.n12]();    d2_gamma11_lll = new double[grid.n12]();
//    d1_gamma12_lll = new double[grid.n12]();    d2_gamma12_lll = new double[grid.n12]();
//    d1_gamma22_lll = new double[grid.n12]();    d2_gamma22_lll = new double[grid.n12]();
//    d1_gamma11_luu = new double[grid.n12]();    d2_gamma11_luu = new double[grid.n12]();
//    d1_gamma12_luu = new double[grid.n12]();    d2_gamma12_luu = new double[grid.n12]();
//    d1_gamma22_luu = new double[grid.n12]();    d2_gamma22_luu = new double[grid.n12]();
//    K11_ll = new double[grid.n12];
//    K12_ll = new double[grid.n12];
//    K22_ll = new double[grid.n12];
//    tetrad00_ul = new double[grid.n12]();    tetrad01_ul = new double[grid.n12]();    tetrad02_ul = new double[grid.n12]();
//    tetrad10_ul = new double[grid.n12]();    tetrad11_ul = new double[grid.n12]();    tetrad12_ul = new double[grid.n12]();
//    tetrad20_ul = new double[grid.n12]();    tetrad21_ul = new double[grid.n12]();    tetrad22_ul = new double[grid.n12]();
//
//    for(int ij=0; ij<grid.n12; ij++)
//    {
//        g00_ll[ij] = metric.g00_ll[ij];    g00_uu[ij] = metric.g00_uu[ij];
//        g01_ll[ij] = metric.g01_ll[ij];    g01_uu[ij] = metric.g01_uu[ij];
//        g02_ll[ij] = metric.g02_ll[ij];    g02_uu[ij] = metric.g02_uu[ij];
//        g11_ll[ij] = metric.g11_ll[ij];    g11_uu[ij] = metric.g11_uu[ij];
//        g12_ll[ij] = metric.g12_ll[ij];    g12_uu[ij] = metric.g12_uu[ij];
//        g22_ll[ij] = metric.g22_ll[ij];    g22_uu[ij] = metric.g22_uu[ij];
//        d0_g00_lll[ij] = metric.d0_g00_lll[ij];    d1_g00_lll[ij] = metric.d1_g00_lll[ij];    d2_g00_lll[ij] = metric.d2_g00_lll[ij];
//        d0_g01_lll[ij] = metric.d0_g01_lll[ij];    d1_g01_lll[ij] = metric.d1_g01_lll[ij];    d2_g01_lll[ij] = metric.d2_g01_lll[ij];
//        d0_g02_lll[ij] = metric.d0_g02_lll[ij];    d1_g02_lll[ij] = metric.d1_g02_lll[ij];    d2_g02_lll[ij] = metric.d2_g02_lll[ij];
//        d0_g11_lll[ij] = metric.d0_g11_lll[ij];    d1_g11_lll[ij] = metric.d1_g11_lll[ij];    d2_g11_lll[ij] = metric.d2_g11_lll[ij];
//        d0_g12_lll[ij] = metric.d0_g12_lll[ij];    d1_g12_lll[ij] = metric.d1_g12_lll[ij];    d2_g12_lll[ij] = metric.d2_g12_lll[ij];
//        d0_g22_lll[ij] = metric.d0_g22_lll[ij];    d1_g22_lll[ij] = metric.d1_g22_lll[ij];    d2_g22_lll[ij] = metric.d2_g22_lll[ij];
//        d0_g00_luu[ij] = metric.d0_g00_luu[ij];    d1_g00_luu[ij] = metric.d1_g00_luu[ij];    d2_g00_luu[ij] = metric.d2_g00_luu[ij];
//        d0_g01_luu[ij] = metric.d0_g01_luu[ij];    d1_g01_luu[ij] = metric.d1_g01_luu[ij];    d2_g01_luu[ij] = metric.d2_g01_luu[ij];
//        d0_g02_luu[ij] = metric.d0_g02_luu[ij];    d1_g02_luu[ij] = metric.d1_g02_luu[ij];    d2_g02_luu[ij] = metric.d2_g02_luu[ij];
//        d0_g11_luu[ij] = metric.d0_g11_luu[ij];    d1_g11_luu[ij] = metric.d1_g11_luu[ij];    d2_g11_luu[ij] = metric.d2_g11_luu[ij];
//        d0_g12_luu[ij] = metric.d0_g12_luu[ij];    d1_g12_luu[ij] = metric.d1_g12_luu[ij];    d2_g12_luu[ij] = metric.d2_g12_luu[ij];
//        d0_g22_luu[ij] = metric.d0_g22_luu[ij];    d1_g22_luu[ij] = metric.d1_g22_luu[ij];    d2_g22_luu[ij] = metric.d2_g22_luu[ij];
//        alpha[ij] = metric.alpha[ij];
//        beta1_u[ij] = metric.beta1_u[ij];    beta1_l[ij] = metric.beta1_l[ij];
//        beta2_u[ij] = metric.beta2_u[ij];    beta2_l[ij] = metric.beta2_l[ij];
//        gamma11_ll[ij] = metric.gamma11_ll[ij];    gamma11_uu[ij] = metric.gamma11_uu[ij];
//        gamma12_ll[ij] = metric.gamma12_ll[ij];    gamma12_uu[ij] = metric.gamma12_uu[ij];
//        gamma22_ll[ij] = metric.gamma22_ll[ij];    gamma22_uu[ij] = metric.gamma22_uu[ij];
//        d1_alpha_l[ij] = metric.d1_alpha_l[ij];
//        d2_alpha_l[ij] = metric.d2_alpha_l[ij];
//        d1_beta1_lu[ij] = metric.d1_beta1_lu[ij];    d1_beta2_lu[ij] = metric.d1_beta2_lu[ij];
//        d2_beta1_lu[ij] = metric.d2_beta1_lu[ij];    d2_beta2_lu[ij] = metric.d2_beta2_lu[ij];
//        d1_beta1_ll[ij] = metric.d1_beta1_ll[ij];    d1_beta2_ll[ij] = metric.d1_beta2_ll[ij];
//        d2_beta1_ll[ij] = metric.d2_beta1_ll[ij];    d2_beta2_ll[ij] = metric.d2_beta2_ll[ij];
//        d1_gamma11_lll[ij] = metric.d1_gamma11_lll[ij];    d2_gamma11_lll[ij] = metric.d2_gamma11_lll[ij];
//        d1_gamma12_lll[ij] = metric.d1_gamma12_lll[ij];    d2_gamma12_lll[ij] = metric.d2_gamma12_lll[ij];
//        d1_gamma22_lll[ij] = metric.d1_gamma22_lll[ij];    d2_gamma22_lll[ij] = metric.d2_gamma22_lll[ij];
//        d1_gamma11_luu[ij] = metric.d1_gamma11_luu[ij];    d2_gamma11_luu[ij] = metric.d2_gamma11_luu[ij];
//        d1_gamma12_luu[ij] = metric.d1_gamma12_luu[ij];    d2_gamma12_luu[ij] = metric.d2_gamma12_luu[ij];
//        d1_gamma22_luu[ij] = metric.d1_gamma22_luu[ij];    d2_gamma22_luu[ij] = metric.d2_gamma22_luu[ij];
//        K11_ll[ij] = metric.K11_ll[ij];
//        K12_ll[ij] = metric.K12_ll[ij];
//        K22_ll[ij] = metric.K22_ll[ij];
//        tetrad00_ul[ij] = metric.tetrad00_ul[ij];    tetrad01_ul[ij] = metric.tetrad01_ul[ij];    tetrad02_ul[ij] = metric.tetrad02_ul[ij];
//        tetrad10_ul[ij] = metric.tetrad10_ul[ij];    tetrad11_ul[ij] = metric.tetrad11_ul[ij];    tetrad12_ul[ij] = metric.tetrad12_ul[ij];
//        tetrad20_ul[ij] = metric.tetrad20_ul[ij];    tetrad21_ul[ij] = metric.tetrad21_ul[ij];    tetrad22_ul[ij] = metric.tetrad22_ul[ij];
//    }
//}
template<class Coord>
Metric2D<Coord>::~Metric2D()
{
    delete[] g00_ll;    delete[] g00_uu;
    delete[] g01_ll;    delete[] g01_uu;
    delete[] g02_ll;    delete[] g02_uu;
    delete[] g11_ll;    delete[] g11_uu;
    delete[] g12_ll;    delete[] g12_uu;
    delete[] g22_ll;    delete[] g22_uu;
    delete[] d0_g00_lll;    delete[] d1_g00_lll;    delete[] d2_g00_lll;
    delete[] d0_g01_lll;    delete[] d1_g01_lll;    delete[] d2_g01_lll;
    delete[] d0_g02_lll;    delete[] d1_g02_lll;    delete[] d2_g02_lll;
    delete[] d0_g11_lll;    delete[] d1_g11_lll;    delete[] d2_g11_lll;
    delete[] d0_g12_lll;    delete[] d1_g12_lll;    delete[] d2_g12_lll;
    delete[] d0_g22_lll;    delete[] d1_g22_lll;    delete[] d2_g22_lll;
    delete[] d0_g00_luu;    delete[] d1_g00_luu;    delete[] d2_g00_luu;
    delete[] d0_g01_luu;    delete[] d1_g01_luu;    delete[] d2_g01_luu;
    delete[] d0_g02_luu;    delete[] d1_g02_luu;    delete[] d2_g02_luu;
    delete[] d0_g11_luu;    delete[] d1_g11_luu;    delete[] d2_g11_luu;
    delete[] d0_g12_luu;    delete[] d1_g12_luu;    delete[] d2_g12_luu;
    delete[] d0_g22_luu;    delete[] d1_g22_luu;    delete[] d2_g22_luu;
    delete[] alpha;
    delete[] beta1_u;    delete[] beta1_l;
    delete[] beta2_u;    delete[] beta2_l;
    delete[] gamma11_ll;    delete[] gamma11_uu;
    delete[] gamma12_ll;    delete[] gamma12_uu;
    delete[] gamma22_ll;    delete[] gamma22_uu;
    delete[] d1_alpha_l;
    delete[] d2_alpha_l;
    delete[] d1_beta1_lu;    delete[] d1_beta2_lu;
    delete[] d2_beta1_lu;    delete[] d2_beta2_lu;
    delete[] d1_beta1_ll;    delete[] d1_beta2_ll;
    delete[] d2_beta1_ll;    delete[] d2_beta2_ll;
    delete[] d1_gamma11_lll;    delete[] d2_gamma11_lll;
    delete[] d1_gamma12_lll;    delete[] d2_gamma12_lll;
    delete[] d1_gamma22_lll;    delete[] d2_gamma22_lll;
    delete[] d1_gamma11_luu;    delete[] d2_gamma11_luu;
    delete[] d1_gamma12_luu;    delete[] d2_gamma12_luu;
    delete[] d1_gamma22_luu;    delete[] d2_gamma22_luu;
    delete[] K11_ll;
    delete[] K12_ll;
    delete[] K22_ll;
    delete[] tetrad00_ul;    delete[] tetrad10_ul;    delete[] tetrad20_ul;
    delete[] tetrad01_ul;    delete[] tetrad11_ul;    delete[] tetrad21_ul;
    delete[] tetrad02_ul;    delete[] tetrad12_ul;    delete[] tetrad22_ul;
}


template<class Coord>
std::string Metric2D<Coord>::Name()
{
    exit_on_error("Metric2D virtual Method (Name) has been called!");
    return "";
}

// Initialization:
template<class Coord>
Tensor3x3<Coord,LF> Metric2D<Coord>::MetricFunction(Coordinate2<Coord> x)
{
    exit_on_error("Metric2D virtual Method (MetricFunction) has been called!");
    return Tensor3x3<Coord,LF>(1,0,0, 0,1,0, 0,0,1);
}
template<class Coord>
void Metric2D<Coord>::InitializeMetricOnGrid()
{
    for(int j=0; j<grid.n2; j++)
    for(int i=0; i<grid.n1; i++)
    {
        Coordinate2<Coord> x = grid.x12Coord(i,j);
        Tensor3x3<Coord,LF> g_ll = MetricFunction(x);
        Tensor3x3<Coord,LF> g_uu = g_ll.Invert();
        int ij = grid.Index(i,j);

        // Kerr-Schild:
        g00_ll[ij] = g_ll[{0,0}];    g00_uu[ij] = g_uu[{0,0}];
        g01_ll[ij] = g_ll[{0,1}];    g01_uu[ij] = g_uu[{0,1}];
        g02_ll[ij] = g_ll[{0,2}];    g02_uu[ij] = g_uu[{0,2}];
        g11_ll[ij] = g_ll[{1,1}];    g11_uu[ij] = g_uu[{1,1}];
        g12_ll[ij] = g_ll[{1,2}];    g12_uu[ij] = g_uu[{1,2}];
        g22_ll[ij] = g_ll[{2,2}];    g22_uu[ij] = g_uu[{2,2}]; 
    }
}

template<class Coord>
void Metric2D<Coord>::InitializeBoostedTetradOnGrid()
{
    for(int j=0; j<grid.n2; j++)
    for(int i=0; i<grid.n1; i++)
    {
        int ij = grid.Index(i,j);
        Tensor2x2<Coord,LF> g1 = GetGamma_ll(ij);
        Tensor2x2<Coord,LF> g1Inv = GetGamma_uu(ij);
        double a = 1.0/sqrt(g1Inv[{1,1}]);
        double b = a*a*g1Inv[{1,2}];
        Tensor2x2<Coord,LF> m1(1/a,0, b/a,1);
        double g2 = m1[{2,2}]*g1[{2,2}]*m1[{2,2}];

        Tensor3 n1 = uEulObs(ij);
        if constexpr(std::is_same<Coord,xy>::value)
        {
            tetrad00_ul[ij]=n1[0];    tetrad01_ul[ij]=0;            tetrad02_ul[ij] = 0;
            tetrad10_ul[ij]=n1[1];    tetrad11_ul[ij]=m1[{1,1}];    tetrad12_ul[ij] = 0;
            tetrad20_ul[ij]=n1[2];    tetrad21_ul[ij]=m1[{2,1}];    tetrad22_ul[ij] = 1/sqrt(g2);
        }
        if constexpr(std::is_same<Coord,rph>::value)
        {
            tetrad00_ul[ij]=n1[0];    tetrad01_ul[ij]=0;            tetrad02_ul[ij] = 0;
            tetrad10_ul[ij]=n1[1];    tetrad11_ul[ij]=m1[{1,1}];    tetrad12_ul[ij] = 0;
            tetrad20_ul[ij]=n1[2];    tetrad21_ul[ij]=m1[{2,1}];    tetrad22_ul[ij] = grid.rCoord(i,j)/sqrt(g2);
        }
    }
}

template<class Coord>
template<int k>
Tensor3x3<Coord,LF> Metric2D<Coord>::MetricDeriv(Coordinate2<Coord> x)
{
    double dk = 1e-8;
    Tensor3x3<Coord,LF> g_ll_m2 = MetricFunction(Coordinate2<Coord>{x[1]-2*dk*(k==1), x[2]-2*dk*(k==2)});
    Tensor3x3<Coord,LF> g_ll_m1 = MetricFunction(Coordinate2<Coord>{x[1]-1*dk*(k==1), x[2]-1*dk*(k==2)});
    Tensor3x3<Coord,LF> g_ll_p1 = MetricFunction(Coordinate2<Coord>{x[1]+1*dk*(k==1), x[2]+1*dk*(k==2)});
    Tensor3x3<Coord,LF> g_ll_p2 = MetricFunction(Coordinate2<Coord>{x[1]+2*dk*(k==1), x[2]+2*dk*(k==2)});

    Tensor3x3<Coord,LF> dk_g_ll;
    for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
        dk_g_ll[{i,j}] = (1.0/12.0*g_ll_m2[{i,j}] - 1/12.0*g_ll_p2[{i,j}] - 2.0/3.0*g_ll_m1[{i,j}] + 2.0/3.0*g_ll_p1[{i,j}]) / dk;
    return dk_g_ll;
}


template<class Coord>
template<int k>
Tensor3x3<Coord,LF> Metric2D<Coord>::InverseMetricDeriv(Coordinate2<Coord> x)
{
    double dk = 1e-8;
    Tensor3x3<Coord,LF> g_uu_m2 = MetricFunction(Coordinate2<Coord>{x[1]-2*dk*(k==1), x[2]-2*dk*(k==2)}).Invert();
    Tensor3x3<Coord,LF> g_uu_m1 = MetricFunction(Coordinate2<Coord>{x[1]-1*dk*(k==1), x[2]-1*dk*(k==2)}).Invert();
    Tensor3x3<Coord,LF> g_uu_p1 = MetricFunction(Coordinate2<Coord>{x[1]+1*dk*(k==1), x[2]+1*dk*(k==2)}).Invert();
    Tensor3x3<Coord,LF> g_uu_p2 = MetricFunction(Coordinate2<Coord>{x[1]+2*dk*(k==1), x[2]+2*dk*(k==2)}).Invert();

    Tensor3x3<Coord,LF> dk_g_uu;
    for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
        dk_g_uu[{i,j}] = (1.0/12.0*g_uu_m2[{i,j}] - 1/12.0*g_uu_p2[{i,j}] - 2.0/3.0*g_uu_m1[{i,j}] + 2/3.0*g_uu_p1[{i,j}]) / dk;
    return dk_g_uu;
}

template<class Coord>
void Metric2D<Coord>::InitializeMetricDerivativesOnGrid()
{
    for(int j=0; j<grid.n2; j++)
    for(int i=0; i<grid.n1; i++)
    {
        Coordinate2<Coord> x = grid.x12Coord(i,j);
        Tensor3x3<Coord,LF> dt_gll(0.0);
        Tensor3x3<Coord,LF> dx_gll = MetricDeriv<1>(x);
        Tensor3x3<Coord,LF> dy_gll = MetricDeriv<2>(x);
        Tensor3x3<Coord,LF> dt_guu(0.0);
        Tensor3x3<Coord,LF> dx_guu = InverseMetricDeriv<1>(x);
        Tensor3x3<Coord,LF> dy_guu = InverseMetricDeriv<2>(x);

        int ij = grid.Index(i,j);
        d0_g00_lll[ij] = dt_gll[{0,0}];    d1_g00_lll[ij] = dx_gll[{0,0}];     d2_g00_lll[ij] = dy_gll[{0,0}];
        d0_g01_lll[ij] = dt_gll[{0,1}];    d1_g01_lll[ij] = dx_gll[{0,1}];     d2_g01_lll[ij] = dy_gll[{0,1}];
        d0_g02_lll[ij] = dt_gll[{0,2}];    d1_g02_lll[ij] = dx_gll[{0,2}];     d2_g02_lll[ij] = dy_gll[{0,2}];
        d0_g11_lll[ij] = dt_gll[{1,1}];    d1_g11_lll[ij] = dx_gll[{1,1}];     d2_g11_lll[ij] = dy_gll[{1,1}];
        d0_g12_lll[ij] = dt_gll[{1,2}];    d1_g12_lll[ij] = dx_gll[{1,2}];     d2_g12_lll[ij] = dy_gll[{1,2}];
        d0_g22_lll[ij] = dt_gll[{2,2}];    d1_g22_lll[ij] = dx_gll[{2,2}];     d2_g22_lll[ij] = dy_gll[{2,2}];

        d0_g00_luu[ij] = dt_guu[{0,0}];    d1_g00_luu[ij] = dx_guu[{0,0}];     d2_g00_luu[ij] = dy_guu[{0,0}];
        d0_g01_luu[ij] = dt_guu[{0,1}];    d1_g01_luu[ij] = dx_guu[{0,1}];     d2_g01_luu[ij] = dy_guu[{0,1}];
        d0_g02_luu[ij] = dt_guu[{0,2}];    d1_g02_luu[ij] = dx_guu[{0,2}];     d2_g02_luu[ij] = dy_guu[{0,2}];
        d0_g11_luu[ij] = dt_guu[{1,1}];    d1_g11_luu[ij] = dx_guu[{1,1}];     d2_g11_luu[ij] = dy_guu[{1,1}];
        d0_g12_luu[ij] = dt_guu[{1,2}];    d1_g12_luu[ij] = dx_guu[{1,2}];     d2_g12_luu[ij] = dy_guu[{1,2}];
        d0_g22_luu[ij] = dt_guu[{2,2}];    d1_g22_luu[ij] = dx_guu[{2,2}];     d2_g22_luu[ij] = dy_guu[{2,2}];
    }
}

template<class Coord>
void Metric2D<Coord>::InitializeAdmComponentsOnGrid()
{
    for(int j=0; j<grid.n2; j++)
    for(int i=0; i<grid.n1; i++)
    {
        int ij = grid.Index(i,j);
        alpha[ij] = 1.0/sqrt(-g00_uu[ij]);
        beta1_u[ij] = alpha[ij] * alpha[ij] * g01_uu[ij];
        beta2_u[ij] = alpha[ij] * alpha[ij] * g02_uu[ij];
        beta1_l[ij] = g01_ll[ij];
        beta2_l[ij] = g02_ll[ij];
        gamma11_ll[ij] = g11_ll[ij];    gamma11_uu[ij] = g11_uu[ij] + beta1_u[ij] * beta1_u[ij] / (alpha[ij] * alpha[ij]);
        gamma12_ll[ij] = g12_ll[ij];    gamma12_uu[ij] = g12_uu[ij] + beta1_u[ij] * beta2_u[ij] / (alpha[ij] * alpha[ij]);
        gamma22_ll[ij] = g22_ll[ij];    gamma22_uu[ij] = g22_uu[ij] + beta2_u[ij] * beta2_u[ij] / (alpha[ij] * alpha[ij]);

        d1_alpha_l[ij] = d1_g00_luu[ij] / (2.0 * pow(-g00_uu[ij],3.0/2.0));
        d2_alpha_l[ij] = d2_g00_luu[ij] / (2.0 * pow(-g00_uu[ij],3.0/2.0));
        d1_beta1_lu[ij] = alpha[ij] * alpha[ij] * d1_g01_luu[ij] + 2.0 * g01_uu[ij] * alpha[ij] * d1_alpha_l[ij];
        d2_beta1_lu[ij] = alpha[ij] * alpha[ij] * d2_g01_luu[ij] + 2.0 * g01_uu[ij] * alpha[ij] * d2_alpha_l[ij];
        d1_beta2_lu[ij] = alpha[ij] * alpha[ij] * d1_g02_luu[ij] + 2.0 * g02_uu[ij] * alpha[ij] * d1_alpha_l[ij];
        d2_beta2_lu[ij] = alpha[ij] * alpha[ij] * d2_g02_luu[ij] + 2.0 * g02_uu[ij] * alpha[ij] * d2_alpha_l[ij];
        d1_beta1_ll[ij] = d1_g01_lll[ij];    d1_beta2_ll[ij] = d1_g02_lll[ij];
        d2_beta1_ll[ij] = d2_g01_lll[ij];    d2_beta2_ll[ij] = d2_g02_lll[ij];
        d1_gamma11_lll[ij] = d1_g11_lll[ij];    d2_gamma11_lll[ij] = d2_g11_lll[ij];
        d1_gamma12_lll[ij] = d1_g12_lll[ij];    d2_gamma12_lll[ij] = d2_g12_lll[ij];
        d1_gamma22_lll[ij] = d1_g22_lll[ij];    d2_gamma22_lll[ij] = d2_g22_lll[ij];
        d1_gamma11_luu[ij] = d1_g11_luu[ij] + (d1_beta1_lu[ij] * beta1_u[ij] + beta1_u[ij] * d1_beta1_lu[ij]) / (alpha[ij] * alpha[ij]) - 2.0 * beta1_u[ij] * beta1_u[ij] * d1_alpha_l[ij]/pow(alpha[ij],3);
        d1_gamma12_luu[ij] = d1_g12_luu[ij] + (d1_beta1_lu[ij] * beta2_u[ij] + beta1_u[ij] * d1_beta2_lu[ij]) / (alpha[ij] * alpha[ij]) - 2.0 * beta1_u[ij] * beta2_u[ij] * d1_alpha_l[ij]/pow(alpha[ij],3);
        d1_gamma22_luu[ij] = d1_g22_luu[ij] + (d1_beta2_lu[ij] * beta2_u[ij] + beta2_u[ij] * d1_beta2_lu[ij]) / (alpha[ij] * alpha[ij]) - 2.0 * beta2_u[ij] * beta2_u[ij] * d1_alpha_l[ij]/pow(alpha[ij],3);
        d2_gamma11_luu[ij] = d2_g11_luu[ij] + (d2_beta1_lu[ij] * beta1_u[ij] + beta1_u[ij] * d2_beta1_lu[ij]) / (alpha[ij] * alpha[ij]) - 2.0 * beta1_u[ij] * beta1_u[ij] * d2_alpha_l[ij]/pow(alpha[ij],3);
        d2_gamma12_luu[ij] = d2_g12_luu[ij] + (d2_beta1_lu[ij] * beta2_u[ij] + beta1_u[ij] * d2_beta2_lu[ij]) / (alpha[ij] * alpha[ij]) - 2.0 * beta1_u[ij] * beta2_u[ij] * d2_alpha_l[ij]/pow(alpha[ij],3);
        d2_gamma22_luu[ij] = d2_g22_luu[ij] + (d2_beta2_lu[ij] * beta2_u[ij] + beta2_u[ij] * d2_beta2_lu[ij]) / (alpha[ij] * alpha[ij]) - 2.0 * beta2_u[ij] * beta2_u[ij] * d2_alpha_l[ij]/pow(alpha[ij],3);

        Tensor2<Coord,LF> beta_u = GetBeta_u(ij);
        Tensor2x2<Coord,LF> dl_beta_u = GetDerivBeta_lu(ij);
        Tensor2x2<Coord,LF> gamma_ll = GetGamma_ll(ij);
        Tensor2x2<Coord,LF> dt_gamma_ll(0.0);
        Tensor2x2x2<Coord,LF> dl_gamma_ll = GetDerivGamma_lll(ij);
        Tensor2x2<Coord,LF> K_ll(0.0);
        for(int a=1; a<3; a++)
        for(int b=1; b<3; b++)
        {
            K_ll[{a,b}] = -dt_gamma_ll[{a,b}];
            for(int c=1; c<3; c++)
                K_ll[{a,b}] += 2.0 * gamma_ll[{a,c}] * dl_beta_u[{b,c}] + dl_gamma_ll[{c,a,b}] * beta_u[c];
            K_ll[{a,b}] /= 2.0 * alpha[ij];
        }
        K11_ll[ij] = K_ll[{1,1}];
        K12_ll[ij] = K_ll[{1,2}];
        K22_ll[ij] = K_ll[{2,2}];
    }
}

template<class Coord>
double Metric2D<Coord>::InterpolateArrayTo_ij(double* array, double i, double j)
{
    int i0 = std::floor(i);
    int j0 = std::floor(j);
    int i1 = std::ceil(i);
    int j1 = std::ceil(j);

    return BilinearInterpolation
    (i-i0,j-j0,
    array[grid.Index(i0,j0)],array[grid.Index(i0,j1)],
    array[grid.Index(i1,j0)],array[grid.Index(i1,j1)]);
}

template<class Coord>
bool Metric2D<Coord>::InsideBH(const int i, const int j)
{
    exit_on_error("Metric2D virtual Method (InsideBH) has been called!");
    return false;
}

template<class Coord>
bool Metric2D<Coord>::InsideBH(const Coordinate2<Coord>& x12)
{
    exit_on_error("Metric2D virtual Method (InsideBH) has been called!");
    return false;
}



// Tensor getters:
template<class Coord>
Tensor3<Coord,LF> Metric2D<Coord>::uEulObs(const int ij)
{
    double alpha = GetAlpha(ij);
    Tensor2<Coord,LF> beta_u = GetBeta_u(ij);
    return Tensor3<Coord,LF>(1.0/alpha, -beta_u[1]/alpha, -beta_u[2]/alpha);
}
template<class Coord>
Tensor3<Coord,LF> Metric2D<Coord>::uEulObs(const Coordinate2<Coord>& x12)
{
    double alpha = GetAlpha(x12);
    Tensor2<Coord,LF> beta_u = GetBeta_u(x12);
    return Tensor3<Coord,LF> (1.0/alpha, -beta_u[1]/alpha, -beta_u[2]/alpha);
}
template<class Coord>
Tensor3x3<Coord,LF> Metric2D<Coord>::GetMetric_ll(const int ij)
{
    return Tensor3x3<Coord,LF>
    (g00_ll[ij], g01_ll[ij], g02_ll[ij],
     g01_ll[ij], g11_ll[ij], g12_ll[ij],
     g02_ll[ij], g12_ll[ij], g22_ll[ij]);
}
template<class Coord>
Tensor3x3<Coord,LF> Metric2D<Coord>::GetMetric_ll(const Coordinate2<Coord>& x12)
{
    Double2 ij = grid.ij(x12);
    double i = (x12[1]-grid.start1)/grid.d1;
    double j = (x12[2]-grid.start2)/grid.d2;
    // TODO: think about polar edge!
    if(ij[0]<0 or ij[1]<0 or ij[0]>grid.n1-1 or ij[1]>grid.n2-1)
    {
        if constexpr(std::is_same<Coord,xy>::value)
        { return Tensor3x3<Coord,LF>(-1,0,0, 0,1,0, 0,0,1); }
        if constexpr(std::is_same<Coord,rph>::value)
        { return Tensor3x3<Coord,LF>(-1,0,0 ,0,1,0, 0,0,x12[1]*x12[1]); }
    }
    else
    {
        Tensor3x3<Coord,LF> g_ll;
        g_ll[{0,0}] = InterpolateArrayTo_ij(g00_ll,i,j);
        g_ll[{1,1}] = InterpolateArrayTo_ij(g11_ll,i,j);
        g_ll[{2,2}] = InterpolateArrayTo_ij(g22_ll,i,j);
        g_ll[{0,1}] = g_ll[{1,0}] = InterpolateArrayTo_ij(g01_ll,i,j);
        g_ll[{0,2}] = g_ll[{2,0}] = InterpolateArrayTo_ij(g02_ll,i,j);
        g_ll[{1,2}] = g_ll[{2,1}] = InterpolateArrayTo_ij(g12_ll,i,j);
        return g_ll;
    }
}

template<class Coord>
Tensor3x3<Coord,LF> Metric2D<Coord>::GetMetric_uu(const int ij)
{
    return Tensor3x3<Coord,LF>
    (g00_uu[ij], g01_uu[ij], g02_uu[ij],
     g01_uu[ij], g11_uu[ij], g12_uu[ij],
     g02_uu[ij], g12_uu[ij], g22_uu[ij]);
}
template<class Coord>
Tensor3x3<Coord,LF> Metric2D<Coord>::GetMetric_uu(const Coordinate2<Coord>& x12)
{
    double i = (x12[1]-grid.start1)/grid.d1;
    double j = (x12[2]-grid.start2)/grid.d2;
    if(i<0 or j<0 or i>grid.n1-1 or j>grid.n2-1)
    {
        if constexpr(std::is_same<Coord,xy>::value)
        { return Tensor3x3<Coord,LF>(-1,0,0, 0,1,0, 0,0,1); }
        if constexpr(std::is_same<Coord,rph>::value)
        { return Tensor3x3<Coord,LF>(-1,0,0 ,0,1,0, 0,0,x12[1]*x12[1]); }
    }
    else
    {
        Tensor3x3<Coord,LF> g_uu;
        g_uu[{0,0}] = InterpolateArrayTo_ij(g00_uu,i,j);
        g_uu[{1,1}] = InterpolateArrayTo_ij(g11_uu,i,j);
        g_uu[{2,2}] = InterpolateArrayTo_ij(g22_uu,i,j);
        g_uu[{0,1}] = g_uu[{1,0}] = InterpolateArrayTo_ij(g01_uu,i,j);
        g_uu[{0,2}] = g_uu[{2,0}] = InterpolateArrayTo_ij(g02_uu,i,j);
        g_uu[{1,2}] = g_uu[{2,1}] = InterpolateArrayTo_ij(g12_uu,i,j);
        return g_uu;
    }
}
template<class Coord>
Tensor3x3<Coord,IF> Metric2D<Coord>::GetMinkowskiMetric_ll(const int ij)
{
    if constexpr(std::is_same<Coord,xy>::value)
    { return Tensor3x3<Coord,IF>(-1,0,0, 0,1,0, 0,0,1); }
    if constexpr(std::is_same<Coord,rph>::value)
    {
        double r = grid.rCoord(ij);
        return Tensor3x3<Coord,IF>(-1,0,0 ,0,1,0, 0,0,r*r);
    }
}
template<class Coord>
Tensor3x3<Coord,IF> Metric2D<Coord>::GetMinkowskiMetric_ll(const Coordinate2<Coord>& x12)
{
    if constexpr(std::is_same<Coord,xy>::value)
    { return Tensor3x3<Coord,IF>(-1,0,0, 0,1,0, 0,0,1); }
    if constexpr(std::is_same<Coord,rph>::value)
    { return Tensor3x3<Coord,IF>(-1,0,0 ,0,1,0, 0,0,x12[1]*x12[1]); }
}
template<class Coord>
Tensor3x3<Coord,IF> Metric2D<Coord>::GetMinkowskiMetric_uu(const int ij)
{
    if constexpr(std::is_same<Coord,xy>::value)
    { return Tensor3x3<Coord,IF>(-1,0,0, 0,1,0, 0,0,1); }
    if constexpr(std::is_same<Coord,rph>::value)
    {
        double r = grid.rCoord(ij);
        return Tensor3x3<Coord,IF>(-1,0,0 ,0,1,0, 0,0,r*r);
    }
}
template<class Coord>
Tensor3x3<Coord,IF> Metric2D<Coord>::GetMinkowskiMetric_uu(const Coordinate2<Coord>& x12)
{
    if constexpr(std::is_same<Coord,xy>::value)
    { return Tensor3x3<Coord,IF>(-1,0,0, 0,1,0, 0,0,1); }
    if constexpr(std::is_same<Coord,rph>::value)
    { return Tensor3x3<Coord,IF>(-1,0,0 ,0,1,0, 0,0,x12[1]*x12[1]); }
}

template<class Coord>
Tensor3x3x3<Coord,LF> Metric2D<Coord>::GetDerivMetric_lll(const int ij)
{
    return Tensor3x3x3<Coord,LF>
    (d0_g00_lll[ij], d0_g01_lll[ij], d0_g02_lll[ij],
     d0_g01_lll[ij], d0_g11_lll[ij], d0_g12_lll[ij],
     d0_g02_lll[ij], d0_g12_lll[ij], d0_g22_lll[ij],

     d1_g00_lll[ij], d1_g01_lll[ij], d1_g02_lll[ij],
     d1_g01_lll[ij], d1_g11_lll[ij], d1_g12_lll[ij],
     d1_g02_lll[ij], d1_g12_lll[ij], d1_g22_lll[ij],

     d2_g00_lll[ij], d2_g01_lll[ij], d2_g02_lll[ij],
     d2_g01_lll[ij], d2_g11_lll[ij], d2_g12_lll[ij],
     d2_g02_lll[ij], d2_g12_lll[ij], d2_g22_lll[ij]);
}
template<class Coord>
Tensor3x3x3<Coord,LF> Metric2D<Coord>::GetDerivMetric_lll(const Coordinate2<Coord>& x12)
{
    double i = (x12[1]-grid.start1)/grid.d1;
    double j = (x12[2]-grid.start2)/grid.d2;
    if(i<0 or j<0 or i>grid.n1-1 or j>grid.n2-1)
        return Tensor3x3x3<Coord,LF>(0.0);
    else
    {
        Tensor3x3x3<Coord,LF> dl_g_ll;
        dl_g_ll[{0,0,0}] = InterpolateArrayTo_ij(d0_g00_lll,i,j);
        dl_g_ll[{0,1,1}] = InterpolateArrayTo_ij(d0_g11_lll,i,j);
        dl_g_ll[{0,2,2}] = InterpolateArrayTo_ij(d0_g22_lll,i,j);
        dl_g_ll[{0,0,1}] = dl_g_ll[{0,1,0}] = InterpolateArrayTo_ij(d0_g01_lll,i,j);
        dl_g_ll[{0,0,2}] = dl_g_ll[{0,2,0}] = InterpolateArrayTo_ij(d0_g02_lll,i,j);
        dl_g_ll[{0,1,2}] = dl_g_ll[{0,2,1}] = InterpolateArrayTo_ij(d0_g12_lll,i,j);
        
        dl_g_ll[{1,0,0}] = InterpolateArrayTo_ij(d1_g00_lll,i,j);
        dl_g_ll[{1,1,1}] = InterpolateArrayTo_ij(d1_g11_lll,i,j);
        dl_g_ll[{1,2,2}] = InterpolateArrayTo_ij(d1_g22_lll,i,j);
        dl_g_ll[{1,0,1}] = dl_g_ll[{1,1,0}] = InterpolateArrayTo_ij(d1_g01_lll,i,j);
        dl_g_ll[{1,0,2}] = dl_g_ll[{1,2,0}] = InterpolateArrayTo_ij(d1_g02_lll,i,j);
        dl_g_ll[{1,1,2}] = dl_g_ll[{1,2,1}] = InterpolateArrayTo_ij(d1_g12_lll,i,j);
        
        dl_g_ll[{2,0,0}] = InterpolateArrayTo_ij(d2_g00_lll,i,j);
        dl_g_ll[{2,1,1}] = InterpolateArrayTo_ij(d2_g11_lll,i,j);
        dl_g_ll[{2,2,2}] = InterpolateArrayTo_ij(d2_g22_lll,i,j);
        dl_g_ll[{2,0,1}] = dl_g_ll[{2,1,0}] = InterpolateArrayTo_ij(d2_g01_lll,i,j);
        dl_g_ll[{2,0,2}] = dl_g_ll[{2,2,0}] = InterpolateArrayTo_ij(d2_g02_lll,i,j);
        dl_g_ll[{2,1,2}] = dl_g_ll[{2,2,1}] = InterpolateArrayTo_ij(d2_g12_lll,i,j);
        return dl_g_ll;
    }
}


template<class Coord>
Tensor3x3x3<Coord,LF> Metric2D<Coord>::GetDerivMetric_luu(const int ij)
{
    return Tensor3x3x3<Coord,LF>
    (d0_g00_luu[ij], d0_g01_luu[ij], d0_g02_luu[ij],
     d0_g01_luu[ij], d0_g11_luu[ij], d0_g12_luu[ij],
     d0_g02_luu[ij], d0_g12_luu[ij], d0_g22_luu[ij],

     d1_g00_luu[ij], d1_g01_luu[ij], d1_g02_luu[ij],
     d1_g01_luu[ij], d1_g11_luu[ij], d1_g12_luu[ij],
     d1_g02_luu[ij], d1_g12_luu[ij], d1_g22_luu[ij],

     d2_g00_luu[ij], d2_g01_luu[ij], d2_g02_luu[ij],
     d2_g01_luu[ij], d2_g11_luu[ij], d2_g12_luu[ij],
     d2_g02_luu[ij], d2_g12_luu[ij], d2_g22_luu[ij]);
}
template<class Coord>
Tensor3x3x3<Coord,LF> Metric2D<Coord>::GetDerivMetric_luu(const Coordinate2<Coord>& x12)
{
    double i = (x12[1]-grid.start1)/grid.d1;
    double j = (x12[2]-grid.start2)/grid.d2;
    if(i<0 or j<0 or i>grid.n1-1 or j>grid.n2-1)
        return Tensor3x3x3<Coord,LF>(0.0);
    else
    {
        Tensor3x3x3<Coord,LF> dl_g_uu;
        dl_g_uu[{0,0,0}] = InterpolateArrayTo_ij(d0_g00_luu,i,j);
        dl_g_uu[{0,1,1}] = InterpolateArrayTo_ij(d0_g11_luu,i,j);
        dl_g_uu[{0,2,2}] = InterpolateArrayTo_ij(d0_g22_luu,i,j);
        dl_g_uu[{0,0,1}] = dl_g_uu[{0,1,0}] = InterpolateArrayTo_ij(d0_g01_luu,i,j);
        dl_g_uu[{0,0,2}] = dl_g_uu[{0,2,0}] = InterpolateArrayTo_ij(d0_g02_luu,i,j);
        dl_g_uu[{0,1,2}] = dl_g_uu[{0,2,1}] = InterpolateArrayTo_ij(d0_g12_luu,i,j);
        
        dl_g_uu[{1,0,0}] = InterpolateArrayTo_ij(d1_g00_luu,i,j);
        dl_g_uu[{1,1,1}] = InterpolateArrayTo_ij(d1_g11_luu,i,j);
        dl_g_uu[{1,2,2}] = InterpolateArrayTo_ij(d1_g22_luu,i,j);
        dl_g_uu[{1,0,1}] = dl_g_uu[{1,1,0}] = InterpolateArrayTo_ij(d1_g01_luu,i,j);
        dl_g_uu[{1,0,2}] = dl_g_uu[{1,2,0}] = InterpolateArrayTo_ij(d1_g02_luu,i,j);
        dl_g_uu[{1,1,2}] = dl_g_uu[{1,2,1}] = InterpolateArrayTo_ij(d1_g12_luu,i,j);
        
        dl_g_uu[{2,0,0}] = InterpolateArrayTo_ij(d2_g00_luu,i,j);
        dl_g_uu[{2,1,1}] = InterpolateArrayTo_ij(d2_g11_luu,i,j);
        dl_g_uu[{2,2,2}] = InterpolateArrayTo_ij(d2_g22_luu,i,j);
        dl_g_uu[{2,0,1}] = dl_g_uu[{2,1,0}] = InterpolateArrayTo_ij(d2_g01_luu,i,j);
        dl_g_uu[{2,0,2}] = dl_g_uu[{2,2,0}] = InterpolateArrayTo_ij(d2_g02_luu,i,j);
        dl_g_uu[{2,1,2}] = dl_g_uu[{2,2,1}] = InterpolateArrayTo_ij(d2_g12_luu,i,j);
        return dl_g_uu;
    }
}

template<class Coord>
Tensor3x3<Coord,Tetrad> Metric2D<Coord>::GetTetrad(const int ij)
{
    return Tensor3x3<Coord,Tetrad>
    (tetrad00_ul[ij], tetrad01_ul[ij], tetrad02_ul[ij],
     tetrad10_ul[ij], tetrad11_ul[ij], tetrad12_ul[ij],
     tetrad20_ul[ij], tetrad21_ul[ij], tetrad22_ul[ij]);
}
template<class Coord>
Tensor3x3<Coord,Tetrad> Metric2D<Coord>::GetTetrad(const Coordinate2<Coord>& x12)
{
    double i = (x12[1]-grid.start1)/grid.d1;
    double j = (x12[2]-grid.start2)/grid.d2;
    if(i<0 or j<0 or i>grid.n1-1 or j>grid.n2-1)
    {
        return Tensor3x3<Coord,Tetrad>
        (1, 0, 0,
         0, 1, 0,
         0, 0, 1);
    }
    else
    {
        return Tensor3x3<Coord,Tetrad>
        (InterpolateArrayTo_ij(tetrad00_ul,i,j), InterpolateArrayTo_ij(tetrad01_ul,i,j), InterpolateArrayTo_ij(tetrad02_ul,i,j),
         InterpolateArrayTo_ij(tetrad10_ul,i,j), InterpolateArrayTo_ij(tetrad11_ul,i,j), InterpolateArrayTo_ij(tetrad12_ul,i,j),
         InterpolateArrayTo_ij(tetrad20_ul,i,j), InterpolateArrayTo_ij(tetrad21_ul,i,j), InterpolateArrayTo_ij(tetrad22_ul,i,j));
    }
}

// ADM getters:
template<class Coord>
double Metric2D<Coord>::GetAlpha(const int ij)
{
    return alpha[ij];
}
template<class Coord>
double Metric2D<Coord>::GetAlpha(const Coordinate2<Coord>& x12)
{
    double i = (x12[1]-grid.start1)/grid.d1;
    double j = (x12[2]-grid.start2)/grid.d2;
    if(i<0 or j<0 or i>grid.n1-1 or j>grid.n2-1)
        return 1.0;
    else
        return InterpolateArrayTo_ij(alpha,i,j);
}

template<class Coord>
Tensor2<Coord,LF> Metric2D<Coord>::GetBeta_u(const int ij)
{
    return Tensor2<Coord,LF>(beta1_u[ij],beta2_u[ij]);
}
template<class Coord>
Tensor2<Coord,LF> Metric2D<Coord>::GetBeta_u(const Coordinate2<Coord>& x12)
{
    double i = (x12[1]-grid.start1)/grid.d1;
    double j = (x12[2]-grid.start2)/grid.d2;
    if(i<0 or j<0 or i>grid.n1-1 or j>grid.n2-1)
        return Tensor2<Coord,LF>(0.0);
    else
        return Tensor2<Coord,LF>(InterpolateArrayTo_ij(beta1_u,i,j), InterpolateArrayTo_ij(beta2_u,i,j));
}

template<class Coord>
Tensor2<Coord,LF> Metric2D<Coord>::GetBeta_l(const int ij)
{
    return Tensor2<Coord,LF>(beta1_l[ij],beta2_l[ij]);
}
template<class Coord>
Tensor2<Coord,LF> Metric2D<Coord>::GetBeta_l(const Coordinate2<Coord>& x12)
{
    double i = (x12[1]-grid.start1)/grid.d1;
    double j = (x12[2]-grid.start2)/grid.d2;
    if(i<0 or j<0 or i>grid.n1-1 or j>grid.n2-1)
        return Tensor2<Coord,LF>(0.0);
    else
        return Tensor2<Coord,LF>(InterpolateArrayTo_ij(beta1_l,i,j), InterpolateArrayTo_ij(beta2_l,i,j));
}

template<class Coord>
Tensor2x2<Coord,LF> Metric2D<Coord>::GetGamma_ll(const int ij)
{
    return Tensor2x2<Coord,LF>
    (gamma11_ll[ij],gamma12_ll[ij],
     gamma12_ll[ij],gamma22_ll[ij]);
}
template<class Coord>
Tensor2x2<Coord,LF> Metric2D<Coord>::GetGamma_ll(const Coordinate2<Coord>& x12)
{
    double i = (x12[1]-grid.start1)/grid.d1;
    double j = (x12[2]-grid.start2)/grid.d2;
    if(i<0 or j<0 or i>grid.n1-1 or j>grid.n2-1)
        return Tensor2x2<Coord,LF>(1.0, 0.0, 1.0, 0.0);
    else
    {
        Tensor2x2<Coord,LF> gamma_ll;
        gamma_ll[{1,1}] = InterpolateArrayTo_ij(gamma11_ll,i,j);
        gamma_ll[{2,2}] = InterpolateArrayTo_ij(gamma22_ll,i,j);
        gamma_ll[{1,2}] = gamma_ll[{2,1}] = InterpolateArrayTo_ij(gamma12_ll,i,j);
        return gamma_ll;
    }
}

template<class Coord>
Tensor2x2<Coord,LF> Metric2D<Coord>::GetGamma_uu(const int ij)
{
    return Tensor2x2<Coord,LF>
    (gamma11_uu[ij],gamma12_uu[ij],
     gamma12_uu[ij],gamma22_uu[ij]);
}
template<class Coord>
Tensor2x2<Coord,LF> Metric2D<Coord>::GetGamma_uu(const Coordinate2<Coord>& x12)
{
    double i = (x12[1]-grid.start1)/grid.d1;
    double j = (x12[2]-grid.start2)/grid.d2;
    if(i<0 or j<0 or i>grid.n1-1 or j>grid.n2-1)
        return Tensor2x2<Coord,LF>(1.0, 0.0, 1.0, 0.0);
    else
    {
        Tensor2x2<Coord,LF> gamma_uu;
        gamma_uu[{1,1}] = InterpolateArrayTo_ij(gamma11_uu,i,j);
        gamma_uu[{2,2}] = InterpolateArrayTo_ij(gamma22_uu,i,j);
        gamma_uu[{1,2}] = gamma_uu[{2,1}] = InterpolateArrayTo_ij(gamma12_uu,i,j);
        return gamma_uu;
    }
}

template<class Coord>
Tensor2x2<Coord,IF> Metric2D<Coord>::GetMinkowskiGamma_ll(const int ij)
{
    if constexpr(std::is_same<Coord,xy>::value)
    { return Tensor2x2<Coord,IF>(1,0, 0,1); }
    if constexpr(std::is_same<Coord,rph>::value)
    {
        double r = grid.rCoord(ij);
        return Tensor2x2<Coord,IF>(1,0, 0,r*r);
    }
}
template<class Coord>
Tensor2x2<Coord,IF> Metric2D<Coord>::GetMinkowskiGamma_ll(const Coordinate2<Coord>& x12)
{
    if constexpr(std::is_same<Coord,xy>::value)
    { return Tensor2x2<Coord,IF>(1,0, 0,1); }
    if constexpr(std::is_same<Coord,rph>::value)
    { return Tensor2x2<Coord,IF>(1,0, 0,x12[1]*x12[1]); }
}

template<class Coord>
Tensor2x2<Coord,IF> Metric2D<Coord>::GetMinkowskiGamma_uu(const int ij)
{
    if constexpr(std::is_same<Coord,xy>::value)
    { return Tensor2x2<Coord,IF>(1,0, 0,1); }
    if constexpr(std::is_same<Coord,rph>::value)
    {
        double r = grid.rCoord(ij);
        return Tensor2x2<Coord,IF>(1,0, 0,r*r);
    }
}
template<class Coord>
Tensor2x2<Coord,IF> Metric2D<Coord>::GetMinkowskiGamma_uu(const Coordinate2<Coord>& x12)
{
    if constexpr(std::is_same<Coord,xy>::value)
    { return Tensor2x2<Coord,IF>(1,0, 0,1); }
    if constexpr(std::is_same<Coord,rph>::value)
    { return Tensor2x2<Coord,IF>(1,0, 0,x12[1]*x12[1]); }
}

template<class Coord>
Tensor2<Coord,LF> Metric2D<Coord>::GetDerivAlpha_l(const int ij)
{
    return Tensor2<Coord,LF>(d1_alpha_l[ij],d2_alpha_l[ij]);
}
template<class Coord>
Tensor2<Coord,LF> Metric2D<Coord>::GetDerivAlpha_l(const Coordinate2<Coord>& x12)
{
    double i = (x12[1]-grid.start1)/grid.d1;
    double j = (x12[2]-grid.start2)/grid.d2;
    if(i<0 or j<0 or i>grid.n1-1 or j>grid.n2-1)
        return Tensor2<Coord,LF>(0.0);
    else
        return Tensor2<Coord,LF>(InterpolateArrayTo_ij(d1_alpha_l,i,j), InterpolateArrayTo_ij(d2_alpha_l,i,j));
}

template<class Coord>
Tensor2x2<Coord,LF> Metric2D<Coord>::GetDerivBeta_lu(const int ij)
{
    return Tensor2x2<Coord,LF>
    (d1_beta1_lu[ij],d1_beta2_lu[ij],
     d2_beta1_lu[ij],d2_beta2_lu[ij]);
}
template<class Coord>
Tensor2x2<Coord,LF> Metric2D<Coord>::GetDerivBeta_lu(const Coordinate2<Coord>& x12)
{
    double i = (x12[1]-grid.start1)/grid.d1;
    double j = (x12[2]-grid.start2)/grid.d2;
    if(i<0 or j<0 or i>grid.n1-1 or j>grid.n2-1)
        return Tensor2x2<Coord,LF>(0.0);
    else
    {
        return Tensor2x2<Coord,LF>
        (InterpolateArrayTo_ij(d1_beta1_lu,i,j), InterpolateArrayTo_ij(d1_beta2_lu,i,j),
         InterpolateArrayTo_ij(d2_beta1_lu,i,j), InterpolateArrayTo_ij(d2_beta2_lu,i,j));
    }
}

template<class Coord>
Tensor2x2<Coord,LF> Metric2D<Coord>::GetDerivBeta_ll(const int ij)
{
    return Tensor2x2<Coord,LF>
    (d1_beta1_ll[ij],d1_beta2_ll[ij],
     d2_beta1_ll[ij],d2_beta2_ll[ij]);
}
template<class Coord>
Tensor2x2<Coord,LF> Metric2D<Coord>::GetDerivBeta_ll(const Coordinate2<Coord>& x12)
{
    double i = (x12[1]-grid.start1)/grid.d1;
    double j = (x12[2]-grid.start2)/grid.d2;
    if(i<0 or j<0 or i>grid.n1-1 or j>grid.n2-1)
        return Tensor2x2<Coord,LF>(0.0);
    else
    {
        return Tensor2x2<Coord,LF>
        (InterpolateArrayTo_ij(d1_beta1_ll,i,j), InterpolateArrayTo_ij(d1_beta2_ll,i,j),
         InterpolateArrayTo_ij(d2_beta1_ll,i,j), InterpolateArrayTo_ij(d2_beta2_ll,i,j));
    }
}

template<class Coord>
Tensor2x2x2<Coord,LF> Metric2D<Coord>::GetDerivGamma_lll(const int ij)
{
    return Tensor2x2x2<Coord,LF>
    (d1_gamma11_lll[ij], d1_gamma12_lll[ij],
     d1_gamma12_lll[ij], d1_gamma22_lll[ij],
     d2_gamma11_lll[ij], d2_gamma12_lll[ij],
     d2_gamma12_lll[ij], d2_gamma22_lll[ij]);
}
template<class Coord>
Tensor2x2x2<Coord,LF> Metric2D<Coord>::GetDerivGamma_lll(const Coordinate2<Coord>& x12)
{
    double i = (x12[1]-grid.start1)/grid.d1;
    double j = (x12[2]-grid.start2)/grid.d2;
    if(i<0 or j<0 or i>grid.n1-1 or j>grid.n2-1)
        return Tensor2x2x2<Coord,LF>(0.0);
    else
    {
        Tensor2x2x2<Coord,LF> dGamma_lll;
        dGamma_lll[{1,1,1}] = InterpolateArrayTo_ij(d1_gamma11_lll,i,j);
        dGamma_lll[{1,2,2}] = InterpolateArrayTo_ij(d1_gamma22_lll,i,j);
        dGamma_lll[{1,1,2}] = dGamma_lll[{1,2,1}] = InterpolateArrayTo_ij(d1_gamma12_lll,i,j);
        dGamma_lll[{2,1,1}] = InterpolateArrayTo_ij(d2_gamma11_lll,i,j);
        dGamma_lll[{2,2,2}] = InterpolateArrayTo_ij(d2_gamma22_lll,i,j);
        dGamma_lll[{2,1,2}] = dGamma_lll[{2,2,1}] = InterpolateArrayTo_ij(d2_gamma12_lll,i,j);
        return dGamma_lll;
    }
}

template<class Coord>
Tensor2x2x2<Coord,LF> Metric2D<Coord>::GetDerivGamma_luu(const int ij)
{
    return Tensor2x2x2<Coord,LF>
    (d1_gamma11_luu[ij], d1_gamma12_luu[ij],
     d1_gamma12_luu[ij], d1_gamma22_luu[ij],
     d2_gamma11_luu[ij], d2_gamma12_luu[ij],
     d2_gamma12_luu[ij], d2_gamma22_luu[ij]);
}
template<class Coord>
Tensor2x2x2<Coord,LF> Metric2D<Coord>::GetDerivGamma_luu(const Coordinate2<Coord>& x12)
{
    double i = (x12[1]-grid.start1)/grid.d1;
    double j = (x12[2]-grid.start2)/grid.d2;
    if(i<0 or j<0 or i>grid.n1-1 or j>grid.n2-1)
        return Tensor2x2x2<Coord,LF>(0.0);
    else
    {
        Tensor2x2x2<Coord,LF> dGamma_luu;
        dGamma_luu[{1,1,1}] = InterpolateArrayTo_ij(d1_gamma11_luu,i,j);
        dGamma_luu[{1,2,2}] = InterpolateArrayTo_ij(d1_gamma22_luu,i,j);
        dGamma_luu[{1,1,2}] = dGamma_luu[{1,2,1}] = InterpolateArrayTo_ij(d1_gamma12_luu,i,j);
        dGamma_luu[{2,1,1}] = InterpolateArrayTo_ij(d2_gamma11_luu,i,j);
        dGamma_luu[{2,2,2}] = InterpolateArrayTo_ij(d2_gamma22_luu,i,j);
        dGamma_luu[{2,1,2}] = dGamma_luu[{2,2,1}] = InterpolateArrayTo_ij(d2_gamma12_luu,i,j);
        return dGamma_luu;
    }
}

template<class Coord>
Tensor2x2<Coord,LF> Metric2D<Coord>::GetExtrCurv_ll(const int ij)
{
    return Tensor2x2<Coord,LF>
    (K11_ll[ij],K12_ll[ij],
     K12_ll[ij],K22_ll[ij]);
}
template<class Coord>
Tensor2x2<Coord,LF> Metric2D<Coord>::GetExtrCurv_ll(const Coordinate2<Coord>& x12)
{
    double i = (x12[1]-grid.start1)/grid.d1;
    double j = (x12[2]-grid.start2)/grid.d2;
    if(i<0 or j<0 or i>grid.n1-1 or j>grid.n2-1)
        return Tensor2x2<Coord,LF>(0.0);
    else
    {
        Tensor2x2<Coord,LF> K_ll;
        K_ll[{1,1}] = InterpolateArrayTo_ij(K11_ll,i,j);
        K_ll[{2,2}] = InterpolateArrayTo_ij(K22_ll,i,j);
        K_ll[{1,2}] = K_ll[{2,1}] = InterpolateArrayTo_ij(K12_ll,i,j);
        return K_ll;
    }
}





template class Metric2D<xy>;
template class Metric2D<rph>;