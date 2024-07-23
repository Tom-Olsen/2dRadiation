#include "Metric.h"

// Constructor/Desctructor
Metric::Metric(Grid &grid_, double m_, double a_) : grid(grid_), m(m_), a(a_)
{
    // Metrik Data:
    g00_ll.resize(grid.nxy);
    g00_uu.resize(grid.nxy);
    g01_ll.resize(grid.nxy);
    g01_uu.resize(grid.nxy);
    g02_ll.resize(grid.nxy);
    g02_uu.resize(grid.nxy);
    g11_ll.resize(grid.nxy);
    g11_uu.resize(grid.nxy);
    g12_ll.resize(grid.nxy);
    g12_uu.resize(grid.nxy);
    g22_ll.resize(grid.nxy);
    g22_uu.resize(grid.nxy);
    d0_g00_lll.resize(grid.nxy);
    d1_g00_lll.resize(grid.nxy);
    d2_g00_lll.resize(grid.nxy);
    d0_g01_lll.resize(grid.nxy);
    d1_g01_lll.resize(grid.nxy);
    d2_g01_lll.resize(grid.nxy);
    d0_g02_lll.resize(grid.nxy);
    d1_g02_lll.resize(grid.nxy);
    d2_g02_lll.resize(grid.nxy);
    d0_g11_lll.resize(grid.nxy);
    d1_g11_lll.resize(grid.nxy);
    d2_g11_lll.resize(grid.nxy);
    d0_g12_lll.resize(grid.nxy);
    d1_g12_lll.resize(grid.nxy);
    d2_g12_lll.resize(grid.nxy);
    d0_g22_lll.resize(grid.nxy);
    d1_g22_lll.resize(grid.nxy);
    d2_g22_lll.resize(grid.nxy);
    d0_g00_luu.resize(grid.nxy);
    d1_g00_luu.resize(grid.nxy);
    d2_g00_luu.resize(grid.nxy);
    d0_g01_luu.resize(grid.nxy);
    d1_g01_luu.resize(grid.nxy);
    d2_g01_luu.resize(grid.nxy);
    d0_g02_luu.resize(grid.nxy);
    d1_g02_luu.resize(grid.nxy);
    d2_g02_luu.resize(grid.nxy);
    d0_g11_luu.resize(grid.nxy);
    d1_g11_luu.resize(grid.nxy);
    d2_g11_luu.resize(grid.nxy);
    d0_g12_luu.resize(grid.nxy);
    d1_g12_luu.resize(grid.nxy);
    d2_g12_luu.resize(grid.nxy);
    d0_g22_luu.resize(grid.nxy);
    d1_g22_luu.resize(grid.nxy);
    d2_g22_luu.resize(grid.nxy);
    // ADM Data:
    alpha.resize(grid.nxy);
    beta1_u.resize(grid.nxy);
    beta1_l.resize(grid.nxy);
    beta2_u.resize(grid.nxy);
    beta2_l.resize(grid.nxy);
    gamma11_ll.resize(grid.nxy);
    gamma11_uu.resize(grid.nxy);
    gamma12_ll.resize(grid.nxy);
    gamma12_uu.resize(grid.nxy);
    gamma22_ll.resize(grid.nxy);
    gamma22_uu.resize(grid.nxy);
    d1_alpha_l.resize(grid.nxy);
    d2_alpha_l.resize(grid.nxy);
    d1_beta1_lu.resize(grid.nxy);
    d1_beta2_lu.resize(grid.nxy);
    d2_beta1_lu.resize(grid.nxy);
    d2_beta2_lu.resize(grid.nxy);
    d1_beta1_ll.resize(grid.nxy);
    d1_beta2_ll.resize(grid.nxy);
    d2_beta1_ll.resize(grid.nxy);
    d2_beta2_ll.resize(grid.nxy);
    d1_gamma11_lll.resize(grid.nxy);
    d2_gamma11_lll.resize(grid.nxy);
    d1_gamma12_lll.resize(grid.nxy);
    d2_gamma12_lll.resize(grid.nxy);
    d1_gamma22_lll.resize(grid.nxy);
    d2_gamma22_lll.resize(grid.nxy);
    d1_gamma11_luu.resize(grid.nxy);
    d2_gamma11_luu.resize(grid.nxy);
    d1_gamma12_luu.resize(grid.nxy);
    d2_gamma12_luu.resize(grid.nxy);
    d1_gamma22_luu.resize(grid.nxy);
    d2_gamma22_luu.resize(grid.nxy);
    K11_ll.resize(grid.nxy);
    K12_ll.resize(grid.nxy);
    K22_ll.resize(grid.nxy);
    // Tetrad Data:
    tetrad00_ul.resize(grid.nxy);
    tetrad01_ul.resize(grid.nxy);
    tetrad02_ul.resize(grid.nxy);
    tetrad10_ul.resize(grid.nxy);
    tetrad11_ul.resize(grid.nxy);
    tetrad12_ul.resize(grid.nxy);
    tetrad20_ul.resize(grid.nxy);
    tetrad21_ul.resize(grid.nxy);
    tetrad22_ul.resize(grid.nxy);
    tetrad00_lu.resize(grid.nxy);
    tetrad01_lu.resize(grid.nxy);
    tetrad02_lu.resize(grid.nxy);
    tetrad10_lu.resize(grid.nxy);
    tetrad11_lu.resize(grid.nxy);
    tetrad12_lu.resize(grid.nxy);
    tetrad20_lu.resize(grid.nxy);
    tetrad21_lu.resize(grid.nxy);
    tetrad22_lu.resize(grid.nxy);
}

std::string Metric::Name()
{
    ExitOnError("Metric virtual Method (Name) has been called!");
    return "";
}

// Initialization:
Tensor3x3 Metric::MetricFunction(const Coord &xy)
{
    ExitOnError("Metric virtual Method (MetricFunction) has been called!");
    return Tensor3x3(0);
}

void Metric::InitializeMetricOnGrid()
{
    PARALLEL_FOR(2)
    for (size_t j = 0; j < grid.ny; j++)
        for (size_t i = 0; i < grid.nx; i++)
        {
            Coord xy = grid.xy(i, j);
            Tensor3x3 g_ll = MetricFunction(xy);
            Tensor3x3 g_uu = g_ll.Invert();
            size_t ij = grid.Index(i, j);

            g00_ll[ij] = g_ll[{0, 0}];
            g00_uu[ij] = g_uu[{0, 0}];
            g01_ll[ij] = g_ll[{0, 1}];
            g01_uu[ij] = g_uu[{0, 1}];
            g02_ll[ij] = g_ll[{0, 2}];
            g02_uu[ij] = g_uu[{0, 2}];
            g11_ll[ij] = g_ll[{1, 1}];
            g11_uu[ij] = g_uu[{1, 1}];
            g12_ll[ij] = g_ll[{1, 2}];
            g12_uu[ij] = g_uu[{1, 2}];
            g22_ll[ij] = g_ll[{2, 2}];
            g22_uu[ij] = g_uu[{2, 2}];
        }
}

void Metric::InitializeBoostedTetradOnGrid()
{
    for (size_t j = 0; j < grid.ny; j++)
        for (size_t i = 0; i < grid.nx; i++)
        {
            size_t ij = grid.Index(i, j);
            Tensor2x2 g1 = GetGamma_ll(ij);
            Tensor2x2 g1Inv = GetGamma_uu(ij);
            double a = 1.0 / sqrt(g1Inv[{1, 1}]);
            double b = a * a * g1Inv[{1, 2}];
            Tensor2x2 m1(1 / a, 0, b / a, 1);
            double g2 = m1[{2, 2}] * g1[{2, 2}] * m1[{2, 2}];

            Tensor3 n1 = uEulObs(ij);
            tetrad00_ul[ij] = n1[0];
            tetrad01_ul[ij] = 0;
            tetrad02_ul[ij] = 0;
            tetrad10_ul[ij] = n1[1];
            tetrad11_ul[ij] = m1[{1, 1}];
            tetrad12_ul[ij] = 0;
            tetrad20_ul[ij] = n1[2];
            tetrad21_ul[ij] = m1[{2, 1}];
            tetrad22_ul[ij] = 1 / sqrt(g2);

            Tensor3x3 inverseTetrad = GetTetrad(ij).Invert();
            tetrad00_lu[ij] = inverseTetrad[{0, 0}];
            tetrad01_lu[ij] = inverseTetrad[{0, 1}];
            tetrad02_lu[ij] = inverseTetrad[{0, 2}];
            tetrad10_lu[ij] = inverseTetrad[{1, 0}];
            tetrad11_lu[ij] = inverseTetrad[{1, 1}];
            tetrad12_lu[ij] = inverseTetrad[{1, 2}];
            tetrad20_lu[ij] = inverseTetrad[{2, 0}];
            tetrad21_lu[ij] = inverseTetrad[{2, 1}];
            tetrad22_lu[ij] = inverseTetrad[{2, 2}];
        }
}

template <int k>
Tensor3x3 Metric::MetricDeriv(const Coord &xy)
{
    double dk = 1e-8;
    Tensor3x3 g_ll_m2 = MetricFunction(Coord{xy[1] - 2 * dk * (k == 1), xy[2] - 2 * dk * (k == 2)});
    Tensor3x3 g_ll_m1 = MetricFunction(Coord{xy[1] - 1 * dk * (k == 1), xy[2] - 1 * dk * (k == 2)});
    Tensor3x3 g_ll_p1 = MetricFunction(Coord{xy[1] + 1 * dk * (k == 1), xy[2] + 1 * dk * (k == 2)});
    Tensor3x3 g_ll_p2 = MetricFunction(Coord{xy[1] + 2 * dk * (k == 1), xy[2] + 2 * dk * (k == 2)});

    Tensor3x3 dk_g_ll;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            dk_g_ll[{i, j}] = (1.0 / 12.0 * g_ll_m2[{i, j}] - 1 / 12.0 * g_ll_p2[{i, j}] - 2.0 / 3.0 * g_ll_m1[{i, j}] + 2.0 / 3.0 * g_ll_p1[{i, j}]) / dk;
    return dk_g_ll;
}

template <int k>
Tensor3x3 Metric::InverseMetricDeriv(const Coord &xy)
{
    double dk = 1e-8;
    Tensor3x3 g_uu_m2 = MetricFunction(Coord{xy[1] - 2 * dk * (k == 1), xy[2] - 2 * dk * (k == 2)}).Invert();
    Tensor3x3 g_uu_m1 = MetricFunction(Coord{xy[1] - 1 * dk * (k == 1), xy[2] - 1 * dk * (k == 2)}).Invert();
    Tensor3x3 g_uu_p1 = MetricFunction(Coord{xy[1] + 1 * dk * (k == 1), xy[2] + 1 * dk * (k == 2)}).Invert();
    Tensor3x3 g_uu_p2 = MetricFunction(Coord{xy[1] + 2 * dk * (k == 1), xy[2] + 2 * dk * (k == 2)}).Invert();

    Tensor3x3 dk_g_uu;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            dk_g_uu[{i, j}] = (1.0 / 12.0 * g_uu_m2[{i, j}] - 1 / 12.0 * g_uu_p2[{i, j}] - 2.0 / 3.0 * g_uu_m1[{i, j}] + 2 / 3.0 * g_uu_p1[{i, j}]) / dk;
    return dk_g_uu;
}

void Metric::InitializeMetricDerivativesOnGrid()
{
    PARALLEL_FOR(2)
    for (size_t j = 0; j < grid.ny; j++)
        for (size_t i = 0; i < grid.nx; i++)
        {
            Coord xy = grid.xy(i, j);
            Tensor3x3 dt_gll(0.0);
            Tensor3x3 dx_gll = MetricDeriv<1>(xy);
            Tensor3x3 dy_gll = MetricDeriv<2>(xy);
            Tensor3x3 dt_guu(0.0);
            Tensor3x3 dx_guu = InverseMetricDeriv<1>(xy);
            Tensor3x3 dy_guu = InverseMetricDeriv<2>(xy);

            size_t ij = grid.Index(i, j);
            d0_g00_lll[ij] = dt_gll[{0, 0}];
            d1_g00_lll[ij] = dx_gll[{0, 0}];
            d2_g00_lll[ij] = dy_gll[{0, 0}];
            d0_g01_lll[ij] = dt_gll[{0, 1}];
            d1_g01_lll[ij] = dx_gll[{0, 1}];
            d2_g01_lll[ij] = dy_gll[{0, 1}];
            d0_g02_lll[ij] = dt_gll[{0, 2}];
            d1_g02_lll[ij] = dx_gll[{0, 2}];
            d2_g02_lll[ij] = dy_gll[{0, 2}];
            d0_g11_lll[ij] = dt_gll[{1, 1}];
            d1_g11_lll[ij] = dx_gll[{1, 1}];
            d2_g11_lll[ij] = dy_gll[{1, 1}];
            d0_g12_lll[ij] = dt_gll[{1, 2}];
            d1_g12_lll[ij] = dx_gll[{1, 2}];
            d2_g12_lll[ij] = dy_gll[{1, 2}];
            d0_g22_lll[ij] = dt_gll[{2, 2}];
            d1_g22_lll[ij] = dx_gll[{2, 2}];
            d2_g22_lll[ij] = dy_gll[{2, 2}];

            d0_g00_luu[ij] = dt_guu[{0, 0}];
            d1_g00_luu[ij] = dx_guu[{0, 0}];
            d2_g00_luu[ij] = dy_guu[{0, 0}];
            d0_g01_luu[ij] = dt_guu[{0, 1}];
            d1_g01_luu[ij] = dx_guu[{0, 1}];
            d2_g01_luu[ij] = dy_guu[{0, 1}];
            d0_g02_luu[ij] = dt_guu[{0, 2}];
            d1_g02_luu[ij] = dx_guu[{0, 2}];
            d2_g02_luu[ij] = dy_guu[{0, 2}];
            d0_g11_luu[ij] = dt_guu[{1, 1}];
            d1_g11_luu[ij] = dx_guu[{1, 1}];
            d2_g11_luu[ij] = dy_guu[{1, 1}];
            d0_g12_luu[ij] = dt_guu[{1, 2}];
            d1_g12_luu[ij] = dx_guu[{1, 2}];
            d2_g12_luu[ij] = dy_guu[{1, 2}];
            d0_g22_luu[ij] = dt_guu[{2, 2}];
            d1_g22_luu[ij] = dx_guu[{2, 2}];
            d2_g22_luu[ij] = dy_guu[{2, 2}];
        }
}

void Metric::InitializeAdmComponentsOnGrid()
{
    PARALLEL_FOR(2)
    for (size_t j = 0; j < grid.ny; j++)
        for (size_t i = 0; i < grid.nx; i++)
        {
            size_t ij = grid.Index(i, j);
            alpha[ij] = 1.0 / sqrt(-g00_uu[ij]);
            beta1_u[ij] = alpha[ij] * alpha[ij] * g01_uu[ij];
            beta2_u[ij] = alpha[ij] * alpha[ij] * g02_uu[ij];
            beta1_l[ij] = g01_ll[ij];
            beta2_l[ij] = g02_ll[ij];
            gamma11_ll[ij] = g11_ll[ij];
            gamma11_uu[ij] = g11_uu[ij] + beta1_u[ij] * beta1_u[ij] / (alpha[ij] * alpha[ij]);
            gamma12_ll[ij] = g12_ll[ij];
            gamma12_uu[ij] = g12_uu[ij] + beta1_u[ij] * beta2_u[ij] / (alpha[ij] * alpha[ij]);
            gamma22_ll[ij] = g22_ll[ij];
            gamma22_uu[ij] = g22_uu[ij] + beta2_u[ij] * beta2_u[ij] / (alpha[ij] * alpha[ij]);

            d1_alpha_l[ij] = d1_g00_luu[ij] / (2.0 * pow(-g00_uu[ij], 3.0 / 2.0));
            d2_alpha_l[ij] = d2_g00_luu[ij] / (2.0 * pow(-g00_uu[ij], 3.0 / 2.0));
            d1_beta1_lu[ij] = alpha[ij] * alpha[ij] * d1_g01_luu[ij] + 2.0 * g01_uu[ij] * alpha[ij] * d1_alpha_l[ij];
            d2_beta1_lu[ij] = alpha[ij] * alpha[ij] * d2_g01_luu[ij] + 2.0 * g01_uu[ij] * alpha[ij] * d2_alpha_l[ij];
            d1_beta2_lu[ij] = alpha[ij] * alpha[ij] * d1_g02_luu[ij] + 2.0 * g02_uu[ij] * alpha[ij] * d1_alpha_l[ij];
            d2_beta2_lu[ij] = alpha[ij] * alpha[ij] * d2_g02_luu[ij] + 2.0 * g02_uu[ij] * alpha[ij] * d2_alpha_l[ij];
            d1_beta1_ll[ij] = d1_g01_lll[ij];
            d1_beta2_ll[ij] = d1_g02_lll[ij];
            d2_beta1_ll[ij] = d2_g01_lll[ij];
            d2_beta2_ll[ij] = d2_g02_lll[ij];
            d1_gamma11_lll[ij] = d1_g11_lll[ij];
            d2_gamma11_lll[ij] = d2_g11_lll[ij];
            d1_gamma12_lll[ij] = d1_g12_lll[ij];
            d2_gamma12_lll[ij] = d2_g12_lll[ij];
            d1_gamma22_lll[ij] = d1_g22_lll[ij];
            d2_gamma22_lll[ij] = d2_g22_lll[ij];
            d1_gamma11_luu[ij] = d1_g11_luu[ij] + (d1_beta1_lu[ij] * beta1_u[ij] + beta1_u[ij] * d1_beta1_lu[ij]) / (alpha[ij] * alpha[ij]) - 2.0 * beta1_u[ij] * beta1_u[ij] * d1_alpha_l[ij] / pow(alpha[ij], 3);
            d1_gamma12_luu[ij] = d1_g12_luu[ij] + (d1_beta1_lu[ij] * beta2_u[ij] + beta1_u[ij] * d1_beta2_lu[ij]) / (alpha[ij] * alpha[ij]) - 2.0 * beta1_u[ij] * beta2_u[ij] * d1_alpha_l[ij] / pow(alpha[ij], 3);
            d1_gamma22_luu[ij] = d1_g22_luu[ij] + (d1_beta2_lu[ij] * beta2_u[ij] + beta2_u[ij] * d1_beta2_lu[ij]) / (alpha[ij] * alpha[ij]) - 2.0 * beta2_u[ij] * beta2_u[ij] * d1_alpha_l[ij] / pow(alpha[ij], 3);
            d2_gamma11_luu[ij] = d2_g11_luu[ij] + (d2_beta1_lu[ij] * beta1_u[ij] + beta1_u[ij] * d2_beta1_lu[ij]) / (alpha[ij] * alpha[ij]) - 2.0 * beta1_u[ij] * beta1_u[ij] * d2_alpha_l[ij] / pow(alpha[ij], 3);
            d2_gamma12_luu[ij] = d2_g12_luu[ij] + (d2_beta1_lu[ij] * beta2_u[ij] + beta1_u[ij] * d2_beta2_lu[ij]) / (alpha[ij] * alpha[ij]) - 2.0 * beta1_u[ij] * beta2_u[ij] * d2_alpha_l[ij] / pow(alpha[ij], 3);
            d2_gamma22_luu[ij] = d2_g22_luu[ij] + (d2_beta2_lu[ij] * beta2_u[ij] + beta2_u[ij] * d2_beta2_lu[ij]) / (alpha[ij] * alpha[ij]) - 2.0 * beta2_u[ij] * beta2_u[ij] * d2_alpha_l[ij] / pow(alpha[ij], 3);

            Tensor2 beta_u = GetBeta_u(ij);
            Tensor2x2 dl_beta_u = GetDerivBeta_lu(ij);
            Tensor2x2 gamma_ll = GetGamma_ll(ij);
            Tensor2x2 dt_gamma_ll(0.0);
            Tensor2x2x2 dl_gamma_ll = GetDerivGamma_lll(ij);
            Tensor2x2 K_ll(0.0);
            for (int a = 1; a < 3; a++)
                for (int b = 1; b < 3; b++)
                {
                    K_ll[{a, b}] = -dt_gamma_ll[{a, b}];
                    for (int c = 1; c < 3; c++)
                        K_ll[{a, b}] += 2.0 * gamma_ll[{a, c}] * dl_beta_u[{b, c}] + dl_gamma_ll[{c, a, b}] * beta_u[c];
                    K_ll[{a, b}] /= 2.0 * alpha[ij];
                }
            K11_ll[ij] = K_ll[{1, 1}];
            K12_ll[ij] = K_ll[{1, 2}];
            K22_ll[ij] = K_ll[{2, 2}];
        }
}

double Metric::InterpolateArrayTo_ij(const DoubleBuffer &array, const Coord &ij)
{
    size_t i0 = std::floor(ij[1]);
    size_t j0 = std::floor(ij[2]);
    size_t i1 = i0 + 1;
    size_t j1 = j0 + 1;

    return BilinearInterpolation(ij[1] - i0, ij[2] - j0,
                                 array[grid.Index(i0, j0)], array[grid.Index(i0, j1)],
                                 array[grid.Index(i1, j0)], array[grid.Index(i1, j1)]);
}
double Metric::InterpolateArrayTo_ij(const DoubleBuffer &array, double i, double j)
{
    size_t i0 = std::floor(i);
    size_t j0 = std::floor(j);
    size_t i1 = i0 + 1;
    size_t j1 = j0 + 1;

    return BilinearInterpolation(i - i0, j - j0,
                                 array[grid.Index(i0, j0)], array[grid.Index(i0, j1)],
                                 array[grid.Index(i1, j0)], array[grid.Index(i1, j1)]);
}

bool Metric::InsideBH(const Coord &xy)
{
    ExitOnError("Metric virtual Method (InsideBH(xy)) has been called!");
    return false;
}

// Tensor getters:
Tensor3 Metric::uEulObs(size_t ij)
{
    double alpha = GetAlpha(ij);
    Tensor2 beta_u = GetBeta_u(ij);
    return Tensor3(1.0 / alpha, -beta_u[1] / alpha, -beta_u[2] / alpha);
}
Tensor3 Metric::uEulObs(const Coord &xy)
{
    double alpha = GetAlpha(xy);
    Tensor2 beta_u = GetBeta_u(xy);
    if (grid.OutsideDomain(xy))
        return Tensor3(1, 0, 0);
    else
        return Tensor3(1.0 / alpha, -beta_u[1] / alpha, -beta_u[2] / alpha);
}
Tensor3x3 Metric::GetMetric_ll(size_t ij)
{
    return Tensor3x3(g00_ll[ij], g01_ll[ij], g02_ll[ij],
                     g01_ll[ij], g11_ll[ij], g12_ll[ij],
                     g02_ll[ij], g12_ll[ij], g22_ll[ij]);
}
Tensor3x3 Metric::GetMetric_ll(const Coord &xy)
{
    if (grid.OutsideDomain(xy))
        return Tensor3x3(-1, 0, 0, 0, 1, 0, 0, 0, 1);
    else
    {
        Coord ij = grid.ij(xy);
        Tensor3x3 g_ll;
        g_ll[{0, 0}] = InterpolateArrayTo_ij(g00_ll, ij);
        g_ll[{1, 1}] = InterpolateArrayTo_ij(g11_ll, ij);
        g_ll[{2, 2}] = InterpolateArrayTo_ij(g22_ll, ij);
        g_ll[{0, 1}] = g_ll[{1, 0}] = InterpolateArrayTo_ij(g01_ll, ij);
        g_ll[{0, 2}] = g_ll[{2, 0}] = InterpolateArrayTo_ij(g02_ll, ij);
        g_ll[{1, 2}] = g_ll[{2, 1}] = InterpolateArrayTo_ij(g12_ll, ij);
        return g_ll;
    }
}
Tensor3x3 Metric::GetMetric_uu(size_t ij)
{
    return Tensor3x3(g00_uu[ij], g01_uu[ij], g02_uu[ij],
                     g01_uu[ij], g11_uu[ij], g12_uu[ij],
                     g02_uu[ij], g12_uu[ij], g22_uu[ij]);
}
Tensor3x3 Metric::GetMetric_uu(const Coord &xy)
{
    if (grid.OutsideDomain(xy))
        return Tensor3x3(-1, 0, 0, 0, 1, 0, 0, 0, 1);
    else
    {
        Coord ij = grid.ij(xy);
        Tensor3x3 g_uu;
        g_uu[{0, 0}] = InterpolateArrayTo_ij(g00_uu, ij);
        g_uu[{1, 1}] = InterpolateArrayTo_ij(g11_uu, ij);
        g_uu[{2, 2}] = InterpolateArrayTo_ij(g22_uu, ij);
        g_uu[{0, 1}] = g_uu[{1, 0}] = InterpolateArrayTo_ij(g01_uu, ij);
        g_uu[{0, 2}] = g_uu[{2, 0}] = InterpolateArrayTo_ij(g02_uu, ij);
        g_uu[{1, 2}] = g_uu[{2, 1}] = InterpolateArrayTo_ij(g12_uu, ij);
        return g_uu;
    }
}

Tensor3x3 Metric::GetMinkowskiMetric_ll(size_t ij)
{
    return Tensor3x3(-1, 0, 0, 0, 1, 0, 0, 0, 1);
}
Tensor3x3 Metric::GetMinkowskiMetric_ll(const Coord &xy)
{
    return Tensor3x3(-1, 0, 0, 0, 1, 0, 0, 0, 1);
}
Tensor3x3 Metric::GetMinkowskiMetric_uu(size_t ij)
{
    return Tensor3x3(-1, 0, 0, 0, 1, 0, 0, 0, 1);
}
Tensor3x3 Metric::GetMinkowskiMetric_uu(const Coord &xy)
{
    return Tensor3x3(-1, 0, 0, 0, 1, 0, 0, 0, 1);
}

Tensor3x3x3 Metric::GetDerivMetric_lll(size_t ij)
{
    return Tensor3x3x3(d0_g00_lll[ij], d0_g01_lll[ij], d0_g02_lll[ij],
                       d0_g01_lll[ij], d0_g11_lll[ij], d0_g12_lll[ij],
                       d0_g02_lll[ij], d0_g12_lll[ij], d0_g22_lll[ij],

                       d1_g00_lll[ij], d1_g01_lll[ij], d1_g02_lll[ij],
                       d1_g01_lll[ij], d1_g11_lll[ij], d1_g12_lll[ij],
                       d1_g02_lll[ij], d1_g12_lll[ij], d1_g22_lll[ij],

                       d2_g00_lll[ij], d2_g01_lll[ij], d2_g02_lll[ij],
                       d2_g01_lll[ij], d2_g11_lll[ij], d2_g12_lll[ij],
                       d2_g02_lll[ij], d2_g12_lll[ij], d2_g22_lll[ij]);
}

Tensor3x3x3 Metric::GetDerivMetric_lll(const Coord &xy)
{
    if (grid.OutsideDomain(xy))
        return Tensor3x3x3(0);
    else
    {
        Coord ij = grid.ij(xy);
        Tensor3x3x3 dl_g_ll;
        dl_g_ll[{0, 0, 0}] = InterpolateArrayTo_ij(d0_g00_lll, ij);
        dl_g_ll[{0, 1, 1}] = InterpolateArrayTo_ij(d0_g11_lll, ij);
        dl_g_ll[{0, 2, 2}] = InterpolateArrayTo_ij(d0_g22_lll, ij);
        dl_g_ll[{0, 0, 1}] = dl_g_ll[{0, 1, 0}] = InterpolateArrayTo_ij(d0_g01_lll, ij);
        dl_g_ll[{0, 0, 2}] = dl_g_ll[{0, 2, 0}] = InterpolateArrayTo_ij(d0_g02_lll, ij);
        dl_g_ll[{0, 1, 2}] = dl_g_ll[{0, 2, 1}] = InterpolateArrayTo_ij(d0_g12_lll, ij);

        dl_g_ll[{1, 0, 0}] = InterpolateArrayTo_ij(d1_g00_lll, ij);
        dl_g_ll[{1, 1, 1}] = InterpolateArrayTo_ij(d1_g11_lll, ij);
        dl_g_ll[{1, 2, 2}] = InterpolateArrayTo_ij(d1_g22_lll, ij);
        dl_g_ll[{1, 0, 1}] = dl_g_ll[{1, 1, 0}] = InterpolateArrayTo_ij(d1_g01_lll, ij);
        dl_g_ll[{1, 0, 2}] = dl_g_ll[{1, 2, 0}] = InterpolateArrayTo_ij(d1_g02_lll, ij);
        dl_g_ll[{1, 1, 2}] = dl_g_ll[{1, 2, 1}] = InterpolateArrayTo_ij(d1_g12_lll, ij);

        dl_g_ll[{2, 0, 0}] = InterpolateArrayTo_ij(d2_g00_lll, ij);
        dl_g_ll[{2, 1, 1}] = InterpolateArrayTo_ij(d2_g11_lll, ij);
        dl_g_ll[{2, 2, 2}] = InterpolateArrayTo_ij(d2_g22_lll, ij);
        dl_g_ll[{2, 0, 1}] = dl_g_ll[{2, 1, 0}] = InterpolateArrayTo_ij(d2_g01_lll, ij);
        dl_g_ll[{2, 0, 2}] = dl_g_ll[{2, 2, 0}] = InterpolateArrayTo_ij(d2_g02_lll, ij);
        dl_g_ll[{2, 1, 2}] = dl_g_ll[{2, 2, 1}] = InterpolateArrayTo_ij(d2_g12_lll, ij);
        return dl_g_ll;
    }
}

Tensor3x3x3 Metric::GetDerivMetric_luu(size_t ij)
{
    return Tensor3x3x3(d0_g00_luu[ij], d0_g01_luu[ij], d0_g02_luu[ij],
                       d0_g01_luu[ij], d0_g11_luu[ij], d0_g12_luu[ij],
                       d0_g02_luu[ij], d0_g12_luu[ij], d0_g22_luu[ij],

                       d1_g00_luu[ij], d1_g01_luu[ij], d1_g02_luu[ij],
                       d1_g01_luu[ij], d1_g11_luu[ij], d1_g12_luu[ij],
                       d1_g02_luu[ij], d1_g12_luu[ij], d1_g22_luu[ij],

                       d2_g00_luu[ij], d2_g01_luu[ij], d2_g02_luu[ij],
                       d2_g01_luu[ij], d2_g11_luu[ij], d2_g12_luu[ij],
                       d2_g02_luu[ij], d2_g12_luu[ij], d2_g22_luu[ij]);
}

Tensor3x3x3 Metric::GetDerivMetric_luu(const Coord &xy)
{
    if (grid.OutsideDomain(xy))
        return Tensor3x3x3(0);
    else
    {
        Coord ij = grid.ij(xy);
        Tensor3x3x3 dl_g_uu;
        dl_g_uu[{0, 0, 0}] = InterpolateArrayTo_ij(d0_g00_luu, ij);
        dl_g_uu[{0, 1, 1}] = InterpolateArrayTo_ij(d0_g11_luu, ij);
        dl_g_uu[{0, 2, 2}] = InterpolateArrayTo_ij(d0_g22_luu, ij);
        dl_g_uu[{0, 0, 1}] = dl_g_uu[{0, 1, 0}] = InterpolateArrayTo_ij(d0_g01_luu, ij);
        dl_g_uu[{0, 0, 2}] = dl_g_uu[{0, 2, 0}] = InterpolateArrayTo_ij(d0_g02_luu, ij);
        dl_g_uu[{0, 1, 2}] = dl_g_uu[{0, 2, 1}] = InterpolateArrayTo_ij(d0_g12_luu, ij);

        dl_g_uu[{1, 0, 0}] = InterpolateArrayTo_ij(d1_g00_luu, ij);
        dl_g_uu[{1, 1, 1}] = InterpolateArrayTo_ij(d1_g11_luu, ij);
        dl_g_uu[{1, 2, 2}] = InterpolateArrayTo_ij(d1_g22_luu, ij);
        dl_g_uu[{1, 0, 1}] = dl_g_uu[{1, 1, 0}] = InterpolateArrayTo_ij(d1_g01_luu, ij);
        dl_g_uu[{1, 0, 2}] = dl_g_uu[{1, 2, 0}] = InterpolateArrayTo_ij(d1_g02_luu, ij);
        dl_g_uu[{1, 1, 2}] = dl_g_uu[{1, 2, 1}] = InterpolateArrayTo_ij(d1_g12_luu, ij);

        dl_g_uu[{2, 0, 0}] = InterpolateArrayTo_ij(d2_g00_luu, ij);
        dl_g_uu[{2, 1, 1}] = InterpolateArrayTo_ij(d2_g11_luu, ij);
        dl_g_uu[{2, 2, 2}] = InterpolateArrayTo_ij(d2_g22_luu, ij);
        dl_g_uu[{2, 0, 1}] = dl_g_uu[{2, 1, 0}] = InterpolateArrayTo_ij(d2_g01_luu, ij);
        dl_g_uu[{2, 0, 2}] = dl_g_uu[{2, 2, 0}] = InterpolateArrayTo_ij(d2_g02_luu, ij);
        dl_g_uu[{2, 1, 2}] = dl_g_uu[{2, 2, 1}] = InterpolateArrayTo_ij(d2_g12_luu, ij);
        return dl_g_uu;
    }
}

Tensor3x3 Metric::GetTetrad(size_t ij)
{
    return Tensor3x3(tetrad00_ul[ij], tetrad01_ul[ij], tetrad02_ul[ij],
                     tetrad10_ul[ij], tetrad11_ul[ij], tetrad12_ul[ij],
                     tetrad20_ul[ij], tetrad21_ul[ij], tetrad22_ul[ij]);
}
Tensor3x3 Metric::GetTetrad(const Coord &xy)
{
    if (grid.OutsideDomain(xy))
        return Tensor3x3(1, 0, 0, 0, 1, 0, 0, 0, 1);
    else
    {
        Coord ij = grid.ij(xy);
        return Tensor3x3(InterpolateArrayTo_ij(tetrad00_ul, ij), InterpolateArrayTo_ij(tetrad01_ul, ij), InterpolateArrayTo_ij(tetrad02_ul, ij),
                         InterpolateArrayTo_ij(tetrad10_ul, ij), InterpolateArrayTo_ij(tetrad11_ul, ij), InterpolateArrayTo_ij(tetrad12_ul, ij),
                         InterpolateArrayTo_ij(tetrad20_ul, ij), InterpolateArrayTo_ij(tetrad21_ul, ij), InterpolateArrayTo_ij(tetrad22_ul, ij));
    }
}

Tensor3x3 Metric::GetTetradInverse(size_t ij)
{
    return Tensor3x3(tetrad00_lu[ij], tetrad01_lu[ij], tetrad02_lu[ij],
                     tetrad10_lu[ij], tetrad11_lu[ij], tetrad12_lu[ij],
                     tetrad20_lu[ij], tetrad21_lu[ij], tetrad22_lu[ij]);
}
Tensor3x3 Metric::GetTetradInverse(const Coord &xy)
{
    if (grid.OutsideDomain(xy))
        return Tensor3x3(1, 0, 0, 0, 1, 0, 0, 0, 1);
    else
    {
        Coord ij = grid.ij(xy);
        return Tensor3x3(InterpolateArrayTo_ij(tetrad00_lu, ij), InterpolateArrayTo_ij(tetrad01_lu, ij), InterpolateArrayTo_ij(tetrad02_lu, ij),
                         InterpolateArrayTo_ij(tetrad10_lu, ij), InterpolateArrayTo_ij(tetrad11_lu, ij), InterpolateArrayTo_ij(tetrad12_lu, ij),
                         InterpolateArrayTo_ij(tetrad20_lu, ij), InterpolateArrayTo_ij(tetrad21_lu, ij), InterpolateArrayTo_ij(tetrad22_lu, ij));
    }
}

// ADM getters:
double Metric::GetAlpha(size_t ij)
{
    return alpha[ij];
}
double Metric::GetAlpha(const Coord &xy)
{
    if (grid.OutsideDomain(xy))
        return 1.0;
    else
    {
        Coord ij = grid.ij(xy);
        return InterpolateArrayTo_ij(alpha, ij);
    }
}

Tensor2 Metric::GetBeta_u(size_t ij)
{
    return Tensor2(beta1_u[ij], beta2_u[ij]);
}
Tensor2 Metric::GetBeta_u(const Coord &xy)
{
    if (grid.OutsideDomain(xy))
        return Tensor2(0.0);
    else
    {
        Coord ij = grid.ij(xy);
        return Tensor2(InterpolateArrayTo_ij(beta1_u, ij), InterpolateArrayTo_ij(beta2_u, ij));
    }
}

Tensor2 Metric::GetBeta_l(size_t ij)
{
    return Tensor2(beta1_l[ij], beta2_l[ij]);
}
Tensor2 Metric::GetBeta_l(const Coord &xy)
{
    if (grid.OutsideDomain(xy))
        return Tensor2(0.0);
    else
    {
        Coord ij = grid.ij(xy);
        return Tensor2(InterpolateArrayTo_ij(beta1_l, ij), InterpolateArrayTo_ij(beta2_l, ij));
    }
}

Tensor2x2 Metric::GetGamma_ll(size_t ij)
{
    return Tensor2x2(gamma11_ll[ij], gamma12_ll[ij],
                     gamma12_ll[ij], gamma22_ll[ij]);
}
Tensor2x2 Metric::GetGamma_ll(const Coord &xy)
{
    if (grid.OutsideDomain(xy))
        return Tensor2x2(1, 0, 0, 1);
    else
    {
        Coord ij = grid.ij(xy);
        Tensor2x2 gamma_ll;
        gamma_ll[{1, 1}] = InterpolateArrayTo_ij(gamma11_ll, ij);
        gamma_ll[{2, 2}] = InterpolateArrayTo_ij(gamma22_ll, ij);
        gamma_ll[{1, 2}] = gamma_ll[{2, 1}] = InterpolateArrayTo_ij(gamma12_ll, ij);
        return gamma_ll;
    }
}

Tensor2x2 Metric::GetGamma_uu(size_t ij)
{
    return Tensor2x2(gamma11_uu[ij], gamma12_uu[ij],
                     gamma12_uu[ij], gamma22_uu[ij]);
}
Tensor2x2 Metric::GetGamma_uu(const Coord &xy)
{
    if (grid.OutsideDomain(xy))
        return Tensor2x2(1, 0, 0, 1);
    else
    {
        Coord ij = grid.ij(xy);
        Tensor2x2 gamma_uu;
        gamma_uu[{1, 1}] = InterpolateArrayTo_ij(gamma11_uu, ij);
        gamma_uu[{2, 2}] = InterpolateArrayTo_ij(gamma22_uu, ij);
        gamma_uu[{1, 2}] = gamma_uu[{2, 1}] = InterpolateArrayTo_ij(gamma12_uu, ij);
        return gamma_uu;
    }
}

Tensor2x2 Metric::GetMinkowskiGamma_ll(size_t ij)
{
    return Tensor2x2(1, 0, 0, 1);
}
Tensor2x2 Metric::GetMinkowskiGamma_ll(const Coord &xy)
{
    return Tensor2x2(1, 0, 0, 1);
}

Tensor2x2 Metric::GetMinkowskiGamma_uu(size_t ij)
{
    return Tensor2x2(1, 0, 0, 1);
}
Tensor2x2 Metric::GetMinkowskiGamma_uu(const Coord &xy)
{
    return Tensor2x2(1, 0, 0, 1);
}

Tensor2 Metric::GetDerivAlpha_l(size_t ij)
{
    return Tensor2(d1_alpha_l[ij], d2_alpha_l[ij]);
}
Tensor2 Metric::GetDerivAlpha_l(const Coord &xy)
{
    if (grid.OutsideDomain(xy))
        return Tensor2(0.0);
    else
    {
        Coord ij = grid.ij(xy);
        return Tensor2(InterpolateArrayTo_ij(d1_alpha_l, ij), InterpolateArrayTo_ij(d2_alpha_l, ij));
    }
}

Tensor2x2 Metric::GetDerivBeta_lu(size_t ij)
{
    return Tensor2x2(d1_beta1_lu[ij], d1_beta2_lu[ij],
                     d2_beta1_lu[ij], d2_beta2_lu[ij]);
}
Tensor2x2 Metric::GetDerivBeta_lu(const Coord &xy)
{
    if (grid.OutsideDomain(xy))
        return Tensor2x2(0.0);
    else
    {
        Coord ij = grid.ij(xy);
        return Tensor2x2(InterpolateArrayTo_ij(d1_beta1_lu, ij), InterpolateArrayTo_ij(d1_beta2_lu, ij),
                         InterpolateArrayTo_ij(d2_beta1_lu, ij), InterpolateArrayTo_ij(d2_beta2_lu, ij));
    }
}

Tensor2x2 Metric::GetDerivBeta_ll(size_t ij)
{
    return Tensor2x2(d1_beta1_ll[ij], d1_beta2_ll[ij],
                     d2_beta1_ll[ij], d2_beta2_ll[ij]);
}
Tensor2x2 Metric::GetDerivBeta_ll(const Coord &xy)
{
    if (grid.OutsideDomain(xy))
        return Tensor2x2(0.0);
    else
    {
        Coord ij = grid.ij(xy);
        return Tensor2x2(InterpolateArrayTo_ij(d1_beta1_ll, ij), InterpolateArrayTo_ij(d1_beta2_ll, ij),
                         InterpolateArrayTo_ij(d2_beta1_ll, ij), InterpolateArrayTo_ij(d2_beta2_ll, ij));
    }
}

Tensor2x2x2 Metric::GetDerivGamma_lll(size_t ij)
{
    return Tensor2x2x2(d1_gamma11_lll[ij], d1_gamma12_lll[ij],
                       d1_gamma12_lll[ij], d1_gamma22_lll[ij],
                       d2_gamma11_lll[ij], d2_gamma12_lll[ij],
                       d2_gamma12_lll[ij], d2_gamma22_lll[ij]);
}
Tensor2x2x2 Metric::GetDerivGamma_lll(const Coord &xy)
{
    if (grid.OutsideDomain(xy))
        return Tensor2x2x2(0.0);
    else
    {
        Coord ij = grid.ij(xy);
        Tensor2x2x2 dGamma_lll;
        dGamma_lll[{1, 1, 1}] = InterpolateArrayTo_ij(d1_gamma11_lll, ij);
        dGamma_lll[{1, 2, 2}] = InterpolateArrayTo_ij(d1_gamma22_lll, ij);
        dGamma_lll[{1, 1, 2}] = dGamma_lll[{1, 2, 1}] = InterpolateArrayTo_ij(d1_gamma12_lll, ij);
        dGamma_lll[{2, 1, 1}] = InterpolateArrayTo_ij(d2_gamma11_lll, ij);
        dGamma_lll[{2, 2, 2}] = InterpolateArrayTo_ij(d2_gamma22_lll, ij);
        dGamma_lll[{2, 1, 2}] = dGamma_lll[{2, 2, 1}] = InterpolateArrayTo_ij(d2_gamma12_lll, ij);
        return dGamma_lll;
    }
}

Tensor2x2x2 Metric::GetDerivGamma_luu(size_t ij)
{
    return Tensor2x2x2(d1_gamma11_luu[ij], d1_gamma12_luu[ij],
                       d1_gamma12_luu[ij], d1_gamma22_luu[ij],
                       d2_gamma11_luu[ij], d2_gamma12_luu[ij],
                       d2_gamma12_luu[ij], d2_gamma22_luu[ij]);
}
Tensor2x2x2 Metric::GetDerivGamma_luu(const Coord &xy)
{
    if (grid.OutsideDomain(xy))
        return Tensor2x2x2(0.0);
    else
    {
        Coord ij = grid.ij(xy);
        Tensor2x2x2 dGamma_luu;
        dGamma_luu[{1, 1, 1}] = InterpolateArrayTo_ij(d1_gamma11_luu, ij);
        dGamma_luu[{1, 2, 2}] = InterpolateArrayTo_ij(d1_gamma22_luu, ij);
        dGamma_luu[{1, 1, 2}] = dGamma_luu[{1, 2, 1}] = InterpolateArrayTo_ij(d1_gamma12_luu, ij);
        dGamma_luu[{2, 1, 1}] = InterpolateArrayTo_ij(d2_gamma11_luu, ij);
        dGamma_luu[{2, 2, 2}] = InterpolateArrayTo_ij(d2_gamma22_luu, ij);
        dGamma_luu[{2, 1, 2}] = dGamma_luu[{2, 2, 1}] = InterpolateArrayTo_ij(d2_gamma12_luu, ij);
        return dGamma_luu;
    }
}

Tensor2x2 Metric::GetExtrCurv_ll(size_t ij)
{
    return Tensor2x2(K11_ll[ij], K12_ll[ij],
                     K12_ll[ij], K22_ll[ij]);
}
Tensor2x2 Metric::GetExtrCurv_ll(const Coord &xy)
{
    if (grid.OutsideDomain(xy))
        return Tensor2x2(0.0);
    else
    {
        Coord ij = grid.ij(xy);
        Tensor2x2 K_ll;
        K_ll[{1, 1}] = InterpolateArrayTo_ij(K11_ll, ij);
        K_ll[{2, 2}] = InterpolateArrayTo_ij(K22_ll, ij);
        K_ll[{1, 2}] = K_ll[{2, 1}] = InterpolateArrayTo_ij(K12_ll, ij);
        return K_ll;
    }
}