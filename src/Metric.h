#ifndef __INCLUDE_GUARD_Metric_h__
#define __INCLUDE_GUARD_Metric_h__
#include "Grid.h"           // Numerical Grid and mapping to physical domain.
#include "Interpolation.h"  // Basic interpolation schemes.



class Metric
{
public:
    // Grid Data:
    Grid& grid;
    double m = 1.0;
    double a = 0.0;

// protected:
public:
    // Metrik Data:
    DoubleBuffer g00_ll;    DoubleBuffer g00_uu;
    DoubleBuffer g01_ll;    DoubleBuffer g01_uu;
    DoubleBuffer g02_ll;    DoubleBuffer g02_uu;
    DoubleBuffer g11_ll;    DoubleBuffer g11_uu;
    DoubleBuffer g12_ll;    DoubleBuffer g12_uu;
    DoubleBuffer g22_ll;    DoubleBuffer g22_uu;
    DoubleBuffer d0_g00_lll;    DoubleBuffer d1_g00_lll;    DoubleBuffer d2_g00_lll;
    DoubleBuffer d0_g01_lll;    DoubleBuffer d1_g01_lll;    DoubleBuffer d2_g01_lll;
    DoubleBuffer d0_g02_lll;    DoubleBuffer d1_g02_lll;    DoubleBuffer d2_g02_lll;
    DoubleBuffer d0_g11_lll;    DoubleBuffer d1_g11_lll;    DoubleBuffer d2_g11_lll;
    DoubleBuffer d0_g12_lll;    DoubleBuffer d1_g12_lll;    DoubleBuffer d2_g12_lll;
    DoubleBuffer d0_g22_lll;    DoubleBuffer d1_g22_lll;    DoubleBuffer d2_g22_lll;
    DoubleBuffer d0_g00_luu;    DoubleBuffer d1_g00_luu;    DoubleBuffer d2_g00_luu;
    DoubleBuffer d0_g01_luu;    DoubleBuffer d1_g01_luu;    DoubleBuffer d2_g01_luu;
    DoubleBuffer d0_g02_luu;    DoubleBuffer d1_g02_luu;    DoubleBuffer d2_g02_luu;
    DoubleBuffer d0_g11_luu;    DoubleBuffer d1_g11_luu;    DoubleBuffer d2_g11_luu;
    DoubleBuffer d0_g12_luu;    DoubleBuffer d1_g12_luu;    DoubleBuffer d2_g12_luu;
    DoubleBuffer d0_g22_luu;    DoubleBuffer d1_g22_luu;    DoubleBuffer d2_g22_luu;
    // ADM Data:
    DoubleBuffer alpha;
    DoubleBuffer beta1_u;    DoubleBuffer beta1_l;
    DoubleBuffer beta2_u;    DoubleBuffer beta2_l;
    DoubleBuffer gamma11_ll;    DoubleBuffer gamma11_uu;
    DoubleBuffer gamma12_ll;    DoubleBuffer gamma12_uu;
    DoubleBuffer gamma22_ll;    DoubleBuffer gamma22_uu;
    DoubleBuffer d1_alpha_l;
    DoubleBuffer d2_alpha_l;
    DoubleBuffer d1_beta1_lu;    DoubleBuffer d1_beta2_lu;
    DoubleBuffer d2_beta1_lu;    DoubleBuffer d2_beta2_lu;
    DoubleBuffer d1_beta1_ll;    DoubleBuffer d1_beta2_ll;
    DoubleBuffer d2_beta1_ll;    DoubleBuffer d2_beta2_ll;
    DoubleBuffer d1_gamma11_lll;    DoubleBuffer d2_gamma11_lll;
    DoubleBuffer d1_gamma12_lll;    DoubleBuffer d2_gamma12_lll;
    DoubleBuffer d1_gamma22_lll;    DoubleBuffer d2_gamma22_lll;
    DoubleBuffer d1_gamma11_luu;    DoubleBuffer d2_gamma11_luu;
    DoubleBuffer d1_gamma12_luu;    DoubleBuffer d2_gamma12_luu;
    DoubleBuffer d1_gamma22_luu;    DoubleBuffer d2_gamma22_luu;
    DoubleBuffer K11_ll;
    DoubleBuffer K12_ll;
    DoubleBuffer K22_ll;
    // Tetrad Data:
    DoubleBuffer tetrad00_ul;    DoubleBuffer tetrad01_ul;    DoubleBuffer tetrad02_ul;
    DoubleBuffer tetrad10_ul;    DoubleBuffer tetrad11_ul;    DoubleBuffer tetrad12_ul;
    DoubleBuffer tetrad20_ul;    DoubleBuffer tetrad21_ul;    DoubleBuffer tetrad22_ul;
    DoubleBuffer tetrad00_lu;    DoubleBuffer tetrad01_lu;    DoubleBuffer tetrad02_lu;
    DoubleBuffer tetrad10_lu;    DoubleBuffer tetrad11_lu;    DoubleBuffer tetrad12_lu;
    DoubleBuffer tetrad20_lu;    DoubleBuffer tetrad21_lu;    DoubleBuffer tetrad22_lu;

public:
    // Constructors/Destructor:
    Metric(Grid& grid_, double m_, double a_);

    virtual std::string Name();

    // Initialization:
    virtual Tensor3x3 MetricFunction(const Coord& xy);
    void InitializeMetricOnGrid();
    void InitializeBoostedTetradOnGrid();
    template<int k>
    Tensor3x3 MetricDeriv(const Coord& xy);
    template<int k>
    Tensor3x3 InverseMetricDeriv(const Coord& xy);
    void InitializeMetricDerivativesOnGrid();
    void InitializeAdmComponentsOnGrid();
    double InterpolateArrayTo_ij(const DoubleBuffer& array, const Coord& ij);
    double InterpolateArrayTo_ij(const DoubleBuffer& array, double i, double j);

public:
    // Boolean checks:
    virtual bool InsideBH(const Coord& xy);

    // Tensor getters:
    Tensor3 uEulObs(size_t ij);
    Tensor3 uEulObs(const Coord& xy);
    Tensor3x3 GetMetric_ll(size_t ij);
    Tensor3x3 GetMetric_ll(const Coord& xy);
    Tensor3x3 GetMetric_uu(size_t ij);
    Tensor3x3 GetMetric_uu(const Coord& xy);
    Tensor3x3 GetMinkowskiMetric_ll(size_t ij);
    Tensor3x3 GetMinkowskiMetric_ll(const Coord& xy);
    Tensor3x3 GetMinkowskiMetric_uu(size_t ij);
    Tensor3x3 GetMinkowskiMetric_uu(const Coord& xy);
    Tensor3x3x3 GetDerivMetric_lll(size_t ij);
    Tensor3x3x3 GetDerivMetric_lll(const Coord& xy);
    Tensor3x3x3 GetDerivMetric_luu(size_t ij);
    Tensor3x3x3 GetDerivMetric_luu(const Coord& xy);
    Tensor3x3 GetTetrad(size_t ij);
    Tensor3x3 GetTetrad(const Coord& xy);
    Tensor3x3 GetTetradInverse(size_t ij);
    Tensor3x3 GetTetradInverse(const Coord& xy);

    // ADM getters:
    double GetAlpha(size_t ij);
    double GetAlpha(const Coord& xy);
    Tensor2 GetBeta_u(size_t ij);
    Tensor2 GetBeta_u(const Coord& xy);
    Tensor2 GetBeta_l(size_t ij);
    Tensor2 GetBeta_l(const Coord& xy);
    Tensor2x2 GetGamma_ll(size_t ij);
    Tensor2x2 GetGamma_ll(const Coord& xy);
    Tensor2x2 GetGamma_uu(size_t ij);
    Tensor2x2 GetGamma_uu(const Coord& xy);
    Tensor2x2 GetMinkowskiGamma_ll(size_t ij);
    Tensor2x2 GetMinkowskiGamma_ll(const Coord& xy);
    Tensor2x2 GetMinkowskiGamma_uu(size_t ij);
    Tensor2x2 GetMinkowskiGamma_uu(const Coord& xy);

    Tensor2 GetDerivAlpha_l(size_t ij);
    Tensor2 GetDerivAlpha_l(const Coord& xy);
    Tensor2x2 GetDerivBeta_lu(size_t ij);
    Tensor2x2 GetDerivBeta_lu(const Coord& xy);
    Tensor2x2 GetDerivBeta_ll(size_t ij);
    Tensor2x2 GetDerivBeta_ll(const Coord& xy);
    Tensor2x2x2 GetDerivGamma_lll(size_t ij);
    Tensor2x2x2 GetDerivGamma_lll(const Coord& xy);
    Tensor2x2x2 GetDerivGamma_luu(size_t ij);
    Tensor2x2x2 GetDerivGamma_luu(const Coord& xy);
    Tensor2x2 GetExtrCurv_ll(size_t ij);
    Tensor2x2 GetExtrCurv_ll(const Coord& xy);
};
#endif //__INCLUDE_GUARD_Metric_h__