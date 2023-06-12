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
    RealBuffer g00_ll;    RealBuffer g00_uu;
    RealBuffer g01_ll;    RealBuffer g01_uu;
    RealBuffer g02_ll;    RealBuffer g02_uu;
    RealBuffer g11_ll;    RealBuffer g11_uu;
    RealBuffer g12_ll;    RealBuffer g12_uu;
    RealBuffer g22_ll;    RealBuffer g22_uu;
    RealBuffer d0_g00_lll;    RealBuffer d1_g00_lll;    RealBuffer d2_g00_lll;
    RealBuffer d0_g01_lll;    RealBuffer d1_g01_lll;    RealBuffer d2_g01_lll;
    RealBuffer d0_g02_lll;    RealBuffer d1_g02_lll;    RealBuffer d2_g02_lll;
    RealBuffer d0_g11_lll;    RealBuffer d1_g11_lll;    RealBuffer d2_g11_lll;
    RealBuffer d0_g12_lll;    RealBuffer d1_g12_lll;    RealBuffer d2_g12_lll;
    RealBuffer d0_g22_lll;    RealBuffer d1_g22_lll;    RealBuffer d2_g22_lll;
    RealBuffer d0_g00_luu;    RealBuffer d1_g00_luu;    RealBuffer d2_g00_luu;
    RealBuffer d0_g01_luu;    RealBuffer d1_g01_luu;    RealBuffer d2_g01_luu;
    RealBuffer d0_g02_luu;    RealBuffer d1_g02_luu;    RealBuffer d2_g02_luu;
    RealBuffer d0_g11_luu;    RealBuffer d1_g11_luu;    RealBuffer d2_g11_luu;
    RealBuffer d0_g12_luu;    RealBuffer d1_g12_luu;    RealBuffer d2_g12_luu;
    RealBuffer d0_g22_luu;    RealBuffer d1_g22_luu;    RealBuffer d2_g22_luu;
    // ADM Data:
    RealBuffer alpha;
    RealBuffer beta1_u;    RealBuffer beta1_l;
    RealBuffer beta2_u;    RealBuffer beta2_l;
    RealBuffer gamma11_ll;    RealBuffer gamma11_uu;
    RealBuffer gamma12_ll;    RealBuffer gamma12_uu;
    RealBuffer gamma22_ll;    RealBuffer gamma22_uu;
    RealBuffer d1_alpha_l;
    RealBuffer d2_alpha_l;
    RealBuffer d1_beta1_lu;    RealBuffer d1_beta2_lu;
    RealBuffer d2_beta1_lu;    RealBuffer d2_beta2_lu;
    RealBuffer d1_beta1_ll;    RealBuffer d1_beta2_ll;
    RealBuffer d2_beta1_ll;    RealBuffer d2_beta2_ll;
    RealBuffer d1_gamma11_lll;    RealBuffer d2_gamma11_lll;
    RealBuffer d1_gamma12_lll;    RealBuffer d2_gamma12_lll;
    RealBuffer d1_gamma22_lll;    RealBuffer d2_gamma22_lll;
    RealBuffer d1_gamma11_luu;    RealBuffer d2_gamma11_luu;
    RealBuffer d1_gamma12_luu;    RealBuffer d2_gamma12_luu;
    RealBuffer d1_gamma22_luu;    RealBuffer d2_gamma22_luu;
    RealBuffer K11_ll;
    RealBuffer K12_ll;
    RealBuffer K22_ll;
    // Tetrad Data:
    RealBuffer tetrad00_ul;    RealBuffer tetrad01_ul;    RealBuffer tetrad02_ul;
    RealBuffer tetrad10_ul;    RealBuffer tetrad11_ul;    RealBuffer tetrad12_ul;
    RealBuffer tetrad20_ul;    RealBuffer tetrad21_ul;    RealBuffer tetrad22_ul;
    RealBuffer tetrad00_lu;    RealBuffer tetrad01_lu;    RealBuffer tetrad02_lu;
    RealBuffer tetrad10_lu;    RealBuffer tetrad11_lu;    RealBuffer tetrad12_lu;
    RealBuffer tetrad20_lu;    RealBuffer tetrad21_lu;    RealBuffer tetrad22_lu;

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
    double InterpolateArrayTo_ij(const RealBuffer& array, const Coord& ij);
    double InterpolateArrayTo_ij(const RealBuffer& array, double i, double j);

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