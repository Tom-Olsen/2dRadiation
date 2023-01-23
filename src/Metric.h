#ifndef __INCLUDE_GUARD_Metric_h__
#define __INCLUDE_GUARD_Metric_h__
#include <string>
#include "ControlFlow.hh"       // used for template arguments
#include "Utility.hh"           // basic utility functions
#include "TensorTypes.hh"       // simple containers for rank 1-3 tensors
#include "Grid.h"             // underlying numerical Grid
#include "Interpolation.hh"     // biilinear interpolation



template<class Coord>
class Metric
{
public:
    // Grid Data:
    Grid<Coord>& grid;
    double m = 1.0;
    double a = 0.0;

// protected:
public:
    // Metrik Data:
    double* g00_ll;    double* g00_uu;
    double* g01_ll;    double* g01_uu;
    double* g02_ll;    double* g02_uu;
    double* g11_ll;    double* g11_uu;
    double* g12_ll;    double* g12_uu;
    double* g22_ll;    double* g22_uu;
    double* d0_g00_lll;    double* d1_g00_lll;    double* d2_g00_lll;
    double* d0_g01_lll;    double* d1_g01_lll;    double* d2_g01_lll;
    double* d0_g02_lll;    double* d1_g02_lll;    double* d2_g02_lll;
    double* d0_g11_lll;    double* d1_g11_lll;    double* d2_g11_lll;
    double* d0_g12_lll;    double* d1_g12_lll;    double* d2_g12_lll;
    double* d0_g22_lll;    double* d1_g22_lll;    double* d2_g22_lll;
    double* d0_g00_luu;    double* d1_g00_luu;    double* d2_g00_luu;
    double* d0_g01_luu;    double* d1_g01_luu;    double* d2_g01_luu;
    double* d0_g02_luu;    double* d1_g02_luu;    double* d2_g02_luu;
    double* d0_g11_luu;    double* d1_g11_luu;    double* d2_g11_luu;
    double* d0_g12_luu;    double* d1_g12_luu;    double* d2_g12_luu;
    double* d0_g22_luu;    double* d1_g22_luu;    double* d2_g22_luu;
    // ADM Data:
    double* alpha;
    double* beta1_u;    double* beta1_l;
    double* beta2_u;    double* beta2_l;
    double* gamma11_ll;    double* gamma11_uu;
    double* gamma12_ll;    double* gamma12_uu;
    double* gamma22_ll;    double* gamma22_uu;
    double* d1_alpha_l;
    double* d2_alpha_l;
    double* d1_beta1_lu;    double* d1_beta2_lu;
    double* d2_beta1_lu;    double* d2_beta2_lu;
    double* d1_beta1_ll;    double* d1_beta2_ll;
    double* d2_beta1_ll;    double* d2_beta2_ll;
    double* d1_gamma11_lll;    double* d2_gamma11_lll;
    double* d1_gamma12_lll;    double* d2_gamma12_lll;
    double* d1_gamma22_lll;    double* d2_gamma22_lll;
    double* d1_gamma11_luu;    double* d2_gamma11_luu;
    double* d1_gamma12_luu;    double* d2_gamma12_luu;
    double* d1_gamma22_luu;    double* d2_gamma22_luu;
    double* K11_ll;
    double* K12_ll;
    double* K22_ll;
    // Tetrad Data:
    double* tetrad00_ul;    double* tetrad01_ul;    double* tetrad02_ul;
    double* tetrad10_ul;    double* tetrad11_ul;    double* tetrad12_ul;
    double* tetrad20_ul;    double* tetrad21_ul;    double* tetrad22_ul;

public:
    // Constructors/Destructor:
    Metric(Grid<Coord>& grid_, double m_, double a_);
    //Metric(const Metric& metric);
    ~Metric();

    virtual std::string Name();

    // Initialization:
    virtual Tensor3x3<Coord,LF> MetricFunction(Coordinate2<Coord> x);
    void InitializeMetricOnGrid();
    void InitializeBoostedTetradOnGrid();
    template<int k>
    Tensor3x3<Coord,LF> MetricDeriv(Coordinate2<Coord> x);
    template<int k>
    Tensor3x3<Coord,LF> InverseMetricDeriv(Coordinate2<Coord> x);
    void InitializeMetricDerivativesOnGrid();
    void InitializeAdmComponentsOnGrid();
    double InterpolateArrayTo_ij(double* array, double i, double j);

public:
    // Boolean checks:
    virtual bool InsideBH(const int i, const int j);
    virtual bool InsideBH(const Coordinate2<Coord>& x12);

    // Tensor getters:
    Tensor3<Coord,LF> uEulObs(const int ij);
    Tensor3<Coord,LF> uEulObs(const Coordinate2<Coord>& x12);
    Tensor3x3<Coord,LF> GetMetric_ll(const int ij);
    Tensor3x3<Coord,LF> GetMetric_ll(const Coordinate2<Coord>& x12);
    Tensor3x3<Coord,LF> GetMetric_uu(const int ij);
    Tensor3x3<Coord,LF> GetMetric_uu(const Coordinate2<Coord>& x12);
    Tensor3x3<Coord,IF> GetMinkowskiMetric_ll(const int ij);
    Tensor3x3<Coord,IF> GetMinkowskiMetric_ll(const Coordinate2<Coord>& x12);
    Tensor3x3<Coord,IF> GetMinkowskiMetric_uu(const int ij);
    Tensor3x3<Coord,IF> GetMinkowskiMetric_uu(const Coordinate2<Coord>& x12);
    Tensor3x3x3<Coord,LF> GetDerivMetric_lll(const int ij);
    Tensor3x3x3<Coord,LF> GetDerivMetric_lll(const Coordinate2<Coord>& x12);
    Tensor3x3x3<Coord,LF> GetDerivMetric_luu(const int ij);
    Tensor3x3x3<Coord,LF> GetDerivMetric_luu(const Coordinate2<Coord>& x12);
    Tensor3x3<Coord,Tetrad> GetTetrad(const int ij);
    Tensor3x3<Coord,Tetrad> GetTetrad(const Coordinate2<Coord>& x12);
    // ADM getters:
    double GetAlpha(const int ij);
    double GetAlpha(const Coordinate2<Coord>& x12);
    Tensor2<Coord,LF> GetBeta_u(const int ij);
    Tensor2<Coord,LF> GetBeta_u(const Coordinate2<Coord>& x12);
    Tensor2<Coord,LF> GetBeta_l(const int ij);
    Tensor2<Coord,LF> GetBeta_l(const Coordinate2<Coord>& x12);
    Tensor2x2<Coord,LF> GetGamma_ll(const int ij);
    Tensor2x2<Coord,LF> GetGamma_ll(const Coordinate2<Coord>& x12);
    Tensor2x2<Coord,LF> GetGamma_uu(const int ij);
    Tensor2x2<Coord,LF> GetGamma_uu(const Coordinate2<Coord>& x12);
    Tensor2x2<Coord,IF> GetMinkowskiGamma_ll(const int ij);
    Tensor2x2<Coord,IF> GetMinkowskiGamma_ll(const Coordinate2<Coord>& x12);
    Tensor2x2<Coord,IF> GetMinkowskiGamma_uu(const int ij);
    Tensor2x2<Coord,IF> GetMinkowskiGamma_uu(const Coordinate2<Coord>& x12);

    Tensor2<Coord,LF> GetDerivAlpha_l(const int ij);
    Tensor2<Coord,LF> GetDerivAlpha_l(const Coordinate2<Coord>& x12);
    Tensor2x2<Coord,LF> GetDerivBeta_lu(const int ij);
    Tensor2x2<Coord,LF> GetDerivBeta_lu(const Coordinate2<Coord>& x12);
    Tensor2x2<Coord,LF> GetDerivBeta_ll(const int ij);
    Tensor2x2<Coord,LF> GetDerivBeta_ll(const Coordinate2<Coord>& x12);
    Tensor2x2x2<Coord,LF> GetDerivGamma_lll(const int ij);
    Tensor2x2x2<Coord,LF> GetDerivGamma_lll(const Coordinate2<Coord>& x12);
    Tensor2x2x2<Coord,LF> GetDerivGamma_luu(const int ij);
    Tensor2x2x2<Coord,LF> GetDerivGamma_luu(const Coordinate2<Coord>& x12);
    Tensor2x2<Coord,LF> GetExtrCurv_ll(const int ij);
    Tensor2x2<Coord,LF> GetExtrCurv_ll(const Coordinate2<Coord>& x12);
};
#endif //__INCLUDE_GUARD_Metric_h__