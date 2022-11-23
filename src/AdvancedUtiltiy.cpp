#include "AdvancedUtiltiy.h"


template<class Coord, class FrameIn, class FrameOut>
Tensor2<Coord,FrameOut> Vec2ObservedByEulObs(const Tensor3<Coord,FrameIn>& u, const Coordinate2<Coord>& x, Metric2D<Coord>& metric)
{
    if constexpr(std::is_same<FrameIn,IF>::value && std::is_same<FrameOut,IF>::value)
    {
        double alpha = metric.GetAlpha(x);
        return Tensor2<Coord,IF>(u[1]/alpha,u[2]/alpha);
    }
    if constexpr(std::is_same<FrameIn,IF>::value && std::is_same<FrameOut,LF>::value)
    {
        double alpha = metric.GetAlpha(x);
        Tensor2<Coord,IF> v(u[1]/alpha,u[2]/alpha);
        return v.template Transform<LF>(metric.GetTetrad(x));
    }
    if constexpr(std::is_same<FrameIn,LF>::value && std::is_same<FrameOut,IF>::value)
    {
        double alpha = metric.GetAlpha(x);
        Tensor2<Coord,LF> beta_u = metric.GetBeta_u(x);
        Tensor2<Coord,LF> v((u[1]+beta_u[1])/alpha,(u[2]+beta_u[2])/alpha);
        return v.template Transform<IF>(metric.GetTetrad(x));
    }
    if constexpr(std::is_same<FrameIn,LF>::value && std::is_same<FrameOut,LF>::value)
    {
        double alpha = metric.GetAlpha(x);
        Tensor2<Coord,LF> beta_u = metric.GetBeta_u(x);
        return Tensor2<Coord,LF>((u[1]+beta_u[1])/alpha,(u[2]+beta_u[2])/alpha);
    }
}

template Tensor2<xy,IF> Vec2ObservedByEulObs(const Tensor3<xy,IF>& u, const Coordinate2<xy>& x, Metric2D<xy>& metric);
template Tensor2<xy,IF> Vec2ObservedByEulObs(const Tensor3<xy,LF>& u, const Coordinate2<xy>& x, Metric2D<xy>& metric);
template Tensor2<xy,LF> Vec2ObservedByEulObs(const Tensor3<xy,IF>& u, const Coordinate2<xy>& x, Metric2D<xy>& metric);
template Tensor2<xy,LF> Vec2ObservedByEulObs(const Tensor3<xy,LF>& u, const Coordinate2<xy>& x, Metric2D<xy>& metric);
template Tensor2<rph,IF> Vec2ObservedByEulObs(const Tensor3<rph,IF>& u, const Coordinate2<rph>& x, Metric2D<rph>& metric);
template Tensor2<rph,IF> Vec2ObservedByEulObs(const Tensor3<rph,LF>& u, const Coordinate2<rph>& x, Metric2D<rph>& metric);
template Tensor2<rph,LF> Vec2ObservedByEulObs(const Tensor3<rph,IF>& u, const Coordinate2<rph>& x, Metric2D<rph>& metric);
template Tensor2<rph,LF> Vec2ObservedByEulObs(const Tensor3<rph,LF>& u, const Coordinate2<rph>& x, Metric2D<rph>& metric);