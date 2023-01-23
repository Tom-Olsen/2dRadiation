#ifndef __INCLUDE_GUARD_Spacetimes_h__
#define __INCLUDE_GUARD_Spacetimes_h__
#include <string>
#include "ControlFlow.hh"       // used for template arguments
#include "TensorTypes.hh"       // simple containers for rank 1-3 tensors
#include "Grid.h"             // underlying numerical Grid
#include "Metric.h"           // metric base class



template<class Coord>
class Minkowski : public Metric<Coord>
{
public:
    Minkowski(Grid<Coord>& grid_, double m_=0, double a_=0);
    Minkowski(const Minkowski& minkowski) = delete;

    bool InsideBH(const int i, const int j) override;
    bool InsideBH(const Coordinate2<Coord>& x12) override;

    std::string Name() override;
private:
    Tensor3x3<Coord,LF> MetricFunction(Coordinate2<Coord> x) override;
};



template<class Coord>
class SchwarzSchild : public Metric<Coord>
{
public:
    SchwarzSchild(Grid<Coord>& grid_, double m_=1, double a_=0);
    SchwarzSchild(const SchwarzSchild& schwarzSchild) = delete;

    bool InsideBH(const int i, const int j) override;
    bool InsideBH(const Coordinate2<Coord>& x12) override;

    std::string Name() override;
private:
    Tensor3x3<Coord,LF> MetricFunction(Coordinate2<Coord> x) override;
};



template<class Coord>
class KerrSchild : public Metric<Coord>
{
public:
    KerrSchild(Grid<Coord>& grid_, double m_=1, double a_=0);
    KerrSchild(const KerrSchild& kerrSchild) = delete;

    bool InsideBH(const int i, const int j) override;
    bool InsideBH(const Coordinate2<Coord>& x12) override;

    std::string Name() override;
private:
    Tensor3x3<Coord,LF> MetricFunction(Coordinate2<Coord> x) override;
};
#endif //__INCLUDE_GUARD_Spacetimes_h__