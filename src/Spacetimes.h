#ifndef __INCLUDE_GUARD_Spacetimes_h__
#define __INCLUDE_GUARD_Spacetimes_h__
#include "Metric.h" // Metric parent class.



class Minkowski : public Metric
{
public:
    Minkowski(Grid& grid_, double m_=0, double a_=0);
    Minkowski(const Minkowski& minkowski) = delete;

    bool InsideBH(const Coord& xy) override;
    std::string Name() override;
private:
    Tensor3x3 MetricFunction(const Coord& xy) override;
};



class SchwarzSchild : public Metric
{
public:
    SchwarzSchild(Grid& grid_, double m_=1, double a_=0);
    SchwarzSchild(const SchwarzSchild& schwarzSchild) = delete;

    bool InsideBH(const Coord& xy) override;
    std::string Name() override;
private:
    Tensor3x3 MetricFunction(const Coord& xy) override;
};



class KerrSchild : public Metric
{
public:
    KerrSchild(Grid& grid_, double m_=1, double a_=0);
    KerrSchild(const KerrSchild& kerrSchild) = delete;

    bool InsideBH(const Coord& xy) override;
    std::string Name() override;
private:
    Tensor3x3 MetricFunction(const Coord& xy) override;
};
#endif //__INCLUDE_GUARD_Spacetimes_h__