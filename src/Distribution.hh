#ifndef __INCLUDE_GUARD_Distribution_hh__
#define __INCLUDE_GUARD_Distribution_hh__

#include <math.h>
#include "Utility.h"



// Distribution base class.
// Contains name of the distribution, the distribution function (Value) and its inverse (InvValue).
// Both these functions Map smoothly from [0,1] to [0,1] and are monotonic increasing => bijective.
// Feel free to add more distributions!
struct Distribution
{
    virtual std::string Name()
    {
        exit_on_error("Distribution virtual Method (Name) has been called!");
        return "";
    }
    inline __attribute__((always_inline)) virtual double Value(double x)
    {
        exit_on_error("Distribution virtual Method (Value) has been called!");
        return 0;
    }
    inline __attribute__((always_inline)) virtual double InvValue(double x)
    {
        exit_on_error("Distribution virtual Method (InvValue) has been called!");
        return 0;
    }
};



struct X2 : Distribution
{
    std::string Name() override
    { return "X^2"; }
    inline __attribute__((always_inline)) double Value(double x) override
    {
        if (x < 0.5)
            return 2.0 * x * x;
        else
            return -2.0 * (x - 1.0) * (x - 1.0) + 1.0;
    }
    inline __attribute__((always_inline)) double InvValue(double x) override
    {
        if (x < 0.5)
            return sqrt(0.5 * x);
        else
            return -sqrt(0.5 * (1.0 - x)) + 1.0;
    }
};



struct X4 : Distribution
{
    std::string Name() override
    { return "X^4"; }
    inline __attribute__((always_inline)) double Value(double x) override
    {
        if (x < 0.5f)
            return 8.0 * pow(x,4);
        else
            return -8.0 * pow(x - 1.0,4) + 1.0;
    }
    inline __attribute__((always_inline)) double InvValue(double x) override
    {
        if (x < 0.5)
            return pow(0.125 * x, 0.25);
        else
            return -pow(0.125 * (1.0 - x), 0.25) + 1.0;
    }
};



struct Cosine : Distribution
{
    std::string Name() override
    { return "Cosine"; }
    inline __attribute__((always_inline)) double Value(double x) override
    {
        return 0.5 * (1 - cos(M_PI * x));
    }
    inline __attribute__((always_inline)) double InvValue(double x) override
    {
        return acos(1 - 2 * x) / M_PI;
    }
};



struct PieceWiseLinear0 : Distribution
{
    std::string Name() override
    { return "PieceWiseLinear0"; }
    inline __attribute__((always_inline)) double Value(double x) override
    {
        if (x < 0.3)
            return (1.0 / 3.0) * x;
        else if (x < 0.7)
            return 2.0 * x - 0.5;
        else
            return (1.0 / 3.0) * x + (2.0 / 3.0);
    }
    inline __attribute__((always_inline)) double InvValue(double x) override
    {
        if (x < 0.1)
            return 3.0 * x;
        else if (x < 0.9)
            return 0.5 * (x + 0.5);
        else
            return 3.0 * x - 2.0;
    }
};



struct PieceWiseLinear1 : Distribution
{
    std::string Name() override
    { return "PieceWiseLinear1"; }
    inline __attribute__((always_inline)) double Value(double x) override
    {
        if (x < 0.15)
            return 1.0 / 6.0 * x;
        else if (x < 0.3)
            return 1.0 / 2.0 * x - 1.0 / 20.0;
        else if (x < 0.7)
            return 2.0 * x - 0.5;
        else if (x < 0.85)
            return 1.0 / 2.0 * x + 11.0 / 20.0;
        else
            return 1.0 / 6.0 * x + 5.0 / 6.0;
    }
    inline __attribute__((always_inline)) double InvValue(double x) override
    {
        if (x < 0.025)
            return 6.0 * x;
        else if (x < 0.1)
            return 2.0 * x + 1.0 / 10.0;
        else if (x < 0.9)
            return 0.5 * (x + 0.5);
        else if (x < 0.975)
            return 2.0 * x - 11.0 / 10.0;
        else
            return 6.0 * x - 5.0;
    }
};

#endif //__INCLUDE_GUARD_Distribution_hh__