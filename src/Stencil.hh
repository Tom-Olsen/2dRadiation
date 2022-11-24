#ifndef __INCLUDE_GUARD_Stencil_hh__
#define __INCLUDE_GUARD_Stencil_hh__
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include "Utility.hh"
#include "TensorTypes.hh"
#include "Distribution.hh"
#include "eigen/Eigen/Dense"



// Singleton Design Pattern Example:
//https://stackoverflow.com/questions/1008019/c-singleton-design-pattern

// TODO:
// -test heap vs stack allocation of Stencil arrays.

// This class holds the discretized velocities and the corresponding weights.
// The direction vectors, c=(cx,cy), are always given in cartesian coordinates.
struct Stencil
{
    int nDir;
    double sigma;

    // Inverse of Phi(d). Returns index d of angle phi.
    inline __attribute__((always_inline)) virtual double Index(double phi, double rotation=0) const
    { exit_on_error("Stencil Index(phi) override is missing!"); return 0; }

    // d = direction € [-0.5,nDir-0.5], continous values possible.
    // Returns angle of direction 'd' in global space, meaning:
    // 0->right, pi/2->up, pi->left, 3pi/2->down, 2pi->right.
    inline __attribute__((always_inline)) virtual double Phi(double d, double rotation=0) const
    { exit_on_error("Stencil Phi(d) override is missing!"); return 0; }

    // Weights are only defined for discrete directions, int d.
    inline __attribute__((always_inline)) virtual double W(int d) const
    { exit_on_error("Stencil W(d) override is missing!"); return 0; }

    // Velocity vector.
    inline __attribute__((always_inline)) double Cx(double d, double rotation=0) const
    { return cos(Phi(d,rotation)); }
    inline __attribute__((always_inline)) double Cy(double d, double rotation=0) const
    { return sin(Phi(d,rotation)); }
    inline __attribute__((always_inline)) Tensor2<xy,IF> Cxy(double d, double rotation=0) const
    { return Tensor2<xy,IF>(Cx(d,rotation), Cy(d,rotation)); }
    
    // Maps the interval [-0.5,nDir-0.5] -> [0,1]
    inline __attribute__((always_inline)) double X(double d) const
    { return (d) / nDir; }
    // { return (d + 0.5) / nDir; }

    void Print()
    {
        std::cout << "phi:        w:\n";
        for(int i=0; i<nDir; i++)
            std::cout << "(" << Format(Phi(i),3) << ", " << Format(W(i),3) << ")\n";
        std::cout << std::endl;
    }
};



struct UniformStencil : Stencil
{
    UniformStencil(int nDir_, double sigma_=1)
    {
        nDir = nDir_;
        sigma = sigma_;
    }

    inline __attribute__((always_inline)) double Index(double phi, double rotation=0) const override
    { return 0.5 * phi * nDir / M_PI - 0.5; }
    inline __attribute__((always_inline)) double Phi(double d, double rotation=0) const override
    { return 2.0 * M_PI * X(d); }
    inline __attribute__((always_inline)) double W(int d) const override
    { return 2.0 / nDir; }
};



struct RotStencil : Stencil
{
    RotStencil(int nDir_, double sigma_=1)
    {
        nDir = nDir_;
        sigma = sigma_;
    }

    inline __attribute__((always_inline)) double Index(double phi, double rotation=0) const override
    {
        double angle = fmod(phi - rotation + 2.0 * M_PI, 2.0 * M_PI);
        return nDir * angle / (2.0 * M_PI) - 0.5;
    }
    inline __attribute__((always_inline)) double Phi(double d, double rotation=0) const override
    { return fmod(2.0 * M_PI * X(d) + rotation, 2.0 * M_PI); }
    inline __attribute__((always_inline)) double W(int d) const override
    { return 2.0 / nDir; }
};



// Directed stencil is more dense around phi=0.
struct DirectedStencil : Stencil
{
    Distribution& distribution;

    DirectedStencil(int nDir_, Distribution& distribution_, double sigma_=1)
    : distribution(distribution_)
    {
        if(nDir_%2 != 0)
            exit_on_error("number of directions in stencil must be even!");
        nDir   = nDir_;
        sigma = sigma_;
    }

    inline __attribute__((always_inline)) double Index(double phi, double rotation=0) const override
    {
        double angle = fmod(phi - rotation + 2.0 * M_PI, 2.0 * M_PI);
        return distribution.InvValue(angle / (2.0 * M_PI)) * nDir - 0.5;
    }
    // Returns angle of direction 'd' in global space, meaning:
    // 0->right, pi/2->up, pi->left, 3pi/2->down, 2pi->right.
    inline __attribute__((always_inline)) double Phi(double d, double rotation=0) const override
    { return fmod(2.0 * M_PI * distribution.Value(X(d)) + rotation, 2.0 * M_PI); }
    inline __attribute__((always_inline)) double W(int d) const override
    {
        exit_on_error("Weights for directed stencil do not exist!");
        return 0;
    }
};


/*
struct MomentStencil : Stencil
{
    int nDir0;
    int nDir1;
    int nDir2;
    int nDir3;
    MomentStencil(int nDir0_, int nDir1_, int nDir2_, int nDir3_)
    {
        nDir0 = nDir0_;
        nDir1 = nDir1_;
        nDir2 = nDir2_;
        nDir3 = nDir3_;
        nDir = nDir0 + nDir1 + nDir2 + nDir3;
    }

    double Index(double phi, double rotation=0) const override
    {// implementation missing!
        return 0;
    }
    // Returns angle of direction 'd' in global space, meaning:
    // 0->right, pi/2->up, pi->left, 3pi/2->down, 2pi->right.
    double Y(double d) const
    {
        if (-0.5 <= d && d <= nDir0 - 0.5)
            return (d + 0.5) / (4.0 * nDir0);
        else if (nDir0 - 0.5 <= d && d <= nDir0 + nDir1 - 0.5)
            return (d + 0.5 - nDir0) / (4.0 * nDir1) + 0.25;
        else if (nDir0 + nDir1 - 0.5 <= d && d <= nDir0 + nDir1 + nDir2 - 0.5)
            return (d + 0.5 - nDir0 - nDir1) / (4.0 * nDir2) + 0.50;
        else if (nDir0 + nDir1 + nDir2 - 0.5 <= d && d <= nDir0 + nDir1 + nDir2 + nDir3 - 0.5)
            return (d + 0.5 - nDir0 - nDir1 - nDir2) / (4.0f * nDir3) + 0.75;
        else
        {
            exit_on_error("Y(d) in stencil is fucked.");
            return -1;
        }
    }
    double Phi(double d, double rotation=0) const override
    { return fmod(2.0 * M_PI * Y(d) + rotation, 2.0 * M_PI); }
    double W(int d) const override
    {
        if(d<nDir0)
            return 1.0 / (2.0 * nDir0);
        else if(nDir0<=d && d<nDir0+nDir1)
            return 1.0 / (2.0 * nDir1);
        else if(nDir0+nDir1<=d && d<nDir0+nDir1+nDir2)
            return 1.0 / (2.0 * nDir2);
        else if(nDir0+nDir1+nDir2<=d && d<nDir0+nDir1+nDir2+nDir3)
            return 1.0 / (2.0 * nDir3);
        else
            return -1;
    }
    double Sigma() const override
    {
        return 1;
    }
};



struct HalfStencil0 : Stencil
{
    HalfStencil0(int nDir_)
    { nDir = nDir_; }

    double Index(double phi, double rotation=0) const override
    { return phi * nDir / M_PI - 0.5; }
    double Phi(double d, double rotation=0) const override
    { return M_PI * X(d); }
    double W(int d) const override
    { return 1.0 / nDir; }
};



struct HalfStencil1 : Stencil
{
    HalfStencil1(int nDir_)
    { nDir = nDir_; }

    double Index(double phi, double rotation=0) const override
    { return (phi - M_PI) * nDir / M_PI - 0.5; }
    double Phi(double d, double rotation=0) const override
    { return M_PI * X(d) + M_PI; }
    double W(int d) const override
    { return 1.0 / nDir; }
};
*/

#endif //__INCLUDE_GUARD_Stencil_hh__