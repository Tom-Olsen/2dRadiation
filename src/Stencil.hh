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

    // d = direction € [-0.5,nDir-0.5], continous values possible.
    // Returns angle of direction 'd' in global space, meaning:
    // 0->right, pi/2->up, pi->left, 3pi/2->down, 2pi->right.
    inline __attribute__((always_inline)) virtual double Phi(double d, double rotation=0) const
    { exit_on_error("Stencil Phi(d) override is missing!"); return 0; }

    // Weights are only defined for discrete directions, int d.
    inline __attribute__((always_inline)) virtual double W(int d) const
    { exit_on_error("Stencil W(d) override is missing!"); return 0; }
    
    // Maps the interval [-0.5,nDir-0.5] -> [0,1]
    inline __attribute__((always_inline)) virtual double X(double d) const
    { exit_on_error("Stencil X(d) override is missing!"); return 0; }

    // Inverse of Phi(d). Returns index d of angle phi.
    inline __attribute__((always_inline)) virtual double Index(double phi, double rotation=0) const
    { exit_on_error("Stencil Index(phi) override is missing!"); return 0; }

    // Velocity vector.
    inline __attribute__((always_inline)) double Cx(double d, double rotation=0) const
    { return cos(Phi(d,rotation)); }
    inline __attribute__((always_inline)) double Cy(double d, double rotation=0) const
    { return sin(Phi(d,rotation)); }
    inline __attribute__((always_inline)) Tensor2<xy,IF> Cxy(double d, double rotation=0) const
    { return Tensor2<xy,IF>(Cx(d,rotation), Cy(d,rotation)); }
    

    void Print(double rotation=0)
    {
        std::cout << "phi:        w:\n";
        for(int i=0; i<nDir; i++)
            std::cout << "(" << Format(Phi(i,rotation),3) << ", " << Format(W(i),3) << ")\n";
        std::cout << std::endl;
    }
};



struct StaticStencil : Stencil
{
    StaticStencil(int nDir_, double sigma_=1)
    {
        nDir = nDir_;
        sigma = sigma_;
    }

    inline __attribute__((always_inline)) double Phi(double d, double rotation=0) const override
    { return 2.0 * M_PI * X(d); }
    inline __attribute__((always_inline)) double W(int d) const override
    { return 2.0 / nDir; }
    inline __attribute__((always_inline)) double X(double d) const override
    { return (d + 0.5) / nDir; }
    inline __attribute__((always_inline)) double Index(double phi, double rotation=0) const override
    { return nDir * phi / (2.0 * M_PI) - 0.5; }
};



struct RotatingStencil : Stencil
{
    RotatingStencil(int nDir_, double sigma_=1)
    {
        nDir = nDir_;
        sigma = sigma_;
    }

    inline __attribute__((always_inline)) double Phi(double d, double rotation=0) const override
    { return fmod(2.0 * M_PI * X(d) + rotation, 2.0 * M_PI); }
    inline __attribute__((always_inline)) double W(int d) const override
    { return 2.0 / nDir; }

    // New:
    inline __attribute__((always_inline)) double X(double d) const override
    { return d / nDir; }
    inline __attribute__((always_inline)) double Index(double phi, double rotation=0) const override
    {
        double angle = fmod(phi - rotation + 2.0 * M_PI, 2.0 * M_PI);
        return nDir * angle / (2.0 * M_PI);
    }
};



#endif //__INCLUDE_GUARD_Stencil_hh__