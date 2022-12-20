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
    INLINE virtual double Phi(double d, double rotation=0) const
    { exit_on_error("Stencil Phi(d) override is missing!"); return 0; }

    // Weights are only defined for discrete directions, int d.
    INLINE virtual double W(int d) const
    { exit_on_error("Stencil W(d) override is missing!"); return 0; }
    
    // Maps the interval [-0.5,nDir-0.5] -> [0,1]
    INLINE virtual double X(double d) const
    { exit_on_error("Stencil X(d) override is missing!"); return 0; }

    // Inverse of Phi(d). Returns index d of angle phi.
    INLINE virtual double Index(double phi, double rotation=0) const
    { exit_on_error("Stencil Index(phi) override is missing!"); return 0; }

    // Velocity vector.
    INLINE double Cx(double d, double rotation=0) const
    { return cos(Phi(d,rotation)); }
    INLINE double Cy(double d, double rotation=0) const
    { return sin(Phi(d,rotation)); }
    INLINE Tensor2<xy,IF> Cxy(double d, double rotation=0) const
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

    INLINE double Phi(double d, double rotation=0) const override
    { return 2.0 * M_PI * X(d); }
    INLINE double W(int d) const override
    { return 2.0 / nDir; }
    INLINE double X(double d) const override
    { return (d + 0.5) / nDir; }
    INLINE double Index(double phi, double rotation=0) const override
    { return nDir * phi / (2.0 * M_PI) - 0.5; }
};



struct DynamicStencil : Stencil
{
    DynamicStencil(int nDir_, double sigma_=1)
    {
        nDir = nDir_;
        sigma = sigma_;
    }

    INLINE double Phi(double d, double rotation=0) const override
    { return fmod(2.0 * M_PI * X(d) + rotation, 2.0 * M_PI); }
    INLINE double W(int d) const override
    { return 2.0 / nDir; }

    // New:
    INLINE double X(double d) const override
    { return d / nDir; }
    INLINE double Index(double phi, double rotation=0) const override
    {
        double angle = fmod(phi - rotation + 2.0 * M_PI, 2.0 * M_PI);
        return nDir * angle / (2.0 * M_PI);
    }
};



#endif //__INCLUDE_GUARD_Stencil_hh__