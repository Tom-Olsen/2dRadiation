#ifndef __INCLUDE_GUARD_TesorTypes_hh__
#define __INCLUDE_GUARD_TesorTypes_hh__
#include <iomanip>                  // std::setprecision()
#include <iostream>                 // cout
#include "ControlFlow.hh"           // used for template arguments
#include "eigen/Eigen/Dense"        // Eigen library for solving linear systems
#include "Utility.hh"               // utility functions and constants



/*
 * Access pattern for all tensors follows the mathematical access pattern:
 * A[{i,j}] = A_ij
 * 
 * This allows for mathematical equations to be translated one to one:
 * C = A*B <=> C_ij = A_ik B_kj         sum over k
 * => C[{i,j}] = A[{i,k}] * B[{k,j}]    sum over k
 * 
 * This leads to the last index being the fastest, => for(i...)for(j...)for(k...)
 *     (B00,B01,B02) 
 * B = (B10,B11,B12) <=> Memory of B = (B00, B01, B02, B10, B11, B12, B20, B21, B22)
 *     (B20,B21,B22) 
 * 
 * Initialization is row major:
 *     (0,1,2) 
 * B = (3,4,5) <=> Tensor3x3 B(0,1,2, 3,4,5, 6,7,8);
 *     (6,7,8) 
 */



template<class Coord, class Frame>
class Tensor2x2;
template<class Coord, class Frame>
class Tensor3x3;



template<class Coord>
class Coordinate2
{
private:
    double data[2];
public:
    Coordinate2(double value=0)
    { data[0] = data[1] = value; }
    Coordinate2(double data1, double data2)
    { data[0] = data1; data[1] = data2; }

    // Transforms Coord to CoordB
    template<class CoordB>
    INLINE Coordinate2<CoordB> Transform() const
    {
        if constexpr(std::is_same<Coord,xy>::value && std::is_same<CoordB,rph>::value)
            return Coordinate2<CoordB>(sqrt(data[0]*data[0] + data[1]*data[1]), fmod(atan2(data[1],data[0]) + 2.0*M_PI, 2.0*M_PI));
        if constexpr(std::is_same<Coord,rph>::value && std::is_same<CoordB,xy>::value)
            return Coordinate2<CoordB>(data[0]*cos(data[1]), data[0]*sin(data[1]));
        if constexpr(std::is_same<Coord,CoordB>::value)
            return Coordinate2<CoordB>(data[0], data[1]);
    }

    INLINE double& operator[](const int index)
    { return data[index-1]; }
    INLINE const double& operator[](const int index) const
    { return data[index-1]; }

    void Print(std::string name, bool newline=false, int precision=6) const
    {
        std::cout << name << " = ("
    	<< Format(data[0],precision) << ","
    	<< Format(data[1],precision) << ")\n";
        if(newline) std::cout << "\n";
    }
};



template<class Coord, class Frame>
class Tensor2
{
private:
    double data[2];
public:
    Tensor2(double value=0)
    { data[0] = data[1] = value; }
    Tensor2(double data1, double data2)
    { data[0] = data1; data[1] = data2; }

    // Transform between Coordinate representations.
    // x12 must be in same coordinates as initial vector.
    template<class CoordB>
    INLINE Tensor2<CoordB,Frame> Transform(const Coordinate2<Coord>& x12) const
    {
        if constexpr(std::is_same<Coord,xy>::value && std::is_same<CoordB,rph>::value)
        {// x12 = xy;
            double r = sqrt(x12[1]*x12[1] + x12[2]*x12[2]);
            return Tensor2<CoordB,Frame>(x12[1]/r*data[0] + x12[2]/r*data[1], -x12[2]/(r*r)*data[0] + x12[1]/(r*r)*data[1]);
        }
        if constexpr(std::is_same<Coord,rph>::value && std::is_same<CoordB,xy>::value)
        {// x12 = rph;
            return Tensor2<CoordB,Frame>(cos(x12[2])*data[0] - x12[1]*sin(x12[2])*data[1], sin(x12[2])*data[0] + x12[1]*cos(x12[2])*data[1]);
        }
        if constexpr(std::is_same<Coord,CoordB>::value)
            return Tensor2<CoordB,Frame>(data[0], data[1]);
    }
    // Transform between observer frames.
    // tetrad must be in same coordinates as initial vector.
    template<class FrameB>
    INLINE Tensor2<Coord,FrameB> Transform(Tensor3x3<Coord,Tetrad>& tetrad) const
    {
        if constexpr(std::is_same<Frame,IF>::value && std::is_same<FrameB,LF>::value)
        {// IF -> LF
            return Tensor2<Coord,FrameB>(tetrad[{1,1}]*data[0] + tetrad[{1,2}]*data[1], tetrad[{2,1}]*data[0] + tetrad[{2,2}]*data[1]);
        }
        if constexpr(std::is_same<Frame,LF>::value && std::is_same<FrameB,IF>::value)
        {// LF -> IF
            Tensor3x3<Coord,Tetrad> tetradInverse = tetrad.Invert();
            return Tensor2<Coord,FrameB>(tetradInverse[{1,1}]*data[0] + tetradInverse[{1,2}]*data[1], tetradInverse[{2,1}]*data[0] + tetradInverse[{2,2}]*data[1]);
        }
    }
    // Transform between observer frames.
    // tetrad must be in same coordinates as initial vector.
    template<class FrameB>
    INLINE Tensor2<Coord,FrameB> Transform(Tensor3x3<Coord,Tetrad>&& tetrad) const
    {
        if constexpr(std::is_same<Frame,IF>::value && std::is_same<FrameB,LF>::value)
        {// IF -> LF
            return Tensor2<Coord,FrameB>(tetrad[{1,1}]*data[0] + tetrad[{1,2}]*data[1], tetrad[{2,1}]*data[0] + tetrad[{2,2}]*data[1]);
        }
        if constexpr(std::is_same<Frame,LF>::value && std::is_same<FrameB,IF>::value)
        {// LF -> IF
            Tensor3x3<Coord,Tetrad> tetradInverse = tetrad.Invert();
            return Tensor2<Coord,FrameB>(tetradInverse[{1,1}]*data[0] + tetradInverse[{1,2}]*data[1], tetradInverse[{2,1}]*data[0] + tetradInverse[{2,2}]*data[1]);
        }
    }

    INLINE double Norm(const Tensor2x2<Coord,Frame>& metric2_ll) const
    {
        double norm = 0;
        for(int j=1; j<3; j++)
        for(int i=1; i<3; i++)
            norm += metric2_ll[{i,j}]*(*this)[i]*(*this)[j];
        return sqrt(abs(norm));
    }
    INLINE double EuklNorm() const
    {
        double norm = data[0]*data[0] + data[1]*data[1];
        return sqrt(norm);
    }
    INLINE double Angle() const
    {
        // transform [-pi,pi] output of atan2 to [0,2pi]:
        double angle = atan2(data[1],data[0])+2*M_PI;
        if(angle > 2*M_PI)
            return angle - 2*M_PI;
        else
            return angle;
        //return std::fmod(atan2(data[1],data[0])+2*M_PI, 2*M_PI);
    }
    INLINE double static Dot(Tensor2<Coord,Frame> v, Tensor2<Coord,Frame> w)
    {
        return v[1] * w[1] + v[2] * w[2];
    }

    INLINE double& operator[](const int index)
    { return data[index-1]; }
    INLINE const double& operator[](const int index) const
    { return data[index-1]; }

    void Print(std::string name, bool newline=false, int precision=6) const
    {
        std::cout << name;
        if constexpr(std::is_same<Frame,IF>::value)
        {
            std::cout << "(IF)" << " = ("
    	    << Format(data[0],precision) << ","
    	    << Format(data[1],precision) << ")\n";
        }
        if constexpr(std::is_same<Frame,LF>::value)
        {
            std::cout << "(LF)" << " = ("
    	    << Format(data[0],precision) << ","
    	    << Format(data[1],precision) << ")\n";
        }
        if(newline) std::cout << "\n";
    }
};



template<class Coord, class Frame>
class Tensor3
{
private:
    double data[3];
public:
    Tensor3(double value=0)
    { data[0] = data[1] = data[2] = value; }
    Tensor3(double data0, double data1, double data2)
    { data[0] = data0; data[1] = data1; data[2] = data2; }

    // Transform between Coordinate representations.
    // x12 must be in same coordinates as initial vector.
    template<class CoordB>
    INLINE Tensor3<CoordB,Frame> Transform(const Coordinate2<Coord>& x12) const
    {
        if constexpr(std::is_same<Coord,xy>::value && std::is_same<CoordB,rph>::value)
        {// x12 = xy;
            double r = sqrt(x12[1]*x12[1] + x12[2]*x12[2]);
            return Tensor3<CoordB,Frame>(data[0], x12[1]/r*data[1] + x12[2]/r*data[2], -x12[2]/(r*r)*data[1] + x12[1]/(r*r)*data[2]);
        }
        if constexpr(std::is_same<Coord,rph>::value && std::is_same<CoordB,xy>::value)
        {// x12 = rph;
            return Tensor3<CoordB,Frame>(data[0], cos(x12[2])*data[1] - x12[1]*sin(x12[2])*data[2], sin(x12[2])*data[1] + x12[1]*cos(x12[2])*data[2]);
        }
        if constexpr(std::is_same<Coord,CoordB>::value)
            return Tensor3<CoordB,Frame>(data[0], data[1], data[2]);
    }
    // Transform between observer frames.
    // tetrad must be in same coordinates as initial vector.
    template<class FrameB>
    INLINE Tensor3<Coord,FrameB> Transform(Tensor3x3<Coord,Tetrad>& tetrad) const
    {
        if constexpr(std::is_same<Frame,IF>::value && std::is_same<FrameB,LF>::value)
        {// IF -> LF
            return Tensor3<Coord,FrameB>
            (tetrad[{0,0}]*data[0] + tetrad[{0,1}]*data[1] + tetrad[{0,2}]*data[2],
             tetrad[{1,0}]*data[0] + tetrad[{1,1}]*data[1] + tetrad[{1,2}]*data[2],
             tetrad[{2,0}]*data[0] + tetrad[{2,1}]*data[1] + tetrad[{2,2}]*data[2]);
        }
        if constexpr(std::is_same<Frame,LF>::value && std::is_same<FrameB,IF>::value)
        {// LF -> IF
            Tensor3x3<Coord,Tetrad> tetradInverse = tetrad.Invert();
            return Tensor3<Coord,FrameB>
            (tetradInverse[{0,0}]*data[0] + tetradInverse[{0,1}]*data[1] + tetradInverse[{0,2}]*data[2],
             tetradInverse[{1,0}]*data[0] + tetradInverse[{1,1}]*data[1] + tetradInverse[{1,2}]*data[2],
             tetradInverse[{2,0}]*data[0] + tetradInverse[{2,1}]*data[1] + tetradInverse[{2,2}]*data[2]);
        }
    }
    // Transform between observer frames.
    // tetrad must be in same coordinates as initial vector.
    template<class FrameB>
    INLINE Tensor3<Coord,FrameB> Transform(Tensor3x3<Coord,Tetrad>&& tetrad) const
    {
        if constexpr(std::is_same<Frame,IF>::value && std::is_same<FrameB,LF>::value)
        {// IF -> LF
            return Tensor3<Coord,FrameB>
            (tetrad[{0,0}]*data[0] + tetrad[{0,1}]*data[1] + tetrad[{0,2}]*data[2],
             tetrad[{1,0}]*data[0] + tetrad[{1,1}]*data[1] + tetrad[{1,2}]*data[2],
             tetrad[{2,0}]*data[0] + tetrad[{2,1}]*data[1] + tetrad[{2,2}]*data[2]);
        }
        if constexpr(std::is_same<Frame,LF>::value && std::is_same<FrameB,IF>::value)
        {// LF -> IF
            Tensor3x3<Coord,Tetrad> tetradInverse = tetrad.Invert();
            return Tensor3<Coord,FrameB>
            (tetradInverse[{0,0}]*data[0] + tetradInverse[{0,1}]*data[1] + tetradInverse[{0,2}]*data[2],
             tetradInverse[{1,0}]*data[0] + tetradInverse[{1,1}]*data[1] + tetradInverse[{1,2}]*data[2],
             tetradInverse[{2,0}]*data[0] + tetradInverse[{2,1}]*data[1] + tetradInverse[{2,2}]*data[2]);
        }
    }

    INLINE void NullNormalize(const Tensor3x3<Coord,Frame>& metric3_ll)
    {
        double a = 0;
        for(int i=1; i<3; i++)
        for(int j=1; j<3; j++)
            a += metric3_ll[{i,j}]*data[i]*data[j];
        double b = 0;
        for(int i=1; i<3; i++)
            b += metric3_ll[{0,i}]*data[0]*data[i];
        double c = metric3_ll[{0,0}]*data[0]*data[0];
        double d = -b/a + sqrt( (b*b)/(a*a) - c/a );
        data[1] *= d;
        data[2] *= d;
    }
    INLINE double Norm(const Tensor3x3<Coord,Frame>& metric3_ll) const
    {
        double norm = 0;
        for(int j=0; j<3; j++)
        for(int i=0; i<3; i++)
            norm += metric3_ll[{i,j}]*data[i]*data[j];
        return sqrt(abs(norm));
    }

    INLINE double& operator[](const int index)
    { return data[index]; }
    INLINE const double& operator[](const int index) const
    { return data[index]; }

    void Print(std::string name, bool newline=false, int precision=6) const
    {
        std::cout << name;
        if constexpr(std::is_same<Frame,IF>::value)
        {
            std::cout << "(IF)" << " = ("
    	    << Format(data[0],precision) << ","
    	    << Format(data[1],precision) << ","
    	    << Format(data[2],precision) << ")\n";
        }
        if constexpr(std::is_same<Frame,LF>::value)
        {
            std::cout << "(LF)" << " = ("
    	    << Format(data[0],precision) << ","
    	    << Format(data[1],precision) << ","
    	    << Format(data[2],precision) << ")\n";
        }
        if(newline) std::cout << "\n";
    }
};



class Int2
{
private:
    int data[2];
public:
    Int2(int value=0)
    { data[0] = data[1] = value; }
    Int2(int data0, int data1)
    { data[0] = data0; data[1] = data1; }

    INLINE int& operator[](const int index)
    { return data[index]; }
    INLINE const int& operator[](const int index) const
    { return data[index]; }

    void Print(std::string name, bool newline=false, int precision=6) const
    {
        std::cout << name << " = ("
    	<< Format(data[0],precision) << ","
    	<< Format(data[1],precision) << ")\n";
        if(newline) std::cout << "\n";
    }
};



class Int3
{
private:
    int data[3];
public:
    Int3(int value=0)
    { data[0] = data[1] = data[2] = value; }
    Int3(int data0, int data1, int data2)
    { data[0] = data0; data[1] = data1; data[2] = data2; }

    INLINE int& operator[](const int index)
    { return data[index]; }
    INLINE const int& operator[](const int index) const
    { return data[index]; }

    void Print(std::string name, bool newline=false, int precision=6) const
    {
        std::cout << name << " = ("
    	<< Format(data[0],precision) << ","
    	<< Format(data[1],precision) << ","
    	<< Format(data[2],precision) << ")\n";
        if(newline) std::cout << "\n";
    }
};



class Double2
{
private:
    double data[2];
public:
    Double2(double value=0)
    { data[0] = data[1] = value; }
    Double2(double data0, double data1)
    { data[0] = data0; data[1] = data1; }

    INLINE double& operator[](const int index)
    { return data[index]; }
    INLINE const double& operator[](const int index) const
    { return data[index]; }

    void Print(std::string name, bool newline=false, int precision=6) const
    {
        std::cout << name << " = ("
    	<< Format(data[0],precision) << ","
    	<< Format(data[1],precision) << ")\n";
        if(newline) std::cout << "\n";
    }
};
class Double3
{
private:
    double data[3];
public:
    Double3(double value=0)
    { data[0] = data[1] = data[2] = value; }
    Double3(double data0, double data1, double data2)
    { data[0] = data0; data[1] = data1; data[2] = data2; }

    INLINE double& operator[](const int index)
    { return data[index]; }
    INLINE const double& operator[](const int index) const
    { return data[index]; }

    void Print(std::string name, bool newline=false, int precision=6) const
    {
        std::cout << name << " = ("
    	<< Format(data[0],precision) << ","
    	<< Format(data[1],precision) << ","
    	<< Format(data[2],precision) << ")\n";
        if(newline) std::cout << "\n";
    }
};



class Double4
{
private:
    double data[4];
public:
    Double4(double value=0)
    { data[0] = data[1] = data[2] = data[3] = value; }
    Double4(double data0, double data1, double data2, double data3)
    { data[0] = data0; data[1] = data1; data[2] = data2; data[3] = data3; }

    INLINE double& operator[](const int index)
    { return data[index]; }
    INLINE const double& operator[](const int index) const
    { return data[index]; }

    void Print(std::string name, bool newline=false, int precision=6) const
    {
        std::cout << name << " = ("
    	<< Format(data[0],precision) << ","
    	<< Format(data[1],precision) << ","
    	<< Format(data[2],precision) << ","
    	<< Format(data[3],precision) << ")\n";
        if(newline) std::cout << "\n";
    }
};



struct rank2Indices
{int i, j;};



template<class Coord, class Frame>
class Tensor2x2
{
private:
    double data[4];
public:
    Tensor2x2(double value=0)
    { data[0] = data[1] = data[2] = data[3] = value; }
    Tensor2x2(double data11, double data12,
              double data21, double data22)
    {
        data[0*2 + 0] = data11; data[0*2 + 1] = data12;
        data[1*2 + 0] = data21; data[1*2 + 1] = data22;
    }
    
    INLINE Tensor2x2<Coord,Frame> Invert()
    {
        Tensor2x2<Coord,Frame> invers;
        using namespace Eigen;
        Map<Matrix<double,2,2,RowMajor>> matrix(data);
        Map<Matrix<double,2,2,RowMajor>> matrixInvers(invers.data);
        matrixInvers = matrix.inverse().eval();
        return invers;
    }

    INLINE double& operator[](const rank2Indices& index)
    { return data[(index.i-1)*2 + (index.j-1)]; }
    INLINE const double& operator[](const rank2Indices& index) const
    { return data[(index.i-1)*2 + (index.j-1)]; }

    void Print(std::string name, bool newline=false, int precision=6) const
    {
        int size = name.size();
        std::string space(size,' ');
        if constexpr(std::is_same<Frame,IF>::value)
        {
            std::cout << name  << "(IF) = (" << Format(data[0*2 + 0],precision) << "," << Format(data[0*2 + 1],precision) << ")" << std::endl;
            std::cout << space << "       (" << Format(data[1*2 + 0],precision) << "," << Format(data[1*2 + 1],precision) << ")" << std::endl;
        }
        if constexpr(std::is_same<Frame,LF>::value)
        {
            std::cout << name  << "(LF) = (" << Format(data[0*2 + 0],precision) << "," << Format(data[0*2 + 1],precision) << ")" << std::endl;
            std::cout << space << "       (" << Format(data[1*2 + 0],precision) << "," << Format(data[1*2 + 1],precision) << ")" << std::endl;
        }
        if(newline)
            std::cout << std::endl;
    }
};
template<class Coord, class Frame>
class Tensor3x3
{
private:
    double data[9];
public:
    Tensor3x3(double value=0)
    {
        for(int ij=0; ij<9; ij++)
            data[ij] = value;
    }
    Tensor3x3(double data00, double data01, double data02,
              double data10, double data11, double data12,
              double data20, double data21, double data22)
    {
        data[0*3 + 0] = data00; data[0*3 + 1] = data01; data[0*3 + 2] = data02;
        data[1*3 + 0] = data10; data[1*3 + 1] = data11; data[1*3 + 2] = data12;
        data[2*3 + 0] = data20; data[2*3 + 1] = data21; data[2*3 + 2] = data22;
    }

    INLINE Tensor3x3<Coord,Frame> Invert()
    {
        Tensor3x3<Coord,Frame> invers;
        using namespace Eigen;
        Map<Matrix<double,3,3,RowMajor>> matrix(data);
        Map<Matrix<double,3,3,RowMajor>> matrixInvers(invers.data);
        matrixInvers = matrix.inverse().eval();
        return invers;
    }

    // Transform between observer frames.
    // tetrad must be in same coordinates as initial vector.
    template<class FrameB>
    INLINE Tensor3x3<Coord,FrameB> Transform(Tensor3x3<Coord,Tetrad>& tetrad) const
    {
        if constexpr(std::is_same<Frame,IF>::value && std::is_same<FrameB,LF>::value)
        {// IF -> LF
            Tensor3x3<Coord,FrameB> result(0.0);
            for(int a=0; a<3; a++)
                for(int b=0; b<3; b++)
                    for(int A=0; A<3; A++)
                        for(int B=0; B<3; B++)
                            result[{a,b}] += (*this)[{A,B}] * tetrad[{a,A}] * tetrad[{b,B}];
            return result;
        }
        if constexpr(std::is_same<Frame,LF>::value && std::is_same<FrameB,IF>::value)
        {// LF -> IF
            Tensor3x3<Coord,Tetrad> tetradInverse = tetrad.Invert();
            Tensor3x3<Coord,FrameB> result(0.0);
            for(int a=0; a<3; a++)
                for(int b=0; b<3; b++)
                    for(int A=0; A<3; A++)
                        for(int B=0; B<3; B++)
                            result[{a,b}] += (*this)[{A,B}] * tetradInverse[{a,A}] * tetradInverse[{b,B}];
            return result;
        }
    }
    // Transform between observer frames.
    // tetrad must be in same coordinates as initial vector.
    template<class FrameB>
    INLINE Tensor3x3<Coord,FrameB> Transform(Tensor3x3<Coord,Tetrad>&& tetrad) const
    {
        if constexpr(std::is_same<Frame,IF>::value && std::is_same<FrameB,LF>::value)
        {// IF -> LF
            Tensor3x3<Coord,FrameB> result(0.0);
            for(int a=0; a<3; a++)
                for(int b=0; b<3; b++)
                    for(int A=0; A<3; A++)
                        for(int B=0; B<3; B++)
                            result[{a,b}] += (*this)[{A,B}] * tetrad[{a,A}] * tetrad[{b,B}];
            return result;
        }
        if constexpr(std::is_same<Frame,LF>::value && std::is_same<FrameB,IF>::value)
        {// LF -> IF
            Tensor3x3<Coord,Tetrad> tetradInverse = tetrad.Invert();
            Tensor3x3<Coord,FrameB> result(0.0);
            for(int a=0; a<3; a++)
                for(int b=0; b<3; b++)
                    for(int A=0; A<3; A++)
                        for(int B=0; B<3; B++)
                            result[{a,b}] += (*this)[{A,B}] * tetradInverse[{a,A}] * tetradInverse[{b,B}];
            return result;
        }
    }

    INLINE double& operator[](const rank2Indices& index)
    { return data[index.i*3 + index.j]; }
    INLINE const double& operator[](const rank2Indices& index) const
    { return data[index.i*3 + index.j]; }

    void Print(std::string name, bool newline=false, int precision=6) const
    {
        int size = name.size();
        std::string space(size,' ');
        if constexpr(std::is_same<Frame,IF>::value)
        {
            std::cout << space << "       (" << Format(data[0*3 + 0],precision) << "," << Format(data[0*3 + 1],precision) << "," << Format(data[0*3 + 2],precision) << ")" << "\n";
            std::cout << name  << "(IF) = (" << Format(data[1*3 + 0],precision) << "," << Format(data[1*3 + 1],precision) << "," << Format(data[1*3 + 2],precision) << ")" << "\n";
            std::cout << space << "       (" << Format(data[2*3 + 0],precision) << "," << Format(data[2*3 + 1],precision) << "," << Format(data[2*3 + 2],precision) << ")" << "\n";
        }
        if constexpr(std::is_same<Frame,LF>::value)
        {
            std::cout << space << "       (" << Format(data[0*3 + 0],precision) << "," << Format(data[0*3 + 1],precision) << "," << Format(data[0*3 + 2],precision) << ")" << "\n";
            std::cout << name  << "(LF) = (" << Format(data[1*3 + 0],precision) << "," << Format(data[1*3 + 1],precision) << "," << Format(data[1*3 + 2],precision) << ")" << "\n";
            std::cout << space << "       (" << Format(data[2*3 + 0],precision) << "," << Format(data[2*3 + 1],precision) << "," << Format(data[2*3 + 2],precision) << ")" << "\n";
        }
        if constexpr(std::is_same<Frame,Tetrad>::value)
        {
            std::cout << space << "   (" << Format(data[0*3 + 0],precision) << "," << Format(data[0*3 + 1],precision) << "," << Format(data[0*3 + 2],precision) << ")" << "\n";
            std::cout << name  << " = (" << Format(data[1*3 + 0],precision) << "," << Format(data[1*3 + 1],precision) << "," << Format(data[1*3 + 2],precision) << ")" << "\n";
            std::cout << space << "   (" << Format(data[2*3 + 0],precision) << "," << Format(data[2*3 + 1],precision) << "," << Format(data[2*3 + 2],precision) << ")" << "\n";
        }
        if(newline) std::cout << "\n";
    }
};



struct rank3Indices
{int i, j, k;};



template<class Coord, class Frame>
class Tensor2x2x2
{
private:
    double data[8];
public:
    Tensor2x2x2(double value=0)
    {
        for(int ijk=0; ijk<8; ijk++)
            data[ijk] = value;
    }
    Tensor2x2x2
    (double data111, double data112,
     double data121, double data122,
        double data211, double data212,
        double data221, double data222)
    {
        data[0*4 + 0*2 + 0] = data111; data[0*4 + 0*2 + 1] = data112;
        data[0*4 + 1*2 + 0] = data121; data[0*4 + 1*2 + 1] = data122;
        data[1*4 + 0*2 + 0] = data211; data[1*4 + 0*2 + 1] = data212;
        data[1*4 + 1*2 + 0] = data221; data[1*4 + 1*2 + 1] = data222;
    }

    INLINE double& operator[](const rank3Indices& index)
    { return data[(index.i-1)*4 + (index.j-1)*2 + (index.k-1)]; }
    INLINE const double& operator[](const rank3Indices& index) const
    { return data[(index.i-1)*4 + (index.j-1)*2 + (index.k-1)]; }

    void Print(std::string name, bool newline=false, int precision=6) const
    {
        int size = name.size();
        std::string space(size,' ');
        if constexpr(std::is_same<Frame,IF>::value)
        {
            std::cout << space << "       (" << Format(data[0*4 + 0*2 + 0],precision) << "," << Format(data[0*4 + 0*2 + 1],precision) << ")" << "\n";
            std::cout << space << "       (" << Format(data[0*4 + 1*2 + 0],precision) << "," << Format(data[0*4 + 1*2 + 1],precision) << ")" << "\n";
            std::cout << name  << "(IF) = (---------------------)" << std::endl;
            std::cout << space << "       (" << Format(data[1*4 + 0*2 + 0],precision) << "," << Format(data[1*4 + 0*2 + 1],precision) << ")" << "\n";
            std::cout << space << "       (" << Format(data[1*4 + 1*2 + 0],precision) << "," << Format(data[1*4 + 1*2 + 1],precision) << ")" << "\n";
        }
        if constexpr(std::is_same<Frame,LF>::value)
        {
            std::cout << space << "       (" << Format(data[0*4 + 0*2 + 0],precision) << "," << Format(data[0*4 + 0*2 + 1],precision) << ")" << "\n";
            std::cout << space << "       (" << Format(data[0*4 + 1*2 + 0],precision) << "," << Format(data[0*4 + 1*2 + 1],precision) << ")" << "\n";
            std::cout << name  << "(LF) = (---------------------)" << std::endl;
            std::cout << space << "       (" << Format(data[1*4 + 0*2 + 0],precision) << "," << Format(data[1*4 + 0*2 + 1],precision) << ")" << "\n";
            std::cout << space << "       (" << Format(data[1*4 + 1*2 + 0],precision) << "," << Format(data[1*4 + 1*2 + 1],precision) << ")" << "\n";
        }
        if(newline) std::cout << "\n";
    }
};



template<class Coord, class Frame>
class Tensor3x3x3
{
private:
    double data[27];
public:
    Tensor3x3x3(double value=0)
    {
        for(int ijk=0; ijk<27; ijk++)
            data[ijk] = value;
    }
    Tensor3x3x3
    (double data000, double data001, double data002,
     double data010, double data011, double data012,
     double data020, double data021, double data022,
        double data100, double data101, double data102,
        double data110, double data111, double data112,
        double data120, double data121, double data122,
            double data200, double data201, double data202,
            double data210, double data211, double data212,
            double data220, double data221, double data222)
    {
        data[0*9 + 0*3 + 0] = data000; data[0*9 + 0*3 + 1] = data001; data[0*9 + 0*3 + 2] = data002;
        data[0*9 + 1*3 + 0] = data010; data[0*9 + 1*3 + 1] = data011; data[0*9 + 1*3 + 2] = data012;
        data[0*9 + 2*3 + 0] = data020; data[0*9 + 2*3 + 1] = data021; data[0*9 + 2*3 + 2] = data022;
        data[1*9 + 0*3 + 0] = data100; data[1*9 + 0*3 + 1] = data101; data[1*9 + 0*3 + 2] = data102;
        data[1*9 + 1*3 + 0] = data110; data[1*9 + 1*3 + 1] = data111; data[1*9 + 1*3 + 2] = data112;
        data[1*9 + 2*3 + 0] = data120; data[1*9 + 2*3 + 1] = data121; data[1*9 + 2*3 + 2] = data122;
        data[2*9 + 0*3 + 0] = data200; data[2*9 + 0*3 + 1] = data201; data[2*9 + 0*3 + 2] = data202;
        data[2*9 + 1*3 + 0] = data210; data[2*9 + 1*3 + 1] = data211; data[2*9 + 1*3 + 2] = data212;
        data[2*9 + 2*3 + 0] = data220; data[2*9 + 2*3 + 1] = data221; data[2*9 + 2*3 + 2] = data222;
    }
    
    INLINE double& operator[](const rank3Indices& index)
    { return data[index.i*9 + index.j*3 + index.k]; }
    INLINE const double& operator[](const rank3Indices& index) const
    { return data[index.i*9 + index.j*3 + index.k]; }

    void Print(std::string name, bool newline=false, int precision=6) const
    {
        int size = name.size();
        std::string space(size,' ');
        if constexpr(std::is_same<Frame,IF>::value)
        {
            std::cout << space << "       (" << Format(data[0*9 + 0*3 + 0],precision) << "," << Format(data[0*9 + 0*3 + 1],precision) << "," << Format(data[0*9 + 0*3 + 2],precision) << ")" << "\n";
            std::cout << space << "       (" << Format(data[0*9 + 1*3 + 0],precision) << "," << Format(data[0*9 + 1*3 + 1],precision) << "," << Format(data[0*9 + 1*3 + 2],precision) << ")" << "\n";
            std::cout << space << "       (" << Format(data[0*9 + 2*3 + 0],precision) << "," << Format(data[0*9 + 2*3 + 1],precision) << "," << Format(data[0*9 + 2*3 + 2],precision) << ")" << "\n";
            std::cout << space << "       (--------------------------------)" << std::endl;
            std::cout << space << "       (" << Format(data[1*9 + 0*3 + 0],precision) << "," << Format(data[1*9 + 0*3 + 1],precision) << "," << Format(data[1*9 + 0*3 + 2],precision) << ")" << "\n";
            std::cout << name  << "(IF) = (" << Format(data[1*9 + 1*3 + 0],precision) << "," << Format(data[1*9 + 1*3 + 1],precision) << "," << Format(data[1*9 + 1*3 + 2],precision) << ")" << "\n";
            std::cout << space << "       (" << Format(data[1*9 + 2*3 + 0],precision) << "," << Format(data[1*9 + 2*3 + 1],precision) << "," << Format(data[1*9 + 2*3 + 2],precision) << ")" << "\n";
            std::cout << space << "       (--------------------------------)" << std::endl;
            std::cout << space << "       (" << Format(data[2*9 + 0*3 + 0],precision) << "," << Format(data[2*9 + 0*3 + 1],precision) << "," << Format(data[2*9 + 0*3 + 2],precision) << ")" << "\n";
            std::cout << space << "       (" << Format(data[2*9 + 1*3 + 0],precision) << "," << Format(data[2*9 + 1*3 + 1],precision) << "," << Format(data[2*9 + 1*3 + 2],precision) << ")" << "\n";
            std::cout << space << "       (" << Format(data[2*9 + 2*3 + 0],precision) << "," << Format(data[2*9 + 2*3 + 1],precision) << "," << Format(data[2*9 + 2*3 + 2],precision) << ")" << "\n";
        }
        if constexpr(std::is_same<Frame,LF>::value)
        {
            std::cout << space << "       (" << Format(data[0*9 + 0*3 + 0],precision) << "," << Format(data[0*9 + 0*3 + 1],precision) << "," << Format(data[0*9 + 0*3 + 2],precision) << ")" << "\n";
            std::cout << space << "       (" << Format(data[0*9 + 1*3 + 0],precision) << "," << Format(data[0*9 + 1*3 + 1],precision) << "," << Format(data[0*9 + 1*3 + 2],precision) << ")" << "\n";
            std::cout << space << "       (" << Format(data[0*9 + 2*3 + 0],precision) << "," << Format(data[0*9 + 2*3 + 1],precision) << "," << Format(data[0*9 + 2*3 + 2],precision) << ")" << "\n";
            std::cout << space << "       (--------------------------------)" << std::endl;
            std::cout << space << "       (" << Format(data[1*9 + 0*3 + 0],precision) << "," << Format(data[1*9 + 0*3 + 1],precision) << "," << Format(data[1*9 + 0*3 + 2],precision) << ")" << "\n";
            std::cout << name  << "(LF) = (" << Format(data[1*9 + 1*3 + 0],precision) << "," << Format(data[1*9 + 1*3 + 1],precision) << "," << Format(data[1*9 + 1*3 + 2],precision) << ")" << "\n";
            std::cout << space << "       (" << Format(data[1*9 + 2*3 + 0],precision) << "," << Format(data[1*9 + 2*3 + 1],precision) << "," << Format(data[1*9 + 2*3 + 2],precision) << ")" << "\n";
            std::cout << space << "       (--------------------------------)" << std::endl;
            std::cout << space << "       (" << Format(data[2*9 + 0*3 + 0],precision) << "," << Format(data[2*9 + 0*3 + 1],precision) << "," << Format(data[2*9 + 0*3 + 2],precision) << ")" << "\n";
            std::cout << space << "       (" << Format(data[2*9 + 1*3 + 0],precision) << "," << Format(data[2*9 + 1*3 + 1],precision) << "," << Format(data[2*9 + 1*3 + 2],precision) << ")" << "\n";
            std::cout << space << "       (" << Format(data[2*9 + 2*3 + 0],precision) << "," << Format(data[2*9 + 2*3 + 1],precision) << "," << Format(data[2*9 + 2*3 + 2],precision) << ")" << "\n";
        }
        if(newline) std::cout << "\n";
    }
};


#endif //__INCLUDE_GUARD_TesorTypes_hh__




//template<int N>
//class Double
//{
//private:
//    double data[N];
//public:
//    Double(double value=0);
//    Double(double data0, double data1, double data2);
//    // Parameter Pack Constructor:	Double<N> d(d1,d2,...)
//	template<typename... Types>
//	Double(const double& Input0, const Types&... Inputs);
//    // base case for recursions
//	void fill(const int& Index, const double& Input);
//    // recursive function for parameter pack
//	template<typename... Types>
//	void fill(int Index, const double& Input0, const Types&... Inputs);
//    void Print(std::string name, bool newline=false, int precision=6) const;
//    double& operator[](const int index);
//    const double& operator[](const int index) const;
//};
//void InstantiateDoubleTemplates();
//template<int N>
//Double<N>::Double(double value)
//{
//    for(int i=0; i<N; i++)
//        data[i] = value;
//}
//template<int N>
//template<typename... Types>
//Double<N>::Double(const double& Input0, const Types&... Inputs)
//{
//    fill(0, Input0, Inputs...);
//}
//template<int N>
//void Double<N>::fill(const int& Index, const double& Input)
//{
//	if(Index<N)
//	    data[Index] = Input;
//}
//template<int N>
//template<typename... Types>
//void Double<N>::fill(int Index, const double& Input0, const Types&... Inputs)
//{
//    fill(Index, Input0);
//    Index++;
//    fill(Index, Inputs...);
//}
//template<int N>
//void Double<N>::Print(std::string name, bool newline, int precision) const
//{
//    std::cout << name << " = (" << Format(data[0],precision);
//    for(int i=1; i<N; i++)
//        std::cout << "," << Format(data[i],precision);
//	std::cout << ")" << std::endl;
//    if(newline)
//        std::cout << std::endl;
//}
//template<int N>
//double& Double<N>::operator[](const int index)
//{
//    return data[index];
//}
//template<int N>
//const double& Double<N>::operator[](const int index) const
//{
//    return data[index];
//}