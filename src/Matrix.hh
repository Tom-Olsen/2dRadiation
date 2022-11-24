#ifndef __INCLUDE_GUARD_Matrix_hh__
#define __INCLUDE_GUARD_Matrix_hh__
#include <iostream>             // cout etc.
#include "Utility.hh"            // basic utility
#include "eigen/Eigen/Dense"    // linear algebra


/*
template<int N>
class VectorN
{
public:
    double data[N];
    VectorN()
    {
        for(int i=0; i<N; i++)
            data[i] = 0;
    }
    VectorN(const std::initializer_list<double>& list)
    {
		int i=0;
		for(auto& element : list)
		{ data[i] = element; i++; }
    }
    void Print(std::string name, bool newline=false, int precision=6)
    {
        std::cout << name + ": ";
        std::cout << "(" << Format(data[0],precision);
        for(int i=1; i<N; i++)
            std::cout << ", " << Format(data[i],precision);
        std::cout << ")" << std::endl;;
    }
    double& operator[](const int& index)
    { return data[index]; }
    const double& operator[](const int& index) const
    { return data[index]; }
};



struct MatrixIndices
{int i, j;};
template<int N>
class MatrixNxN
{
public:
    double data[N*N];
    MatrixNxN()
    {
        for(int i=0; i<N*N; i++)
            data[i] = 0;
    }
    MatrixNxN(const std::initializer_list<double>& list)
    {
		int i=0;
		for(auto& element : list)
		{ data[i] = element; i++; }
    }
    void Print(std::string name, bool newline=false, int precision=6)
    {
        std::cout << name + ":\n";
        for(int i=0; i<N; i++)
        {
            std::cout << "(" << Format(data[i*N + 0],precision);
            for(int j=1; j<N; j++)
                std::cout << ", " << Format(data[i*N + j],precision);
            std::cout << ")" << std::endl;
        }
        if(newline)
            std::cout << "\n";
    }
    double& operator[](const MatrixIndices& index)
    { return data[index.i*N + index.j]; }
    const double& operator[](const MatrixIndices& index) const
    { return data[index.i*N + index.j]; }
    MatrixNxN Invert()
    {
        MatrixNxN<N> invers;
        using namespace Eigen;
        Map<Matrix<double,N,N,RowMajor>> matrix(data);
        Map<Matrix<double,N,N,RowMajor>> matrixInvers(invers.data);
        matrixInvers = matrix.inverse().eval();
        return invers;
    }
    VectorN<N> SolveLS(VectorN<N> b)
    {
        using namespace Eigen;
        Map<Matrix<double,N,N,RowMajor>> matM(data);
        Map<Vector<double,N>> vecB(b.data);
        VectorN<N> x;
        Map<Vector<double,N>> vecX(x.data);

        vecX = matM.partialPivLu().solve(vecB);

        return x;
    }
};*/
#endif //__INCLUDE_GUARD_Matrix_hh__