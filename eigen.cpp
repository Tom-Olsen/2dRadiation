#include <iostream>
#include <iomanip>
#include <vector>
#include "src/eigen/Eigen/Dense"
#include "src/Utility.h"

using namespace Eigen;
using namespace std;


// returns number between 0 and 1 for all d € [-0.5,nDir-0.5].
double X(double d, int N)
{ return 0.5 * (d + 0.5) / N; }
double Distribution(double x)
{
    //return x;    // => Gauß, UniformStencil
    if (x < 0.5)
        return 2.0 * x * x;
    else
        return -2.0 * (x - 1.0) * (x - 1.0) + 1.0;
}
double Phi(double d, int N)
{ return 2.0 * M_PI * Distribution(X(d,N)); }



void Print(VectorXd vec, string name)
{
    for(int d=0; d<vec.size(); d++)
        cout << name << "[" << to_string(d) << "]: " << Format(vec[d]) << endl;
    cout << endl;
}
void Print(MatrixXd mat, string name)
{
    cout << name << ":" << endl;
    for(int i=0; i<mat.rows(); i++)
    {
        for(int d=0; d<mat.cols(); d++)
            cout << Format(mat(i,d)) << "\t";
        cout << endl;
    }
    cout << endl;
}
VectorXd Mul(MatrixXd mat, VectorXd vec)
{
    VectorXd result(mat.rows());
    for(int i=0; i<mat.rows(); i++)
    {
        result[i] = 0;
        for(int j=0; j<mat.cols(); j++)
            result[i] += mat(i,j) * vec(j);
    }
    return result;
}



VectorXd GetWeightsLagrange(int numberOfWeights, int numberOfEquations)
{
    int N = numberOfWeights;
    int M = numberOfEquations;

    // All vectors and matrices:
    double angles[N];
    VectorXd I(M);
    VectorXd lambda(N);
    VectorXd x(M+N);
    VectorXd y(M+N);
    VectorXd w(N);
    MatrixXd A(M,N);
    MatrixXd AT(N,M);
    MatrixXd Zero(M,M);
    MatrixXd Unit(N,N);
    
    // Set angles:
    for(int d=0; d<N; d++)
        angles[d] = Phi(d,N);

    // Set vectors:
    I(0) = 1.0;
    for(int d=1; d<M; d++)
        I(d) = (pow(-1,d) + 1) * tgamma((d+1.0)/2.0) / (sqrt(M_PI)*d*tgamma(d/2.0));
    for(int d=0; d<N; d++)
        lambda(d) = 1.0/N;
    y << I, lambda;

    // Set matrices:
    for(int i=0; i<M; i++)
        for(int j=0; j<N; j++)
            A(i,j) = pow(cos(angles[j]),i);
    AT = A.transpose();
    Zero = MatrixXd::Constant(M, M, 0);
    Unit.setIdentity(N,N);

    // Concatinate Matrices:
    MatrixXd UpperPart(A.rows(), A.cols() + Zero.cols());
    UpperPart << A, Zero;
    MatrixXd LowerPart(AT.rows(), AT.cols() + Unit.cols());
    LowerPart << Unit, AT;
    MatrixXd FullMatrix(M+N,M+N);
    FullMatrix << UpperPart, LowerPart;

    // Solve linear system:
    x = FullMatrix.colPivHouseholderQr().solve(y);
    for(int d=0; d<N; d++)
        w(d) = x(d);

    return w;
}


VectorXd GetPositiveWeights(int numberOfWeights)
{
    int N = numberOfWeights;
    
    for(int M=N; M>0; M--)
    {
        VectorXd w = GetWeightsLagrange(N, M);
        for(int k=0; k<M; k++)
        {
            if(w(k) < 0)
                break;
            else
                return w;
        }
    }

    exit_on_error("no solution for positive weights found!");
    return GetWeightsLagrange(0, 0); // no negate warning
}


int main()
{
    VectorXd w = GetPositiveWeights(50);

    Print(w,"w");
}