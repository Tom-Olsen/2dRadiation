#ifndef __INCLUDE_GUARD_DataTypes_hh__
#define __INCLUDE_GUARD_DataTypes_hh__
#include <iomanip>           // std::setprecision()
#include <iostream>          // cout
#include <fstream>           // File input/output.
#include "ControlFlow.hh"    // used for template arguments
#include "eigen/Eigen/Dense" // Eigen library for solving linear systems
#include "Utility.hh"        // utility functions and constants

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

struct rank2Indices
{
    int i, j;
};
struct rank3Indices
{
    int i, j, k;
};

struct Coord
{
private:
    double data[2];

public:
    // Constructors:
    Coord(double value = 0)
    {
        data[0] = data[1] = value;
    }
    Coord(double x, double y)
    {
        data[0] = x;
        data[1] = y;
    }

    // Accessors:
    double &operator[](const int index)
    {
        return data[index - 1];
    }
    const double &operator[](const int index) const
    {
        return data[index - 1];
    }

    // Lengths and Angles:
    double EuklNormSquared() const
    {
        return Dot((*this), (*this));
    }
    double EuklNorm() const
    {
        return sqrt(EuklNormSquared());
    }
    double Phi()
    {
        double phi = MyAtan2(data[1], data[0]);
        return (phi < 0) ? phi + 2.0 * M_PI : phi;
    }

    // Operator overloading:
    Coord operator+(const Coord &other) const
    {
        return Coord((*this)[1] + other[1], (*this)[2] + other[2]);
    }
    Coord operator-(const Coord &other) const
    {
        return Coord((*this)[1] - other[1], (*this)[2] - other[2]);
    }

    // Basic Math:
    static double Dot(const Coord &a, const Coord &b)
    {
        return a[1] * b[1] + a[2] * b[2];
    }
    static double Cross(const Coord &a, const Coord &b)
    {
        return a[1] * b[2] - b[1] * a[2];
    }

    // Output:
    void Print(std::string name, bool newline = false, int precision = 6) const
    {
        std::cout << name << " = ("
                  << Format(data[0], precision) << ","
                  << Format(data[1], precision) << ")\n";
        if (newline)
            std::cout << "\n";
    }
    friend std::ostream &operator<<(std::ostream &os, const Coord &x);
};
inline std::ostream &operator<<(std::ostream &os, const Coord &x)
{
    os << "(" << x[1] << "," << x[2] << ")";
    return os;
}
inline Coord operator*(const double &scalar, const Coord &coord)
{
    return Coord(scalar * coord[1], scalar * coord[2]);
}
inline Coord operator*(const Coord &coord, const double &scalar)
{
    return Coord(scalar * coord[1], scalar * coord[2]);
}
inline Coord operator/(const Coord &coord, const double &scalar)
{
    return Coord(coord[1] / scalar, coord[2] / scalar);
}

struct Tensor2
{
private:
    double data[2];

public:
    // Constructors:
    Tensor2(double value = 0)
    {
        data[0] = data[1] = value;
    }
    Tensor2(double x, double y)
    {
        data[0] = x;
        data[1] = y;
    }

    // Accessors:
    double &operator[](const int index)
    {
        return data[index - 1];
    }
    const double &operator[](const int index) const
    {
        return data[index - 1];
    }

    // Lengths and Angles:
    double EuklNorm() const
    {
        return sqrt(Dot((*this), (*this)));
    }
    Tensor2 EuklNormalized() const
    {
        double norm = EuklNorm();
        return Tensor2(data[0] / norm, data[1] / norm);
    }
    double Phi() const
    {
        double phi = MyAtan2(data[1], data[0]);
        return (phi < 0) ? phi + 2.0 * M_PI : phi;
    }

    // Operator overloading:
    Tensor2 operator+(const Tensor2 &other) const
    {
        return Tensor2((*this)[1] + other[1], (*this)[2] + other[2]);
    }
    Tensor2 operator-(const Tensor2 &other) const
    {
        return Tensor2((*this)[1] - other[1], (*this)[2] - other[2]);
    }

    // Basic Math:
    static double Dot(const Tensor2 &a, const Tensor2 &b)
    {
        return a[1] * b[1] + a[2] * b[2];
    }
    static double Cross(const Tensor2 &a, const Tensor2 &b)
    {
        return a[1] * b[2] - b[1] * a[2];
    }

    // Output:
    void Print(std::string name, bool newline = false, int precision = 6) const
    {
        std::cout << name << " = ("
                  << Format(data[0], precision) << ","
                  << Format(data[1], precision) << ")\n";
        if (newline)
            std::cout << "\n";
    }
    friend std::ostream &operator<<(std::ostream &os, const Tensor2 &v);
};
inline std::ostream &operator<<(std::ostream &os, const Tensor2 &v)
{
    os << "(" << v[1] << "," << v[2] << ")";
    return os;
}
inline Tensor2 operator*(const double &scalar, const Tensor2 &tensor)
{
    return Tensor2(scalar * tensor[1], scalar * tensor[2]);
}
inline Tensor2 operator*(const Tensor2 &tensor, const double &scalar)
{
    return Tensor2(scalar * tensor[1], scalar * tensor[2]);
}
inline Tensor2 operator/(const Tensor2 &tensor, const double &scalar)
{
    return Tensor2(tensor[1] / scalar, tensor[2] / scalar);
}

struct Tensor3
{
private:
    double data[3];

public:
    // Constructors:
    Tensor3(double value = 0)
    {
        data[0] = data[1] = data[2] = value;
    }
    Tensor3(double t, double x, double y)
    {
        data[0] = t;
        data[1] = x;
        data[2] = y;
    }

    // Accessors:
    double &operator[](const int index)
    {
        return data[index];
    }
    const double &operator[](const int index) const
    {
        return data[index];
    }

    // Operator overloading:
    Tensor3 operator+(const Tensor3 &other) const
    {
        return Tensor3((*this)[0] + other[0], (*this)[1] + other[1], (*this)[2] + other[2]);
    }
    Tensor3 operator-(const Tensor3 &other) const
    {
        return Tensor3((*this)[0] - other[0], (*this)[1] - other[1], (*this)[2] - other[2]);
    }

    // Basic Math:
    static double Dot(const Tensor3 &a, const Tensor3 &b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }


    // Output:
    void Print(std::string name, bool newline = false, int precision = 6) const
    {
        std::cout << name << " = ("
                  << Format(data[0], precision) << ","
                  << Format(data[1], precision) << ","
                  << Format(data[2], precision) << ")\n";
        if (newline)
            std::cout << "\n";
    }
    friend std::ostream &operator<<(std::ostream &os, const Tensor3 &v);
};
inline std::ostream &operator<<(std::ostream &os, const Tensor3 &v)
{
    os << "(" << v[0] << "," << v[1] << "," << v[2] << ")";
    return os;
}
inline Tensor3 operator*(const double &scalar, const Tensor3 &tensor)
{
    return Tensor3(scalar * tensor[0], scalar * tensor[1], scalar * tensor[2]);
}
inline Tensor3 operator*(const Tensor3 &tensor, const double &scalar)
{
    return Tensor3(scalar * tensor[0], scalar * tensor[1], scalar * tensor[2]);
}
inline Tensor3 operator/(const Tensor3 &tensor, const double &scalar)
{
    return Tensor3(tensor[0] / scalar, tensor[1] / scalar, tensor[2] / scalar);
}

struct Tensor2x2
{
private:
    double data[4];

public:
    // Constructors:
    Tensor2x2(double value = 0)
    {
        data[0] = data[1] = data[2] = data[3] = value;
    }
    Tensor2x2(double xx, double xy,
              double yx, double yy)
    {
        data[0] = xx;
        data[1] = xy;
        data[2] = yx;
        data[3] = yy;
    }

    // Accessors:
    double &operator[](const rank2Indices &index)
    {
        return data[(index.i - 1) * 2 + (index.j - 1)];
    }
    const double &operator[](const rank2Indices &index) const
    {
        return data[(index.i - 1) * 2 + (index.j - 1)];
    }

    // Operator overloading:
    Tensor2x2 operator*(const Tensor2x2 &other) const
    {
        Tensor2x2 result(0);
        for(int i = 1; i < 3; i++)
            for(int j = 1; j < 3; j++)
                for(int k = 1; k < 3; k++)
                    result[{i,j}] += (*this)[{i,k}] * other[{k,j}];
        return result;
    }


    // Basic Math:
    Tensor2x2 Invert()
    {
        Tensor2x2 invers;
        using namespace Eigen;
        Map<Matrix<double, 2, 2, RowMajor>> matrix(data);
        Map<Matrix<double, 2, 2, RowMajor>> matrixInvers(invers.data);
        matrixInvers = matrix.inverse().eval();
        return invers;
    }

    // Output
    void Print(std::string name, bool newline = false, int precision = 6) const
    {
        int size = name.size();
        std::string space(size, ' ');
        std::cout << name << " = (" << Format(data[0 * 2 + 0], precision) << "," << Format(data[0 * 2 + 1], precision) << ")"
                  << "\n";
        std::cout << space << " = (" << Format(data[1 * 2 + 0], precision) << "," << Format(data[1 * 2 + 1], precision) << ")"
                  << "\n";
        if (newline)
            std::cout << "\n";
    }
    friend std::ostream &operator<<(std::ostream &os, const Tensor2x2 &v);
};
inline std::ostream &operator<<(std::ostream &os, const Tensor2x2 &v)
{
    os << "(" << v[{1, 1}] << "," << v[{1, 2}]
       << "|" << v[{2, 1}] << "," << v[{2, 2}] << ")";
    return os;
}
inline Tensor2 operator*(const Tensor2x2 &matrix, const Tensor2 &vector)
{
    return Tensor2(matrix[{1, 1}] * vector[1] + matrix[{1, 2}] * vector[2],
                   matrix[{2, 1}] * vector[1] + matrix[{2, 2}] * vector[2]);
}

struct Tensor3x3
{
private:
    double data[9];

public:
    // Constructors:
    Tensor3x3(double value = 0)
    {
        for (int ij = 0; ij < 9; ij++)
            data[ij] = value;
    }
    Tensor3x3(double tt, double tx, double ty,
              double xt, double xx, double xy,
              double yt, double yx, double yy)
    {
        data[0] = tt;
        data[1] = tx;
        data[2] = ty;
        data[3] = xt;
        data[4] = xx;
        data[5] = xy;
        data[6] = yt;
        data[7] = yx;
        data[8] = yy;
    }

    // Accessors:
    double &operator[](const rank2Indices &index)
    {
        return data[index.i * 3 + index.j];
    }
    const double &operator[](const rank2Indices &index) const
    {
        return data[index.i * 3 + index.j];
    }
    const Tensor3 GetColumn(int index) const
    {
        return Tensor3(data[index], data[index + 3], data[index + 6]);
    }

    // Basic Math:
    Tensor3x3 Invert()
    {
        Tensor3x3 invers;
        using namespace Eigen;
        Map<Matrix<double, 3, 3, RowMajor>> matrix(data);
        Map<Matrix<double, 3, 3, RowMajor>> matrixInvers(invers.data);
        matrixInvers = matrix.inverse().eval();
        return invers;
    }

    void Print(std::string name, bool newline = false, int precision = 6) const
    {
        int size = name.size();
        std::string space(size, ' ');
        std::cout << space << "   (" << Format(data[0], precision) << "," << Format(data[1], precision) << "," << Format(data[2], precision) << ")"
                  << "\n";
        std::cout << name << " = (" << Format(data[3], precision) << "," << Format(data[4], precision) << "," << Format(data[5], precision) << ")"
                  << "\n";
        std::cout << space << "   (" << Format(data[6], precision) << "," << Format(data[7], precision) << "," << Format(data[8], precision) << ")"
                  << "\n";
        if (newline)
            std::cout << "\n";
    }
    friend std::ostream &operator<<(std::ostream &os, const Tensor3x3 &v);
};
inline std::ostream &operator<<(std::ostream &os, const Tensor3x3 &v)
{
    os << "(" << v[{0, 0}] << "," << v[{0, 1}] << "," << v[{0, 2}]
       << "|" << v[{1, 0}] << "," << v[{1, 1}] << "," << v[{1, 2}]
       << "|" << v[{2, 0}] << "," << v[{2, 1}] << "," << v[{2, 2}] << ")";
    return os;
}

struct Tensor2x2x2
{
private:
    double data[8];

public:
    // Constructors:
    Tensor2x2x2(double value = 0)
    {
        for (int ijk = 0; ijk < 8; ijk++)
            data[ijk] = value;
    }
    Tensor2x2x2(double xxx, double xxy,
                double xyx, double xyy,
                double yxx, double yxy,
                double yyx, double yyy)
    {
        data[0] = xxx;
        data[1] = xxy;
        data[2] = xyx;
        data[3] = xyy;
        data[4] = yxx;
        data[5] = yxy;
        data[6] = yyx;
        data[7] = yyy;
    }

    // Accessors:
    double &operator[](const rank3Indices &index)
    {
        return data[(index.i - 1) * 4 + (index.j - 1) * 2 + (index.k - 1)];
    }
    const double &operator[](const rank3Indices &index) const
    {
        return data[(index.i - 1) * 4 + (index.j - 1) * 2 + (index.k - 1)];
    }

    // Output
    void Print(std::string name, bool newline = false, int precision = 6) const
    {
        int size = name.size();
        std::string space(size, ' ');
        std::cout << space << "   (" << Format(data[0], precision) << "," << Format(data[1], precision) << ")"
                  << "\n";
        std::cout << space << "   (" << Format(data[2], precision) << "," << Format(data[3], precision) << ")"
                  << "\n";
        std::cout << name << " = (---------------------)" << std::endl;
        std::cout << space << "   (" << Format(data[4], precision) << "," << Format(data[5], precision) << ")"
                  << "\n";
        std::cout << space << "   (" << Format(data[6], precision) << "," << Format(data[7], precision) << ")"
                  << "\n";
        if (newline)
            std::cout << "\n";
    }
    friend std::ostream &operator<<(std::ostream &os, const Tensor2x2x2 &v);
};
inline std::ostream &operator<<(std::ostream &os, const Tensor2x2x2 &v)
{
    os << "(" << v[{1, 1, 1}] << "," << v[{1, 1, 2}]
       << "|" << v[{1, 2, 1}] << "," << v[{1, 2, 2}] << ")";
    os << "(" << v[{2, 1, 1}] << "," << v[{2, 1, 2}]
       << "|" << v[{2, 2, 1}] << "," << v[{2, 2, 2}] << ")";
    return os;
}

struct Tensor3x3x3
{
private:
    double data[27];

public:
    Tensor3x3x3(double value = 0)
    {
        for (int ijk = 0; ijk < 27; ijk++)
            data[ijk] = value;
    }
    // Constructors:
    Tensor3x3x3(double ttt, double ttx, double tty,
                double txt, double txx, double txy,
                double tyt, double tyx, double tyy,
                double xtt, double xtx, double xty,
                double xxt, double xxx, double xxy,
                double xyt, double xyx, double xyy,
                double ytt, double ytx, double yty,
                double yxt, double yxx, double yxy,
                double yyt, double yyx, double yyy)
    {
        data[0] = ttt;
        data[1] = ttx;
        data[2] = tty;
        data[3] = txt;
        data[4] = txx;
        data[5] = txy;
        data[6] = tyt;
        data[7] = tyx;
        data[8] = tyy;
        data[9] = xtt;
        data[10] = xtx;
        data[11] = xty;
        data[12] = xxt;
        data[13] = xxx;
        data[14] = xxy;
        data[15] = xyt;
        data[16] = xyx;
        data[17] = xyy;
        data[18] = ytt;
        data[19] = ytx;
        data[20] = yty;
        data[21] = yxt;
        data[22] = yxx;
        data[23] = yxy;
        data[24] = yyt;
        data[25] = yyx;
        data[26] = yyy;
    }

    // Accessors:
    double &operator[](const rank3Indices &index)
    {
        return data[index.i * 9 + index.j * 3 + index.k];
    }
    const double &operator[](const rank3Indices &index) const
    {
        return data[index.i * 9 + index.j * 3 + index.k];
    }

    // Output
    void Print(std::string name, bool newline = false, int precision = 6) const
    {
        int size = name.size();
        std::string space(size, ' ');
        std::cout << space << "   (" << Format(data[0], precision) << "," << Format(data[1], precision) << "," << Format(data[2], precision) << ")"
                  << "\n";
        std::cout << space << "   (" << Format(data[3], precision) << "," << Format(data[4], precision) << "," << Format(data[5], precision) << ")"
                  << "\n";
        std::cout << space << "   (" << Format(data[6], precision) << "," << Format(data[7], precision) << "," << Format(data[8], precision) << ")"
                  << "\n";
        std::cout << space << "   (--------------------------------)" << std::endl;
        std::cout << space << "   (" << Format(data[9], precision) << "," << Format(data[10], precision) << "," << Format(data[11], precision) << ")"
                  << "\n";
        std::cout << name << " = (" << Format(data[12], precision) << "," << Format(data[13], precision) << "," << Format(data[14], precision) << ")"
                  << "\n";
        std::cout << space << "   (" << Format(data[15], precision) << "," << Format(data[16], precision) << "," << Format(data[17], precision) << ")"
                  << "\n";
        std::cout << space << "   (--------------------------------)" << std::endl;
        std::cout << space << "   (" << Format(data[18], precision) << "," << Format(data[19], precision) << "," << Format(data[20], precision) << ")"
                  << "\n";
        std::cout << space << "   (" << Format(data[21], precision) << "," << Format(data[22], precision) << "," << Format(data[23], precision) << ")"
                  << "\n";
        std::cout << space << "   (" << Format(data[24], precision) << "," << Format(data[25], precision) << "," << Format(data[26], precision) << ")"
                  << "\n";
        if (newline)
            std::cout << "\n";
    }
    friend std::ostream &operator<<(std::ostream &os, const Tensor3x3x3 &v);
};
inline std::ostream &operator<<(std::ostream &os, const Tensor3x3x3 &v)
{
    os << "(" << v[{0, 0, 0}] << "," << v[{0, 0, 1}] << "," << v[{0, 0, 2}]
       << "|" << v[{0, 1, 0}] << "," << v[{0, 1, 1}] << "," << v[{0, 1, 2}]
       << "|" << v[{0, 2, 0}] << "," << v[{0, 2, 1}] << "," << v[{0, 2, 2}] << ")";
    os << "(" << v[{1, 0, 0}] << "," << v[{1, 0, 1}] << "," << v[{1, 0, 2}]
       << "|" << v[{1, 1, 0}] << "," << v[{1, 1, 1}] << "," << v[{1, 1, 2}]
       << "|" << v[{1, 2, 0}] << "," << v[{1, 2, 1}] << "," << v[{1, 2, 2}] << ")";
    os << "(" << v[{2, 0, 0}] << "," << v[{2, 0, 1}] << "," << v[{2, 0, 2}]
       << "|" << v[{2, 1, 0}] << "," << v[{2, 1, 1}] << "," << v[{2, 1, 2}]
       << "|" << v[{2, 2, 0}] << "," << v[{2, 2, 1}] << "," << v[{2, 2, 2}] << ")";
    return os;
}

struct RotationMatrix
{
public:
    double angle;
    RotationMatrix(double angle) : angle(angle) {}
    RotationMatrix Inverse()
    {
        return RotationMatrix(fmod(-angle + 2.0 * M_PI, 2.0 * M_PI));
    }
};
inline Tensor2 operator*(const RotationMatrix &matrix, const Tensor2 &tensor)
{
    double s = MySin(matrix.angle);
    double c = MyCos(matrix.angle);
    return Tensor2(tensor[1] * c - tensor[2] * s, tensor[1] * s + tensor[2] * c);
}

template <typename T>
struct LookUpTable
{
    int size = 0;
    std::vector<T> inputs;
    std::vector<T> outputs;

    void Add(T input, T output)
    {
        if (inputs.size() == 0 || input < inputs[0])
        {
            inputs.insert(inputs.begin(), input);
            outputs.insert(outputs.begin(), output);
        }
        else if (input > inputs.back())
        {
            inputs.push_back(input);
            outputs.push_back(output);
        }
        else
        {
            // Find the correct spot where new values must be inserted:
            auto iterator = std::lower_bound(inputs.begin(), inputs.end(), input);

            // Only add new values if they are not in the vectors yet:
            if (iterator == inputs.end() || input < *iterator)
            {
                int index = std::distance(inputs.begin(), iterator);
                inputs.insert(inputs.begin() + index, input);
                outputs.insert(outputs.begin() + index, output);
            }
        }
        size = inputs.size();
    }

    T Evaluate(T x)
    {
        if (size == 0)
            ExitOnError("Trying to read from empty look up table.");
        if (x <= inputs[0] || size == 1)
            return outputs[0];
        if (x >= inputs[size - 1])
            return outputs[size - 1];

        auto lowerIterator = std::lower_bound(inputs.begin(), inputs.end(), x) - 1;
        auto upperIterator = lowerIterator + 1;

        // Handle the case where x is exactly one of the input values
        if (lowerIterator != inputs.end() && *lowerIterator == x)
        {
            int index = std::distance(inputs.begin(), lowerIterator);
            return outputs[index];
        }

        int index0 = std::distance(inputs.begin(), lowerIterator);
        int index1 = std::distance(inputs.begin(), upperIterator);

        T x0 = inputs[index0];
        T x1 = inputs[index1];
        T y0 = outputs[index0];
        T y1 = outputs[index1];

        T t = (x - x0) / (x1 - x0);
        return y0 + t * (y1 - y0);
    }

    void WriteToCsv(std::string name, int resolution)
    {
        if (inputs.size() == 0)
            ExitOnError("Trying to read from empty look up table.");

        std::ofstream fileOut(name + ".csv");

        T x0 = inputs[0];
        T x1 = inputs[inputs.size() - 1];

        fileOut << "#x,y\n";
        if (resolution < 0)
        {
            for (size_t i = 0; i < size; i++)
                fileOut << Format(inputs[i], 8, false, 3) << "," << Format(outputs[i], 8, false, 3) << "\n";
        }
        else
        {
            for (size_t i = 0; i < resolution; i++)
            {
                T t = i / (resolution - 1.0);
                T x = x0 + t * (x1 - x0);
                fileOut << Format(x, 8, false, 3) << "," << Format(Evaluate(x), 8, false, 3) << "\n";
            }
        }

        fileOut.close();
    }

    void Print()
    {
        PrintList(inputs, "inputs");
        PrintList(outputs, "outputs");
    }
};
#endif //__INCLUDE_GUARD_DataTypes_hh__