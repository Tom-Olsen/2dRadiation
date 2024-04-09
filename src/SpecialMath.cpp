#include "SpecialMath.h"

double Dot(const Tensor2 &vector0, const Tensor2 &vector1, const Tensor2x2 &gamma_ll)
{
    double dot = 0;
    for (int i = 1; i < 3; i++)
        for (int j = 1; j < 3; j++)
            dot += gamma_ll[{i, j}] * vector0[i] * vector1[j];
    return dot;
}
double Dot(const Tensor3 &vector0, const Tensor3 &vector1, const Tensor3x3 &g_ll)
{
    double dot = 0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            dot += g_ll[{i, j}] * vector0[i] * vector1[j];
    return dot;
}

double Norm2(const Tensor2 &vector, const Tensor2x2 &gamma_ll)
{
    double norm2 = 0;
    for (int i = 1; i < 3; i++)
        for (int j = 1; j < 3; j++)
            norm2 += gamma_ll[{i, j}] * vector[i] * vector[j];
    return norm2;
}

double Norm2(const Tensor3 &vector, const Tensor3x3 &g_ll)
{
    double norm2 = 0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            norm2 += g_ll[{i, j}] * vector[i] * vector[j];
    return norm2;
}

Tensor3 NullNormalize(const Tensor3 &vector, const Tensor3x3 &g_ll)
{
    double a = 0;
    for (int i = 1; i < 3; i++)
        for (int j = 1; j < 3; j++)
            a += g_ll[rank2Indices{i, j}] * vector[i] * vector[j];
    double b = 0;
    for (int i = 1; i < 3; i++)
        b += g_ll[rank2Indices{0, i}] * vector[0] * vector[i];
    double c = g_ll[rank2Indices{0, 0}] * vector[0] * vector[0];
    double d = -b / a + sqrt((b * b) / (a * a) - c / a);

    return Tensor3(vector[0], vector[1] * d, vector[2] * d);
}

Tensor2 TransformIFtoLF(const Tensor2 &vector, const Tensor3x3 &tetrad)
{
    return Tensor2(tetrad[{1, 1}] * vector[1] + tetrad[{1, 2}] * vector[2],
                   tetrad[{2, 1}] * vector[1] + tetrad[{2, 2}] * vector[2]);
}
Tensor2 TransformLFtoIF(const Tensor2 &vector, const Tensor3x3 &tetradInverse)
{
    return Tensor2(tetradInverse[{1, 1}] * vector[1] + tetradInverse[{1, 2}] * vector[2],
                   tetradInverse[{2, 1}] * vector[1] + tetradInverse[{2, 2}] * vector[2]);
}

Tensor3 TransformIFtoLF(const Tensor3 &vector, const Tensor3x3 &tetrad)
{
    return Tensor3(tetrad[{0, 0}] * vector[0] + tetrad[{0, 1}] * vector[1] + tetrad[{0, 2}] * vector[2],
                   tetrad[{1, 0}] * vector[0] + tetrad[{1, 1}] * vector[1] + tetrad[{1, 2}] * vector[2],
                   tetrad[{2, 0}] * vector[0] + tetrad[{2, 1}] * vector[1] + tetrad[{2, 2}] * vector[2]);
}
Tensor3 TransformLFtoIF(const Tensor3 &vector, const Tensor3x3 &tetradInverse)
{
    return Tensor3(tetradInverse[{0, 0}] * vector[0] + tetradInverse[{0, 1}] * vector[1] + tetradInverse[{0, 2}] * vector[2],
                   tetradInverse[{1, 0}] * vector[0] + tetradInverse[{1, 1}] * vector[1] + tetradInverse[{1, 2}] * vector[2],
                   tetradInverse[{2, 0}] * vector[0] + tetradInverse[{2, 1}] * vector[1] + tetradInverse[{2, 2}] * vector[2]);
}

Tensor3x3 TransformIFtoLF(const Tensor3x3 &tensor, const Tensor3x3 &tetrad)
{
    Tensor3x3 result(0.0);
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            for (int A = 0; A < 3; A++)
                for (int B = 0; B < 3; B++)
                    result[{a, b}] += tensor[{A, B}] * tetrad[{a, A}] * tetrad[{b, B}];
    return result;
}
Tensor3x3 TransformLFtoIF(const Tensor3x3 &tensor, const Tensor3x3 &tetradInverse)
{
    Tensor3x3 result(0.0);
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            for (int A = 0; A < 3; A++)
                for (int B = 0; B < 3; B++)
                    result[{a, b}] += tensor[{A, B}] * tetradInverse[{a, A}] * tetradInverse[{b, B}];
    return result;
}

template <class FrameIn, class FrameOut>
Tensor2 Vec2ObservedByEulObs(const Tensor3 &u, const Coord &xy, Metric &metric)
{
    if constexpr (std::is_same<FrameIn, IF>::value && std::is_same<FrameOut, IF>::value)
    { // IF -> IF
        double alpha = metric.GetAlpha(xy);
        return Tensor2(u[1] / alpha, u[2] / alpha);
    }
    if constexpr (std::is_same<FrameIn, IF>::value && std::is_same<FrameOut, LF>::value)
    { // IF -> LF
        double alpha = metric.GetAlpha(xy);
        Tensor2 v(u[1] / alpha, u[2] / alpha);
        return TransformIFtoLF(v, metric.GetTetrad(xy));
    }
    if constexpr (std::is_same<FrameIn, LF>::value && std::is_same<FrameOut, IF>::value)
    { // LF -> IF
        double alpha = metric.GetAlpha(xy);
        Tensor2 beta_u = metric.GetBeta_u(xy);
        Tensor2 v((u[1] + beta_u[1]) / alpha, (u[2] + beta_u[2]) / alpha);
        return TransformLFtoIF(v, metric.GetTetradInverse(xy));
    }
    if constexpr (std::is_same<FrameIn, LF>::value && std::is_same<FrameOut, LF>::value)
    { // LF -> LF
        double alpha = metric.GetAlpha(xy);
        Tensor2 beta_u = metric.GetBeta_u(xy);
        return Tensor2((u[1] + beta_u[1]) / alpha, (u[2] + beta_u[2]) / alpha);
    }
}

template Tensor2 Vec2ObservedByEulObs<IF, IF>(const Tensor3 &u, const Coord &xy, Metric &metric);
template Tensor2 Vec2ObservedByEulObs<IF, LF>(const Tensor3 &u, const Coord &xy, Metric &metric);
template Tensor2 Vec2ObservedByEulObs<LF, IF>(const Tensor3 &u, const Coord &xy, Metric &metric);
template Tensor2 Vec2ObservedByEulObs<LF, LF>(const Tensor3 &u, const Coord &xy, Metric &metric);