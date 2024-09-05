#include "Radiation.h"

Radiation::Radiation(Metric &metric, Stencil &stencil, Stencil &streamingStencil, Config config)
    : grid(metric.grid), metric(metric), stencil(stencil), streamingStencil(streamingStencil), config(config), logger(stencil, streamingStencil, metric)
{
    isInitialGridPoint = new bool[grid.nxy]();
    initialE_LF.resize(grid.nxy);
    initialFx_LF.resize(grid.nxy);
    initialFy_LF.resize(grid.nxy);
    initialPxx_LF.resize(grid.nxy);
    initialPxy_LF.resize(grid.nxy);
    initialPyy_LF.resize(grid.nxy);
    initialI.resize(grid.nxy * stencil.nDir);
    initialFluxAngle_IF.resize(grid.nxy);

    rotationAngle.resize(grid.nxy);
    rotationAngleNew.resize(grid.nxy);

    E.resize(grid.nxy);
    Fx.resize(grid.nxy);
    Fy.resize(grid.nxy);
    Pxx.resize(grid.nxy);
    Pxy.resize(grid.nxy);
    Pyy.resize(grid.nxy);
    E_LF.resize(grid.nxy);
    Fx_LF.resize(grid.nxy);
    Fy_LF.resize(grid.nxy);
    Pxx_LF.resize(grid.nxy);
    Pxy_LF.resize(grid.nxy);
    Pyy_LF.resize(grid.nxy);
    F_LF.resize(grid.nxy);

    kappa0.resize(grid.nxy);
    kappa1.resize(grid.nxy);
    kappaA.resize(grid.nxy);
    eta.resize(grid.nxy);
    ux.resize(grid.nxy);
    uy.resize(grid.nxy);

    I.resize(grid.nxy * stencil.nDir);
    Inew.resize(grid.nxy * stencil.nDir);

    coefficientsS.resize(grid.nxy * streamingStencil.nCoefficients);
    coefficientsX.resize(grid.nxy * streamingStencil.nCoefficients);
    coefficientsY.resize(grid.nxy * streamingStencil.nCoefficients);
    coefficientsCx.resize(grid.nxy * streamingStencil.nCoefficients);
    coefficientsCy.resize(grid.nxy * streamingStencil.nCoefficients);

    itterationCount.resize(grid.nxy);

    // Initialize all rotations to identity:
    PARALLEL_FOR(1)
    for (size_t ij = 0; ij < grid.nxy; ij++)
    {
        rotationAngle[ij] = rotationAngleNew[ij] = 0.0;
        itterationCount[ij] = 0;
    }

    // Initialize all Fourier Coefficinets to identity:
    PARALLEL_FOR(1)
    for (size_t ij = 0; ij < grid.nxy; ij++)
    {
            double dataS[streamingStencil.nDir];
            double dataX[streamingStencil.nDir];
            double dataY[streamingStencil.nDir];
            double dataCx[streamingStencil.nDir];
            double dataCy[streamingStencil.nDir];
            for (size_t d = 0; d < streamingStencil.nDir; d++)
            {
                Coord xy = grid.xy(ij);
                dataS[d] = 1.0;
                dataX[d] = xy[1];
                dataY[d] = xy[2];
                dataCx[d] = stencil.Cx(d);
                dataCy[d] = stencil.Cy(d);
            }
            Fourier::GetCoefficients(streamingStencil, dataS, &coefficientsS[HarmonicIndex(0, ij)]);
            Fourier::GetCoefficients(streamingStencil, dataX, &coefficientsX[HarmonicIndex(0, ij)]);
            Fourier::GetCoefficients(streamingStencil, dataY, &coefficientsY[HarmonicIndex(0, ij)]);
            Fourier::GetCoefficients(streamingStencil, dataCx, &coefficientsCx[HarmonicIndex(0, ij)]);
            Fourier::GetCoefficients(streamingStencil, dataCy, &coefficientsCy[HarmonicIndex(0, ij)]);
    }
}
Radiation::~Radiation()
{
    delete[] isInitialGridPoint;
}

size_t Radiation::Index(size_t ij, size_t d)
{
    return d + ij * stencil.nDir;
}
size_t Radiation::Index(size_t i, size_t j, size_t d)
{
    size_t ij = grid.Index(i, j);
    return d + ij * stencil.nDir;
}
size_t Radiation::HarmonicIndex(size_t f, size_t ij)
{
    return f + ij * streamingStencil.nCoefficients;
}

Tensor3 Radiation::InitialDataLFtoIF(size_t ij)
{
    // Transform given initial E_LF, Fx_LF, Fy_LF into initial E_IF, Fx_IF, Fy_IF:
    Tensor3x3 EnergyMomentumTensorLF(initialE_LF[ij], initialFx_LF[ij], initialFy_LF[ij],
                                     initialFx_LF[ij], initialPxx_LF[ij], initialPxy_LF[ij],
                                     initialFy_LF[ij], initialPxy_LF[ij], initialPyy_LF[ij]);
    Tensor3x3 EnergyMomentumTensorIF = TransformLFtoIF(EnergyMomentumTensorLF, metric.GetTetrad(ij).Invert());

    double initialE_IF = EnergyMomentumTensorIF[{0, 0}];
    double initialFx_IF = EnergyMomentumTensorLF[{0, 1}];
    double initialFy_IF = EnergyMomentumTensorLF[{0, 2}];
    double initialPxx_IF = EnergyMomentumTensorLF[{1, 1}];
    double initialPxy_IF = EnergyMomentumTensorLF[{1, 2}];
    double initialPyy_IF = EnergyMomentumTensorLF[{2, 2}];

    // |F| > E is unphysical:
    double initialF_IF = Tensor2(initialFx_IF, initialFy_IF).EuklNorm();
    if (initialF_IF > initialE_IF)
    {
        initialFx_IF *= initialE_IF / initialF_IF;
        initialFy_IF *= initialE_IF / initialF_IF;
    }

    return Tensor3(initialE_IF, initialFx_IF, initialFy_IF);
}
void Radiation::LoadInitialData()
{
    PROFILE_FUNCTION();
    bool isAdaptiveStreaming = (config.streamingType == StreamingType::FlatAdaptive || config.streamingType == StreamingType::CurvedAdaptive);


    double normalizationOverwrite;
    if (0 <= sigmaOverwrite)
        normalizationOverwrite = NormalizationOverwrite();

    PARALLEL_FOR(1)
    for (size_t ij = 0; ij < grid.nxy; ij++)
    {
        if (!isInitialGridPoint[ij])
            continue;

        // Convert given LF initial data to IF:
        Tensor3 initialDataIF = InitialDataLFtoIF(ij);
        double initialE_IF = initialDataIF[0];
        double initialFx_IF = initialDataIF[1];
        double initialFy_IF = initialDataIF[2];

        // Flux direction and magnitude in IF:
        Tensor2 initialFxy_IF(initialFx_IF, initialFy_IF);
        double initialF_IF = initialFxy_IF.EuklNorm();
        double relativeF_IF = initialF_IF / initialE_IF;
        Tensor2 dirInitialF = (isAdaptiveStreaming) ? Tensor2(1, 0) : Tensor2(initialFx_IF / initialF_IF, initialFy_IF / initialF_IF);
        if (initialF_IF < MIN_FLUX_NORM)
            dirInitialF = Tensor2(1, 0);

        double sigma = stencil.fluxToSigmaTable.Evaluate(relativeF_IF);
        double normalization = stencil.fluxToNormalizationTable.Evaluate(relativeF_IF);
        if (0 <= sigmaOverwrite)
        {
            sigma = sigmaOverwrite;
            normalization = normalizationOverwrite;
        }

        if (config.initialDataType == InitialDataType::Moments)
        {
            // Von Mises distribution:
            // https://en.wikipedia.org/wiki/Von_Mises_distribution
            if (!(initialF_IF < MIN_FLUX_NORM))
                rotationAngle[ij] = (isAdaptiveStreaming) ? fmod(MyAtan2(initialFy_IF, initialFx_IF) + 2.0 * M_PI, 2.0 * M_PI) : 0;
            for (size_t d = 0; d < stencil.nDir; d++)
                I[Index(ij, d)] = initialE_IF * exp(sigma * Tensor2::Dot(dirInitialF, stencil.C(d)) - normalization);
        }
        else if (config.initialDataType == InitialDataType::Intensities)
        {
            rotationAngle[ij] = (isAdaptiveStreaming) ? initialFluxAngle_IF[ij] : 0;
            for (size_t d = 0; d < stencil.nDir; d++)
                I[Index(ij, d)] = initialI[Index(ij, d)];
        }
    }
}
double Radiation::NormalizationOverwrite()
{
    double normalization = 0;
    for (int d = 0; d < stencil.nDir; d++)
        normalization += stencil.W(d) * exp(sigmaOverwrite * Tensor2::Dot(Tensor2(1, 0), stencil.C(d)));
    return log(normalization);
}

void Radiation::UpdateFourierCoefficients()
{
    PROFILE_FUNCTION();
    PARALLEL_FOR(2)
    for (size_t j = HALO; j < grid.ny - HALO; j++)
        for (size_t i = HALO; i < grid.nx - HALO; i++)
        {
            size_t ij = grid.Index(i, j);
            double dataS[streamingStencil.nDir];
            double dataX[streamingStencil.nDir];
            double dataY[streamingStencil.nDir];
            double dataCx[streamingStencil.nDir];
            double dataCy[streamingStencil.nDir];
            Coord xy0 = grid.xy(i, j);
            double alpha = metric.GetAlpha(ij);

            for (size_t d = 0; d < streamingStencil.nDir; d++)
            {
                // Initial data for geodesic equation:
                double s = 1;
                Coord xy = xy0;
                Tensor2 c = streamingStencil.C(d);
                Tensor3 uIF(alpha, c[1] * alpha, c[2] * alpha);
                Tensor2 vLF = Vec2ObservedByEulObs<IF, LF>(uIF, xy, metric);

                // Solve geodesic equation backward:
                if (!metric.InsideBH(xy))
                    s *= RK45_GeodesicEquation<-1>(grid.dt, xy, vLF, metric);
                else // inside BH tetrad destroys the velocity stencil. Thus set it to 0.
                    vLF = Tensor2(0.0);

                Tensor2 vIF = TransformLFtoIF(vLF, metric.GetTetradInverse(xy));

                // Final data points for fourier expansion:
                dataS[d] = 1.0 / s; // due to backward integration s must be inverted.
                dataX[d] = xy[1];
                dataY[d] = xy[2];
                dataCx[d] = vIF[1];
                dataCy[d] = vIF[2];
            }
            Fourier::GetCoefficients(streamingStencil, dataS, &coefficientsS[HarmonicIndex(0, ij)]);
            Fourier::GetCoefficients(streamingStencil, dataX, &coefficientsX[HarmonicIndex(0, ij)]);
            Fourier::GetCoefficients(streamingStencil, dataY, &coefficientsY[HarmonicIndex(0, ij)]);
            Fourier::GetCoefficients(streamingStencil, dataCx, &coefficientsCx[HarmonicIndex(0, ij)]);
            Fourier::GetCoefficients(streamingStencil, dataCy, &coefficientsCy[HarmonicIndex(0, ij)]);
        }
}

void Radiation::ComputeMomentsIF()
{
    PROFILE_FUNCTION();
    PARALLEL_FOR(1)
    for (size_t ij = 0; ij < grid.nxy; ij++)
    {
        E[ij] = 0.0;
        Fx[ij] = 0.0;
        Fy[ij] = 0.0;
        Pxx[ij] = 0.0;
        Pxy[ij] = 0.0;
        Pyy[ij] = 0.0;
        for (size_t d = 0; d < stencil.nDir; d++)
        {
            // Skip ghost directions:
            if (stencil.W(d) == 0.0)
                continue;

            Tensor2 dir = RotationMatrix(rotationAngle[ij]) * stencil.C(d);
            size_t index = Index(ij, d);
            double c = stencil.W(d) * I[index];

            E[ij] += c;
            Fx[ij] += c * dir[1];
            Fy[ij] += c * dir[2];
            Pxx[ij] += c * dir[1] * dir[1];
            Pxy[ij] += c * dir[1] * dir[2];
            Pyy[ij] += c * dir[2] * dir[2];
        }
    }
}
void Radiation::ComputeMomentsLF()
{
    PROFILE_FUNCTION();
    PARALLEL_FOR(2)
    for (size_t j = 0; j < grid.ny; j++)
        for (size_t i = 0; i < grid.nx; i++)
        {
            size_t ij = grid.Index(i, j);
            if (metric.InsideBH(grid.xy(i, j)))
            {
                E_LF[ij] = 0.0;
                Fx_LF[ij] = 0.0;
                Fy_LF[ij] = 0.0;
                Pxx_LF[ij] = 0.0;
                Pxy_LF[ij] = 0.0;
                Pyy_LF[ij] = 0.0;
                continue;
            }
            Tensor3x3 EnergyMomentumTensorIF(E[ij], Fx[ij], Fy[ij],
                                             Fx[ij], Pxx[ij], Pxy[ij],
                                             Fy[ij], Pxy[ij], Pyy[ij]);
            Tensor3x3 EnergyMomentumTensorLF = TransformIFtoLF(EnergyMomentumTensorIF, metric.GetTetrad(ij));

            E_LF[ij] = EnergyMomentumTensorLF[{0, 0}];
            Fx_LF[ij] = EnergyMomentumTensorLF[{0, 1}];
            Fy_LF[ij] = EnergyMomentumTensorLF[{0, 2}];
            Pxx_LF[ij] = EnergyMomentumTensorLF[{1, 1}];
            Pxy_LF[ij] = EnergyMomentumTensorLF[{1, 2}];
            Pyy_LF[ij] = EnergyMomentumTensorLF[{2, 2}];
            F_LF[ij] = sqrt(abs(Norm2(Tensor2(Fx_LF[ij], Fy_LF[ij]), metric.GetGamma_ll(ij))));
        }
}

Coord Radiation::GetTempCoordinate(size_t ij, double angle)
{
    Coord xyTemp;
    xyTemp[1] = Fourier::GetValue(angle, &coefficientsX[HarmonicIndex(0, ij)], streamingStencil.nCoefficients);
    xyTemp[2] = Fourier::GetValue(angle, &coefficientsY[HarmonicIndex(0, ij)], streamingStencil.nCoefficients);
    return xyTemp;
}
Tensor2 Radiation::GetTemp2VelocityIF(size_t ij, double angle)
{
    Tensor2 vTempIF;
    vTempIF[1] = Fourier::GetValue(angle, &coefficientsCx[HarmonicIndex(0, ij)], streamingStencil.nCoefficients);
    vTempIF[2] = Fourier::GetValue(angle, &coefficientsCy[HarmonicIndex(0, ij)], streamingStencil.nCoefficients);
    return vTempIF;
}
double Radiation::GetFrequencyShift(size_t ij, double angle)
{
    return Fourier::GetValue(angle, &coefficientsS[HarmonicIndex(0, ij)], streamingStencil.nCoefficients);
}

double Radiation::IntensityAt(size_t ij, Tensor2 vTempIF)
{
    vTempIF = RotationMatrix(rotationAngle[ij]).Inverse() * vTempIF;
    double angle = vTempIF.Phi();
    double d = stencil.interpolationGrid.d(angle);
    size_t d0 = (int)std::floor(d) % stencil.interpolationGrid.nGrid;
    size_t d1 = (d0 + 1) % stencil.interpolationGrid.nGrid;

    std::array<size_t, 4> neighbourIndexesCubic0 = stencil.interpolationGrid.neighbourIndexesCubic[d0];
    std::array<size_t, 4> neighbourIndexesCubic1 = stencil.interpolationGrid.neighbourIndexesCubic[d1];
    std::array<double, 4> neighbourWeightsCubic0 = stencil.interpolationGrid.neighbourWeightsCubic[d0];
    std::array<double, 4> neighbourWeightsCubic1 = stencil.interpolationGrid.neighbourWeightsCubic[d1];

    double value0 = 0;
    double value1 = 0;
    for (size_t k = 0; k < 4; k++)
    {
        value0 += neighbourWeightsCubic0[k] * I[Index(ij, neighbourIndexesCubic0[k])];
        value1 += neighbourWeightsCubic1[k] * I[Index(ij, neighbourIndexesCubic1[k])];
    }
    return std::max(0.0, LinearInterpolation(d - d0, value0, value1));
}

Tensor2 Radiation::AverageF(size_t i, size_t j)
{
    // i,j >= 1, thus i+a etc will never be negative.
    Tensor2 averageF(0.0);
    for (int b = -1; b <= 1; b++)
        for (int a = -1; a <= 1; a++)
        {
            size_t index = grid.Index(i + a, j + b);
            averageF[1] += Fx[index];
            averageF[2] += Fy[index];
        }
    averageF[1] /= 9.0;
    averageF[2] /= 9.0;
    return averageF;
}

void Radiation::UpdateRotationMatrizes()
{
    PROFILE_FUNCTION();
    PARALLEL_FOR(2)
    for (size_t j = HALO; j < grid.ny - HALO; j++)
        for (size_t i = HALO; i < grid.nx - HALO; i++)
        {
            size_t ij = grid.Index(i, j);
            Tensor2 averageF = AverageF(i, j);
            double norm = averageF.EuklNorm();

            // At least 1% of the light points in the direction of the first momentum
            if (E[ij] > MIN_ENERGY_DENSITY && norm / E[ij] > 0.01)
                rotationAngleNew[ij] = averageF.Phi();
            else
                rotationAngleNew[ij] = rotationAngle[ij];
        }
}

void Radiation::StreamFlatFixed()
{
    PROFILE_FUNCTION();
    PARALLEL_FOR(2)
    for (size_t j = HALO; j < grid.ny - HALO; j++)
        for (size_t i = HALO; i < grid.nx - HALO; i++)
        {
            // Index of lattice point ij:
            size_t ij = grid.Index(i, j);
            for (size_t d = 0; d < stencil.nDir; d++)
            {
                // Index of population d at lattice point ij:
                size_t index = Index(ij, d);

                // Get temp velocity:
                Tensor2 direction = stencil.C(d);

                // Get temp lattice point:
                Coord xyTemp = grid.xy(i, j);
                xyTemp[1] -= direction[1] * grid.dt;
                xyTemp[2] -= direction[2] * grid.dt;

                // Get 4 nearest Grid Points:
                double iTemp = grid.i(xyTemp[1]);
                double jTemp = grid.j(xyTemp[2]);
                size_t i0 = std::floor(iTemp);
                size_t i1 = i0 + 1;
                size_t j0 = std::floor(jTemp);
                size_t j1 = j0 + 1;

                Inew[index] = BilinearInterpolation(iTemp - i0, jTemp - j0, I[Index(i0, j0, d)], I[Index(i0, j1, d)], I[Index(i1, j0, d)], I[Index(i1, j1, d)]);
            }
        }
    std::swap(I, Inew);
}
void Radiation::StreamFlatAdaptive()
{
    PROFILE_FUNCTION();
    PARALLEL_FOR(2)
    for (size_t j = HALO; j < grid.ny - HALO; j++)
        for (size_t i = HALO; i < grid.nx - HALO; i++)
        {
            // Index of lattice point ij:
            size_t ij = grid.Index(i, j);
            for (size_t d = 0; d < stencil.nDir; d++)
            {
                // Index of population d at lattice point ij:
                size_t index = Index(ij, d);

                // Get temp velocity:
                Tensor2 direction = RotationMatrix(rotationAngleNew[ij]) * stencil.C(d);

                // Get temp lattice point:
                Coord xyTemp = grid.xy(i, j);
                xyTemp[1] -= direction[1] * grid.dt;
                xyTemp[2] -= direction[2] * grid.dt;

                // Get 4 nearest Grid Points:
                double iTemp = grid.i(xyTemp[1]);
                double jTemp = grid.j(xyTemp[2]);
                size_t i0 = std::floor(iTemp);
                size_t i1 = i0 + 1;
                size_t j0 = std::floor(jTemp);
                size_t j1 = j0 + 1;
                double intensityAt_i0j0 = IntensityAt(grid.Index(i0, j0), direction);
                double intensityAt_i0j1 = IntensityAt(grid.Index(i0, j1), direction);
                double intensityAt_i1j0 = IntensityAt(grid.Index(i1, j0), direction);
                double intensityAt_i1j1 = IntensityAt(grid.Index(i1, j1), direction);
                Inew[index] = BilinearInterpolation(iTemp - i0, jTemp - j0, intensityAt_i0j0, intensityAt_i0j1, intensityAt_i1j0, intensityAt_i1j1);
            }
        }
    std::swap(I, Inew);
    std::swap(rotationAngle, rotationAngleNew);
}
void Radiation::StreamCurvedFixed()
{
    PROFILE_FUNCTION();
    PARALLEL_FOR(2)
    for (size_t j = HALO; j < grid.ny - HALO; j++)
        for (size_t i = HALO; i < grid.nx - HALO; i++)
        {
            // Index of lattice point ij:
            size_t ij = grid.Index(i, j);

            // Skip LPs which are inside BH:
            if (metric.InsideBH(grid.xy(i, j)))
            {
                for (size_t d = 0; d < stencil.nDir; d++)
                    Inew[Index(ij, d)] = 0;
                continue;
            }

            double alpha = metric.GetAlpha(ij);
            for (size_t d = 0; d < stencil.nDir; d++)
            {
                // Index of population d at lattice point ij:
                size_t index = Index(ij, d);

                // Get velocity direction in IF:
                double angle = stencil.Phi(d);

                // Get quantities at emission point:
                double s = GetFrequencyShift(ij, angle);
                Coord xyTemp = GetTempCoordinate(ij, angle);
                Tensor2 vTempIF = GetTemp2VelocityIF(ij, angle);

                // Skip temporary Grid Points inside BH:
                if (metric.InsideBH(xyTemp))
                {
                    Inew[index] = 0;
                    continue;
                }

                // Get 4 nearest Grid Points:
                double iTemp = grid.i(xyTemp[1]);
                double jTemp = grid.j(xyTemp[2]);
                size_t i0 = std::floor(iTemp);
                size_t i1 = i0 + 1;
                size_t j0 = std::floor(jTemp);
                size_t j1 = j0 + 1;

                // Intensity interpolation:
                double intensityAt_i0j0 = IntegerPow<3>(alpha / metric.GetAlpha(grid.Index(i0, j0))) * IntensityAt(grid.Index(i0, j0), vTempIF);
                double intensityAt_i0j1 = IntegerPow<3>(alpha / metric.GetAlpha(grid.Index(i0, j1))) * IntensityAt(grid.Index(i0, j1), vTempIF);
                double intensityAt_i1j0 = IntegerPow<3>(alpha / metric.GetAlpha(grid.Index(i1, j0))) * IntensityAt(grid.Index(i1, j0), vTempIF);
                double intensityAt_i1j1 = IntegerPow<3>(alpha / metric.GetAlpha(grid.Index(i1, j1))) * IntensityAt(grid.Index(i1, j1), vTempIF);

                // Interpolate intensity from neighbouring 4 lattice points to temporary point:
                Inew[index] = IntegerPow<3>(s) * BilinearInterpolation(iTemp - i0, jTemp - j0, intensityAt_i0j0, intensityAt_i0j1, intensityAt_i1j0, intensityAt_i1j1);
            }
        }
    std::swap(I, Inew);
}
void Radiation::StreamCurvedAdaptive()
{
    PROFILE_FUNCTION();
    PARALLEL_FOR(2)
    for (size_t j = HALO; j < grid.ny - HALO; j++)
        for (size_t i = HALO; i < grid.nx - HALO; i++)
        {
            // Index of lattice point ij:
            size_t ij = grid.Index(i, j);

            // Skip LPs which are inside BH:
            if (metric.InsideBH(grid.xy(i, j)))
            {
                for (size_t d = 0; d < stencil.nDir; d++)
                    Inew[Index(ij, d)] = 0;
                continue;
            }

            double alpha = metric.GetAlpha(ij);
            Tensor3x3 tetrad = metric.GetTetrad(ij);
            Tensor2x2 subTetrad = (tetrad[{1,1}], tetrad[{1,2}], tetrad[{2,1}], tetrad[{2,2}]);
            for (size_t d = 0; d < stencil.nDir; d++)
            {
                // Index of population d at lattice point ij:
                size_t index = Index(ij, d);

                // Get velocity direction in IF:
                Tensor2 direction = RotationMatrix(rotationAngleNew[ij]) * stencil.C(d);
                double angle = direction.Phi();

                // Get quantities at emission point:
                double s = GetFrequencyShift(ij, angle);
                Coord xyTemp = GetTempCoordinate(ij, angle);
                Tensor2 vTempIF = GetTemp2VelocityIF(ij, angle);

                // Skip temporary Grid Points inside BH:
                if (metric.InsideBH(xyTemp))
                {
                    Inew[index] = 0;
                    continue;
                }

                // Get 4 nearest Grid Points:
                double iTemp = grid.i(xyTemp[1]);
                double jTemp = grid.j(xyTemp[2]);
                size_t i0 = std::floor(iTemp);
                size_t i1 = i0 + 1;
                size_t j0 = std::floor(jTemp);
                size_t j1 = j0 + 1;

                // Intensity interpolation:
                double intensityAt_i0j0 = IntegerPow<3>(alpha / metric.GetAlpha(grid.Index(i0, j0))) * IntensityAt(grid.Index(i0, j0), vTempIF);
                double intensityAt_i0j1 = IntegerPow<3>(alpha / metric.GetAlpha(grid.Index(i0, j1))) * IntensityAt(grid.Index(i0, j1), vTempIF);
                double intensityAt_i1j0 = IntegerPow<3>(alpha / metric.GetAlpha(grid.Index(i1, j0))) * IntensityAt(grid.Index(i1, j0), vTempIF);
                double intensityAt_i1j1 = IntegerPow<3>(alpha / metric.GetAlpha(grid.Index(i1, j1))) * IntensityAt(grid.Index(i1, j1), vTempIF);

                // Interpolate intensity from neighbouring 4 lattice points to temporary point:
                Inew[index] = IntegerPow<3>(s) * BilinearInterpolation(iTemp - i0, jTemp - j0, intensityAt_i0j0, intensityAt_i0j1, intensityAt_i1j0, intensityAt_i1j1);
            }
        }
    std::swap(I, Inew);
    std::swap(rotationAngle, rotationAngleNew);
}
void Radiation::StreamGeodesicFixed()
{
    // Never use this, it is only a performance test for the paper.
    PROFILE_FUNCTION();
    PARALLEL_FOR(2)
    for (size_t j = HALO; j < grid.ny - HALO; j++)
        for (size_t i = HALO; i < grid.nx - HALO; i++)
        {
            // Index of lattice point ij:
            size_t ij = grid.Index(i, j);

            // Skip LPs which are inside BH:
            if (metric.InsideBH(grid.xy(i, j)))
            {
                for (size_t d = 0; d < stencil.nDir; d++)
                    Inew[Index(ij, d)] = 0;
                continue;
            }

            Coord xy0 = grid.xy(i, j);
            double alpha = metric.GetAlpha(ij);
            for (size_t d = 0; d < stencil.nDir; d++)
            {
                // Index of population d at lattice point ij:
                size_t index = Index(ij, d);

                // Get quantities at emission point:
                double s = 1;
                Coord xyTemp = xy0;
                Tensor2 c = streamingStencil.C(d);
                Tensor3 uIF(alpha, c[1] * alpha, c[2] * alpha);
                Tensor2 vLF = Vec2ObservedByEulObs<IF, LF>(uIF, xyTemp, metric);

                s /= RK45_GeodesicEquation<-1>(grid.dt, xyTemp, vLF, metric);
                Tensor2 vTempIF = TransformLFtoIF(vLF, metric.GetTetradInverse(xyTemp));

                // Skip temporary Grid Points inside BH:
                if (metric.InsideBH(xyTemp))
                {
                    Inew[index] = 0;
                    continue;
                }

                // Get 4 nearest Grid Points:
                double iTemp = grid.i(xyTemp[1]);
                double jTemp = grid.j(xyTemp[2]);
                size_t i0 = std::floor(iTemp);
                size_t i1 = i0 + 1;
                size_t j0 = std::floor(jTemp);
                size_t j1 = j0 + 1;

                // Intensity interpolation:
                double alpha = metric.GetAlpha(ij);
                double intensityAt_i0j0 = IntegerPow<3>(alpha / metric.GetAlpha(grid.Index(i0, j0))) * IntensityAt(grid.Index(i0, j0), vTempIF);
                double intensityAt_i0j1 = IntegerPow<3>(alpha / metric.GetAlpha(grid.Index(i0, j1))) * IntensityAt(grid.Index(i0, j1), vTempIF);
                double intensityAt_i1j0 = IntegerPow<3>(alpha / metric.GetAlpha(grid.Index(i1, j0))) * IntensityAt(grid.Index(i1, j0), vTempIF);
                double intensityAt_i1j1 = IntegerPow<3>(alpha / metric.GetAlpha(grid.Index(i1, j1))) * IntensityAt(grid.Index(i1, j1), vTempIF);

                // Interpolate intensity from neighbouring 4 lattice points to temporary point:
                Inew[index] = s * s * s * BilinearInterpolation(iTemp - i0, jTemp - j0, intensityAt_i0j0, intensityAt_i0j1, intensityAt_i1j0, intensityAt_i1j1);
            }
        }
    std::swap(I, Inew);
}

void Radiation::Collide()
{
    PROFILE_FUNCTION();
    PARALLEL_FOR(2)
    for (size_t j = HALO; j < grid.ny - HALO; j++)
        for (size_t i = HALO; i < grid.nx - HALO; i++)
        {
            if (metric.InsideBH(grid.xy(i, j)))
                continue;
            size_t ij = grid.Index(i, j);
            
            Tensor2 uLF(ux[ij], uy[ij]);
            Tensor2 uIF = TransformLFtoIF(uLF, metric.GetTetrad(ij));

            double alpha = metric.GetAlpha(ij);
            double uu = uIF[1] * uIF[1] + uIF[2] * uIF[2];
            double lorentz = 1.0 / sqrt(1.0 - uu);

            int cycles = 0;
            double diff = 1;
            while (abs(diff) > LAMBDA_ITTERATION_TOLERENCE && (cycles < MAX_LAMBDA_ITERATIONS))
            {
                cycles++;
                // Note that all (previous) moments are in the IF.
                double prevE = E[ij];
                double prevFx = Fx[ij];
                double prevFy = Fy[ij];
                double prevPxx = Pxx[ij];
                double prevPxy = Pxy[ij];
                double prevPyy = Pyy[ij];

                double uF  = uIF[1] *  prevFx + uIF[2] *  prevFy; // u_i F^i
                double uxP = uIF[1] * prevPxx + uIF[2] * prevPxy; // x component of u_i P^ij
                double uyP = uIF[1] * prevPxy + uIF[2] * prevPyy; // y component of u_i P^ij
                double uuP = uIF[1] * uxP + uIF[2] * uyP; // u_i u_j P^ij
                double prevFluidE  = lorentz * lorentz * (prevE - 2.0 * uF + uuP);
                double prevFluidFx = lorentz * lorentz * (prevFx - lorentz * prevE * uIF[1] - uxP + uIF[1] * lorentz / (1.0 + lorentz) * ((2.0 * lorentz + 1.0) * uF - lorentz * uuP));
                double prevFluidFy = lorentz * lorentz * (prevFy - lorentz * prevE * uIF[2] - uyP + uIF[2] * lorentz / (1.0 + lorentz) * ((2.0 * lorentz + 1.0) * uF - lorentz * uuP));

                E[ij] = 0.0;
                Fx[ij] = 0.0;
                Fy[ij] = 0.0;
                Pxx[ij] = 0.0;
                Pxy[ij] = 0.0;
                Pyy[ij] = 0.0;
                for (size_t d = 0; d < stencil.nDir; d++)
                {
                    size_t index = Index(ij, d);

                    // Doppler factor:
                    Tensor2 nIF = RotationMatrix(rotationAngle[ij]) * stencil.C(d);
                    double un = uIF[1] * nIF[1] + uIF[2] * nIF[2];
                    double A = lorentz * (1.0 - un);

                    // nFF contracted with fluxFF:
                    Tensor2 nFF = (nIF - (1 - (lorentz * un) / (lorentz + 1)) * lorentz * uIF) / A;
                    double nF = nFF[1] * prevFluidFx + nFF[2] * prevFluidFy;

                    // Moment collision term:
                    double M = kappa0[ij] * prevFluidE + 3.0 * kappa1[ij] * nF;

                    // Collision:
                    Inew[index] = (I[index] + alpha * grid.dt * (eta[ij] + M) / (A * A)) / (1.0 + alpha * grid.dt * (kappaA[ij] + kappa0[ij]) * A);

                    // Compute new moments:
                    double c = stencil.W(d) * Inew[index];
                    E[ij] += c;
                    Fx[ij] += c * nIF[1];
                    Fy[ij] += c * nIF[2];
                    Pxx[ij] += c * nIF[1] * nIF[1];
                    Pxy[ij] += c * nIF[1] * nIF[2];
                    Pyy[ij] += c * nIF[2] * nIF[2];
                }

                double diffE = IntegerPow<2>((prevE - E[ij]) / E[ij]);
                double diffFx = IntegerPow<2>((prevFx - Fx[ij]) / Fx[ij]);
                double diffFy = IntegerPow<2>((prevFy - Fy[ij]) / Fy[ij]);
                double diffPxx = IntegerPow<2>((prevPxx - Pxx[ij]) / Pxx[ij]);
                double diffPxy = IntegerPow<2>((prevPxy - Pxy[ij]) / Pxy[ij]);
                double diffPyy = IntegerPow<2>((prevPyy - Pyy[ij]) / Pyy[ij]);
                diff = sqrt(diffE + diffFx + diffFy + diffPxx + diffPxy + diffPyy);
            }
            itterationCount[ij] = cycles;
        }

    int maxItteration = 0;
    for(int ij = 0; ij < grid.nxy; ij++)
    {
        maxItteration = std::max(maxItteration, itterationCount[ij]);
        averageItterationCount += itterationCount[ij];
        itterationCount[ij] = 0;
    }
    maxItterationCount = std::max(maxItterationCount, maxItteration);

    std::swap(I, Inew);
}

void Radiation::RunSimulation()
{
    // Initialize Profiler:
    Profiler::Session &session = Profiler::Session::Get();
    session.Start(config.name, "output/" + config.name + "/profileResults.json");

    // -------------------- Initialization --------------------
    LoadInitialData();
    UpdateFourierCoefficients();

    int timeSteps = ceil(config.simTime / grid.dt);
    config.simTime = timeSteps * grid.dt;
    logger.SetValues(config.name, config.simTime);

    // Initial data output:
    if (config.printSetup)
    {
        std::cout << " nx           = " << grid.nx << "\n";
        std::cout << " ny           = " << grid.ny << "\n";
        std::cout << " nDir         = " << stencil.nDir << "\n";
        std::cout << " nFourier     = " << streamingStencil.nDir << "\n";
        std::cout << " simTime      = " << config.simTime << "\n";
        std::cout << " writePeriod  = " << config.writePeriod << "\n";
        std::cout << " dx           = " << grid.dx << "\n";
        std::cout << " dy           = " << grid.dy << "\n";
        std::cout << " dt           = " << grid.dt << "\n";
        std::cout << " timeSteps    = " << logger.timeSteps << "\n";
        std::cout << " filesToWrite = " << std::floor(config.simTime / config.writePeriod) << std::endl;
    }
    // --------------------------------------------------------

    std::chrono::time_point<std::chrono::steady_clock> startTime = std::chrono::steady_clock::now();
    // ----------------- Main simulation Loop -----------------
    {
        PROFILE_SCOPE("Total Time");
        double currentTime = config.t0;
        double timeSinceLastFrame = 0;

        // Save initial data:
        if (config.writeData && config.saveInitialData)
        {
            ComputeMomentsIF();
            ComputeMomentsLF();
            grid.WriteFrametoCsv(currentTime, E_LF, Fx_LF, Fy_LF, F_LF, logger.directoryPath + "/Moments/");
            timeSinceLastFrame = 0;
        }
        for (int n = 0; n < logger.timeSteps; n++)
        {
            if (config.printProgress)
                std::cout << "\n" << n << "," << Format(currentTime, 4) << "," << std::flush;

            // Update stuff:
            if (config.updateFourierHarmonics)
                UpdateFourierCoefficients();
            if (config.keepSourceNodesActive)
                LoadInitialData();
            // Stream:
            switch (config.streamingType)
            {
            case (StreamingType::FlatFixed):
                StreamFlatFixed();
                break;
            case (StreamingType::FlatAdaptive):
                UpdateRotationMatrizes();
                StreamFlatAdaptive();
                break;
            case (StreamingType::CurvedFixed):
                StreamCurvedFixed();
                break;
            case (StreamingType::CurvedAdaptive):
                UpdateRotationMatrizes();
                StreamCurvedAdaptive();
                break;
            case (StreamingType::GeodesicFixed):
                UpdateRotationMatrizes();
                StreamGeodesicFixed();
                break;
            }

            // Collide:
            ComputeMomentsIF();
            Collide();

            currentTime += grid.dt;
            timeSinceLastFrame += grid.dt;

            // Save data:
            if (config.writeData && timeSinceLastFrame >= config.writePeriod - 1e-8)
            {
                ComputeMomentsIF();
                ComputeMomentsLF();
                grid.WriteFrametoCsv(currentTime, E_LF, Fx_LF, Fy_LF, F_LF, logger.directoryPath + "/Moments/");
                timeSinceLastFrame = 0;
            }
        }

        // Save final data:
        if (config.writeData && timeSinceLastFrame != 0)
        {
            ComputeMomentsIF();
            ComputeMomentsLF();
            grid.WriteFrametoCsv(currentTime, E_LF, Fx_LF, Fy_LF, F_LF, logger.directoryPath + "/Moments/");
            timeSinceLastFrame = 0;
        }
    }
    // --------------------------------------------------------
    // Terminate profiler:
    session.End();
    std::chrono::time_point<std::chrono::steady_clock> endTime = std::chrono::steady_clock::now();
    long long start = std::chrono::time_point_cast<std::chrono::microseconds>(startTime).time_since_epoch().count();
    long long end = std::chrono::time_point_cast<std::chrono::microseconds>(endTime).time_since_epoch().count();
    std::cout << "total time: " << (end - start) / 1000000.0 << "s" << std::endl;

    // ---------------------- Termination ---------------------
    if (config.printResults)
        std::cout << std::endl;
        
    logger.maxItterationCount = maxItterationCount;
    averageItterationCount /= (logger.timeSteps * (grid.nx - 2) * (grid.ny - 2));
    logger.averageItterationCount = averageItterationCount;
    if (config.printResults)
    {
        std::cout << "Max Lambda Itteration Count     = " << maxItterationCount << std::endl;
        std::cout << "Average Lambda Itteration Count = " << averageItterationCount << std::endl;
    }

    std::vector<std::string> names = session.GetAllFunctionNames();
    double totalTime;
    double writingTime;
    for (int i = 0; i < names.size(); i++)
    {
        if (config.printResults)
            session.PrintFunctionDuration(names[i]);
        if (names[i] == "Total Time")
            totalTime = session.GetTotalTime(names[i]);
        if (names[i] == "void Grid::WriteFrametoCsv(float, const DoubleBuffer&, const DoubleBuffer&, const DoubleBuffer&, const DoubleBuffer&, std::string, std::string)")
            writingTime = session.GetTotalTime(names[i]);
        logger.AddTimeMeasurement(names[i], session.GetTotalTime(names[i]));
    }
    if (config.printResults)
    {
        std::cout << "Computation Time: " << totalTime - writingTime << "s" << std::endl;
        std::cout << "Benchmark: " << grid.nxy * timeSteps / (1e6 * (totalTime - writingTime)) << "MLUPS" << std::endl;
    }
    if (config.writeData)
        logger.Write();
    // --------------------------------------------------------
}