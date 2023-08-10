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
    initialKappa0.resize(grid.nxy);
    initialKappa1.resize(grid.nxy);
    initialKappaA.resize(grid.nxy);
    initialEta.resize(grid.nxy);
    initialI.resize(grid.nxy * stencil.nDir);
    initialFluxAngle_IF.resize(grid.nxy);

    sigma.resize(grid.nxy);
    normalization.resize(grid.nxy);
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
    kappa0.resize(grid.nxy);
    kappa1.resize(grid.nxy);
    kappaA.resize(grid.nxy);
    eta.resize(grid.nxy);
    I.resize(grid.nxy * stencil.nDir);
    Inew.resize(grid.nxy * stencil.nDir);
    coefficientsS.resize(grid.nxy * streamingStencil.nCoefficients);
    coefficientsX.resize(grid.nxy * streamingStencil.nCoefficients);
    coefficientsY.resize(grid.nxy * streamingStencil.nCoefficients);
    coefficientsCx.resize(grid.nxy * streamingStencil.nCoefficients);
    coefficientsCy.resize(grid.nxy * streamingStencil.nCoefficients);

    // Initialize all Quaternions to identity:
    PARALLEL_FOR(1)
    for (size_t ij = 0; ij < grid.nxy; ij++)
        rotationAngle[ij] = rotationAngleNew[ij] = 0.0;
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
void Radiation::InitSigmaAndNormalization()
{
    PROFILE_FUNCTION();
    bool isAdaptiveStreaming = (config.streamingType == StreamingType::FlatAdaptive || config.streamingType == StreamingType::CurvedAdaptive);

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

        // Flux angle and magnitude in IF:
        Tensor2 initialFxy_IF(initialFx_IF, initialFy_IF);
        double initialF_IF = initialFxy_IF.EuklNorm();
        initialFluxAngle_IF[ij] = initialFxy_IF.Phi();
        Tensor2 dirInitialF = (isAdaptiveStreaming) ? Tensor2(1, 0) : Tensor2(initialFx_IF / initialF_IF, initialFy_IF / initialF_IF);

        // Catch uniform distibution:
        if (initialF_IF < MIN_FLUX_NORM)
        {
            initialFluxAngle_IF[ij] = 0;
            sigma[ij] = 0;
            normalization[ij] = 0;
            for (int d = 0; d < stencil.nDir; d++)
                normalization[ij] += stencil.W(d) * 1.0;
            normalization[ij] = log(normalization[ij]);
            continue;
        }

        // Prepare normalization factor and sigma search:
        sigma[ij] = 0;
        double currentF = 0;
        int refinement = 2; // start with 1e2=100 steps in sigma search.
        double lastGoodSigma = -1;
        while (abs(currentF - initialF_IF) / initialF_IF > 0.001) // while difference bigger 0.1%
        {
            // Adjust sigma:
            if (currentF < initialF_IF)
                // Increase sigma:
                sigma[ij] += pow(10, refinement);
            else if (currentF > initialF_IF)
            {
                // Go back to previous sigma and make sigma step one magnitude smaller:
                sigma[ij] -= pow(10, refinement);
                refinement--;
                if (refinement < -10) // 10 digits after decimal point.
                {
                    sigma[ij] = lastGoodSigma;
                    break;
                }
                currentF = 0; // trigger above branch in next iteration.
                continue;
            }
            else
                // found exact sigma (unlickely).
                break;

            // Determine normalization for current sigma value:
            normalization[ij] = 0;
            for (int d = 0; d < stencil.nDir; d++)
                normalization[ij] += stencil.W(d) * exp(sigma[ij] * Tensor2::Dot(dirInitialF, stencil.C(d)));
            normalization[ij] = log(normalization[ij]);

            // Calculate moments with current sigma value:
            double currentE = 0;
            double currentFx = 0;
            double currentFy = 0;
            for (int d = 0; d < stencil.nDir; d++)
            {
                I[Index(ij, d)] = initialE_IF * exp(sigma[ij] * Tensor2::Dot(dirInitialF, stencil.C(d)) - normalization[ij]);
                Tensor2 dir = RotationMatrix((isAdaptiveStreaming) ? initialFluxAngle_IF[ij] : 0) * stencil.C(d);
                double c = stencil.W(d) * I[Index(ij, d)];
                currentE += c;
                currentFx += c * dir[1];
                currentFy += c * dir[2];
            }
            currentF = Tensor2(currentFx, currentFy).EuklNorm();

            // Debugging:
            // int i = grid.i(ij);
            // int j = grid.j(ij);
            // if ((i == 1) && (j == grid.ny / 2))
            // {
            //     std::cout << Format(currentE) << ", " << Format(currentFx) << ", " << Format(currentFy) << ", " << Format(sigma[ij]) << std::endl;
            // }

            // Check if interpolation error is acceptable:
            double averageError = 0;
            int count = 0;
            for (int d = 0; d < stencil.interpolationGrid.nGrid; d++)
            {
                // Skip angles outside the +-deltaPhi/2 range in which most of the ghost directions are arranged:
                if (acos(Tensor2::Dot(dirInitialF, stencil.interpolationGrid.C(d))) > stencil.deltaPhi / 2.0)
                    continue;

                // Analytic intensity:
                double analyticValue = initialE_IF * exp(sigma[ij] * Tensor2::Dot(dirInitialF, stencil.interpolationGrid.C(d)) - normalization[ij]);

                // Cubic interpolated intensity:
                std::array<size_t, 4> neighbourIndexesCubic = stencil.interpolationGrid.neighbourIndexesCubic[d];
                std::array<double, 4> neighbourWeightsCubic = stencil.interpolationGrid.neighbourWeightsCubic[d];
                double interpolatetValue = 0;
                for (size_t p = 0; p < 4; p++)
                    interpolatetValue += neighbourWeightsCubic[p] * I[Index(ij, neighbourIndexesCubic[p])];

                // Error:
                double error = std::abs((analyticValue - interpolatetValue) / analyticValue);
                averageError += error;
                count++;
            }
            averageError /= count;

            if (averageError > MAX_INTERPOLATION_ERROR)
                currentF = initialF_IF + 1e100; // trigger search refinement!
            else
                lastGoodSigma = sigma[ij];
        }
    }
}

double Radiation::SigmaMax()
{
    double sigmaMax = 0;
    for (int ij = 0; ij < grid.nxy; ij++)
        sigmaMax = std::max(sigmaMax, sigma[ij]);
    return sigmaMax;
}
double Radiation::FluxMax()
{
    double fluxMax = 0;
    for (int ij = 0; ij < grid.nxy; ij++)
    {
        if (!isInitialGridPoint[ij])
            continue;

        Tensor3 initialDataIF = InitialDataLFtoIF(ij);
        double initialE_IF = initialDataIF[0];
        double initialFx_IF = initialDataIF[1];
        double initialFy_IF = initialDataIF[2];

        // Intendet initial data flux norm:
        Tensor2 initialFxy_IF(initialFx_IF, initialFy_IF);
        double initialF_IF = initialFxy_IF.EuklNorm();

        // Actual initial data flux norm:
        double currentFx = 0.0;
        double currentFy = 0.0;
        for (size_t d = 0; d < stencil.nDir; d++)
        {
            Tensor2 dir = RotationMatrix(rotationAngle[ij]) * stencil.C(d);
            size_t index = Index(ij, d);
            double c = stencil.W(d) * I[index];
            currentFx += c * dir[1];
            currentFy += c * dir[2];
        }
        double currentF = Tensor2(currentFx, currentFy).EuklNorm();

        if (initialF_IF != 0)
            fluxMax = std::max(fluxMax, currentF / initialF_IF);
    }
    return fluxMax;
}

void Radiation::LoadInitialData()
{
    PROFILE_FUNCTION();
    bool isAdaptiveStreaming = (config.streamingType == StreamingType::FlatAdaptive || config.streamingType == StreamingType::CurvedAdaptive);

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
        Tensor2 dirInitialF = (isAdaptiveStreaming) ? Tensor2(1, 0) : Tensor2(initialFx_IF / initialF_IF, initialFy_IF / initialF_IF);
        if (initialF_IF < MIN_FLUX_NORM)
            dirInitialF = Tensor2(1, 0);

        // CSG to code unit coversion:
        kappa0[ij] = kappaCGStoCode * initialKappa0[ij];
        kappa1[ij] = kappaCGStoCode * initialKappa1[ij];
        kappaA[ij] = kappaCGStoCode * initialKappaA[ij];
        eta[ij] = etaCGStoCode * initialEta[ij];
        if (config.initialDataType == InitialDataType::Moments)
        {
            // Von Mises distribution:
            // https://en.wikipedia.org/wiki/Von_Mises_distribution
            if (isInitialGridPoint[ij])
            {
                rotationAngle[ij] = (isAdaptiveStreaming) ? initialFluxAngle_IF[ij] : 0;
                for (size_t d = 0; d < stencil.nDir; d++)
                    I[Index(ij, d)] = initialE_IF * exp(sigma[ij] * Tensor2::Dot(dirInitialF, stencil.C(d)) - normalization[ij]);
            }
        }
        else if (config.initialDataType == InitialDataType::Intensities)
        {
            if (isInitialGridPoint[ij])
            {
                rotationAngle[ij] = (isAdaptiveStreaming) ? initialFluxAngle_IF[ij] : 0;
                for (size_t d = 0; d < stencil.nDir; d++)
                    I[Index(ij, d)] = initialE_LF[ij] * initialI[Index(ij, d)];
            }
        }
    }
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

                // Solve geodesic equation backwards:
                if (!metric.InsideBH(xy))
                    s *= RK45_GeodesicEquation<-1>(grid.dt, xy, vLF, metric);
                else // inside BH tetrad destroys the velocity stencil. Thus set it to 0.
                    vLF = Tensor2(0.0);

                Tensor2 vIF = TransformLFtoIF(vLF, metric.GetTetradInverse(xy));

                // Final data points for fourier expansion:
                dataS[d] = 1.0 / s;
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
                double alpha = metric.GetAlpha(ij);
                double intensityAt_i0j0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0, j0))) * IntensityAt(grid.Index(i0, j0), vTempIF);
                double intensityAt_i0j1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0, j1))) * IntensityAt(grid.Index(i0, j1), vTempIF);
                double intensityAt_i1j0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1, j0))) * IntensityAt(grid.Index(i1, j0), vTempIF);
                double intensityAt_i1j1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1, j1))) * IntensityAt(grid.Index(i1, j1), vTempIF);

                // Interpolate intensity from neighbouring 4 lattice points to temporary point:
                Inew[index] = IntegerPow<4>(s) * BilinearInterpolation(iTemp - i0, jTemp - j0, intensityAt_i0j0, intensityAt_i0j1, intensityAt_i1j0, intensityAt_i1j1);
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
                double alpha = metric.GetAlpha(ij);
                double intensityAt_i0j0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0, j0))) * IntensityAt(grid.Index(i0, j0), vTempIF);
                double intensityAt_i0j1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0, j1))) * IntensityAt(grid.Index(i0, j1), vTempIF);
                double intensityAt_i1j0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1, j0))) * IntensityAt(grid.Index(i1, j0), vTempIF);
                double intensityAt_i1j1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1, j1))) * IntensityAt(grid.Index(i1, j1), vTempIF);

                // Interpolate intensity from neighbouring 4 lattice points to temporary point:
                Inew[index] = IntegerPow<4>(s) * BilinearInterpolation(iTemp - i0, jTemp - j0, intensityAt_i0j0, intensityAt_i0j1, intensityAt_i1j0, intensityAt_i1j1);
            }
        }
    std::swap(I, Inew);
    std::swap(rotationAngle, rotationAngleNew);
}

void Radiation::CollideStaticFluidForwardEuler()
{
    PROFILE_FUNCTION();

    PARALLEL_FOR(2)
    for (size_t j = HALO; j < grid.ny - HALO; j++)
        for (size_t i = HALO; i < grid.nx - HALO; i++)
        {
            if (metric.InsideBH(grid.xy(i, j)))
                continue;
            size_t ij = grid.Index(i, j);

            double alpha = metric.GetAlpha(ij);

            for (size_t d = 0; d < stencil.nDir; d++)
            {
                Tensor2 dir = RotationMatrix(rotationAngle[ij]) * stencil.C(d);
                size_t index = Index(ij, d);
                double cDotF = dir[1] * Fx[ij] + dir[2] * Fy[ij];

                double Gamma = eta[ij] + kappa0[ij] * E[ij] + 3.0 * kappa1[ij] * cDotF - I[index] * (kappaA[ij] + kappa0[ij]);
                // I[index] = std::max(I[index] + alpha * grid.dt * Gamma, 0.0);
                I[index] = I[index] + alpha * grid.dt * Gamma;
            }
        }
}
void Radiation::CollideStaticFluidBackwardEuler()
{
    PROFILE_FUNCTION();

    // int highestIteration = 0;
    // std::mutex mutex;
    PARALLEL_FOR(2)
    for (size_t j = HALO; j < grid.ny - HALO; j++)
        for (size_t i = HALO; i < grid.nx - HALO; i++)
        {
            if (metric.InsideBH(grid.xy(i, j)))
                continue;
            size_t ij = grid.Index(i, j);

            double alpha = metric.GetAlpha(ij);
            double guessGammaLinear = 1.0 + alpha * grid.dt * (kappaA[ij] + kappa0[ij]);

            int cycles = 0;
            double diff = 1;
            while (abs(diff) > LAMBDA_ITTERATION_TOLERENCE && (cycles < MAX_LAMBDA_ITERATIONS))
            {
                cycles++;
                double guessE = E[ij];
                double guessFx = Fx[ij];
                double guessFy = Fy[ij];
                double guessPxx = Pxx[ij];
                double guessPxy = Pxy[ij];
                double guessPyy = Pyy[ij];
                double partOfGammaNoneLinear = eta[ij] + kappa0[ij] * guessE;

                E[ij] = 0.0;
                Fx[ij] = 0.0;
                Fy[ij] = 0.0;
                Pxx[ij] = 0.0;
                Pxy[ij] = 0.0;
                Pyy[ij] = 0.0;
                for (size_t d = 0; d < stencil.nDir; d++)
                {
                    Tensor2 dir = RotationMatrix(rotationAngle[ij]) * stencil.C(d);
                    size_t index = Index(ij, d);
                    double cDotGuessF = dir[1] * guessFx + dir[2] * guessFy;

                    double guessGammaNoneLinear = partOfGammaNoneLinear + 3.0 * kappa1[ij] * cDotGuessF;
                    Inew[index] = (I[index] + alpha * grid.dt * guessGammaNoneLinear) / guessGammaLinear;

                    double c = stencil.W(d) * Inew[index];
                    E[ij] += c;
                    Fx[ij] += c * dir[1];
                    Fy[ij] += c * dir[2];
                    Pxx[ij] += c * dir[1] * dir[1];
                    Pxy[ij] += c * dir[1] * dir[2];
                    Pyy[ij] += c * dir[2] * dir[2];
                }
                double diffE = IntegerPow<2>((guessE - E[ij]) / E[ij]);
                double diffFx = IntegerPow<2>((guessFx - Fx[ij]) / Fx[ij]);
                double diffFy = IntegerPow<2>((guessFy - Fy[ij]) / Fy[ij]);
                double diffPxx = IntegerPow<2>((guessPxx - Pxx[ij]) / Pxx[ij]);
                double diffPxy = IntegerPow<2>((guessPxy - Pxy[ij]) / Pxy[ij]);
                double diffPyy = IntegerPow<2>((guessPyy - Pyy[ij]) / Pyy[ij]);
                diff = sqrt(diffE + diffFx + diffFy + diffPxx + diffPxy + diffPyy);
            }
            //{
            //    std::lock_guard<std::mutex> lock(mutex);
            //    highestIteration = std::max(highestIteration, cycles);
            //}
        }
    // std::cout << " highestIteration = " << highestIteration << std::endl;
    std::swap(I, Inew);
}
void Radiation::CollideStaticFluidBackwardEuler2()
{
    // Uses inverse of E matrix instead of only linear term. Does not converge at all!
    PROFILE_FUNCTION();

    double wSum = 0;
    for (int d = 0; d < stencil.nDir; d++)
        wSum += stencil.W(d);

    PARALLEL_FOR(2)
    for (size_t j = HALO; j < grid.ny - HALO; j++)
        for (size_t i = HALO; i < grid.nx - HALO; i++)
        {
            if (metric.InsideBH(grid.xy(i, j)))
                continue;
            size_t ij = grid.Index(i, j);

            // Constants:
            double alpha = metric.GetAlpha(ij);
            double a = 1.0 + kappa0[ij] + kappaA[ij];
            double b = a - wSum;
            double D = a * b;
            double GuessGammaLinear = 1.0 + alpha * grid.dt * (kappaA[ij] + kappa0[ij]);

            int cycles = 0;
            double diff = 1;
            while (abs(diff) > LAMBDA_ITTERATION_TOLERENCE && (cycles < MAX_LAMBDA_ITERATIONS))
            {
                cycles++;
                double guessFx = Fx[ij];
                double guessFy = Fy[ij];

                Fx[ij] = 0.0;
                Fy[ij] = 0.0;
                for (size_t d = 0; d < stencil.nDir; d++)
                {
                    size_t index = Index(ij, d);

                    Inew[index] = 0;
                    for (size_t k = 0; k < stencil.nDir; k++)
                    {
                        Tensor2 dir = RotationMatrix(rotationAngle[ij]) * stencil.C(k);
                        double cDotGuessF = dir[1] * guessFx + dir[2] * guessFy;
                        double GuessGammaNoneLinear = eta[ij] + 3.0 * kappa1[ij] * cDotGuessF;
                        Inew[index] += (stencil.W(d) + (d == k) * b) * (I[Index(ij, k)] + alpha * grid.dt * GuessGammaNoneLinear);
                    }
                    Inew[index] /= D;

                    Tensor2 dir = RotationMatrix(rotationAngle[ij]) * stencil.C(d);
                    double c = stencil.W(d) * Inew[index];
                    Fx[ij] += c * dir[1];
                    Fy[ij] += c * dir[2];
                }
                double diffFx = abs((guessFx - Fx[ij]) / Fx[ij]);
                double diffFy = abs((guessFy - Fy[ij]) / Fy[ij]);
                diff = diffFx + diffFy;
            }
        }
    std::swap(I, Inew);
}
void Radiation::CollideForwardEuler()
{
    PROFILE_FUNCTION();

    PARALLEL_FOR(2)
    for (size_t j = HALO; j < grid.ny - HALO; j++)
        for (size_t i = HALO; i < grid.nx - HALO; i++)
        {
            if (metric.InsideBH(grid.xy(i, j)))
                continue;
            size_t ij = grid.Index(i, j);

            double alpha = metric.GetAlpha(ij);
            Tensor2 u(0.0);                                                                              // fluid 2 velocity as seen by Eulerian observer
            double W = 1.0 / sqrt(1.0 - Norm2(u, metric.GetGamma_ll(ij)));                               // Lorentz factor
            double uDotF = u[1] * Fx[ij] + u[2] * Fy[ij];                                                // u_i F^i
            double uuDotP = u[1] * u[1] * Pxx[ij] + u[2] * u[2] * Pyy[ij] + 2.0 * u[1] * u[2] * Pxy[ij]; // u_i u_j P^ij
            double uDotPx = u[1] * Pxx[ij] + u[2] * Pxy[ij];                                             // u_j P^ij, i=1
            double uDotPy = u[1] * Pxy[ij] + u[2] * Pyy[ij];                                             // u_j P^ij, i=2
            double fluidE = W * W * (E[ij] - 2.0 * uDotF + uuDotP);
            double fluidFx = W * W * W * (2.0 * uDotF - E[ij] - uuDotP) * u[1] + W * (Fx[ij] - uDotPx);
            double fluidFy = W * W * W * (2.0 * uDotF - E[ij] - uuDotP) * u[2] + W * (Fy[ij] - uDotPy);

            for (size_t d = 0; d < stencil.nDir; d++)
            {
                Tensor2 dir = RotationMatrix(rotationAngle[ij]) * stencil.C(d);
                size_t index = Index(ij, d);
                double A = W * (1.0 - Tensor2::Dot(dir, u));
                double cDotFluidF = dir[1] * fluidFx + dir[2] * fluidFy;

                double Gamma = (eta[ij] + kappa0[ij] * fluidE + 3.0 * kappa1[ij] * cDotFluidF) / (A * A * A) - A * I[index] * (kappaA[ij] + kappa0[ij]);
                I[index] = std::max(I[index] + alpha * grid.dt * Gamma, 0.0);
            }
        }
}
void Radiation::CollideBackwardEuler()
{
    PROFILE_FUNCTION();

    int highestIteration = 0;
    // std::mutex mutex;
    PARALLEL_FOR(2)
    for (size_t j = HALO; j < grid.ny - HALO; j++)
        for (size_t i = HALO; i < grid.nx - HALO; i++)
        {
            if (metric.InsideBH(grid.xy(i, j)))
                continue;
            size_t ij = grid.Index(i, j);

            double alpha = metric.GetAlpha(ij);
            double guessGammaLinear = 1.0 + alpha * grid.dt * (kappaA[ij] + kappa0[ij]);
            Tensor2 u(0.0);                                                // fluid 2 velocity as seen by Eulerian observer
            double W = 1.0 / sqrt(1.0 - Norm2(u, metric.GetGamma_ll(ij))); // Lorentz factor

            int cycles = 0;
            double diff = 1;
            while (abs(diff) > LAMBDA_ITTERATION_TOLERENCE && (cycles < MAX_LAMBDA_ITERATIONS))
            {
                cycles++;
                double guessE = E[ij];
                double guessFx = Fx[ij];
                double guessFy = Fy[ij];
                double guessPxx = Pxx[ij];
                double guessPxy = Pxy[ij];
                double guessPyy = Pyy[ij];

                double uDotF = u[1] * guessFx + u[2] * guessFy;                                                 // u_i F^i
                double uuDotP = u[1] * u[1] * guessPxx + u[2] * u[2] * guessPyy + 2.0 * u[1] * u[2] * guessPxy; // u_i u_j P^ij
                double uDotPx = u[1] * guessPxx + u[2] * guessPxy;                                              // u_j P^ij, i=1
                double uDotPy = u[1] * guessPxy + u[2] * guessPyy;                                              // u_j P^ij, i=2
                double guessFluidE = W * W * (guessE - 2.0 * uDotF + uuDotP);
                double guessFluidFx = W * W * W * (2.0 * uDotF - guessE - uuDotP) * u[1] + W * (guessFx - uDotPx);
                double guessFluidFy = W * W * W * (2.0 * uDotF - guessE - uuDotP) * u[2] + W * (guessFy - uDotPy);
                double partOfGammaNoneLinear = eta[ij] + kappa0[ij] * guessFluidE;

                E[ij] = 0.0;
                Fx[ij] = 0.0;
                Fy[ij] = 0.0;
                Pxx[ij] = 0.0;
                Pxy[ij] = 0.0;
                Pyy[ij] = 0.0;
                for (size_t d = 0; d < stencil.nDir; d++)
                {
                    Tensor2 dir = RotationMatrix(rotationAngle[ij]) * stencil.C(d);
                    size_t index = Index(ij, d);
                    double A = W * (1.0 - Tensor2::Dot(dir, u));
                    double cDotGuessFluidF = dir[1] * guessFluidFx + dir[2] * guessFluidFy;

                    double guessGammaNoneLinear = (partOfGammaNoneLinear + 3.0 * kappa1[ij] * cDotGuessFluidF) / (A * A * A);
                    Inew[index] = (I[index] + alpha * grid.dt * guessGammaNoneLinear) / guessGammaLinear;

                    double c = stencil.W(d) * Inew[index];
                    E[ij] += c;
                    Fx[ij] += c * dir[1];
                    Fy[ij] += c * dir[2];
                    Pxx[ij] += c * dir[1] * dir[1];
                    Pxy[ij] += c * dir[1] * dir[2];
                    Pyy[ij] += c * dir[2] * dir[2];
                }
                double diffE = abs((guessE - E[ij]) / E[ij]);
                double diffFx = abs((guessFx - Fx[ij]) / Fx[ij]);
                double diffFy = abs((guessFy - Fy[ij]) / Fy[ij]);
                double diffPxx = abs((guessPxx - Pxx[ij]) / Pxx[ij]);
                double diffPxy = abs((guessPxy - Pxy[ij]) / Pxy[ij]);
                double diffPyy = abs((guessPyy - Pyy[ij]) / Pyy[ij]);
                diff = diffE + diffFx + diffFy + diffPxx + diffPxy + diffPyy;
            }
            //{
            //    std::lock_guard<std::mutex> lock(mutex);
            //    highestIteration = std::max(highestIteration, cycles);
            //}
        }
    // std::cout << " highestIteration = " << highestIteration << std::endl;
    std::swap(I, Inew);
}

void Radiation::RunSimulation()
{
    // Initialize Profiler:
    Profiler::Session &session = Profiler::Session::Get();
    session.Start(config.name, "output/" + config.name + "/profileResults.json");

    // -------------------- Initialization --------------------
    InitSigmaAndNormalization();
    LoadInitialData();
    UpdateFourierCoefficients();

    int timeSteps = ceil(config.simTime / grid.dt);
    config.simTime = timeSteps * grid.dt;
    logger.SetValues(config.name, config.simTime, MAX_INTERPOLATION_ERROR, SigmaMax(), FluxMax());

    // Initial data output:
    if (config.printToTerminal)
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

    // ----------------- Main simulation Loop -----------------
    {
        PROFILE_SCOPE("Total Time");
        double currentTime = config.t0;
        double timeSinceLastFrame = 0;

        // Save initial data:
        if (config.writeData)
        {
            ComputeMomentsIF();
            ComputeMomentsLF();
            grid.WriteFrametoCsv(currentTime, E_LF, Fx_LF, Fy_LF, E_LF, logger.directoryPath + "/Moments/");
            timeSinceLastFrame = 0;
        }

        for (int n = 0; n < logger.timeSteps; n++)
        {
            if (config.printToTerminal)
                std::cout << "\n"
                          << n << "," << Format(currentTime, 4) << "," << std::flush;

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
            }

            // Collide:
            ComputeMomentsIF();
            // CollideStaticFluidForwardEuler();
            CollideStaticFluidBackwardEuler();
            // CollideStaticFluidBackwardEuler2();
            // CollideForwardEuler();
            // CollideBackwardEuler();

            currentTime += grid.dt;
            timeSinceLastFrame += grid.dt;

            // Save data:
            if (config.writeData && timeSinceLastFrame >= config.writePeriod)
            {
                ComputeMomentsLF();
                grid.WriteFrametoCsv(currentTime, E_LF, Fx_LF, Fy_LF, E_LF, logger.directoryPath + "/Moments/");
                timeSinceLastFrame = 0;
            }

            // Update other stuff:
            if (config.updateFourierHarmonics)
                UpdateFourierCoefficients();
            if (config.keepSourceNodesActive)
                LoadInitialData();
        }

        // Save final data:
        if (config.writeData)
        {
            ComputeMomentsIF();
            ComputeMomentsLF();
            grid.WriteFrametoCsv(currentTime, E_LF, Fx_LF, Fy_LF, E_LF, logger.directoryPath + "/Moments/");
            timeSinceLastFrame = 0;
        }
    }
    // --------------------------------------------------------
    // Terminate profiler:
    session.End();

    // ---------------------- Termination ---------------------
    std::vector<std::string> names = session.GetAllFunctionNames();
    if (config.printToTerminal)
    {
        std::cout << std::endl;
    }
    for (int i = 0; i < names.size(); i++)
    {
        if (config.printToTerminal)
            session.PrintFunctionDuration(names[i]);
        logger.AddTimeMeasurement(names[i], session.GetTotalTime(names[i]));
    }
    // --------------------------------------------------------
}