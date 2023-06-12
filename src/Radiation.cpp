#include "Radiation.h"

#define HALO 2

Radiation::Radiation(Metric& metric, Stencil& stencil, Stencil& streamingStencil, StreamingType streamingType):
grid(metric.grid), metric(metric), stencil(stencil), streamingStencil(streamingStencil), streamingType(streamingType)
{
	isInitialGridPoint = new bool[grid.nxy]();
	initialE.resize(grid.nxy);
	initialNx.resize(grid.nxy);
	initialNy.resize(grid.nxy);
	initialKappa0.resize(grid.nxy);
	initialKappa1.resize(grid.nxy);
	initialKappaA.resize(grid.nxy);
	initialEta.resize(grid.nxy);
	rotationAngle.resize(grid.nxy);
	rotationAngleNew.resize(grid.nxy);
	E.resize(grid.nxy);
	Fx.resize(grid.nxy);
	Fy.resize(grid.nxy);
	Pxx.resize(grid.nxy);
	Pxy.resize(grid.nxy);
	Pyy.resize(grid.nxy);
	E_LF .resize(grid.nxy);
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
	for(size_t ij=0; ij<grid.nxy; ij++)
        rotationAngle[ij] = rotationAngleNew[ij] = 0.0;
}
Radiation::~Radiation()
{
	delete[] isInitialGridPoint;
}



size_t Radiation::Index(size_t ij, size_t d)
{
	#ifdef ijd
		return ij + d * grid.nxy;
	#endif
	#ifdef dij
		return d + ij * stencil.nDir;
	#endif
}
size_t Radiation::Index(size_t i, size_t j, size_t d)
{
	size_t ij = grid.Index(i,j);
	#ifdef ijd
		return ij + d * grid.nxy;
	#endif
	#ifdef dij
		return d + ij * stencil.nDir;
	#endif
}
size_t Radiation::HarmonicIndex(size_t f, size_t ij)
{
	return f + ij * streamingStencil.nCoefficients;
}



void Radiation::NormalizeInitialDirections()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(1)
	for(size_t ij=0; ij<grid.nxy; ij++)
	{
		Tensor2 n(initialNx[ij], initialNy[ij]);
		double norm = n.EuklNorm();

		// If norm to low we assume no direction.
		if(norm < 0.0001)
		{
			initialNx[ij] = 0;
			initialNy[ij] = 0;
		}
		else
		{
			initialNx[ij] = n[1] / norm;
			initialNy[ij] = n[2] / norm;
		}
	}
}



void Radiation::LoadInitialData()
{
	PROFILE_FUNCTION();
	bool isAdaptiveStreaming = (streamingType == StreamingType::FlatAdaptive || streamingType == StreamingType::CurvedAdaptive);
	srand((unsigned) time(NULL));

	PARALLEL_FOR(1)
	for(size_t ij=0; ij<grid.nxy; ij++)
	{
		kappa0[ij] = initialKappa0[ij];
		kappa1[ij] = initialKappa1[ij];
		kappaA[ij] = initialKappaA[ij];
		eta[ij] = initialEta[ij];
		
		if(isInitialGridPoint[ij])
		{
			if(initialNx[ij] == 0 && initialNy[ij] == 0)
			{// Uniform intensity distribution:
				if (isAdaptiveStreaming)
				{// Random direction for uniformity:
					Tensor2 n;
                    double norm = 2;
                    while(norm > 1)
                    {
					    n[1] = RandomRange(-1,1);
					    n[2] = RandomRange(-1,1);
                        norm = n.EuklNorm();
                    }
					rotationAngle[ij] = n.Phi();
				}
				else	// (0,0) is an invalid direction. Set to identity:
					rotationAngle[ij] = 0.0;

				for(int d=0; d<stencil.nDir; d++)
					I[Index(ij,d)] = initialE[ij];
			}
			else
			{// Von Mises distribution:
			 // https://en.wikipedia.org/wiki/Von_Mises_distribution
				// In case of adaptive streaming the stencil is rotated such that its east side points towards n, whenever stencil.C is used.
				// This means that the Von Mises distribution needs to point east for the adaptive streaming case.
                if (isAdaptiveStreaming)
                {
				    rotationAngle[ij] = Tensor2(initialNx[ij],initialNy[ij]).Phi();
				    for(size_t d=0; d<stencil.nDir; d++)
				    	I[Index(ij,d)] = initialE[ij] * exp(sigma * MyCos(stencil.Phi(d) - 0.0));
                }
                else
                {
				    rotationAngle[ij] = 0.0;
                    double vonMisesAngle = Tensor2(initialNx[ij],initialNy[ij]).Phi();
				    for(size_t d=0; d<stencil.nDir; d++)
                    {
                        double angle = fmod(stencil.Phi(d) - vonMisesAngle + 2.0 * M_PI, 2.0 * M_PI);
				    	I[Index(ij,d)] = initialE[ij] * exp(sigma * MyCos(angle));
                    }
                }
			}
		}
	}
}



void Radiation::NormalizeInitialIntensities()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(1)
	for(size_t ij=0; ij<grid.nxy; ij++)
		if(isInitialGridPoint[ij])
			initialE[ij] *= initialE[ij] / E_LF[ij];
}



void Radiation::UpdateFourierCoefficients()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(2)
	for(size_t j=HALO; j<grid.ny-HALO; j++)
	for(size_t i=HALO; i<grid.nx-HALO; i++)
	{
		size_t ij = grid.Index(i,j);
		double dataS[streamingStencil.nDir];
		double dataX[streamingStencil.nDir];
		double dataY[streamingStencil.nDir];
    	double dataCx[streamingStencil.nDir];
    	double dataCy[streamingStencil.nDir];
		Coord xy0 = grid.xy(i,j);
		double alpha = metric.GetAlpha(ij);

		for(size_t d=0; d<streamingStencil.nDir; d++)
    	{
			// Initial data for geodesic equation:
			double s = 1;
			Coord xy = xy0;
            Tensor2 c = streamingStencil.C(d);
            Tensor3 uIF(alpha, c[1] * alpha, c[2] * alpha);
            Tensor2 vLF = Vec2ObservedByEulObs<IF,LF>(uIF, xy, metric);

			// Solve geodesic equation backwards:
			if(!metric.InsideBH(xy))
    	    	s *= RK45_GeodesicEquation<-1>(grid.dt, xy, vLF, metric);
			else // inside BH tetrad destroys the velocity stencil. Thus set it to 0.
				vLF = Tensor2(0.0);

			Tensor2 vIF = TransformLFtoIF(vLF,metric.GetTetradInverse(xy));

			// Final data points for fourier expansion:
			dataS[d] = 1.0/s;
			dataX[d] = xy[1];
			dataY[d] = xy[2];
			dataCx[d] = vIF[1];
			dataCy[d] = vIF[2];
    	}
		Fourier::GetCoefficients(streamingStencil, dataS , &coefficientsS [HarmonicIndex(0,ij)]);
		Fourier::GetCoefficients(streamingStencil, dataX , &coefficientsX [HarmonicIndex(0,ij)]);
		Fourier::GetCoefficients(streamingStencil, dataY , &coefficientsY [HarmonicIndex(0,ij)]);
		Fourier::GetCoefficients(streamingStencil, dataCx, &coefficientsCx[HarmonicIndex(0,ij)]);
		Fourier::GetCoefficients(streamingStencil, dataCy, &coefficientsCy[HarmonicIndex(0,ij)]);
	}
}



void Radiation::ComputeMomentsIF()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(1)
	for(size_t ij=0; ij<grid.nxy; ij++)
	{
		E[ij]   = 0.0;
		Fx[ij]  = 0.0;
		Fy[ij]  = 0.0;
		Pxx[ij] = 0.0;
		Pxy[ij] = 0.0;
		Pyy[ij] = 0.0;
		for(size_t d=0; d<stencil.nDir; d++)
		{
			Tensor2 dir = RotationMatrix(rotationAngle[ij]) * stencil.C(d);
			size_t index = Index(ij,d);
            double c = stencil.W(d) * I[index];
			E[ij]   += c;
			Fx[ij]  += c * dir[1];
			Fy[ij]  += c * dir[2];
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
	for(size_t j=0; j<grid.ny; j++)
	for(size_t i=0; i<grid.nx; i++)
	{
        int ij = grid.Index(i,j);
        if(metric.InsideBH(grid.xy(i,j)))
        {
		    E_LF[ij]   = 0.0;
		    Fx_LF[ij]  = 0.0;
		    Fy_LF[ij]  = 0.0;
		    Pxx_LF[ij] = 0.0;
		    Pxy_LF[ij] = 0.0;
		    Pyy_LF[ij] = 0.0;
            continue;
        }
		Tensor3x3 EnergyMomentumTensorIF
		( E[ij], Fx[ij], Fy[ij],
		 Fx[ij],Pxx[ij],Pxy[ij],
		 Fy[ij],Pxy[ij],Pyy[ij]);
		Tensor3x3 EnergyMomentumTensorLF = TransformIFtoLF(EnergyMomentumTensorIF, metric.GetTetrad(ij));

		E_LF[ij]   = EnergyMomentumTensorLF[{0,0}];
		Fx_LF[ij]  = EnergyMomentumTensorLF[{0,1}];
		Fy_LF[ij]  = EnergyMomentumTensorLF[{0,2}];
		Pxx_LF[ij] = EnergyMomentumTensorLF[{1,1}];
		Pxy_LF[ij] = EnergyMomentumTensorLF[{1,2}];
		Pyy_LF[ij] = EnergyMomentumTensorLF[{2,2}];
	}
}



Coord Radiation::GetTempCoordinate(size_t ij, double angle)
{
	Coord xyTemp;
	xyTemp[1] = Fourier::GetValue(angle, &coefficientsX[HarmonicIndex(0,ij)], streamingStencil.nCoefficients);
	xyTemp[2] = Fourier::GetValue(angle, &coefficientsY[HarmonicIndex(0,ij)], streamingStencil.nCoefficients);
	return xyTemp;
}
Tensor2 Radiation::GetTemp2VelocityIF(size_t ij, double angle)
{
	Tensor2 vTempIF;
	vTempIF[1] = Fourier::GetValue(angle, &coefficientsCx[HarmonicIndex(0,ij)], streamingStencil.nCoefficients);
	vTempIF[2] = Fourier::GetValue(angle, &coefficientsCy[HarmonicIndex(0,ij)], streamingStencil.nCoefficients);
	return vTempIF;
}
double Radiation::GetFrequencyShift(size_t ij, double angle)
{
	return Fourier::GetValue(angle, &coefficientsS[HarmonicIndex(0,ij)], streamingStencil.nCoefficients);
}



double Radiation::IntensityAt(size_t ij, Tensor2 vTempIF)
{
    vTempIF = RotationMatrix(rotationAngle[ij]).Inverse() * vTempIF;
    double angle = vTempIF.Phi();
    double d = stencil.interpolationGrid.d(angle);
    size_t d0 = (int)std::floor(d) % stencil.interpolationGrid.nGrid;
    size_t d1 = (d0 + 1) % stencil.interpolationGrid.nGrid;

    std::array<size_t,4> neighbourIndexes0 = stencil.interpolationGrid.neighbourIndexes[d0];
    std::array<size_t,4> neighbourIndexes1 = stencil.interpolationGrid.neighbourIndexes[d1];
    std::array<double,4> neighbourWeights0 = stencil.interpolationGrid.neighbourWeights[d0];
    std::array<double,4> neighbourWeights1 = stencil.interpolationGrid.neighbourWeights[d1];

    double value0 = 0;
    double value1 = 0;
    for(size_t k=0; k<4; k++)
    {
        value0 += neighbourWeights0[k] * I[Index(ij,neighbourIndexes0[k])];
        value1 += neighbourWeights1[k] * I[Index(ij,neighbourIndexes1[k])];
    }
    return std::max(0.0, LinearInterpolation(d-d0, value0, value1));
}



Tensor2 Radiation::AverageF(size_t i, size_t j)
{
	// i,j >= 1, thus i+a etc will never be negative.
	Tensor2 averageF(0.0);
	for(int b=-1; b<=1; b++)
	for(int a=-1; a<=1; a++)
	{
		size_t index = grid.Index(i+a, j+b);
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
	for(size_t j=HALO; j<grid.ny-HALO; j++)
	for(size_t i=HALO; i<grid.nx-HALO; i++)
	{
		size_t ij = grid.Index(i,j);
		Tensor2 averageF = AverageF(i,j);
		double norm = averageF.EuklNorm();

		// At least 1% of the light points in the direction of the first momentum
		if (E[ij] > 1e-16 && norm / E[ij] > 0.01)
            rotationAngleNew[ij] = averageF.Phi();
		else
            rotationAngleNew[ij] = rotationAngle[ij];
	}
}



void Radiation::StreamFlatFixed()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(3)
	#ifdef ijd
	for(size_t d=0; d<stencil.nDir; d++)
	for(size_t j=HALO; j<grid.ny-HALO; j++)
	for(size_t i=HALO; i<grid.nx-HALO; i++)
	#endif
	#ifdef dij
	for(size_t j=HALO; j<grid.ny-HALO; j++)
	for(size_t i=HALO; i<grid.nx-HALO; i++)
	for(size_t d=0; d<stencil.nDir; d++)
	#endif
    {
	    size_t ij = grid.Index(i,j);    // Index of lattice point ij
	    size_t index = Index(ij,d);     // Index of population d at lattice point ij

	    // Get temp velocity:
	    Tensor2 direction = stencil.C(d);

	    // Get temp lattice point:
	    Coord xyTemp = grid.xy(i,j);
	    xyTemp[1] -= direction[1] * grid.dt;
	    xyTemp[2] -= direction[2] * grid.dt;

	    // Get 4 nearest Grid Points:
	    double iTemp = grid.i(xyTemp[1]);
	    double jTemp = grid.j(xyTemp[2]);
	    size_t i0 = std::floor(iTemp);	size_t i1 = i0 + 1;
	    size_t j0 = std::floor(jTemp);	size_t j1 = j0 + 1;
        
	    Inew[index] = BilinearInterpolation(iTemp-i0, jTemp-j0, I[Index(i0,j0,d)], I[Index(i0,j1,d)], I[Index(i1,j0,d)], I[Index(i1,j1,d)]);
    }
	std::swap(I,Inew);
}

void Radiation::StreamFlatAdaptive()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(3)
	#ifdef ijd
	for(size_t d=0; d<stencil.nDir; d++)
	for(size_t j=HALO; j<grid.ny-HALO; j++)
	for(size_t i=HALO; i<grid.nx-HALO; i++)
	#endif
	#ifdef dij
	for(size_t j=HALO; j<grid.ny-HALO; j++)
	for(size_t i=HALO; i<grid.nx-HALO; i++)
	for(size_t d=0; d<stencil.nDir; d++)
	#endif
    {
	    size_t ij = grid.Index(i,j);    // Index of lattice point ij
	    size_t index = Index(ij,d);     // Index of population d at lattice point ij

	    // Get temp velocity:
	    Tensor2 direction = RotationMatrix(rotationAngleNew[ij]) * stencil.C(d);

	    // Get temp lattice point:
	    Coord xyTemp = grid.xy(i,j);
	    xyTemp[1] -= direction[1] * grid.dt;
	    xyTemp[2] -= direction[2] * grid.dt;

	    // Get 4 nearest Grid Points:
	    double iTemp = grid.i(xyTemp[1]);
	    double jTemp = grid.j(xyTemp[2]);
	    size_t i0 = std::floor(iTemp);	size_t i1 = i0 + 1;
	    size_t j0 = std::floor(jTemp);	size_t j1 = j0 + 1;

		double intensityAt_i0j0 = IntensityAt(grid.Index(i0,j0), direction);
		double intensityAt_i0j1 = IntensityAt(grid.Index(i0,j1), direction);
		double intensityAt_i1j0 = IntensityAt(grid.Index(i1,j0), direction);
		double intensityAt_i1j1 = IntensityAt(grid.Index(i1,j1), direction);
		Inew[index] = BilinearInterpolation(iTemp-i0, jTemp-j0, intensityAt_i0j0, intensityAt_i0j1, intensityAt_i1j0, intensityAt_i1j1);
    }
	std::swap(I,Inew);
	std::swap(rotationAngle,rotationAngleNew);
}



void Radiation::StreamCurvedFixed()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(3)
	#ifdef ijd
	for(size_t d=0; d<stencil.nDir; d++)
	for(size_t j=HALO; j<grid.ny-HALO; j++)
	for(size_t i=HALO; i<grid.nx-HALO; i++)
	#endif
	#ifdef dij
	for(size_t j=HALO; j<grid.ny-HALO; j++)
	for(size_t i=HALO; i<grid.nx-HALO; i++)
	for(size_t d=0; d<stencil.nDir; d++)
	#endif
    {
        size_t ij = grid.Index(i,j);    // Index of lattice point ij
	    size_t index = Index(ij,d);     // Index of population d at lattice point ij

	    // Skip LPs which are inside BH:
	    if(metric.InsideBH(grid.xy(i,j)))
	    {
            Inew[index] = 0;
	    	continue;
	    }

	    // Get velocity direction in IF:
        double angle = stencil.Phi(d);
    
	    // Get quantities at emission point:
	    double s = GetFrequencyShift(ij, angle);
	    Coord xyTemp = GetTempCoordinate(ij, angle);
	    Tensor2 vTempIF = GetTemp2VelocityIF(ij, angle);

	    // Skip temporary Grid Points inside BH:
	    if(metric.InsideBH(xyTemp))
	    {
            Inew[index] = 0;
	    	continue;
	    }

	    // Get 4 nearest Grid Points:
	    double iTemp = grid.i(xyTemp[1]);
	    double jTemp = grid.j(xyTemp[2]);
	    size_t i0 = std::floor(iTemp);	size_t i1 = i0 + 1;
	    size_t j0 = std::floor(jTemp);	size_t j1 = j0 + 1;

	    // Intensity interpolation:
	    double alpha = metric.GetAlpha(ij);
	    double intensityAt_i0j0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j0))) * IntensityAt(grid.Index(i0,j0), vTempIF);
	    double intensityAt_i0j1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j1))) * IntensityAt(grid.Index(i0,j1), vTempIF);
	    double intensityAt_i1j0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j0))) * IntensityAt(grid.Index(i1,j0), vTempIF);
	    double intensityAt_i1j1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j1))) * IntensityAt(grid.Index(i1,j1), vTempIF);

	    // Interpolate intensity from neighbouring 4 lattice points to temporary point:
	    Inew[index] = IntegerPow<4>(s) * BilinearInterpolation(iTemp-i0, jTemp-j0, intensityAt_i0j0, intensityAt_i0j1, intensityAt_i1j0, intensityAt_i1j1);
    }
	std::swap(I,Inew);
}

void Radiation::StreamCurvedAdaptive()
{
	PROFILE_FUNCTION();
	PARALLEL_FOR(3)
	#ifdef ijd
	for(size_t d=0; d<stencil.nDir; d++)
	for(size_t j=HALO; j<grid.ny-HALO; j++)
	for(size_t i=HALO; i<grid.nx-HALO; i++)
	#endif
	#ifdef dij
	for(size_t j=HALO; j<grid.ny-HALO; j++)
	for(size_t i=HALO; i<grid.nx-HALO; i++)
	for(size_t d=0; d<stencil.nDir; d++)
	#endif
    {
        size_t ij = grid.Index(i,j);    // Index of lattice point ij
	    size_t index = Index(ij,d);     // Index of population d at lattice point ij

	    // Skip LPs which are inside BH:
	    if(metric.InsideBH(grid.xy(i,j)))
	    {
            Inew[index] = 0;
	    	continue;
	    }

	    // Get velocity direction in IF:
	    Tensor2 direction = RotationMatrix(rotationAngleNew[ij]) * stencil.C(d);
        double angle = direction.Phi();
        // double angle = fmod(rotationAngleNew[ij] + stencil.Phi(d) + 2.0 * M_PI, 2.0 * M_PI);
    
	    // Get quantities at emission point:
	    double s = GetFrequencyShift(ij, angle);
	    Coord xyTemp = GetTempCoordinate(ij, angle);
	    Tensor2 vTempIF = GetTemp2VelocityIF(ij, angle);

	    // Skip temporary Grid Points inside BH:
	    if(metric.InsideBH(xyTemp))
	    {
            Inew[index] = 0;
	    	continue;
	    }

	    // Get 4 nearest Grid Points:
	    double iTemp = grid.i(xyTemp[1]);
	    double jTemp = grid.j(xyTemp[2]);
	    size_t i0 = std::floor(iTemp);  size_t i1 = i0 + 1;
	    size_t j0 = std::floor(jTemp);	size_t j1 = j0 + 1;

	    // Intensity interpolation:
	    double alpha = metric.GetAlpha(ij);
	    double intensityAt_i0j0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j0))) * IntensityAt(grid.Index(i0,j0), vTempIF);
	    double intensityAt_i0j1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j1))) * IntensityAt(grid.Index(i0,j1), vTempIF);
	    double intensityAt_i1j0 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j0))) * IntensityAt(grid.Index(i1,j0), vTempIF);
	    double intensityAt_i1j1 = IntegerPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j1))) * IntensityAt(grid.Index(i1,j1), vTempIF);

	    // Interpolate intensity from neighbouring 4 lattice points to temporary point:
	    Inew[index] = IntegerPow<4>(s) * BilinearInterpolation(iTemp-i0, jTemp-j0, intensityAt_i0j0, intensityAt_i0j1, intensityAt_i1j0, intensityAt_i1j1);
    }
	std::swap(I,Inew);
	std::swap(rotationAngle,rotationAngleNew);
}



void Radiation::Collide()
{
	PROFILE_FUNCTION();
	// TODO: Steife DGL?
	
	PARALLEL_FOR(2)
	for(size_t j=HALO; j<metric.grid.ny-HALO; j++)
	for(size_t i=HALO; i<metric.grid.nx-HALO; i++)
	{
		if(metric.InsideBH(grid.xy(i,j)))
			continue;

		size_t ij = grid.Index(i,j);

		// Simulate stationary fluid, u^k=(0,0):
		double alpha = metric.GetAlpha(ij);
		Tensor2 u(0.0);	// fluid 2 velocity as seen by Eulerian observer
		double W = 1.0 / sqrt(1.0 - Norm2(u, metric.GetGamma_ll(ij)));	// Lorentz factor
		double uDotF = Fx_LF[ij] * u[1] + Fy_LF[ij] * u[2];             // F^i u_i
		double uuDotP = Pxx_LF[ij] * u[1] * u[1] + Pyy_LF[ij] * u[2] * u[2] + 2.0 * Pxy_LF[ij] * u[1] * u[2];   // P^ij u_i u_j
		double fluidE = W * W * (E_LF[ij] - 2.0 * uDotF + uuDotP);
		
		for(size_t d = 0; d < stencil.nDir; d++)
		{
			size_t index = Index(ij,d);
			double A = W * (1.0 - Tensor2::Dot(RotationMatrix(rotationAngle[ij]) * stencil.C(d), u));

			double Gamma = (eta[ij] + kappa0[ij]*fluidE) / (A*A*A) - A*I[index] * (kappaA[ij] + kappa0[ij]);
			I[index] += alpha * metric.grid.dt * Gamma;
		}
	}
}



void Radiation::WriteIntensitiesToCsv(double time, const int frameNumber, std::string directory, std::string name)
{
	PROFILE_FUNCTION();
    CreateDirectory(directory);

    name = name + FrameNumber(frameNumber) + ".csv";
    std::ofstream fileOut(directory + "/" + name);

    fileOut << "#nx=" << grid.nx << "\n";
    fileOut << "#ny=" << grid.ny << "\n";
    fileOut << "#startx=" << grid.startx << "\n";
    fileOut << "#starty=" << grid.starty << "\n";
    fileOut << "#endx=" << grid.endx << "\n";
    fileOut << "#endy=" << grid.endy << "\n";
	fileOut << "#x, y, z, color\n";
	for(size_t j = 2; j < metric.grid.ny - 2; j++)
	for(size_t i = 2; i < metric.grid.nx - 2; i++)
	{
		size_t ij = grid.Index(i,j);
		Coord xy = grid.xy(i,j);
		fileOut << xy[1] << ", " << xy[2] << ", " << 0 << ", " << 0 << "\n";
		// Bulk:
		for(size_t d=0; d<stencil.nDir; d++)
		{
			size_t index = Index(ij,d);
			if(I[index] > 1e-8)
			{
				Tensor2 dir = RotationMatrix(rotationAngle[ij]) * stencil.C(d);
                double value = I[index];

				Coord pos = xy;
				pos[1] += dir[1] * grid.dt * 0.5;// * (0.1 + 0.9 * value);
				pos[2] += dir[2] * grid.dt * 0.5;// * (0.1 + 0.9 * value);

				fileOut << pos[1] << ", " << pos[2] << ", " << 0 << ", " << value << "\n";
			}
		}
	}
	
    fileOut.close();
}



void Radiation::RunSimulation(Config config)
{
	// -------------------- Initialization --------------------
	int timeSteps = ceil(config.simTime / grid.dt);
	config.simTime = timeSteps * grid.dt;
	Log logger(config.name, config.simTime, stencil, streamingStencil, metric);
	NormalizeInitialDirections();
	LoadInitialData();
	ComputeMomentsIF();
	ComputeMomentsLF();
	NormalizeInitialIntensities();
	LoadInitialData(); // loads normalized data
	UpdateFourierCoefficients();
	
	// Initial data output:
	if (config.printToTerminal)
	{
		std::cout << " sigma        = " << sigma << "\n";
		std::cout << " nx           = " << grid.nx << "\n";
		std::cout << " ny           = " << grid.ny << "\n";
		std::cout << " nDir         = " << stencil.nDir << "\n";
		std::cout << " nFourier     = " << streamingStencil.nDir << "\n";
		std::cout << " simTime      = " << config.simTime << "\n";
		std::cout << " dx           = " << grid.dx << "\n";
		std::cout << " dy           = " << grid.dy << "\n";
		std::cout << " dt           = " << grid.dt << "\n";
		std::cout << " timeSteps    = " << timeSteps << "\n";
		std::cout << " filesToWrite = " << timeSteps / config.writeFrequency << std::endl;
	}
	// --------------------------------------------------------



	Profiler::Session& session = Profiler::Session::Get();
	session.Start(config.name, "output/" + config.name + "/profileResults.json");
	// ----------------- Main simulation Loop -----------------
	{
		PROFILE_SCOPE("Total Time");
		for(int n=0; n<timeSteps; n++)
		{
			if(config.printToTerminal)
			{ std::cout << n << "," << std::flush; }
			ComputeMomentsIF();
			Collide();

			if(config.writeData && (n % config.writeFrequency) == 0)
            {
				ComputeMomentsLF();
                grid.WriteFrametoCsv(n*grid.dt, E_LF, Fx_LF, Fy_LF, E, n, logger.directoryPath + "/Moments/");
                // WriteIntensitiesToCsv(n*grid.dt, n, logger.directoryPath + "/Intensities", "I");
            }

			if (config.updateFourierHarmonics)
			{ UpdateFourierCoefficients(); }

			// Streaming:
			switch(streamingType)
			{
				case(StreamingType::FlatFixed):         StreamFlatFixed();                                      break;
				case(StreamingType::FlatAdaptive):		UpdateRotationMatrizes();    StreamFlatAdaptive();	    break;
				case(StreamingType::CurvedFixed):		StreamCurvedFixed();                                    break;
				case(StreamingType::CurvedAdaptive):    UpdateRotationMatrizes();    StreamCurvedAdaptive();    break;
			}

			if(config.keepSourceNodesActive)
			{ LoadInitialData(); }
		}
	
		if(config.writeData)
        {
			ComputeMomentsIF();
			ComputeMomentsLF();
        }
		if(config.writeData)
        {
            grid.WriteFrametoCsv(timeSteps*grid.dt, E_LF, Fx_LF, Fy_LF, E, timeSteps, logger.directoryPath + "/Moments/");
		    // WriteIntensitiesToCsv(timeSteps*grid.dt, timeSteps, logger.directoryPath + "/Intensities", "I");
        }
	}
	// --------------------------------------------------------
	session.End();


	// ---------------------- Termination ---------------------
	std::vector<std::string> names = session.GetAllFunctionNames();
	if (config.printToTerminal)
	{ std::cout << std::endl; }
	for(int i=0; i<names.size(); i++)
	{
		if (config.printToTerminal)
			session.PrintFunctionDuration(names[i]);
		logger.AddTimeMeasurement(names[i], session.GetTotalTime(names[i]));
	}
	// --------------------------------------------------------
}