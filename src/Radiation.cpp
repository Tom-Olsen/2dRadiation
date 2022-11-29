#include "Radiation.h"

std::string StreamingName(int n)
{
	std::string s("unknown");
	switch (n)
	{
   		case 0: { s = "StreamCurvedRot";   } break;
   		case 1: { s = "StreamCurvedNoRot"; } break;
   		case 2: { s = "StreamFlatRot";     } break;
   		case 3: { s = "StreamFlatNoRot";   } break;
	}
	return s;
}



template<class Coord>
Radiation<Coord>::Radiation(Grid2D<Coord>& grid_, Metric2D<Coord>& metric_, Stencil& stencil_, Stencil& fourierStencil_, StreamingType streamingType_):
grid(grid_), metric(metric_), stencil(stencil_), fourierStencil(fourierStencil_), nDir(stencil_.nDir), streamingType(streamingType_)
{
	isInitialGridPoint = new bool[grid.n12]();
	initialE           = new double[grid.n12]();
	initialRotation    = new double[grid.n12]();
	initialKappa0      = new double[grid.n12]();
	initialKappa1      = new double[grid.n12]();
	initialKappaA      = new double[grid.n12]();
	initialEta         = new double[grid.n12]();
	rotation    = new double[grid.n12]();
	rotationNew = new double[grid.n12]();
	E           = new double[grid.n12]();
	Fx          = new double[grid.n12]();
	Fy          = new double[grid.n12]();
	Pxx         = new double[grid.n12]();
	Pxy         = new double[grid.n12]();
	Pyy         = new double[grid.n12]();
	E_LF        = new double[grid.n12]();
	Fx_LF       = new double[grid.n12]();
	Fy_LF       = new double[grid.n12]();
	Pxx_LF      = new double[grid.n12]();
	Pxy_LF      = new double[grid.n12]();
	Pyy_LF      = new double[grid.n12]();
	Fnorm_LF    = new double[grid.n12]();
	kappa0      = new double[grid.n12]();
	kappa1      = new double[grid.n12]();
	kappaA      = new double[grid.n12]();
	eta         = new double[grid.n12]();
	I           = new double[grid.n12 * nDir]();
	temp        = new double[grid.n12 * nDir]();
	coefficientsS   = new double[grid.n12 * fourierStencil.nDir];
	coefficientsX1  = new double[grid.n12 * fourierStencil.nDir];
	coefficientsX2  = new double[grid.n12 * fourierStencil.nDir];
	coefficientsCx1 = new double[grid.n12 * fourierStencil.nDir];
	coefficientsCx2 = new double[grid.n12 * fourierStencil.nDir];
}



template<class Coord>
Radiation<Coord>::~Radiation()
{
	delete[] isInitialGridPoint;
	delete[] initialE;
	delete[] initialRotation;
	delete[] initialKappa0;
	delete[] initialKappa1;
	delete[] initialKappaA;
	delete[] initialEta;
	delete[] rotation;
	delete[] rotationNew;
	delete[] E;
	delete[] Fx;
	delete[] Fy;
	delete[] Pxx;
	delete[] Pxy;
	delete[] Pyy;
	delete[] E_LF;
	delete[] Fx_LF;
	delete[] Fy_LF;
	delete[] Pxx_LF;
	delete[] Pxy_LF;
	delete[] Pyy_LF;
	delete[] Fnorm_LF;
	delete[] kappa0;
	delete[] kappa1;
	delete[] kappaA;
	delete[] eta;
	delete[] I;
	delete[] temp;
	delete[] coefficientsS;
	delete[] coefficientsX1;
	delete[] coefficientsX2;
	delete[] coefficientsCx1;
	delete[] coefficientsCx2;
}



template<class Coord>
void Radiation<Coord>::LoadInitialData()
{
	PROFILE_FUNCTION();
	auto IntensityDistribution = [](double phi, double sigma)
	{
		return exp(-0.5 * pow((phi) / sigma, 2))
		     + exp(-0.5 * pow((phi - 2.0*M_PI) / sigma, 2))
		     + exp(-0.5 * pow((phi + 2.0*M_PI) / sigma, 2));
	};

	// Note:
	// For the initial data the static stencil looks like it is being rotated, while the rotating stencil looks like being static.
	// By using phi=stencil.phi(d,-rot), for the static stencil, the gauß peak of the distribution is maximal (phi=0) when phi(d)=rot.
	// On the other hand for the rotating stencil we initially set the peak to be at phi=stencil.phi(d=0).
	// When running the code that direction vector is rotated to point toward rotation[ij].
	if(streamingType == StreamingType::StreamCurvedRot || streamingType == StreamingType::StreamFlatRot)
	{
		#pragma omp parallel for
		for(int ij=0; ij<grid.n12; ij++)
		{
			if(isInitialGridPoint[ij])
			{
				rotation[ij] = initialRotation[ij];
				kappa0[ij]   = initialKappa0[ij];
				kappa1[ij]   = initialKappa1[ij];
				kappaA[ij]   = initialKappaA[ij];
				eta[ij]      = initialEta[ij];

				if(rotation[ij] >= 0)
				{
					for(int d=0; d<nDir; d++)
						I[ij + d*grid.n12] = initialE[ij] * IntensityDistribution(stencil.Phi(d), stencil.sigma);
				}
				else
				{
					for(int d=0; d<nDir; d++)
						I[ij + d*grid.n12] = initialE[ij];
				}
			}
		}
	}
	else
	{
		#pragma omp parallel for
		for(int ij=0; ij<grid.n12; ij++)
		{
			if(isInitialGridPoint[ij])
			{
				rotation[ij] = initialRotation[ij];
				kappa0[ij]   = initialKappa0[ij];
				kappa1[ij]   = initialKappa1[ij];
				kappaA[ij]   = initialKappaA[ij];
				eta[ij]      = initialEta[ij];

				if(rotation[ij] >= 0)
				{
					for(int d=0; d<nDir; d++)
						I[ij + d*grid.n12] = initialE[ij] * IntensityDistribution(stencil.Phi(d, -rotation[ij]), stencil.sigma);
				}
				else
				{
					for(int d=0; d<nDir; d++)
						I[ij + d*grid.n12] = initialE[ij];
				}
			}
		}
	}
}



template<class Coord>
void Radiation<Coord>::NormalizeInitialData()
{
	PROFILE_FUNCTION();
	#pragma omp parallel for
	for(int ij=0; ij<grid.n12; ij++)
		if(isInitialGridPoint[ij])
			if(E_LF[ij] > 0)
				initialE[ij] *= initialE[ij] / E_LF[ij];
}



template<class Coord>
void Radiation<Coord>::UpdateFourierCoefficients()
{
	PROFILE_FUNCTION();
	#pragma omp parallel for
	for(int j=2; j<grid.n2-2; j++)
	for(int i=2; i<grid.n1-2; i++)
	{
		int ij = grid.Index(i,j);
		double dataS  [fourierStencil.nDir];
		double dataX1 [fourierStencil.nDir];
		double dataX2 [fourierStencil.nDir];
    	double dataCx1[fourierStencil.nDir];
    	double dataCx2[fourierStencil.nDir];
		Coordinate2<Coord> xStart = grid.x12Coord(i,j);
		Coordinate2<xy> xyCoord = grid.xyCoord(i,j);
		double alpha = metric.GetAlpha(ij);
		for(int d=0; d<fourierStencil.nDir; d++)
    	{
			// Initial data for geodesic equation:
			double s = 1;
			Coordinate2<Coord> x = xStart;
            Tensor2<xy,IF> cxy = fourierStencil.Cxy(d);
            Tensor2<Coord,IF> c = cxy.template Transform<Coord>(xyCoord);
            Tensor3<Coord,IF> u(alpha, c[1]*alpha, c[2]*alpha);
            Tensor2<Coord,LF> v = Vec2ObservedByEulObs<Coord,IF,LF>(u,x,metric);

			// Solve geodesic equation backwards:
			if(!metric.InsideBH(x))
			{
    	    	// s *= Euler_GeodesicEquation<-1>(grid.dt,x,v,metric);
    	    	s *= RK45_GeodesicEquation<-1>(grid.dt,x,v,metric);
			}
			else // inside BH tetrad destroys the velocity stencil. Thus set it to 0.
				v = Tensor2<Coord,LF>(0.0);

			// Final data points for fourier expansion:
			dataS[d]   = 1.0/s;
			dataX1[d]  = x[1];
			dataX2[d]  = x[2];
			dataCx1[d] = v[1];
			dataCx2[d] = v[2];
    	}
		std::vector<double> cS  = Fourier::Expansion::GetCoefficients(fourierStencil,dataS);
		std::vector<double> cX1 = Fourier::Expansion::GetCoefficients(fourierStencil,dataX1);
		std::vector<double> cX2 = Fourier::Expansion::GetCoefficients(fourierStencil,dataX2);
		std::vector<double> cCx1= Fourier::Expansion::GetCoefficients(fourierStencil,dataCx1);
		std::vector<double> cCx2= Fourier::Expansion::GetCoefficients(fourierStencil,dataCx2);
		for(int f=0; f<fourierStencil.nDir; f++)
		{
			int fij = f + i*fourierStencil.nDir + j*fourierStencil.nDir*grid.n1;
			coefficientsS  [fij] = cS  [f];
			coefficientsX1 [fij] = cX1 [f];
			coefficientsX2 [fij] = cX2 [f];
			coefficientsCx1[fij] = cCx1[f];
			coefficientsCx2[fij] = cCx2[f];
		}
	}
}



template<class Coord>
void Radiation<Coord>::ComputeMomentsIF()
{
	PROFILE_FUNCTION();
	constexpr double twoPiInv = 1.0 / (2.0 * M_PI);
	
	if(streamingType == StreamingType::StreamCurvedRot || streamingType == StreamingType::StreamFlatRot)
	{
		#pragma omp parallel for
		for(int ij=0; ij<grid.n12; ij++)
		{
			E [ij]  = 0.0;
			Fx[ij]  = 0.0;
			Fy[ij]  = 0.0;
			Pxx[ij] = 0.0;
			Pxy[ij] = 0.0;
			Pyy[ij] = 0.0;
			for(int d = 0; d < nDir; d++)
			{
				E [ij]  += stencil.W(d) * I[ij + d*grid.n12];
				Fx[ij]  += stencil.W(d) * I[ij + d*grid.n12] * stencil.Cx(d,rotation[ij]);
				Fy[ij]  += stencil.W(d) * I[ij + d*grid.n12] * stencil.Cy(d,rotation[ij]);
				Pxx[ij] += stencil.W(d) * I[ij + d*grid.n12] * stencil.Cx(d,rotation[ij]) * stencil.Cx(d,rotation[ij]);
				Pxy[ij] += stencil.W(d) * I[ij + d*grid.n12] * stencil.Cx(d,rotation[ij]) * stencil.Cy(d,rotation[ij]);
				Pyy[ij] += stencil.W(d) * I[ij + d*grid.n12] * stencil.Cy(d,rotation[ij]) * stencil.Cy(d,rotation[ij]);
			}
    	    E [ij]  *= twoPiInv;
    	    Fx[ij]  *= twoPiInv;
    	    Fy[ij]  *= twoPiInv;
    	    Pxx[ij] *= twoPiInv;
    	    Pxy[ij] *= twoPiInv;
    	    Pyy[ij] *= twoPiInv;
		}
	}
	else
	{
		#pragma omp parallel for
		for(int ij=0; ij<grid.n12; ij++)
		{
			E [ij]  = 0.0;
			Fx[ij]  = 0.0;
			Fy[ij]  = 0.0;
			Pxx[ij] = 0.0;
			Pxy[ij] = 0.0;
			Pyy[ij] = 0.0;
			for(int d = 0; d < nDir; d++)
			{
				E [ij]  += stencil.W(d) * I[ij + d*grid.n12];
				Fx[ij]  += stencil.W(d) * I[ij + d*grid.n12] * stencil.Cx(d);
				Fy[ij]  += stencil.W(d) * I[ij + d*grid.n12] * stencil.Cy(d);
				Pxx[ij] += stencil.W(d) * I[ij + d*grid.n12] * stencil.Cx(d) * stencil.Cx(d);
				Pxy[ij] += stencil.W(d) * I[ij + d*grid.n12] * stencil.Cx(d) * stencil.Cy(d);
				Pyy[ij] += stencil.W(d) * I[ij + d*grid.n12] * stencil.Cy(d) * stencil.Cy(d);
			}
    	    E [ij]  *= twoPiInv;
    	    Fx[ij]  *= twoPiInv;
    	    Fy[ij]  *= twoPiInv;
    	    Pxx[ij] *= twoPiInv;
    	    Pxy[ij] *= twoPiInv;
    	    Pyy[ij] *= twoPiInv;
		}
	}
}

template<class Coord>
void Radiation<Coord>::ComputeMomentsLF()
{
	PROFILE_FUNCTION();
	#pragma omp parallel for
	for(int ij=0; ij<grid.n12; ij++)
	{
		Tensor3x3<Coord,IF> EnergyMomentumTensorIF(E[ij],Fx[ij],Fy[ij], Fx[ij],Pxx[ij],Pxy[ij], Fy[ij],Pxy[ij],Pyy[ij]);
		Tensor3x3<Coord,LF> EnergyMomentumTensorLF = EnergyMomentumTensorIF.template Transform<LF>(metric.GetTetrad(ij));
		E_LF[ij] = EnergyMomentumTensorLF[{0,0}];
		Fx_LF[ij] = EnergyMomentumTensorLF[{0,1}];
		Fy_LF[ij] = EnergyMomentumTensorLF[{0,2}];
		Pxx_LF[ij] = EnergyMomentumTensorLF[{1,1}];
		Pxy_LF[ij] = EnergyMomentumTensorLF[{1,2}];
		Pyy_LF[ij] = EnergyMomentumTensorLF[{2,2}];
		Fnorm_LF[ij] = Tensor2<Coord,LF>(Fx_LF[ij],Fy_LF[ij]).Norm(metric.GetGamma_ll(ij));
	}
}



template<class Coord>
inline __attribute__((always_inline)) Coordinate2<Coord> Radiation<Coord>::GetTempCoordinate(int i, int j, double phi)
{
	std::vector<double> cX1(fourierStencil.nDir);
	std::vector<double> cX2(fourierStencil.nDir);
	for(int f=0; f<fourierStencil.nDir; f++)
	{
		cX1[f] = coefficientsX1[f + i*fourierStencil.nDir + j*fourierStencil.nDir*grid.n1];
		cX2[f] = coefficientsX2[f + i*fourierStencil.nDir + j*fourierStencil.nDir*grid.n1];
	}
	return Coordinate2<Coord>(Fourier::Expansion::GetValue(phi,cX1), Fourier::Expansion::GetValue(phi,cX2));
}
template<class Coord>
inline __attribute__((always_inline)) Tensor2<Coord,LF> Radiation<Coord>::GetTemp2Velocity(int i, int j, double phi)
{
	std::vector<double> cCx1(fourierStencil.nDir);
	std::vector<double> cCx2(fourierStencil.nDir);
	for(int f=0; f<fourierStencil.nDir; f++)
	{
		cCx1[f] = coefficientsCx1[f + i*fourierStencil.nDir + j*fourierStencil.nDir*grid.n1];
		cCx2[f] = coefficientsCx2[f + i*fourierStencil.nDir + j*fourierStencil.nDir*grid.n1];
	}
	return Tensor2<Coord,LF>(Fourier::Expansion::GetValue(phi,cCx1), Fourier::Expansion::GetValue(phi,cCx2));
}
template<class Coord>
inline __attribute__((always_inline)) double Radiation<Coord>::GetFrequencyShift(int i, int j, double phi)
{
	std::vector<double> cS(fourierStencil.nDir);
	for(int f=0; f<fourierStencil.nDir; f++)
		cS[f] = coefficientsS[f + i*fourierStencil.nDir + j*fourierStencil.nDir*grid.n1];
	return Fourier::Expansion::GetValue(phi,cS);
}
template<class Coord>
inline __attribute__((always_inline)) double Radiation<Coord>::IntensityAt(int ij, double phi)
{
	double k = stencil.Index(phi,rotation[ij]);

	int kFloor = floor(k);
	int kM1 = (kFloor - 1 + nDir) % nDir;
	int kP0 = (kFloor + 0 + nDir) % nDir;
	int kP1 = (kFloor + 1 + nDir) % nDir;
	int kP2 = (kFloor + 2 + nDir) % nDir;

	double value = CubicInterpolation(k - kFloor,I[ij + kM1*grid.n12],I[ij + kP0*grid.n12],I[ij + kP1*grid.n12],I[ij + kP2*grid.n12]);
	// double value = LinearInterpolation(k - kFloor,I[ij + kP0*grid.n12],I[ij + kP1*grid.n12]);
	return std::max(value,0.0);
}
template<class Coord>
inline __attribute__((always_inline)) Tensor2<Coord,IF> Radiation<Coord>::AverageF(int i, int j)
{
	Tensor2<Coord,IF> averageF(0.0);
	for(int a=-1; a<=1; a++)
		for(int b=-1; b<=1; b++)
		{
			int index = grid.Index(i+a,j+b);
			averageF[1] += Fx[index];
			averageF[2] += Fy[index];
		}
	averageF[1] /= 9.0;
	averageF[2] /= 9.0;
	return averageF;
}



template<class Coord>
void Radiation<Coord>::StreamCurvedRot()
{
	PROFILE_FUNCTION();
	#pragma omp parallel for
	for(int j=2; j<grid.n2-2; j++)
	for(int i=2; i<grid.n1-2; i++)
	{
		int ij = grid.Index(i,j);
		Tensor2<Coord,IF> averageF = AverageF(i,j);
		constexpr double rotationThreshold = 1e-3;
		rotationNew[ij] = (averageF.EuklNorm()>rotationThreshold) ? averageF.Angle() : rotation[ij];
	}

	#pragma omp parallel for
	for(int d=0; d<nDir; d++)
	for(int j=2; j<grid.n2-2; j++)
	for(int i=2; i<grid.n1-2; i++)
	{
		int ij = grid.Index(i,j);		// Index of lattice point ij
		int ijd = grid.Index(i,j,d);	// Index of population d at lattice point ij

		// Skip LPs which are inside BH:
		if(metric.InsideBH(i,j))
		{
			temp[ijd] = 0;
			continue;
		}

		// Curved Fourier Streaming:
		double phi = stencil.Phi(d,rotationNew[ij]);
		double s = GetFrequencyShift(i,j,phi);
		Coordinate2<Coord> xTemp = GetTempCoordinate(i,j,phi);
		Tensor2<Coord,LF> vTempLF = GetTemp2Velocity(i,j,phi);
		Tensor2<Coord,IF> vTempIF = vTempLF.template Transform<IF>(metric.GetTetrad(xTemp));
		Tensor2<xy,IF> vTempIFxy = vTempIF.template Transform<xy>(xTemp);

		// Skip temporary Grid Points inside BH:
		if(metric.InsideBH(xTemp))
		{
			temp[ijd] = 0;
			continue;
		}

		// Get 4 nearest Grid Points:
		double iTemp = grid.i(xTemp[1]);
		double jTemp = grid.j(xTemp[2]);
		int i0 = std::floor(iTemp);	int i1 = i0 + 1;
		int j0 = std::floor(jTemp);	int j1 = j0 + 1;

		// Intensity interpolation:
		double angle = vTempIFxy.Angle();
		double alpha = metric.GetAlpha(ij);
		double intensityAt_i0j0 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j0))) * IntensityAt(grid.Index(i0,j0),angle);
		double intensityAt_i0j1 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j1))) * IntensityAt(grid.Index(i0,j1),angle);
		double intensityAt_i1j0 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j0))) * IntensityAt(grid.Index(i1,j0),angle);
		double intensityAt_i1j1 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j1))) * IntensityAt(grid.Index(i1,j1),angle);

		// Interpolate intensity from neighbouring 4 lattice points to temporary point:
		temp[ijd]
		= (s*s*s*s) * BilinearInterpolation(iTemp-i0,jTemp-j0,
		 intensityAt_i0j0,intensityAt_i0j1,
		 intensityAt_i1j0,intensityAt_i1j1);
	}
	// Copy rotationNew to rotation:
	#pragma omp parallel for
	for(int ij=0; ij<grid.n12; ij++)
		rotation[ij] = rotationNew[ij];	
	// Copy temp to intensities:
	#pragma omp parallel for
	for(int ijd=0; ijd<grid.n12*nDir; ijd++)
		I[ijd] = temp[ijd];
}
template<class Coord>
void Radiation<Coord>::StreamCurvedNoRot()
{
	PROFILE_FUNCTION();
	#pragma omp parallel for
	for(int d=0; d<nDir; d++)
	for(int j=2; j<grid.n2-2; j++)
	for(int i=2; i<grid.n1-2; i++)
	{
		int ij = grid.Index(i,j);		// Index of lattice point ij
		int ijd = grid.Index(i,j,d);	// Index of population d at lattice point ij

		// Skip LPs which are inside BH:
		if(metric.InsideBH(i,j))
		{
			temp[ijd] = 0;
			continue;
		}

		// Curved Fourier Streaming:
		double phi = stencil.Phi(d);
		double s = GetFrequencyShift(i,j,phi);
		Coordinate2<Coord> xTemp = GetTempCoordinate(i,j,phi);
		Tensor2<Coord,LF> vTempLF = GetTemp2Velocity(i,j,phi);
		Tensor2<Coord,IF> vTempIF = vTempLF.template Transform<IF>(metric.GetTetrad(xTemp));
		Tensor2<xy,IF> vTempIFxy = vTempIF.template Transform<xy>(xTemp);

		// Skip temporary Grid Points inside BH:
		if(metric.InsideBH(xTemp))
		{
			temp[ijd] = 0;
			continue;
		}

		// Get 4 nearest Grid Points:
		double iTemp = grid.i(xTemp[1]);
		double jTemp = grid.j(xTemp[2]);
		int i0 = std::floor(iTemp);	int i1 = i0 + 1;
		int j0 = std::floor(jTemp);	int j1 = j0 + 1;

		// Intensity interpolation:
		double angle = vTempIFxy.Angle();
		double alpha = metric.GetAlpha(ij);
		double intensityAt_i0j0 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j0))) * IntensityAt(grid.Index(i0,j0),angle);
		double intensityAt_i0j1 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i0,j1))) * IntensityAt(grid.Index(i0,j1),angle);
		double intensityAt_i1j0 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j0))) * IntensityAt(grid.Index(i1,j0),angle);
		double intensityAt_i1j1 = MyPow<4>(alpha / metric.GetAlpha(grid.Index(i1,j1))) * IntensityAt(grid.Index(i1,j1),angle);

		// Interpolate intensity from neighbouring 4 lattice points to temporary point:
		temp[ijd]
		= (s*s*s*s) * BilinearInterpolation(iTemp-i0,jTemp-j0,
		 intensityAt_i0j0,intensityAt_i0j1,
		 intensityAt_i1j0,intensityAt_i1j1);
	}
	// Copy temp to intensities:
	#pragma omp parallel for
	for(int ijd=0; ijd<grid.n12*nDir; ijd++)
		I[ijd] = temp[ijd];
}



template<class Coord>
void Radiation<Coord>::StreamFlatRot()
{
	PROFILE_FUNCTION();
	#pragma omp parallel for
	for(int d=0; d<nDir; d++)
	for(int j=2; j<grid.n2-2; j++)
	for(int i=2; i<grid.n1-2; i++)
	{
		int ij = grid.Index(i,j);		// Index of lattice point ij
		int ijd = grid.Index(i,j,d);	// Index of population d at lattice point ij

		// Find new rotation:
		Tensor2<Coord,IF> averageF = AverageF(i,j);
		constexpr double rotationThreshold = 1e-3;
		rotationNew[ij] = (averageF.EuklNorm()>rotationThreshold) ? averageF.Angle() : rotation[ij];

		// Flat Streaming:
		Tensor2<xy,IF> vTempIFxy = stencil.Cxy(d,rotationNew[ij]);
		Coordinate2<Coord> xTemp = grid.x12Coord(i,j);
		xTemp[1] -= vTempIFxy[1] * grid.dt;
		xTemp[2] -= vTempIFxy[2] * grid.dt;

		// Get 4 nearest Grid Points:
		double iTemp = grid.i(xTemp[1]);
		double jTemp = grid.j(xTemp[2]);
		int i0 = std::floor(iTemp);	int i1 = i0 + 1;
		int j0 = std::floor(jTemp);	int j1 = j0 + 1;

		// Interpolate intensity from neighbouring 4 lattice points to temporary point:
		double angle = stencil.Phi(d,rotationNew[ij]);
		temp[ijd]
		= BilinearInterpolation(iTemp-i0,jTemp-j0,
		 IntensityAt(grid.Index(i0,j0),angle),IntensityAt(grid.Index(i0,j1),angle),
		 IntensityAt(grid.Index(i1,j0),angle),IntensityAt(grid.Index(i1,j1),angle));
	}
	// Copy rotationNew to rotation:
	#pragma omp parallel for
	for(int ij=0; ij<grid.n12; ij++)
		rotation[ij] = rotationNew[ij];	
	// Copy temp to intensities:
	#pragma omp parallel for
	for(int ijd=0; ijd<grid.n12*nDir; ijd++)
		I[ijd] = temp[ijd];
}
template<class Coord>
void Radiation<Coord>::StreamFlatNoRot()
{
	PROFILE_FUNCTION();
	#pragma omp parallel for
	for(int d=0; d<nDir; d++)
	for(int j=2; j<grid.n2-2; j++)
	for(int i=2; i<grid.n1-2; i++)
	{
		int ij = grid.Index(i,j);		// Index of lattice point ij
		int ijd = grid.Index(i,j,d);	// Index of population d at lattice point ij

		// Flat Streaming:
		Tensor2<xy,IF> vTempIFxy = stencil.Cxy(d);
		Coordinate2<Coord> xTemp = grid.x12Coord(i,j);
		xTemp[1] -= vTempIFxy[1] * grid.dt;
		xTemp[2] -= vTempIFxy[2] * grid.dt;

		// Get 4 nearest Grid Points:
		double iTemp = grid.i(xTemp[1]);
		double jTemp = grid.j(xTemp[2]);
		int i0 = std::floor(iTemp);	int i1 = i0 + 1;
		int j0 = std::floor(jTemp);	int j1 = j0 + 1;

		// Interpolate intensity from neighbouring 4 lattice points to temporary point:
		temp[ijd]
		= BilinearInterpolation(iTemp-i0,jTemp-j0,
		 I[grid.Index(i0,j0,d)],I[grid.Index(i0,j1,d)],
		 I[grid.Index(i1,j0,d)],I[grid.Index(i1,j1,d)]);
	}
	// Copy temp to intensities:
	#pragma omp parallel for
	for(int ijd=0; ijd<grid.n12*nDir; ijd++)
		I[ijd] = temp[ijd];
}



template<class Coord>
void Radiation<Coord>::Collide()
{
	PROFILE_FUNCTION();
//	// TODO:
//	// Steife DGL?
//	
//	#pragma omp parallel for
//	for(int j = 0; j < metric.grid.n2; j++)
//	for(int i = 0; i < metric.grid.n1; i++)
//	{
//		Tensor3 x = metric.grid.xyzCoord(i,j,1);
//		if(metric.ToCloseToBH(x))
//			continue;
//
//		int ij = i+j*metric.grid.n1;
//
//		// Select correct stencil:
//		Stencil& stencil = rotation[ij] < 0 ? uniStencil : dirStencil;
//		
//		double Gamma = eta[ij];
//		for(int d = 0; d < velSt.nDir; d++)
//		{
//			int ijd = ij + d*metric.grid.n12;
//			Gamma += -kappaA[ij] * I[ijd] + kappa0[ij]*(E[ij]-I[ijd]) + 3*(velSt.vx[d]*Fx[ij] + velSt.vy[d]*Fy[ij]);
//			I[ijd] += metric.GetAlpha(x) * Gamma * metric.grid.dt;
//		}
//	}
//	
}



template<class Coord>
void Radiation<Coord>::RunSimulation(RunParameters params)
{
	// -------------------- Initialization --------------------
	int timeSteps = ceil(params.simTime/grid.dt);
	params.simTime = timeSteps * grid.dt;
	Log logger(params.name,params.simTime,stencil,fourierStencil,metric);

	LoadInitialData();
	ComputeMomentsIF();
	ComputeMomentsLF();
	NormalizeInitialData();
	LoadInitialData(); // loads normalized data into rotated stencil.
	UpdateFourierCoefficients();

	// Initial data output:
	if (params.printToTerminal)
	{
		std::cout << " sigma=" << stencil.sigma << "\n";
		std::cout << " n1=" << grid.n1 << "\n";
		std::cout << " n2=" << grid.n2 << "\n";
		std::cout << " nDir=" << stencil.nDir << "\n";
		std::cout << " nFourier=" << fourierStencil.nDir << "\n";
		std::cout << " simTime=" << params.simTime << "\n";
		std::cout << " dt=" << grid.dt << "\n";
		std::cout << " timeSteps=" << timeSteps << "\n";
		std::cout << " filesToWrite=" << timeSteps/params.writeFrequency << std::endl;
	}
	// --------------------------------------------------------



	Profiler::Session& session = Profiler::Session::Get();
	session.Start(params.name, "output/" + params.name + "/profileResults.json");
	// ----------------- Main simulation Loop -----------------
	{
		PROFILE_SCOPE("Total Time");
		for(int n=0; n<timeSteps; n++)
		{
			if (params.printToTerminal)
			{ std::cout << n << "," << std::flush; }
			ComputeMomentsIF();
			Collide();

			if (params.writeData && (n % params.writeFrequency) == 0)
			{
				ComputeMomentsLF();
				grid.WriteFrametoJson(n*grid.dt,E_LF,Fx_LF,Fy_LF,Fnorm_LF,n,logger.directoryPath + "/Moments");
			}

			if (params.updateFourierHarmonics)
			{ UpdateFourierCoefficients(); }

			// Streaming:
			switch (streamingType)
			{
				case(StreamingType::StreamCurvedRot):	StreamCurvedRot();		break;
				case(StreamingType::StreamCurvedNoRot):	StreamCurvedNoRot();	break;
				case(StreamingType::StreamFlatRot):		StreamFlatRot();		break;
				case(StreamingType::StreamFlatNoRot):	StreamFlatNoRot();		break;
			}

			if (params.keepSourceNodesActive)
			{ LoadInitialData(); }
		}
	
		if (params.writeData)
		{ grid.WriteFrametoJson(timeSteps*grid.dt,E_LF,Fx_LF,Fy_LF,Fnorm_LF,timeSteps,logger.directoryPath + "/Moments"); }
	}
	// --------------------------------------------------------
	session.End();



	// ---------------------- Termination ---------------------
	std::vector<std::string> names = session.GetAllFunctionNames();
	if (params.printToTerminal)
	{ std::cout << std::endl; }
	for(int i=0; i<names.size(); i++)
	{
		if (params.printToTerminal)
			session.PrintFunctionDuration(names[i]);
		logger.AddTimeMeasurement(names[i], session.GetTotalTime(names[i]));
	}
	// --------------------------------------------------------
}



template class Radiation<xy>;
template class Radiation<rph>;